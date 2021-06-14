from . import connection_nodes as connection_nodes_module
from . import cross_sections as cross_sections_module
from dataclasses import dataclass
from dataclasses import fields
from threedigrid_builder.base import Lines
from threedigrid_builder.base import Nodes
from threedigrid_builder.base import Pumps
from threedigrid_builder.base.settings import GridSettings
from threedigrid_builder.base.settings import TablesSettings
from threedigrid_builder.constants import CalculationType
from threedigrid_builder.constants import ContentType
from threedigrid_builder.constants import LineType
from threedigrid_builder.constants import NodeType
from typing import Optional
from typing import Tuple

import itertools
import numpy as np
import pygeos
import threedigrid_builder


__all__ = ["Grid", "GridMeta", "QuadtreeStats"]


@dataclass
class GridMeta:
    """Metadata that needs to end up in the gridadmin file."""

    epsg_code: int
    model_name: str  # name from sqlite globalsettings.name

    grid_settings: GridSettings
    tables_settings: TablesSettings

    model_slug: Optional[str] = None  # from repository.slug
    revision_hash: Optional[str] = None  # from repository.revision.hash
    revision_nr: Optional[int] = None  # from repository.revision.number
    threedi_version: Optional[str] = None  # threedi-api version
    threedicore_version: Optional[str] = None  # not used anymore
    threedigrid_builder_version: str = ""  # filled in __post_init__
    threedi_tables_version: Optional[str] = None  # filled in threedi-tables

    # TODO what to do with use_1d_flow, use_2d_flow, manhole_storage_area
    has_1d: bool = False
    has_2d: bool = False
    has_breaches: bool = False
    has_groundwater: bool = False
    has_groundwater_flow: bool = False
    has_interception: bool = False
    has_pumpstations: bool = False
    has_simple_infiltration: bool = False
    has_interflow: bool = False

    extent_1d: Optional[Tuple[float, float, float, float]] = None
    extent_2d: Optional[Tuple[float, float, float, float]] = None

    @classmethod
    def from_dict(cls, dct):
        """Construct skipping unknown fields and None values"""
        class_fields = {f.name for f in fields(cls)}
        return cls(
            **{k: v for k, v in dct.items() if k in class_fields and v is not None}
        )

    def __post_init__(self):
        if not self.threedigrid_builder_version:
            self.threedigrid_builder_version = threedigrid_builder.__version__


@dataclass
class QuadtreeStats:
    lgrmin: int
    kmax: int
    mmax: Tuple[int, ...]
    nmax: Tuple[int, ...]
    dx: Tuple[float, ...]
    dxp: float
    x0p: float
    y0p: float


class Grid:
    def __init__(
        self,
        nodes: Nodes,
        lines: Lines,
        pumps: Pumps = None,
        meta=None,
        quadtree_stats=None,
    ):
        if not isinstance(nodes, Nodes):
            raise TypeError(f"Expected Nodes instance, got {type(nodes)}")
        if not isinstance(lines, Lines):
            raise TypeError(f"Expected Lines instance, got {type(lines)}")
        self.nodes = nodes
        self.lines = lines
        self.meta = meta
        self.quadtree_stats = quadtree_stats
        self.pumps = pumps

    def __add__(self, other):
        """Concatenate two grids without renumbering nodes."""
        if self.__class__ is not other.__class__:
            raise TypeError(
                "Cannot concatenate {self} with {other} as they are not of "
                "equal types."
            )
        new_attrs = {}
        for k, v in other.__dict__.items():
            if isinstance(v, Nodes) or isinstance(v, Lines):
                new_attrs[k] = getattr(self, k) + v
            else:
                new_attrs[k] = v if v is not None else getattr(self, k)
        return self.__class__(**new_attrs)

    def __repr__(self):
        return f"<Grid object with {len(self.nodes)} nodes and {len(self.lines)} lines>"

    @classmethod
    def from_meta(cls, **kwargs):
        """Construct a grid with only metadata"""
        meta = GridMeta.from_dict(kwargs)
        # set flags
        s = meta.tables_settings  # shorthand
        meta.has_interception = s.interception_type is not None
        meta.has_groundwater_flow = s.groundwater_hydro_connectivity_type is not None
        meta.has_simple_infiltration = s.infiltration_rate_type is not None
        meta.has_groundwater = s.groundwater_impervious_layer_level_type is not None
        meta.has_interflow = s.interflow_type is not None and s.interflow_type != 0
        return cls(Nodes(id=[]), Lines(id=[]), meta=meta)

    @classmethod
    def from_quadtree(cls, quadtree, area_mask, node_id_counter, line_id_counter):
        """Construct the 2D grid based on the quadtree object.

        Args:
            quadtree (QuadTree): grid of refinement levels and active cell idx
              are used.
            area_mask (numpy.ndarray) : Array with raster of active pixels.
            node_id_counter (iterable): an iterable yielding integers
            line_id_counter (iterable): an iterable yielding integers

        Returns:
            Grid with data in the following columns:
            - nodes.id: integer generated from node_id_counter
            - nodes.coordinates: node coordinates
            - nodes.bounds: node bounds
            - nodes.node_type: Node_type (initially NODE_2D_OPEN_WATER)
            - nodes.nodk: Grid refinement level of node
            - nodes.nodm: horizontal index of node at grid refinement level nodk
            - nodes.nodn: vertical index of node at grid refinement level nodk
            - lines.id: integer generated from line_id_counter
            - lines.line: lines between connecting nodes.
            - lines.lik: Grid refinement level of line (smallest at refinements.)
            - lines.lim: horizontal index of line at grid refimenent level lik
            - lines.lin: vertical index line at grid refimenent level lik
        """

        nodes, lines = quadtree.get_nodes_lines(
            area_mask, node_id_counter, line_id_counter
        )

        # Some general quadtree grid statistics we need in the .h5 later on.
        quadtree_stats = QuadtreeStats(
            lgrmin=quadtree.lgrmin,
            kmax=quadtree.kmax,
            mmax=quadtree.mmax,
            nmax=quadtree.nmax,
            dx=quadtree.dx,
            dxp=quadtree.pixel_size,
            x0p=quadtree.origin[0],
            y0p=quadtree.origin[1],
        )

        return cls(nodes=nodes, lines=lines, quadtree_stats=quadtree_stats)

    @classmethod
    def from_connection_nodes(cls, connection_nodes, node_id_counter):
        """Construct a grid (only nodes) for the connection nodes

        Args:
            connection_nodes (ConnectionNodes): id and the_geom are used
            node_id_counter (iterable): an iterable yielding integers

        Returns:
            Grid with data in the following columns:
            - nodes.id: integer generated from node_id_counter
            - nodes.coordinates: node coordinates (from self.the_geom)
            - nodes.content_type: TYPE_V2_CONNECTION_NODES
            - nodes.content_pk: the user-supplied id
            - nodes.node_type: NODE_1D_NO_STORAGE / NODE_1D_STORAGE
            - nodes.calculation_type: only if set on Manhole
        """
        return cls(connection_nodes.get_nodes(node_id_counter), Lines(id=[]))

    @classmethod
    def from_channels(
        cls,
        connection_nodes,
        channels,
        global_dist_calc_points,
        node_id_counter,
        line_id_counter,
        connection_node_offset=0,
    ):
        """Construct a grid for the channels

        Args:
            connection_nodes (ConnectionNodes): used to map ids to indices
            channels (Channels)
            global_dist_calc_points (float): Default node interdistance.
            node_id_counter (iterable): an iterable yielding integers
            line_id_counter (iterable): an iterable yielding integers
            connection_node_offset (int): offset to give connection node
              indices in the returned lines.line. Default 0.

        Returns:
            Grid with data in the following columns:
            - nodes.id: integer generated from node_id_counter
            - nodes.coordinates
            - nodes.content_type: ContentType.TYPE_V2_CHANNEL
            - nodes.content_pk: the id of the Channel from which this node originates
            - nodes.node_type: NODE_1D_NO_STORAGE
            - nodes.calculation_type: from the channel
            - lines.id: integer generated from line_id_counter
            - lines.line: lines between connection nodes and added channel
              nodes. The indices are offset using the respective parameters.
            - lines.content_type: ContentType.TYPE_V2_CHANNEL
            - lines.content_pk: the id of the Channel from which this line originates
            - lines.kcu: from the channel's calculation_type
            - lines.line_geometries: a segment of the channel geometry
        """
        nodes = channels.interpolate_nodes(node_id_counter, global_dist_calc_points)
        lines = channels.get_lines(
            connection_nodes,
            nodes,
            line_id_counter,
            connection_node_offset=connection_node_offset,
        )
        return cls(nodes, lines)

    def set_channel_weights(self, cross_sections, channels):
        """Set cross section weights to channel nodes and lines.

        The attributes cross1, cross2, cross_weight are changed in place for nodes and
        lines whose content_type equals TYPE_V2_CHANNEL.

        Args:
            cross_sections (CrossSectionLocations)
            channels (Channels): Used to lookup the channel geometry
        """
        # Mask the nodes to only the Channel nodes
        node_mask = self.nodes.content_type == ContentType.TYPE_V2_CHANNEL
        cross1, cross2, cross_weight = cross_sections_module.compute_weights(
            self.nodes.content_pk[node_mask],
            self.nodes.ds1d[node_mask],
            cross_sections,
            channels,
        )
        self.nodes.cross1[node_mask] = cross1
        self.nodes.cross2[node_mask] = cross2
        self.nodes.cross_weight[node_mask] = cross_weight

        # Mask the lines to only the Channel lines
        line_mask = self.lines.content_type == ContentType.TYPE_V2_CHANNEL
        cross1, cross2, cross_weight = cross_sections_module.compute_weights(
            self.lines.content_pk[line_mask],
            self.lines.ds1d[line_mask],
            cross_sections,
            channels,
        )
        self.lines.cross1[line_mask] = cross1
        self.lines.cross2[line_mask] = cross2
        self.lines.cross_weight[line_mask] = cross_weight

    @classmethod
    def from_pipes(
        cls,
        connection_nodes,
        pipes,
        global_dist_calc_points,
        node_id_counter,
        line_id_counter,
        connection_node_offset=0,
    ):
        """Construct a grid for the pipes

        Args:
            connection_nodes (ConnectionNodes): used to map ids to indices
            pipes (Pipes)
            global_dist_calc_points (float): Default node interdistance.
            node_id_counter (iterable): an iterable yielding integers
            line_id_counter (iterable): an iterable yielding integers
            connection_node_offset (int): offset to give connection node
              indices in the returned lines.line. Default 0.

        Returns:
            Grid with data in the following columns:
            - nodes.id: integer generated from node_id_counter
            - nodes.coordinates
            - nodes.content_type: ContentType.TYPE_V2_PIPE
            - nodes.content_pk: the id of the Pipe from which this node originates
            - nodes.node_type: NODE_1D_NO_STORAGE
            - nodes.calculation_type: from the pipe
            - lines.id: integer generated from line_id_counter
            - lines.line: lines between connection nodes and added pipe
              nodes. The indices are offset using the respective parameters.
            - lines.content_type: ContentType.TYPE_V2_PIPE
            - lines.content_pk: the id of the Pipe from which this line originates
            - lines.kcu: from the pipe's calculation_type
            - lines.line_geometries: a segment of the pipe geometry
        """
        pipes.set_geometries(connection_nodes)
        nodes = pipes.interpolate_nodes(node_id_counter, global_dist_calc_points)
        lines = pipes.get_lines(
            connection_nodes,
            nodes,
            line_id_counter,
            connection_node_offset=connection_node_offset,
        )
        return cls(nodes, lines)

    @classmethod
    def from_structures(
        cls,
        connection_nodes,
        culverts,
        weirs,
        orifices,
        global_dist_calc_points,
        node_id_counter,
        line_id_counter,
        connection_node_offset=0,
    ):
        """Construct a grid for the culverts, weirs, and orifices

        Args:
            connection_nodes (ConnectionNodes): used to map ids to indices
            culverts (Culverts)
            weirs (Weirs)
            orifices (Orifices)
            global_dist_calc_points (float): Default node interdistance.
            node_id_counter (iterable): an iterable yielding integers
            line_id_counter (iterable): an iterable yielding integers
            connection_node_offset (int): offset to give connection node
              indices in the returned lines.line. Default 0.

        Returns:
            Grid with data in the following columns:
            - nodes.id: integer generated from node_id_counter
            - nodes.coordinates
            - nodes.content_type: culvert, weir, or orifice
            - nodes.content_pk: the id of the Culvert from which this node originates
            - nodes.node_type: NODE_1D_NO_STORAGE
            - nodes.calculation_type: from the culvert
            - lines.id: integer generated from line_id_counter
            - lines.line: lines between connection nodes and added culvert nodes
            - lines.content_type: culvert, weir, or orifice
            - lines.content_pk: the id of the object from which this line originates
            - lines.kcu: from the object's calculation_type
            - lines.line_geometries: a segment of the culvert's geometry or none
        """
        culverts.set_geometries(connection_nodes)

        nodes = culverts.interpolate_nodes(node_id_counter, global_dist_calc_points)
        lines = culverts.get_lines(
            connection_nodes,
            nodes,
            line_id_counter,
            connection_node_offset=connection_node_offset,
        )
        lines += weirs.get_lines(
            connection_nodes, line_id_counter, connection_node_offset
        )
        lines += orifices.get_lines(
            connection_nodes, line_id_counter, connection_node_offset
        )
        return cls(nodes, lines)

    def set_calculation_types(self):
        """Set the calculation types for connection nodes that do not yet have one.

        ConnectionNodes, Channels and Pipes should be present in the grid.
        """
        connection_nodes_module.set_calculation_types(self.nodes, self.lines)

    def set_bottom_levels(
        self, cross_sections, channels, pipes, weirs, orifices, culverts
    ):
        """Set the bottom levels (dmax and dpumax) for 1D nodes and lines

        This assumes that the channel weights have been computed already.

        The levels are based on:
        1. channel nodes: interpolate between crosssection locations
        2. pipe and culvert nodes: interpolate between invert level start & end
        3. connection nodes: see connection_nodes.set_bottom_levels
        4. lines: dpumax = greatest of the two neighboring nodes dmax
          - except for channels with no interpolated nodes: take reference level, but
            only if that is higher than the two neighboring nodes.
        """
        # Channels, interpolated nodes
        mask = self.nodes.content_type == ContentType.TYPE_V2_CHANNEL
        self.nodes.dmax[mask] = cross_sections_module.interpolate(
            self.nodes.cross1[mask],
            self.nodes.cross2[mask],
            self.nodes.cross_weight[mask],
            cross_sections,
            "reference_level",
        )

        # Pipes, interpolated nodes
        mask = self.nodes.content_type == ContentType.TYPE_V2_PIPE
        self.nodes.dmax[mask] = pipes.compute_bottom_level(
            self.nodes.content_pk[mask], self.nodes.ds1d[mask]
        )

        # Culverts, interpolated nodes
        mask = self.nodes.content_type == ContentType.TYPE_V2_CULVERT
        self.nodes.dmax[mask] = culverts.compute_bottom_level(
            self.nodes.content_pk[mask], self.nodes.ds1d[mask]
        )

        # Connection nodes: complex logic based on the connected objects
        connection_nodes_module.set_bottom_levels(
            self.nodes,
            self.lines,
            cross_sections,
            channels,
            pipes,
            weirs,
            orifices,
            culverts,
        )

        # Lines: based on the nodes
        self.lines.set_bottom_levels(self.nodes, allow_nan=True)

        # Fix channel lines: set dpumax of channel lines that have no interpolated nodes
        cross_sections_module.fix_dpumax(self.lines, self.nodes, cross_sections)

    def set_pumps(self, pumps):
        """Set the pumps on this grid object

        Args:
            pumps (Pumps): the pumps to set, this is changed inplace
        """
        self.pumps = pumps
        self.pumps.renumber()
        self.pumps.set_lines(self.nodes)
        self.pumps.set_node_data(self.nodes)

    def add_1d2d(
        self, connection_nodes, channels, pipes, locations, culverts, line_id_counter
    ):
        """Connect 1D and 2D elements by adding 1D-2D lines.

        Every (double) connected node gets a 1D-2D connection to the cell in which it
        is located. Double connected gives two lines to the same node.

        In addition to id and line attributes, also the kcu (line type) and dpumax
        (bottom level) are computed.
        """
        self.lines += get_1d2d_lines(
            self.nodes,
            connection_nodes,
            channels,
            pipes,
            locations,
            culverts,
            line_id_counter,
        )

    def finalize(self):
        """Finalize the Grid, computing and setting derived attributes"""
        self.lines.set_line_coords(self.nodes)
        self.lines.fix_line_geometries()
        self.lines.set_discharge_coefficients()
        if self.pumps is not None and len(self.pumps) > 0:
            self.meta.has_pumpstations = True
        self.meta.extent_1d = self.nodes.get_extent_1d()
        self.meta.extent_2d = self.nodes.get_extent_2d()
        self.meta.has_1d = self.meta.extent_1d is not None
        self.meta.has_2d = self.meta.extent_2d is not None


def get_1d2d_lines(
    nodes, connection_nodes, channels, pipes, locations, culverts, line_id_counter
):
    """Compute 1D-2D flowlines for (double) connected 1D nodes.

    Currently implemented 1D nodes are: connection nodes and channel nodes.

    The line type and bottom level (kcu and dpumax) are set according to the
    "get_1d2d_properties" on the corresponding objects (ConnectionNodes, Channels, ...).

    If the 2D bottom level will turn out higher, this will be corrected later by
    threedi-tables.

    The line will connect to the cell in which the 1D node is located. For edge cases,
    the line will connect to the cell with the lowest id. Users may want to influence
    which cells are connected to. This is currently unimplemented.

    Args:
        nodes (Nodes): All 1D and 2D nodes to compute 1D-2D lines for. Be sure that
            the 1D nodes already have a dmax and calculation_type set.
        connection_nodes (ConnectionNodes): for looking up manhole_id and drain_level
        channels (Channels)
        pipes (Pipes)
        locations (CrossSectionLocations): for looking up bank_level
        culverts (Culverts)
        line_id_counter (iterable): An iterable yielding integers

    Returns:
        Lines with data in the following columns:
        - id: integer generated from line_id_counter
        - kcu: LINE_1D2D_* type (see above)
        - dpumax: based on drain_level or dmax (see above)
    """
    is_2d_open = nodes.node_type == NodeType.NODE_2D_OPEN_WATER
    cell_ids = nodes.id[is_2d_open]
    cell_tree = pygeos.STRtree(pygeos.box(*nodes.bounds[is_2d_open].T))

    connected_idx = np.where(
        np.isin(
            nodes.content_type,
            [ContentType.TYPE_V2_CONNECTION_NODES, ContentType.TYPE_V2_CHANNEL],
        )
        & np.isin(
            nodes.calculation_type,
            [CalculationType.CONNECTED, CalculationType.DOUBLE_CONNECTED],
        )
    )[0]

    # The query_bulk returns 2 1D arrays: one with indices into the supplied node
    # geometries and one with indices into the tree of cells.
    idx = cell_tree.query_bulk(pygeos.points(nodes.coordinates[connected_idx]))
    # Address edge cases: just take the first line
    _, first_unique_index = np.unique(idx[0], return_index=True)
    idx = idx[:, first_unique_index]
    n_lines = idx.shape[1]
    if n_lines == 0:
        return Lines(id=[])
    node_idx = connected_idx[idx[0]]  # convert to node indexes
    node_id = nodes.index_to_id(node_idx)  # convert to node ids
    cell_id = cell_ids[idx[1]]  # convert to cell ids

    # create a 'duplicator' array that duplicates lines from double connected nodes
    is_double = nodes.calculation_type[node_idx] == CalculationType.DOUBLE_CONNECTED
    duplicator = np.ones(n_lines, dtype=int)
    duplicator[is_double] = 2
    duplicator = np.repeat(np.arange(n_lines), duplicator)

    # Identify different types of objects and dispatch to the associated functions
    is_ch = nodes.content_type[node_idx] == ContentType.TYPE_V2_CHANNEL
    is_cn = nodes.content_type[node_idx] == ContentType.TYPE_V2_CONNECTION_NODES

    is_sewerage = np.zeros(n_lines, dtype=bool)
    dpumax = np.full(n_lines, fill_value=np.nan, dtype=np.float64)

    is_sewerage[is_cn], dpumax[is_cn] = connection_nodes.get_1d2d_properties(
        nodes, node_idx[is_cn]
    )
    is_sewerage[is_ch], dpumax[is_ch] = channels.get_1d2d_properties(
        nodes, node_idx[is_ch], locations
    )

    # map "is_sewerage" to "kcu" (including double/single connected properties)
    # map the two binary arrays on numbers 0, 1, 2, 3
    options = is_double * 2 + is_sewerage
    kcu = np.choose(
        options,
        choices=[
            LineType.LINE_1D2D_SINGLE_CONNECTED_OPEN_WATER,
            LineType.LINE_1D2D_SINGLE_CONNECTED_SEWERAGE,
            LineType.LINE_1D2D_DOUBLE_CONNECTED_OPEN_WATER,
            LineType.LINE_1D2D_DOUBLE_CONNECTED_SEWERAGE,
        ],
    )

    return Lines(
        id=itertools.islice(line_id_counter, len(duplicator)),
        line=np.array([node_id, cell_id]).T[duplicator],
        kcu=kcu[duplicator],
        dpumax=dpumax[duplicator],
    )
