from . import connection_nodes as connection_nodes_module
from . import cross_section_locations as csl_module
from . import embedded as embedded_module
from . import obstacles as obstacles_module
from .cross_section_definitions import CrossSections
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
from threedigrid_builder.exceptions import SchematisationError
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
        pumps: Optional[Pumps] = None,
        cross_sections: Optional[CrossSections] = None,
        nodes_embedded=None,
        meta=None,
        quadtree_stats=None,
    ):
        if not isinstance(nodes, Nodes):
            raise TypeError(f"Expected Nodes instance, got {type(nodes)}")
        if not isinstance(lines, Lines):
            raise TypeError(f"Expected Lines instance, got {type(lines)}")
        if pumps is None:
            pumps = Pumps(id=[])
        elif not isinstance(pumps, Pumps):
            raise TypeError(f"Expected Pumps instance, got {type(pumps)}")
        if cross_sections is None:
            cross_sections = CrossSections(id=[])
        elif not isinstance(cross_sections, CrossSections):
            raise TypeError(
                f"Expected CrossSections instance, got {type(cross_sections)}"
            )
        if nodes_embedded is None:
            nodes_embedded = Nodes(id=[])
        elif not isinstance(nodes_embedded, Nodes):
            raise TypeError(f"Expected Nodes instance, got {type(nodes_embedded)}")
        self.nodes = nodes
        self.lines = lines
        self.meta = meta
        self.quadtree_stats = quadtree_stats
        self.pumps = pumps
        self.cross_sections = cross_sections
        self.nodes_embedded = nodes_embedded
        self._cell_tree = None

    def __add__(self, other):
        """Concatenate two grids without renumbering nodes."""
        if self.__class__ is not other.__class__:
            raise TypeError(
                "Cannot concatenate {self} with {other} as they are not of "
                "equal types."
            )
        new_attrs = {}
        for name in ("nodes", "lines", "nodes_embedded"):
            new_attrs[name] = getattr(self, name) + getattr(other, name)
        for name in ("meta", "quadtree_stats", "pumps", "cross_sections"):
            if getattr(other, name) is None:
                new_attrs[name] = getattr(self, name)
            else:
                new_attrs[name] = getattr(other, name)
        return self.__class__(**new_attrs)

    def __repr__(self):
        return f"<Grid object with {len(self.nodes)} nodes and {len(self.lines)} lines>"

    @property
    def cell_tree(self):
        """A pygeos STRtree of the cells in this grid

        The indices in the tree equal the node indices.
        """
        if self._cell_tree is None:
            is_2d_open = self.nodes.node_type == NodeType.NODE_2D_OPEN_WATER
            geoms = np.empty(len(self.nodes), dtype=object)
            geoms[is_2d_open] = pygeos.box(*self.nodes.bounds[is_2d_open].T)
            self._cell_tree = pygeos.STRtree(geoms)
        return self._cell_tree

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

        lines.set_line_coords(nodes)
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
    def from_linear_objects(
        cls,
        connection_nodes,
        objects,
        definitions,
        cell_tree,
        global_dist_calc_points,
        embedded_cutoff_threshold,
        node_id_counter,
        embedded_node_id_counter,
        line_id_counter,
        connection_node_offset=0,
    ):
        """Construct a grid for linear objects (channels, pipes, culverts)

        Args:
            connection_nodes (ConnectionNodes): used to map ids to indices
            objects (Channels, Pipes, Culverts)
            definitions (CrossSectionDefinitions): to map definition ids
            cell_tree (pygeos.STRtree): strtree of the 2D cells (for embedded channels)
            global_dist_calc_points (float): Default node interdistance.
            embedded_cutoff_threshold (float): The min length of an embedded line.
            node_id_counter (iterable): an iterable yielding integers
            embedded_node_id_counter (iterable): an iterable yielding integers
            line_id_counter (iterable): an iterable yielding integers
            connection_node_offset (int): offset to give connection node
              indices in the returned lines.line. Default 0.

        Returns:
            Grid with data in the following columns:
            - nodes.id: integer generated from node_id_counter
            - nodes.coordinates
            - nodes.content_type: depends on linear object type
            - nodes.content_pk: the id of the object from which this node originates
            - nodes.node_type: NODE_1D_NO_STORAGE
            - nodes.calculation_type: copied from the object
            - nodes.s1d: the node's position along the channel
            - lines.id: integer generated from line_id_counter
            - lines.line: lines between nodes (including connection nodes)
            - lines.content_type: depends on linear object type
            - lines.content_pk: the id of the object from which this line originates
            - lines.kcu: from the object's calculation_type
            - lines.line_geometries: a segment of the channel geometry
            - lines.ds1d: the length of the line
            - lines.s1d: the position of the line's velocity point along the object
            - lines.cross1: the cross section definition id of the line (pipe/culvert)
            - lines.cross_weight: always 1.0 (pipe/culvert)
            - lines.invert_level_start_point: invert level at line start (pipe/culvert)
            - lines.invert_level_end_point: invert level at line end (pipe/culvert)
        """
        objects.set_geometries(connection_nodes)
        nodes = objects.interpolate_nodes(node_id_counter, global_dist_calc_points)
        lines = objects.get_lines(
            connection_nodes,
            definitions,
            nodes,
            line_id_counter,
            connection_node_offset=connection_node_offset,
        )
        nodes_embedded, lines_s1d = objects.get_embedded(
            cell_tree,
            embedded_cutoff_threshold,
            embedded_node_id_counter,
        )
        if len(lines_s1d) > 0:
            # And construct the lines (also connecting to connection nodes)
            lines_embedded = objects.get_lines(
                connection_nodes,
                definitions,
                nodes_embedded,
                line_id_counter,
                connection_node_offset=connection_node_offset,
                embedded_mode=True,
            )
            # Override the velocity point locations (the defaulted to the line midpoint,
            # while for embedded objects we force them to the cell edges)
            lines_embedded.s1d[:] = lines_s1d
            lines += lines_embedded
            # Later gridbuilder functions expect ordering by content_pk
            lines.reorder_by("content_pk")
        return cls(nodes, lines, nodes_embedded=nodes_embedded)

    @classmethod
    def from_structures(
        cls,
        connection_nodes,
        weirs,
        orifices,
        definitions,
        line_id_counter,
        connection_node_offset=0,
    ):
        """Construct a grid for the weirs and orifices

        Args:
            connection_nodes (ConnectionNodes): used to map ids to indices
            weirs (Weirs)
            orifices (Orifices)
            definitions (CrossSectionDefinitions)
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
        lines = weirs.get_lines(
            connection_nodes, definitions, line_id_counter, connection_node_offset
        )
        lines += orifices.get_lines(
            connection_nodes, definitions, line_id_counter, connection_node_offset
        )
        return cls(Nodes(id=[]), lines)

    def set_calculation_types(self):
        """Set the calculation types for connection nodes that do not yet have one.

        ConnectionNodes, Channels and Pipes should be present in the grid.
        """
        connection_nodes_module.set_calculation_types(self.nodes, self.lines)

    def set_bottom_levels(self):
        """Set the bottom levels (dmax and dpumax) for 1D nodes and lines

        This assumes that the channel weights have been computed already.

        The levels are based on:
        1. channel, pipe, culvert nodes: use invert level start & end of the lines
        2. connection nodes: see connection_nodes.set_bottom_levels
        3. lines: dpumax = greatest of the two neighboring nodes dmax
          - except for channels with no interpolated nodes: take reference level, but
            only if that is higher than the two neighboring nodes.
        """
        CH = ContentType.TYPE_V2_CHANNEL
        PI = ContentType.TYPE_V2_PIPE
        CV = ContentType.TYPE_V2_CULVERT

        # Channels, Pipes, Culverts, interpolated nodes which require dmax update:
        mask = np.isin(self.nodes.content_type, [CH, PI, CV]) & ~np.isfinite(
            self.nodes.dmax
        )

        # Get the invert level from the connecting line. Each interpolated node has
        # by definition exactly one line going from it. So we can do this:
        sorter = np.argsort(self.lines.line[:, 0])
        idx = np.searchsorted(self.lines.line[sorter, 0], self.nodes.id[mask])
        assert np.all(self.lines.line[sorter, 0][idx] == self.nodes.id[mask])
        self.nodes.dmax[mask] = self.lines.invert_level_start_point[sorter][idx]

        # Connection nodes: logic based on the connected objects
        connection_nodes_module.set_bottom_levels(self.nodes, self.lines)

        # Lines: based on the nodes
        self.lines.set_bottom_levels(self.nodes, allow_nan=True)

        # Fix channel lines: set dpumax of channel lines that have no interpolated nodes
        csl_module.fix_dpumax(self.lines, self.nodes)

    def set_obstacles(self, obstacles):
        """Set obstacles on 2D lines by determining intersection between
           line_coords (these must be knows at this point) and obstacle geometry.
           Set kcu to 101 and changes flod and flou to crest_level.

        Args:
            obstacles (Obstacles)
        """
        if len(obstacles) > 0:
            obstacles_module.apply_obstacles(self.lines, obstacles)

    def set_pumps(self, pumps):
        """Set the pumps on this grid object

        Args:
            pumps (Pumps): the pumps to set, this is changed inplace
        """
        self.pumps = pumps
        self.pumps.renumber()
        self.pumps.set_lines(self.nodes)
        self.pumps.set_node_data(self.nodes)

    def set_cross_sections(self, definitions):
        """Set the cross sections on this grid object

        Args:
            definitions (CrossSectionDefinitions)
        """
        # TODO Skip definitions that are not used (and remap cross1 and cross2)
        self.cross_sections = definitions.convert()

    def embed_nodes(self, embedded_node_id_counter):
        """Integrate embedded connection nodes into the 2D cells.

        grid.nodes_embedded will contain the removed (embedded) nodes.
        """
        # before connecting lines to 2D nodes, first set the line_coords based on the
        # original coordinate of the embedded nodes:
        self.lines.set_line_coords(self.nodes)

        self.nodes_embedded += embedded_module.embed_nodes(
            self, embedded_node_id_counter
        )

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
            self.cell_tree,
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
        if len(self.pumps) > 0:
            self.meta.has_pumpstations = True
        self.meta.extent_1d = self.nodes.get_extent_1d()
        self.meta.extent_2d = self.nodes.get_extent_2d()
        self.meta.has_1d = self.meta.extent_1d is not None
        self.meta.has_2d = self.meta.extent_2d is not None


def get_1d2d_lines(
    cell_tree,
    nodes,
    connection_nodes,
    channels,
    pipes,
    locations,
    culverts,
    line_id_counter,
):
    """Compute 1D-2D flowlines for (double) connected 1D nodes.

    The line type and bottom level (kcu and dpumax) are set according to the
    "get_1d2d_properties" on the corresponding objects (ConnectionNodes, Channels,
    Pipes, Culverts).

    If the 2D bottom level will turn out higher, this will be corrected later by
    threedi-tables.

    The line will connect to the cell in which the 1D node is located. For edge cases,
    the line will connect to the cell with the lowest id. Users may want to influence
    which cells are connected to. This is currently unimplemented.

    Args:
        cell_tree (pygeos.STRtree): An STRtree containing the cells. The indices into
          the STRtree must be equal to the indices into the nodes.
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
    HAS_1D2D_CONTENT_TYPES = (
        ContentType.TYPE_V2_CONNECTION_NODES,
        ContentType.TYPE_V2_CHANNEL,
        ContentType.TYPE_V2_PIPE,
        ContentType.TYPE_V2_CULVERT,
    )
    HAS_1D2D_CALC_TYPES = (CalculationType.CONNECTED, CalculationType.DOUBLE_CONNECTED)
    connected_idx = np.where(
        np.isin(nodes.content_type, HAS_1D2D_CONTENT_TYPES)
        & np.isin(nodes.calculation_type, HAS_1D2D_CALC_TYPES)
    )[0]
    n_connected_nodes = connected_idx.shape[0]
    if n_connected_nodes == 0:
        return Lines(id=[])

    # The query_bulk returns 2 1D arrays: one with indices into the supplied node
    # geometries and one with indices into the tree of cells.
    idx = cell_tree.query_bulk(pygeos.points(nodes.coordinates[connected_idx]))
    # Address edge cases of multiple 1D-2D lines per node: just take the one
    _, first_unique_index = np.unique(idx[0], return_index=True)
    idx = idx[:, first_unique_index]
    n_lines = idx.shape[1]
    # Error if there is a node without a 1D-2D line
    if n_lines != n_connected_nodes:
        # The code in this if clause is only for pretty error formatting.
        out_of_bounds = np.delete(connected_idx, idx[0])
        types = nodes.content_type[out_of_bounds]
        pks = nodes.content_pk[out_of_bounds]
        pretty_names = ("connection nodes", "channels", "pipes", "culverts")
        object_pk_list = [
            f"{pretty_name} {sorted(set(pks[types == content_type]))}"
            for content_type, pretty_name in zip(HAS_1D2D_CONTENT_TYPES, pretty_names)
            if np.any(types == content_type)
        ]
        raise SchematisationError(
            f"The following object(s) have a connected calculation type but are "
            f"(partially) outside of the 2D calculation cells: "
            f"{', '.join(object_pk_list)}."
        )
    if n_lines == 0:
        return Lines(id=[])
    node_idx = connected_idx[idx[0]]  # convert to node indexes
    node_id = nodes.index_to_id(node_idx)  # convert to node ids
    cell_id = nodes.index_to_id(idx[1])  # convert to cell ids

    # create a 'duplicator' array that duplicates lines from double connected nodes
    is_double = nodes.calculation_type[node_idx] == CalculationType.DOUBLE_CONNECTED
    duplicator = np.ones(n_lines, dtype=int)
    duplicator[is_double] = 2
    duplicator = np.repeat(np.arange(n_lines), duplicator)

    # Identify different types of objects and dispatch to the associated functions
    is_ch = nodes.content_type[node_idx] == ContentType.TYPE_V2_CHANNEL
    is_cn = nodes.content_type[node_idx] == ContentType.TYPE_V2_CONNECTION_NODES
    is_pipe = nodes.content_type[node_idx] == ContentType.TYPE_V2_PIPE
    is_culvert = nodes.content_type[node_idx] == ContentType.TYPE_V2_CULVERT

    is_closed = np.zeros(n_lines, dtype=bool)
    dpumax = np.full(n_lines, fill_value=np.nan, dtype=np.float64)

    is_closed[is_cn], dpumax[is_cn] = connection_nodes.get_1d2d_properties(
        nodes, node_idx[is_cn]
    )
    is_closed[is_ch], dpumax[is_ch] = channels.get_1d2d_properties(
        nodes, node_idx[is_ch], locations
    )
    is_closed[is_pipe], dpumax[is_pipe] = pipes.get_1d2d_properties(
        nodes,
        node_idx[is_pipe],
        connection_nodes,
    )
    is_closed[is_culvert], dpumax[is_culvert] = culverts.get_1d2d_properties(
        nodes,
        node_idx[is_culvert],
        connection_nodes,
    )

    # map "is_closed" to "kcu" (including double/single connected properties)
    # map the two binary arrays on numbers 0, 1, 2, 3
    options = is_double * 2 + is_closed
    kcu = np.choose(
        options,
        choices=[
            LineType.LINE_1D2D_SINGLE_CONNECTED_OPEN_WATER,
            LineType.LINE_1D2D_SINGLE_CONNECTED_CLOSED,
            LineType.LINE_1D2D_DOUBLE_CONNECTED_OPEN_WATER,
            LineType.LINE_1D2D_DOUBLE_CONNECTED_CLOSED,
        ],
    )

    return Lines(
        id=itertools.islice(line_id_counter, len(duplicator)),
        line=np.array([node_id, cell_id]).T[duplicator],
        kcu=kcu[duplicator],
        dpumax=dpumax[duplicator],
    )
