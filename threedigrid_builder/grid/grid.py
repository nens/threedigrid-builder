from . import connection_nodes as connection_nodes_module
from . import cross_section_locations as csl_module
from . import embedded as embedded_module
from . import initial_waterlevels as initial_waterlevels_module
from . import obstacles as obstacles_module
from .cross_section_definitions import CrossSections
from dataclasses import dataclass
from dataclasses import fields
from threedigrid_builder.base import Lines
from threedigrid_builder.base import Nodes
from threedigrid_builder.base import Pumps
from threedigrid_builder.base import replace
from threedigrid_builder.base.settings import GridSettings
from threedigrid_builder.base.settings import TablesSettings
from threedigrid_builder.constants import ContentType
from threedigrid_builder.constants import LineType
from threedigrid_builder.constants import NodeType
from typing import Optional
from typing import Tuple

import numpy as np
import pygeos
import threedigrid_builder


__all__ = ["Grid", "GridMeta", "QuadtreeStats"]

NODE_ORDER = {
    NodeType.NODE_2D_OPEN_WATER: 1,
    NodeType.NODE_2D_GROUNDWATER: 2,
    NodeType.NODE_1D_NO_STORAGE: 3,
    NodeType.NODE_1D_STORAGE: 3,
    NodeType.NODE_2D_BOUNDARIES: 4,
    NodeType.NODE_2D_GROUNDWATER_BOUNDARIES: 5,
    NodeType.NODE_1D_BOUNDARIES: 6,
}
# failsafe for future adding of types:
assert len(NodeType) == len(NODE_ORDER)

LINE_ORDER = {
    LineType.LINE_2D_U: 1,
    LineType.LINE_2D_OBSTACLE_U: 1,
    LineType.LINE_2D_V: 2,
    LineType.LINE_2D_OBSTACLE_V: 2,
    LineType.LINE_2D: 2,
    LineType.LINE_2D_OBSTACLE: 2,
    LineType.LINE_2D_VERTICAL: 3,
    LineType.LINE_2D_GROUNDWATER: 4,
    LineType.LINE_1D_EMBEDDED: 5,
    LineType.LINE_1D_ISOLATED: 5,
    LineType.LINE_1D_CONNECTED: 5,
    LineType.LINE_1D_LONG_CRESTED: 5,
    LineType.LINE_1D_SHORT_CRESTED: 5,
    LineType.LINE_1D_DOUBLE_CONNECTED: 5,
    LineType.LINE_1D2D_SINGLE_CONNECTED_CLOSED: 6,
    LineType.LINE_1D2D_SINGLE_CONNECTED_OPEN_WATER: 6,
    LineType.LINE_1D2D_DOUBLE_CONNECTED_CLOSED: 6,
    LineType.LINE_1D2D_DOUBLE_CONNECTED_OPEN_WATER: 6,
    LineType.LINE_1D2D_POSSIBLE_BREACH: 6,
    LineType.LINE_1D2D_ACTIVE_BREACH: 6,
    LineType.LINE_1D2D_GROUNDWATER_OPEN_WATER: 7,
    LineType.LINE_1D2D_GROUNDWATER_SEWER: 7,
    LineType.LINE_2D_BOUNDARY_WEST: 8,
    LineType.LINE_2D_BOUNDARY_EAST: 8,
    LineType.LINE_2D_BOUNDARY_SOUTH: 8,
    LineType.LINE_2D_BOUNDARY_NORTH: 8,
    # LineType.LINE_2D_GROUNDWATER_BOUNDARY: 9, (to be implemented)
    LineType.LINE_1D_BOUNDARY: 10,
}
# failsafe for future adding of types:
assert len(LineType) == len(LINE_ORDER)


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
    has_initial_waterlevels: bool = True

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

    def set_initial_waterlevels( self, connection_nodes, channels, pipes, culverts):
        """Apply initial waterlevels (global or per connection nodes) to all 1D nodes.

        Bottom levels (dmax) should be set already.

        Args:
            connection_nodes (ConnectionNodes): used to map ids to indices
            channels (Channels)
            pipes (Pipes)
            culverts (Culverts)

        """
        initial_waterlevels_module.compute_initial_waterlevels(
            self.nodes,
            connection_nodes=connection_nodes,
            channels=channels,
            pipes=pipes,
            culverts=culverts,
        )

    def set_obstacles(self, obstacles):
        """Set obstacles on 2D lines by determining intersection between
           line_coords (these must be knows at this point) and obstacle geometry.
           Set kcu to 101 and changes flod and flou to crest_level.

        Args:
            obstacles (Obstacles)
        """
        if len(obstacles) > 0:
            obstacles_module.apply_obstacles(self.lines, obstacles)

    def set_boundary_conditions_1d(self, boundary_conditions_1d):
        boundary_conditions_1d.apply(self)

    def set_boundary_conditions_2d(
        self,
        boundary_conditions_2d,
        quadtree,
        node_id_counter,
        line_id_counter,
    ):
        nodes, lines = boundary_conditions_2d.get_nodes_and_lines(
            self.nodes,
            self.cell_tree,
            quadtree,
            node_id_counter,
            line_id_counter,
        )
        self.nodes += nodes
        self.lines += lines

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
        self,
        connected_points,
        connection_nodes,
        channels,
        pipes,
        locations,
        culverts,
        line_id_counter,
    ):
        """Connect 1D and 2D elements by adding 1D-2D lines.

        Every (double) connected node gets a 1D-2D connection to the cell in which it
        is located. Double connected gives two lines to the same node.

        In addition to id and line attributes, also the kcu (line type) and dpumax
        (bottom level) are computed.
        """
        self.lines += connected_points.get_lines(
            self.cell_tree,
            self.nodes,
            connection_nodes,
            channels,
            pipes,
            locations,
            culverts,
            line_id_counter,
        )

    def sort(self):
        """Sort the nodes and lines into the order required by the calculation core.

        Resets the nodes and lines ids to consecutive indices starting from 0.
        Also maps lines.line and pumps.line according to the new nodes indices.

        See NODE_ORDER and LINE_ORDER for the order.
        """
        node_sorter = np.argsort(replace(self.nodes.node_type, NODE_ORDER))
        line_sorter = np.argsort(replace(self.lines.kcu, LINE_ORDER))

        # now sort the nodes and lines and reset their ids
        old_node_ids = self.nodes.id.copy()
        self.nodes.reorder(node_sorter)
        self.nodes.id[:] = np.arange(len(self.nodes))
        self.lines.reorder(line_sorter)
        self.lines.id[:] = np.arange(len(self.lines))

        # create a mapping with new node ids on the position of the old node ids
        new_ids = np.empty(old_node_ids[-1] + 1, dtype=self.nodes.id.dtype)
        new_ids[old_node_ids[node_sorter]] = self.nodes.id

        # apply the mapping to lines.line and optionally to pumps and embedded nodes
        self.lines.line[:] = np.take(new_ids, self.lines.line)
        if self.pumps is not None:
            mask = self.pumps.line != -9999
            self.pumps.line[mask] = np.take(new_ids, self.pumps.line[mask])
        if self.nodes_embedded is not None:
            self.nodes_embedded.embedded_in = np.take(
                new_ids, self.nodes_embedded.embedded_in
            )

        # sort boundary lines so that they internally match the boundary node order
        BOUNDARY_NODE_TYPES = (
            NodeType.NODE_1D_BOUNDARIES,
            NodeType.NODE_2D_BOUNDARIES,
            NodeType.NODE_2D_GROUNDWATER_BOUNDARIES,
        )
        bc_node_ids = self.nodes.id[np.isin(self.nodes.node_type, BOUNDARY_NODE_TYPES)]
        self.lines.sort_by_nodes(bc_node_ids)

    def finalize(self):
        """Finalize the Grid, computing and setting derived attributes"""
        self.sort()
        self.lines.set_line_coords(self.nodes)
        self.lines.fix_line_geometries()
        self.lines.set_discharge_coefficients()
        if len(self.pumps) > 0:
            self.meta.has_pumpstations = True
        self.meta.has_initial_waterlevels = np.isfinite(self.nodes.initial_waterlevel).any()
        self.meta.extent_1d = self.nodes.get_extent_1d()
        self.meta.extent_2d = self.nodes.get_extent_2d()
        self.meta.has_1d = self.meta.extent_1d is not None
        self.meta.has_2d = self.meta.extent_2d is not None
