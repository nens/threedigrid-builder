from . import connection_nodes as connection_nodes_module
from . import dem_average_area as dem_average_area_module
from . import embedded as embedded_module
from . import initial_waterlevels as initial_waterlevels_module
from . import obstacles as obstacles_module
from .cross_section_definitions import CrossSections
from dataclasses import dataclass
from dataclasses import fields
from threedigrid_builder.base import Breaches
from threedigrid_builder.base import Levees
from threedigrid_builder.base import Lines
from threedigrid_builder.base import Nodes
from threedigrid_builder.base import Pumps
from threedigrid_builder.base import Surfaces
from threedigrid_builder.base import SurfaceMaps
from threedigrid_builder.base import replace
from threedigrid_builder.base.settings import GridSettings
from threedigrid_builder.base.settings import TablesSettings
from threedigrid_builder.constants import ContentType
from threedigrid_builder.constants import LineType
from threedigrid_builder.constants import NodeType
from threedigrid_builder.grid import zero_d

from typing import Optional
from typing import Tuple, Union

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
}
LINE_ORDER_BOUNDARY_1D = 10
# failsafe for future adding of types:
assert len(LineType) == len(LINE_ORDER)


@dataclass
class GridMeta:
    """Metadata that needs to end up in the gridadmin file."""

    model_name: str  # name from sqlite globalsettings.name

    grid_settings: GridSettings
    tables_settings: TablesSettings

    epsg_code: int = 28992  # SQLite EPSG is ignored if DEM file is present
    model_slug: Optional[str] = None  # from repository.slug
    revision_hash: Optional[str] = None  # from repository.revision.hash
    revision_nr: Optional[int] = None  # from repository.revision.number
    threedi_version: Optional[str] = None  # threedi-api version
    threedicore_version: Optional[str] = None  # not used anymore
    threedigrid_builder_version: str = ""  # filled in __post_init__
    threedi_tables_version: Optional[str] = None  # filled in threedi-tables

    # TODO what to do with use_1d_flow, use_2d_flow, manhole_storage_area
    has_0d: bool = False
    has_1d: bool = False
    has_2d: bool = False
    has_embedded: bool = False
    has_breaches: bool = False
    has_groundwater: bool = False
    has_groundwater_flow: bool = False
    has_interception: bool = False
    has_pumpstations: bool = False
    has_simple_infiltration: bool = False
    has_max_infiltration_capacity: bool = False
    has_interflow: bool = False
    has_initial_waterlevels: bool = True

    extent_1d: Optional[Tuple[float, float, float, float]] = None
    extent_2d: Optional[Tuple[float, float, float, float]] = None
    zero_dim_extent: Optional[Tuple[float, float, float, float]] = None

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
        surfaces: Optional[Surfaces] = None,
        surface_maps : Optional[SurfaceMaps] = None,
        nodes_embedded=None,
        levees=None,
        breaches=None,
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

        if surfaces is None:
            surfaces = Surfaces(id=[])
        elif not isinstance(surfaces, Surfaces):
            raise TypeError(f"Expected Surfaces instance, got {type(surfaces)}")

        if surface_maps is None:
            surface_maps = SurfaceMaps(id=[])
        elif not isinstance(surface_maps, SurfaceMaps):
            raise TypeError(f"Expected SurfaceMaps instance, got {type(surface_maps)}")

        if levees is not None and not isinstance(levees, Levees):
            raise TypeError(f"Expected Levees instance, got {type(levees)}")
        if breaches is not None and not isinstance(breaches, Breaches):
            raise TypeError(f"Expected Breaches instance, got {type(breaches)}")
        self.nodes = nodes
        self.lines = lines
        self.meta = meta
        self.surfaces = surfaces
        self.surface_maps = surface_maps
        self.quadtree_stats = quadtree_stats
        self.pumps = pumps
        self.cross_sections = cross_sections
        self.nodes_embedded = nodes_embedded
        self.levees = levees
        self.breaches = breaches
        self._cell_tree = None

    def __add__(self, other):
        """Concatenate two grids without renumbering nodes."""
        if self.__class__ is not other.__class__:
            raise TypeError(
                "Cannot concatenate {self} with {other} as they are not of "
                "equal types."
            )
        new_attrs = {}
        for name in ("nodes", "lines", "nodes_embedded", "surfaces"):
            new_attrs[name] = getattr(self, name) + getattr(other, name)
        for name in (
            "meta",
            "quadtree_stats",
            "pumps",
            "cross_sections",
            "levees",
            "breaches",
        ):
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
        meta.has_max_infiltration_capacity = (
            s.max_infiltration_capacity_type is not None
        )
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
             These attributes are filled for pipes/culverts only:
            - lines.cross_id1 & cross_id2: the id of the cross section definition
            - lines.cross_weight: 1.0 (which means that cross_id2 should be ignored)
            - lines.frict_type1 & frict_type2: the friction type (both are equal)
            - lines.frict_value1 & frict_value2: the friction value (both are equal)
            - lines.invert_level_start_point: copied from pipe/culvert
            - lines.invert_level_end_point: copied from pipe/culvert
            - lines.dpumax: largest of the two invert levels
            - lines.discharge_coefficient_positive & _positive: culverts only
        """
        objects.set_geometries(connection_nodes)
        nodes = objects.interpolate_nodes(node_id_counter, global_dist_calc_points)
        lines = objects.get_lines(
            connection_nodes,
            nodes,
            line_id_counter,
            connection_node_offset=connection_node_offset,
        )
        nodes_embedded, lines_s1d, lines_ds1d_half = objects.get_embedded(
            cell_tree,
            embedded_cutoff_threshold,
            embedded_node_id_counter,
        )
        if len(lines_s1d) > 0:
            # And construct the lines (also connecting to connection nodes)
            lines_embedded = objects.get_lines(
                connection_nodes,
                nodes_embedded,
                line_id_counter,
                connection_node_offset=connection_node_offset,
                embedded_mode=True,
            )
            # Override the velocity point locations (the defaulted to the line midpoint,
            # while for embedded objects we force them to the cell edges)
            lines_embedded.s1d[:] = lines_s1d
            lines_embedded.ds1d_half[:] = lines_ds1d_half
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
        line_id_counter,
        connection_node_offset=0,
    ):
        """Construct a grid for the weirs and orifices

        Args:
            connection_nodes (ConnectionNodes): used to map ids to indices
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
            - lines.dpumax: the crest_level of the structure
            - lines.cross_id1 & cross_id2: the id of the cross section definition
            - lines.cross_weight: 1.0 (which means that cross_id2 should be ignored)
            - lines.frict_type1 & frict_type2: the friction type (both are equal)
            - lines.frict_value1 & frict_value2: the friction value (both are equal)
            - lines.invert_level_start_point: the crest_level of the structure
            - lines.invert_level_end_point: the crest_level of the structure
            - lines.discharge_coefficient_positive and _negative: taken from the structure
        """
        lines = weirs.get_lines(
            connection_nodes, line_id_counter, connection_node_offset
        )
        lines += orifices.get_lines(
            connection_nodes, line_id_counter, connection_node_offset
        )
        return cls(Nodes(id=[]), lines)

    def set_calculation_types(self):
        """Set the calculation types for connection nodes that do not yet have one.

        ConnectionNodes, Channels and Pipes should be present in the grid.
        """
        connection_nodes_module.set_calculation_types(self.nodes, self.lines)

    def set_bottom_levels(self):
        """Set the bottom levels (dmax) for 1D nodes

        This assumes that the channel weights have been computed already.

        The levels are based on:
        1. channel, pipe, culvert nodes: use line invert levels
        2. connection nodes: see connection_nodes.set_bottom_levels
        """
        # All 1D lines should already have a dpumax
        LINE_TYPE_1D = [
            LineType.LINE_1D_EMBEDDED,
            LineType.LINE_1D_ISOLATED,
            LineType.LINE_1D_CONNECTED,
            LineType.LINE_1D_LONG_CRESTED,
            LineType.LINE_1D_SHORT_CRESTED,
            LineType.LINE_1D_DOUBLE_CONNECTED,
        ]
        is_1d_line = np.isin(self.lines.kcu, LINE_TYPE_1D)
        if is_1d_line.any() and not np.isfinite(self.lines.dpumax[is_1d_line]).all():
            raise RuntimeError("Encountered 1D lines without dpumax set.")

        # Channels, Pipes, Culverts, interpolated nodes
        CONTENT_TYPE_LINEAR = [
            ContentType.TYPE_V2_CHANNEL,
            ContentType.TYPE_V2_PIPE,
            ContentType.TYPE_V2_CULVERT,
        ]
        mask = np.isin(self.nodes.content_type, CONTENT_TYPE_LINEAR)
        # Get the invert level from the connecting line. Each interpolated node has
        # by definition exactly one line going from it. So we can do this:

        sorter = np.argsort(self.lines.line[:, 0])
        idx = np.searchsorted(self.lines.line[sorter, 0], self.nodes.id[mask])
        assert np.all(self.lines.line[sorter, 0][idx] == self.nodes.id[mask])
        self.nodes.dmax[mask] = self.lines.invert_level_start_point[sorter][idx]

        # Connection nodes: logic based on the connected objects
        connection_nodes_module.set_bottom_levels(self.nodes, self.lines)

    def set_initial_waterlevels(self, connection_nodes, channels, pipes, culverts):
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

    def set_obstacles(self, obstacles, levees):
        """Set obstacles on 2D lines by determining intersection between
           line_coords (these must be knows at this point) and obstacle geometry.
           Set kcu to 101 and changes flod and flou to crest_level.

        Also store levees on this grid for later output (and Breach determination)

        Args:
            obstacles (Obstacles)
            levees (Levees)
        """
        self.levees = levees
        obstacles_module.apply_obstacles(self.lines, obstacles, levees)

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

        Only the definitions that are actually used will be used.

        Args:
            definitions (CrossSectionDefinitions)
        """
        cs_in_use = np.union1d(self.lines.cross_id1, self.lines.cross_id2)
        cs_in_use = cs_in_use[cs_in_use != -9999]
        self.cross_sections = definitions.convert(cs_in_use)

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

    def add_0d(self, surfaces: Union[zero_d.Surfaces, zero_d.ImperviousSurfaces]):
        """
        Zero dimension admin derived from 'v2_surfaces' and 'v2_impervious_surfaces'.
        """
        self.surfaces = surfaces.as_grid_surfaces()
        self.surface_maps = surfaces.as_surface_maps(self.nodes)

    def add_breaches(self, connected_points):
        """The breaches are derived from the ConnectedPoints: if a ConnectedPoint
        references a Levee, it will result in a Breach. The Breach gets the properties
        from the Levee (max_breach_depth, material) but these may be unset.
        """
        self.lines.set_line_coords(self.nodes)
        self.lines.fix_line_geometries()
        self.breaches = connected_points.get_breaches(self.lines, self.levees)

    def set_dem_averaged_cells(self, dem_average_areas):
        """Determine which nodes need to be dem averaged during tables preprocessing.

        Args:
            dem_average_areas (DemAverageAreas)
        """
        if len(dem_average_areas) > 0:
            dem_average_area_module.apply_dem_average_areas(
                self.nodes, self.cell_tree, dem_average_areas
            )

    def sort(self):
        """Sort the nodes and lines into the order required by the calculation core.

        Resets the nodes and lines ids to consecutive indices starting from 0.
        Also maps lines.line and pumps.line according to the new nodes indices.

        See NODE_ORDER and LINE_ORDER for the order.
        """
        node_sorter = np.argsort(replace(self.nodes.node_type, NODE_ORDER))
        line_sort_groups = replace(self.lines.kcu, LINE_ORDER)
        line_sort_groups[self.lines.is_1d_boundary == 1] = LINE_ORDER_BOUNDARY_1D
        line_sorter = np.argsort(line_sort_groups)

        # now sort the nodes and lines and reset their ids
        old_node_ids = self.nodes.id.copy()
        self.nodes.reorder(node_sorter)
        self.nodes.id[:] = np.arange(len(self.nodes))
        old_line_ids = self.lines.id.copy()
        self.lines.reorder(line_sorter)
        self.lines.id[:] = np.arange(len(self.lines))

        # create mapping with the node new ids on the position of the old node ids
        new_node_ids = np.empty(old_node_ids[-1] + 1, dtype=self.nodes.id.dtype)
        new_node_ids[old_node_ids[node_sorter]] = self.nodes.id

        # apply the node id mappings to lines.line
        self.lines.line[:] = np.take(new_node_ids, self.lines.line)

        # sort boundary lines so that they internally match the boundary node order
        BOUNDARY_NODE_TYPES = (
            NodeType.NODE_1D_BOUNDARIES,
            NodeType.NODE_2D_BOUNDARIES,
            NodeType.NODE_2D_GROUNDWATER_BOUNDARIES,
        )
        bc_node_ids = self.nodes.id[np.isin(self.nodes.node_type, BOUNDARY_NODE_TYPES)]
        self.lines.sort_by_nodes(bc_node_ids)

        # create mapping with the line new ids on the position of the old line ids
        new_line_ids = np.empty(old_line_ids[-1] + 1, dtype=self.lines.id.dtype)
        new_line_ids[old_line_ids[line_sorter]] = self.lines.id

        # apply the mappings to other datasets that contain references to nodes or lines
        if self.pumps is not None:
            mask = self.pumps.line != -9999
            self.pumps.line[mask] = np.take(new_node_ids, self.pumps.line[mask])
        if self.nodes_embedded is not None:
            self.nodes_embedded.embedded_in = np.take(
                new_node_ids, self.nodes_embedded.embedded_in
            )
        if self.breaches is not None:
            self.breaches.levl = np.take(new_line_ids, self.breaches.levl)

    def finalize(self):
        """Finalize the Grid, computing and setting derived attributes"""
        self.sort()
        self.lines.set_line_coords(self.nodes)
        self.lines.fix_line_geometries()
        self.lines.set_discharge_coefficients()
        self.meta: GridMeta
        if len(self.pumps) > 0:
            self.meta.has_pumpstations = True
        self.meta.has_initial_waterlevels = np.isfinite(
            self.nodes.initial_waterlevel
        ).any()
        self.meta.extent_1d = self.nodes.get_extent_1d()
        self.meta.extent_2d = self.nodes.get_extent_2d()

        if self.surfaces:
            self.meta.zero_dim_extent = self.surfaces.get_extent()

        self.meta.has_1d = (
            self.meta.extent_1d is not None or len(self.nodes_embedded) > 0
        )
        self.meta.has_0d = self.meta.zero_dim_extent is not None
        self.meta.has_2d = self.meta.extent_2d is not None
        self.meta.has_breaches = self.breaches is not None and len(self.breaches) > 0
        if len(self.nodes_embedded) > 0:
            self.meta.has_embedded = True
