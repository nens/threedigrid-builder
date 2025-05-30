from dataclasses import dataclass, fields
from typing import Optional, Tuple

import numpy as np
import shapely
from osgeo import osr
from pyproj import CRS
from pyproj.exceptions import CRSError

import threedigrid_builder
from threedigrid_builder.base import (
    Lines,
    Nodes,
    PointsOnLine,
    Pumps,
    Quarters,
    SurfaceMaps,
    Surfaces,
)
from threedigrid_builder.base.settings import GridSettings, TablesSettings
from threedigrid_builder.constants import ContentType, LineType, NodeType, WKT_VERSION
from threedigrid_builder.exceptions import SchematisationError
from threedigrid_builder.grid import ConnectionNodes, zero_d

from . import connection_nodes as connection_nodes_module
from . import dem_average_area as dem_average_area_module
from . import embedded as embedded_module
from . import groundwater as groundwater_module
from . import initial_waterlevels as initial_waterlevels_module
from .cross_section_definitions import CrossSections
from .linear import BaseLinear
from .lines_1d2d import Lines1D2D
from .obstacles import ObstacleAffectsType, Obstacles
from .potential_breaches import PotentialBreaches, PotentialBreachPoints

osr.UseExceptions()


__all__ = ["Grid", "GridMeta", "QuadtreeStats"]

NODE_ORDER = [
    [NodeType.NODE_2D_OPEN_WATER],
    [NodeType.NODE_2D_GROUNDWATER],
    [NodeType.NODE_1D_NO_STORAGE, NodeType.NODE_1D_STORAGE],
    [NodeType.NODE_2D_BOUNDARIES],
    [NodeType.NODE_2D_GROUNDWATER_BOUNDARIES],
    [NodeType.NODE_1D_BOUNDARIES],
]

LINE_ORDER = [
    [
        LineType.LINE_2D_U,
        LineType.LINE_2D_OBSTACLE_U,
        LineType.LINE_2D_V,
        LineType.LINE_2D_OBSTACLE_V,
        LineType.LINE_2D,
        LineType.LINE_2D_OBSTACLE,
    ],
    [LineType.LINE_2D_VERTICAL],
    [LineType.LINE_2D_GROUNDWATER],
    [
        LineType.LINE_1D_EMBEDDED,
        LineType.LINE_1D_ISOLATED,
        LineType.LINE_1D_CONNECTED,
        LineType.LINE_1D_LONG_CRESTED,
        LineType.LINE_1D_SHORT_CRESTED,
        LineType.LINE_1D_DOUBLE_CONNECTED,
    ],
    [
        LineType.LINE_1D2D_SINGLE_CONNECTED_CLOSED,
        LineType.LINE_1D2D_SINGLE_CONNECTED_OPEN_WATER,
        LineType.LINE_1D2D_DOUBLE_CONNECTED_CLOSED,
        LineType.LINE_1D2D_DOUBLE_CONNECTED_OPEN_WATER,
    ],
    [
        LineType.LINE_1D2D_GROUNDWATER,
    ],
    [
        LineType.LINE_2D_BOUNDARY_WEST,
        LineType.LINE_2D_BOUNDARY_EAST,
        LineType.LINE_2D_BOUNDARY_SOUTH,
        LineType.LINE_2D_BOUNDARY_NORTH,
    ],
    [
        LineType.LINE_2D_GROUNDWATER_BOUNDARY_WEST,
        LineType.LINE_2D_GROUNDWATER_BOUNDARY_EAST,
        LineType.LINE_2D_GROUNDWATER_BOUNDARY_SOUTH,
        LineType.LINE_2D_GROUNDWATER_BOUNDARY_NORTH,
    ],
]
LINE_ORDER_BOUNDARY_1D = 10


@dataclass
class GridMeta:
    """Metadata that needs to end up in the gridadmin file."""

    model_name: str  # name from sqlite globalsettings.name

    grid_settings: GridSettings
    tables_settings: TablesSettings

    epsg_code: Optional[int] = 28992  # SQLite EPSG is ignored if DEM file is present
    crs_wkt: Optional[str] = None
    model_slug: Optional[str] = None  # from repository.slug
    revision_hash: Optional[str] = None  # from repository.revision.hash
    revision_nr: Optional[int] = None  # from repository.revision.number
    threedi_version: Optional[str] = None  # threedi-api version
    threedicore_version: Optional[str] = None  # not used anymore
    threedigrid_builder_version: str = ""  # filled in __post_init__
    threedi_tables_version: Optional[str] = None  # filled in threedi-tables

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
    has_vegetation: bool = False

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
        if not self.crs_wkt and (self.epsg_code is not None):
            try:
                self.crs_wkt = CRS.from_epsg(self.epsg_code).to_wkt(WKT_VERSION)
            except CRSError:
                raise SchematisationError(f"Invalid EPSG code: {self.epsg_code}")


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
        surface_maps: Optional[SurfaceMaps] = None,
        nodes_embedded=None,
        obstacles=None,
        breaches=None,
        meta=None,
        quadtree_stats=None,
        quarters: Optional[Quarters] = None,
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

        if obstacles is not None and not isinstance(obstacles, Obstacles):
            raise TypeError(f"Expected Obstacles instance, got {type(obstacles)}")
        if breaches is not None and not isinstance(breaches, PotentialBreaches):
            raise TypeError(f"Expected Breaches instance, got {type(breaches)}")

        if quarters is None:
            quarters = Quarters(id=[])
        elif not isinstance(quarters, Quarters):
            raise TypeError(f"Expected Quarters instance, got {type(quarters)}")

        self.nodes = nodes
        self.lines = lines
        self.quarters = quarters
        self.meta = meta
        self.surfaces = surfaces
        self.surface_maps = surface_maps
        self.quadtree_stats = quadtree_stats
        self.pumps = pumps
        self.cross_sections = cross_sections
        self.nodes_embedded = nodes_embedded
        self.obstacles = obstacles
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
        for name in ("nodes", "lines", "nodes_embedded", "surfaces", "quarters"):
            new_attrs[name] = getattr(self, name) + getattr(other, name)
        for name in (
            "meta",
            "quadtree_stats",
            "pumps",
            "cross_sections",
            "obstacles",
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
        """A shapely STRtree of the cells in this grid

        The indices in the tree equal the node indices.
        """
        if self._cell_tree is None:
            is_2d_open = self.nodes.node_type == NodeType.NODE_2D_OPEN_WATER
            geoms = np.empty(len(self.nodes), dtype=object)
            geoms[is_2d_open] = shapely.box(*self.nodes.bounds[is_2d_open].T)
            self._cell_tree = shapely.STRtree(geoms)
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
        meta.has_vegetation = s.vegetation_height_type is not None
        return cls(Nodes(id=[]), Lines(id=[]), meta=meta)

    def set_crs(self, crs: CRS):
        """Overwrite the epsg_code and crs_wkt attributes of grid.meta"""
        if crs.is_geographic:
            raise SchematisationError(
                f"A calculation grid cannot have geographic projection (supplied: '{crs.name}')"
            )
        self.meta.crs_wkt = crs.srs
        # We currently need the epsg_code for post-processing; use a low confidence to
        # make that happen. It will be better always to use the wkt.
        epsg_code = crs.to_epsg(min_confidence=20)
        if epsg_code is None:
            # Fallback to GDAL, extract the EPSG code from the (original) WKT
            sr = osr.SpatialReference(crs.srs)
            if sr.GetAuthorityName("PROJCS") != "EPSG":
                raise SchematisationError(
                    f"The supplied DEM file has a non-EPSG projection '{sr.GetName()}'"
                )
            epsg_code = int(sr.GetAuthorityCode("PROJCS"))
        self.meta.epsg_code = epsg_code

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
            - lines.line_coords
            - lines.line_geometries
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
        lines.fix_line_geometries()
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
            - nodes.dmax: from bottom_level (which comes from manhole)
            - nodes.is_manhole: node is a manhole
            - nodes.drain_level: drain_level of associated manhole.
            - nodes.storage_area: area of connection_node.
        """
        return cls(connection_nodes.get_nodes(node_id_counter), Lines(id=[]))

    @classmethod
    def from_linear_objects(
        cls,
        connection_nodes: connection_nodes_module.ConnectionNodes,
        objects: BaseLinear,
        fixed_nodes: PointsOnLine,
        cell_tree: shapely.STRtree,
        global_dist_calc_points: float,
        embedded_cutoff_threshold: float,
        node_id_counter,
        embedded_node_id_counter,
        line_id_counter,
        connection_node_offset: int = 0,
    ):
        """Construct a grid for linear objects (channels, pipes, culverts)

        Args:
            connection_nodes (ConnectionNodes): used to map ids to indices
            objects (Channels, Pipes, Culverts)
            fixed_nodes (PointsOnLine): to optionally fix interpolated node positions
            cell_tree (shapely.STRtree): strtree of the 2D cells (for embedded channels)
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
        nodes = objects.interpolate_nodes(
            node_id_counter, global_dist_calc_points, fixed_nodes
        )
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

    def set_obstacles_2d(
        self,
        obstacles: Obstacles,
    ):
        """Set obstacles on 2D lines by determining intersection between
           line_coords (these must be knows at this point) and obstacle geometry.
           Set kcu to LINE_2D_OBSTACLE and changes flod and flou to crest_level.

        Args:
            obstacles (Obstacles)
            affects_type (ObstacleAffectsType): which affects attribute of obstacles are considered
        """
        line_2d = [LineType.LINE_2D_U, LineType.LINE_2D_V]
        selection = np.where(np.isin(self.lines.kcu, line_2d))[0]
        crest_level = obstacles.compute_dpumax(
            self.lines, where=selection, affects_type=ObstacleAffectsType.AFFECTS_2D
        )[0]
        self.lines.set_2d_crest_levels(crest_level, where=selection)
        self.obstacles = obstacles

    def set_obstacles_1d2d(
        self,
        obstacles: Obstacles,
    ):
        """Set obstacles on 2D lines by determining intersection between
           line_coords (these must be knows at this point) and obstacle geometry.
           Set kcu to LINE_2D_OBSTACLE and changes flod and flou to crest_level.

        Args:
            obstacles (Obstacles)
            affects_type (ObstacleAffectsType): which affects attribute of obstacles are considered
        """
        line_2d = [LineType.LINE_2D_U, LineType.LINE_2D_V]
        selection = np.where(np.isin(self.lines.kcu, line_2d))[0]
        crest_level = obstacles.compute_dpumax(
            self.lines, where=selection, affects_type=ObstacleAffectsType.AFFECTS_2D
        )[0]
        self.lines.set_2d_crest_levels(crest_level, where=selection)
        self.obstacles = obstacles

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

    def set_quarter_administration(self, quadtree):
        """Sets the quarter cell administration for lines and neighbouring calculation cells."""
        self.quarters = quadtree.get_quarters_admin(self.nodes, self.lines)

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

    def add_1d2d_legacy(
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

    def add_1d2d_lines(
        self,
        exchange_lines,
        connection_nodes,
        channels,
        pipes,
        locations,
        culverts,
        potential_breaches,
        line_id_counter,
        node_open_water_detection,
    ):
        """Connect 1D and 2D elements by computing 1D-2D lines.

        Every (double) connected node gets a 1D-2D connection to the cell in which it
        is located. Double connected gives two lines to the same node.

        In addition to id and line attributes, also the kcu (line type) and dpumax
        (bottom level) are computed.

        Sets self.breaches and appends self.lines.
        """
        lines_1d2d = Lines1D2D.create(self.nodes, line_id_counter)
        if len(lines_1d2d) == 0:
            return
        lines_1d2d.assign_line_coords(self.nodes)
        lines_1d2d.assign_connection_nodes_to_channels_from_breaches(
            self.nodes, potential_breaches=potential_breaches
        )
        lines_1d2d.assign_connection_nodes_to_channels_from_lines(
            self.nodes, lines=self.lines
        )
        lines_1d2d.assign_exchange_lines(exchange_lines)
        lines_1d2d.assign_2d_side_from_exchange_lines(exchange_lines)
        lines_1d2d.assign_breaches(self.nodes, potential_breaches)
        lines_1d2d.assign_2d_node(self.cell_tree)
        lines_1d2d.check_unassigned(self.nodes)
        lines_1d2d.set_line_coords(self.nodes)
        lines_1d2d.fix_line_geometries()
        # Go through objects and dispatch to get_1d2d_properties
        node_idx = lines_1d2d.get_1d_node_idx(self.nodes)
        for objects in (channels, connection_nodes, pipes, culverts):
            mask = self.nodes.content_type[node_idx] == objects.content_type
            content_pk = self.nodes.content_pk[node_idx[mask]]
            if node_open_water_detection == 0 and isinstance(objects, ConnectionNodes):
                is_closed = objects.is_channel(content_pk, channels)
            else:
                is_closed = objects.is_closed(content_pk)
            lines_1d2d.assign_kcu(mask, is_closed=is_closed)
        lines_1d2d.assign_dpumax_from_breaches(potential_breaches)
        lines_1d2d.assign_dpumax_from_exchange_lines(exchange_lines)
        lines_1d2d.assign_dpumax_from_obstacles_open(self.obstacles)
        lines_1d2d.assign_dpumax_from_obstacles_closed(self.obstacles)
        # Go through objects and dispatch to get_1d2d_properties
        for objects in (channels, connection_nodes, pipes, culverts):
            mask = self.nodes.content_type[node_idx] == objects.content_type
            content_pk = self.nodes.content_pk[node_idx[mask]]
            lines_1d2d.assign_dpumax(
                mask,
                objects.get_1d2d_exchange_levels(
                    content_pk,
                    s1d=self.nodes.s1d[node_idx[mask]],
                    locations=locations,
                    channels=channels,
                    connection_nodes=connection_nodes,
                ),
            )

        lines_1d2d.assign_ds1d(self.nodes)
        lines_1d2d.assign_ds1d_half()

        self.breaches = lines_1d2d.output_breaches(potential_breaches)
        self.lines += lines_1d2d

    def add_1d2d_groundwater_lines(self, line_id_counter, connection_nodes) -> None:
        """Connect 1D and 2D groundwater elements by computing 1D-2D lines.

        Appends self.lines.
        """
        lines_1d2d_gw = Lines1D2D.create_groundwater(self.nodes, line_id_counter)
        if len(lines_1d2d_gw) == 0:
            return
        lines_1d2d_gw.assign_groundwater_exchange(self.nodes, cn=connection_nodes)
        lines_1d2d_gw.assign_line_coords(self.nodes)
        lines_1d2d_gw.assign_2d_node(self.cell_tree)
        lines_1d2d_gw.check_unassigned(self.nodes, is_groundwater=True)
        lines_1d2d_gw.transfer_2d_node_to_groundwater(self.nodes.n_groundwater_cells)
        self.lines += lines_1d2d_gw

    def set_breach_ids(self, breach_points: PotentialBreachPoints):
        breach_points.assign_to_connection_nodes(self.nodes, self.lines)
        breach_points.match_breach_ids_with_calculation_types(self.nodes)

    def add_0d(self, surfaces: zero_d.Surfaces):
        """
        Zero dimension admin derived from 'surfaces'.
        """
        self.surfaces = surfaces.as_grid_surfaces()
        self.surface_maps = surfaces.as_surface_maps(self.nodes, self.nodes_embedded)

    def add_groundwater(self, has_groundwater_flow, node_id_counter, line_id_counter):
        """Add groundwater nodes and lines.

        - groundwater nodes are copied from 2D open water cells
        - groundwater lines are copied from open waterand obstacle lines
        - vertical lines between open water - and groundwater nodes are created
        """
        gw_nodes, offset = groundwater_module.get_nodes(self.nodes, node_id_counter)
        gw_lines = groundwater_module.get_vertical_lines(
            self.nodes, offset, line_id_counter
        )
        if has_groundwater_flow:
            gw_lines += groundwater_module.get_lines(
                self.lines, offset, line_id_counter
            )

        self.nodes += gw_nodes
        self.lines += gw_lines

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
        # if one of the nodes is null, the numpy id reset will silently overflow and assign a
        # different new id, leading to an incorrect result
        if np.any(self.lines.line == -9999):
            raise ValueError(
                "Some lines are not fully connected to nodes, causing a null value to be set for the line node"
            )

        node_sorter = np.concatenate(
            [np.where(np.isin(self.nodes.node_type, group))[0] for group in NODE_ORDER]
        )
        is_1d_boundary = self.lines.is_1d_boundary == 1
        line_sorter = np.concatenate(
            [
                np.where(np.isin(self.lines.kcu, group) & ~is_1d_boundary)[0]
                for group in LINE_ORDER
            ]
            + [np.where(is_1d_boundary)[0]]
        )
        if len(node_sorter) != len(self.nodes):
            raise RuntimeError("Length mismatch while sorting nodes.")
        if len(line_sorter) != len(self.lines):
            raise RuntimeError("Length mismatch while sorting lines.")

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
        if len(old_line_ids) > 0:
            new_line_ids = np.empty(old_line_ids[-1] + 1, dtype=self.lines.id.dtype)
            new_line_ids[old_line_ids[line_sorter]] = self.lines.id
        else:  # Edge case: no lines present
            new_line_ids = old_line_ids

        # apply the mappings to other datasets that contain references to nodes or lines
        if self.pumps is not None:
            mask = self.pumps.line != -9999
            self.pumps.line[mask] = np.take(new_node_ids, self.pumps.line[mask])
        if self.breaches is not None:
            self.breaches.line_id[:] = np.take(new_line_ids, self.breaches.line_id)
        if self.nodes_embedded is not None:
            self.nodes_embedded.embedded_in[:] = np.take(
                new_node_ids, self.nodes_embedded.embedded_in
            )
        if self.surface_maps is not None:
            self.surface_maps.cci[:] = np.take(new_node_ids, self.surface_maps.cci)

    def finalize(self):
        """Finalize the Grid, computing and setting derived attributes"""
        self.sort()
        self.lines.set_line_coords(self.nodes)
        self.lines.fix_line_geometries()
        self.lines.fix_ds1d()
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
