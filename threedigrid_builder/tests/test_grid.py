from unittest import mock

import numpy as np
import pytest
import shapely
from numpy.testing import assert_array_equal

from threedigrid_builder.base import GridSettings, Lines, Nodes, Pumps, TablesSettings
from threedigrid_builder.constants import (
    ContentType,
    InitializationType,
    LineType,
    NodeType,
)
from threedigrid_builder.grid import (
    ConnectionNodes,
    Grid,
    GridMeta,
    PotentialBreaches,
    QuadtreeStats,
)


@pytest.fixture
def connection_nodes():
    return ConnectionNodes(
        the_geom=shapely.points([(0, 0), (10, 0)]),
        id=np.array([1, 3]),
        code=np.array(["one", "two"]),
    )


@pytest.fixture
def grid():
    return Grid(nodes=Nodes(id=[]), lines=Lines(id=[]))


@pytest.fixture
def meta():
    return GridMeta(
        epsg_code=28992,
        model_name="test-name",
        grid_settings=GridSettings(
            use_2d=True,
            use_1d_flow=True,
            use_2d_flow=True,
            use_0d_inflow=0,
            grid_space=20.0,
            dist_calc_points=25.0,
            kmax=4,
            node_open_water_detection=1,
        ),
        tables_settings=TablesSettings(
            table_step_size=0.05,
            frict_coef=0.03,
            frict_coef_type=9,
        ),
    )


@pytest.fixture
def grid2d(meta):
    quadtree_stats = QuadtreeStats(
        **{
            "lgrmin": 2,
            "kmax": 1,
            "mmax": np.array([1]),
            "nmax": np.array([1]),
            "dx": np.array([2.0]),
            "dxp": 0.5,
            "x0p": 10.0,
            "y0p": 10.0,
        }
    )
    return Grid(
        nodes=Nodes(
            id=[0, 1],
            node_type=NodeType.NODE_2D_OPEN_WATER,
            bounds=[(0, 0, 1, 1), (1, 0, 2, 1)],
        ),
        lines=Lines(id=[0]),
        meta=meta,
        quadtree_stats=quadtree_stats,
    )


@pytest.fixture
def grid1d(meta):
    return Grid(nodes=Nodes(id=[2, 3]), lines=Lines(id=[1]), meta=meta)


@pytest.mark.parametrize(
    "setting,expected_true",
    [
        ("interception_type", "has_interception"),
        ("groundwater_hydro_connectivity_type", "has_groundwater_flow"),
        ("infiltration_rate_type", "has_simple_infiltration"),
        ("groundwater_impervious_layer_level_type", "has_groundwater"),
        ("interflow_type", "has_interflow"),
        ("vegetation_height_type", "has_vegetation"),
    ],
)
def test_from_meta(meta, setting, expected_true):
    setattr(meta.tables_settings, setting, InitializationType.GLOBAL)
    grid = Grid.from_meta(
        epsg_code=28992,
        model_name="test",
        grid_settings=meta.grid_settings,
        tables_settings=meta.tables_settings,
    )
    for attr in (
        "has_interception",
        "has_groundwater_flow",
        "has_simple_infiltration",
        "has_groundwater",
        "has_interflow",
        "has_vegetation",
    ):
        assert getattr(grid.meta, attr) is (attr == expected_true)


def test_set_crs_pyproj(grid2d, crs_wkt_28992):
    grid2d.meta.epsg_code = None
    grid2d.meta.crs_wkt = None

    grid2d.set_crs(crs_wkt_28992)
    assert grid2d.meta.epsg_code == 28992
    assert grid2d.meta.crs_wkt == crs_wkt_28992.srs


def test_set_crs_gdal(grid2d, crs_wkt_28992_legacy):
    grid2d.meta.epsg_code = None
    grid2d.meta.crs_wkt = None
    grid2d.set_crs(crs_wkt_28992_legacy)
    assert grid2d.meta.epsg_code == 28992
    assert grid2d.meta.crs_wkt == crs_wkt_28992_legacy.srs


def test_from_quadtree():
    quadtree = mock.Mock()
    area_mask = mock.Mock()
    counter = mock.Mock()
    nodes = Nodes(id=[])
    lines = Lines(id=[])

    quadtree.origin = (0.0, 0.0)
    quadtree.get_nodes_lines.return_value = nodes, lines

    grid = Grid.from_quadtree(
        quadtree=quadtree,
        area_mask=area_mask,
        node_id_counter=counter,
        line_id_counter=counter,
    )

    assert isinstance(grid, Grid)
    assert grid.nodes is nodes
    assert grid.lines is lines


def test_from_connection_nodes():
    connection_nodes = mock.Mock()
    counter = mock.Mock()
    nodes = Nodes(id=[])

    connection_nodes.get_nodes.return_value = nodes
    grid = Grid.from_connection_nodes(connection_nodes, counter)

    connection_nodes.get_nodes.assert_called_with(counter)

    assert grid.nodes is nodes
    assert len(grid.lines) == 0


def test_concatenate_grid(grid2d, grid1d):
    grid = grid2d + grid1d
    assert grid.meta == grid1d.meta
    assert grid.quadtree_stats == grid2d.quadtree_stats
    assert_array_equal(grid.nodes.id[0:2], grid2d.nodes.id)
    assert_array_equal(grid.nodes.id[2:], grid1d.nodes.id)


@mock.patch("threedigrid_builder.grid.connection_nodes.set_calculation_types")
def test_set_calculation_types(set_calculation_types, grid):
    grid.set_calculation_types()
    set_calculation_types.assert_called_with(grid.nodes, grid.lines)


@mock.patch("threedigrid_builder.grid.connection_nodes.set_bottom_levels")
def test_set_bottom_levels(cn_compute):
    CH = ContentType.TYPE_V2_CHANNEL
    CN = ContentType.TYPE_V2_CONNECTION_NODES
    PI = ContentType.TYPE_V2_PIPE
    CV = ContentType.TYPE_V2_CULVERT

    # set an input grid and mock the compute_weights return value
    nodes = Nodes(
        id=[1, 2, 3, 4, 5, 6, 7],
        content_type=[CH, -9999, CH, PI, CV, CN, CN],
        content_pk=[1, 1, 1, 2, 5, 1, 2],
        s1d=[2.0, 12.0, 21.0, 15.0, 0.5, np.nan, np.nan],
        dmax=[np.nan, 12.0, np.nan, np.nan, np.nan, np.nan, np.nan],
    )
    lines = Lines(
        id=range(7),
        content_type=[CH, CH, CH, PI, PI, CV, CV],
        line=[(6, 1), (1, 3), (3, 7), (6, 4), (4, 7), (6, 5), (5, 7)],
        invert_level_start_point=[1.0, 2.0, 3.0, 4.0, 5.0, 6.0, 7.0],
        invert_level_end_point=[2.0, 3.0, 3.0, 5.0, 6.0, 7.0, 8.0],
    )
    grid = Grid(nodes=nodes, lines=lines)
    grid.set_bottom_levels()

    # node attribute dmax is adapted correctly
    assert_array_equal(grid.nodes.dmax, [2.0, 12.0, 3.0, 5.0, 7.0, np.nan, np.nan])

    # connection node set_bottom_levels was called correctly
    cn_compute.assert_called_with(grid.nodes, grid.lines)


@mock.patch("threedigrid_builder.grid.initial_waterlevels.compute_initial_waterlevels")
def test_set_initial_waterlevels(compute_initial_waterlevels, grid):
    connection_nodes = mock.Mock()
    channels = mock.Mock()
    pipes = mock.Mock()
    culverts = mock.Mock()

    grid.set_initial_waterlevels(connection_nodes, channels, pipes, culverts)

    compute_initial_waterlevels.assert_called_with(
        grid.nodes,
        connection_nodes=connection_nodes,
        channels=channels,
        pipes=pipes,
        culverts=culverts,
    )


@pytest.fixture
def grid_for_sorting():
    return Grid(
        Nodes(
            id=[0, 1, 2, 3, 4, 5],
            node_type=[
                NodeType.NODE_1D_BOUNDARIES,
                NodeType.NODE_2D_OPEN_WATER,
                NodeType.NODE_2D_OPEN_WATER,
                NodeType.NODE_2D_OPEN_WATER,
                NodeType.NODE_1D_STORAGE,
                NodeType.NODE_1D_NO_STORAGE,
            ],
            dmax=[0, 1, 2, 3, 4, 5],
        ),
        Lines(
            id=[0, 1, 6],
            kcu=[
                LineType.LINE_1D_SHORT_CRESTED,
                LineType.LINE_1D_ISOLATED,
                LineType.LINE_1D2D_SINGLE_CONNECTED_OPEN_WATER,
            ],
            is_1d_boundary=[1, -9999, -9999],
            line=[(0, 4), (4, 5), (5, 1)],
            dpumax=[0, 1, 6],
        ),
        pumps=Pumps(id=[0], line=[(4, 5)]),
        nodes_embedded=Nodes(id=[0], embedded_in=[1]),
        breaches=PotentialBreaches(id=[2], line_id=[6]),
    )


def test_sort(grid_for_sorting):
    grid = grid_for_sorting
    grid.sort()

    assert_array_equal(grid.nodes.id, [0, 1, 2, 3, 4, 5])
    assert_array_equal(grid.nodes.dmax, [1, 2, 3, 4, 5, 0])
    assert_array_equal(grid.lines.id, [0, 1, 2])
    assert_array_equal(grid.lines.dpumax, [1, 6, 0])
    assert_array_equal(grid.lines.line, [(3, 4), (4, 0), (5, 3)])
    assert_array_equal(grid.pumps.line, [(3, 4)])
    assert_array_equal(grid.nodes_embedded.embedded_in, [0])
    assert_array_equal(grid.breaches.line_id, [1])


def test_sort_no_lines(grid_for_sorting):
    grid_for_sorting.lines = grid_for_sorting.lines[:0]
    grid_for_sorting.pumps = grid_for_sorting.pumps[:0]
    grid_for_sorting.breaches = grid_for_sorting.breaches[:0]

    grid_for_sorting.sort()


def test_sort_null_lines_err(grid_for_sorting):
    grid_for_sorting.lines.line[0] = (-9999, 2)
    with pytest.raises(ValueError):
        grid_for_sorting.sort()


def test_sort_boundary_conditions():
    grid = Grid(
        Nodes(
            id=[0, 1, 2, 3, 4, 5, 6, 7, 8],
            node_type=[
                NodeType.NODE_2D_OPEN_WATER,
                NodeType.NODE_1D_BOUNDARIES,
                NodeType.NODE_2D_BOUNDARIES,
                NodeType.NODE_2D_BOUNDARIES,
                NodeType.NODE_1D_BOUNDARIES,
                NodeType.NODE_2D_BOUNDARIES,
                NodeType.NODE_1D_NO_STORAGE,
                NodeType.NODE_2D_GROUNDWATER_BOUNDARIES,
                NodeType.NODE_2D_GROUNDWATER_BOUNDARIES,
            ],
            dmax=[0, 1, 2, 3, 4, 5, 6, 7, 8],
        ),
        Lines(
            id=[0, 1, 2, 3, 4, 5, 6],
            kcu=[
                LineType.LINE_1D_ISOLATED,
                LineType.LINE_1D_ISOLATED,
                LineType.LINE_2D_BOUNDARY_EAST,
                LineType.LINE_2D_BOUNDARY_WEST,
                LineType.LINE_2D_BOUNDARY_NORTH,
                LineType.LINE_1D2D_SINGLE_CONNECTED_CLOSED,
                LineType.LINE_2D_GROUNDWATER_BOUNDARY_NORTH,
            ],
            is_1d_boundary=[1, 1] + [-9999] * 5,
            line=[(6, 4), (1, 6), (0, 3), (0, 5), (2, 0), (0, 6), (7, 8)],
            dpumax=[0, 1, 2, 3, 4, 5, 6],
        ),
    )
    grid.sort()

    # expected order:
    #   2D (0), 1D (6)
    #   2D boundaries (2, 3, 5)
    #   2D groundwater (7, 8)
    #   1D boundaries (1, 4)
    assert_array_equal(grid.nodes.id, [0, 1, 2, 3, 4, 5, 6, 7, 8])
    assert_array_equal(grid.nodes.dmax, [0, 6, 2, 3, 5, 7, 8, 1, 4])
    # expected order:
    #   1D2D (5)
    #   2D boundaries (4, 2, 3)
    #   2D Groundwater (6)
    #   1D boundaries (1, 0)
    # note that the 1D and 2D boundaries are sorted internally so that the order matches
    # the order in the corresponding nodes
    assert_array_equal(grid.lines.id, [0, 1, 2, 3, 4, 5, 6])
    assert_array_equal(grid.lines.dpumax, [5, 4, 2, 3, 6, 1, 0])
    # the correct line order can also be seen by the increasing node ids:
    assert_array_equal(
        grid.lines.line, [[0, 1], [2, 0], [0, 3], [0, 4], [5, 6], [7, 1], [1, 8]]
    )


def test_set_cross_sections(grid):
    definitions = mock.Mock()
    grid.lines = Lines(
        id=[0, 1, 2, 3], cross_id1=[2, 5, 2, -9999], cross_id2=[2, -9999, 3, 4]
    )

    grid.set_cross_sections(definitions)

    assert_array_equal(definitions.convert.call_args[0][0], [2, 3, 4, 5])
    assert grid.cross_sections is definitions.convert.return_value
