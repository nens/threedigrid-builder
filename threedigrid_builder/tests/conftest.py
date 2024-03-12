import os
import pathlib
import shutil

import numpy as np
import pytest
import shapely
from pyproj import CRS

from threedigrid_builder.base import (
    GridSettings,
    Lines,
    Nodes,
    Pumps,
    Surfaces,
    TablesSettings,
)
from threedigrid_builder.base.surfaces import SurfaceMaps
from threedigrid_builder.constants import ContentType, NodeType
from threedigrid_builder.grid import (
    CrossSections,
    Grid,
    GridMeta,
    Obstacles,
    PotentialBreaches,
    QuadtreeStats,
)
from threedigrid_builder.interface.db import SQLite

data_path = pathlib.Path(__file__).resolve().parents[0] / "data"


@pytest.fixture(scope="session")
def db(tmp_path_factory):
    """Yields a threedigrid_builder.interface.db.SQLite object with access
    to the test v2_bergermeer.sqlite."""
    fn = tmp_path_factory.mktemp("data") / "v2_bergermeer.gpkg"
    sqlite_path = data_path / "v2_bergermeer.gpkg"
    shutil.copyfile(sqlite_path, fn)
    if not os.path.isfile(fn):
        pytest.skip("sample sqlite is not available", allow_module_level=True)
    return SQLite(fn, upgrade=True)


@pytest.fixture
def dem_path():
    """Yields a dem path"""
    return data_path / "dem_test_5m.tif"


@pytest.fixture(scope="session")
def grid_all():
    """A Grid with all features"""
    nodes = Nodes(
        id=[0, 1, 2],
        dmax=[1.2, 2.2, 3.3],
        content_pk=[0, 1, 2],
        coordinates=[(1, 1), (2, 2), (3, 3)],
        bounds=[(0, 0, 1, 1), (0, 0, 0, 0), (0, 0, 0, 0)],
        node_type=[NodeType.NODE_2D_OPEN_WATER] + 2 * [NodeType.NODE_1D_NO_STORAGE],
        display_name=[b"foo", "金蟾", None],
    )
    lines = Lines(
        id=[0, 1, 2, 3, 4],
        dpumax=[1.2, 2.2, 3.3, 4.2, 5.1],
        line=[[0, 1], [1, 2], [2, 0], [0, 2], [2, 1]],
        line_geometries=[
            shapely.linestrings([[1, 1], [2, 2]]),
            shapely.linestrings([[1, 1], [2, 2], [3, 3]]),
            None,
            shapely.linestrings([[0, 0], [10, 0]]),
            shapely.linestrings([[0, 10], [0, 0]]),
        ],
        cross_id1=[4, -9999, 3, 3, 4],
        cross_id2=[-9999, -9999, -9999, -9999, 6],
        content_type=[
            -9999,
            -9999,
            -9999,
            ContentType.TYPE_V2_BREACH,
            ContentType.TYPE_V2_EXCHANGE_LINE,
        ],
        content_pk=[-9999, -9999, -9999, 1, 15],
        ds1d_half=[0.5, 1.0, 1.0, 5.0, 8.0],
    )
    pumps = Pumps(
        id=[0, 1, 2],
        capacity=[0.1, 1.2, 2.5],
        line=[[0, 1], [1, 2], [2, 3]],
    )
    meta = GridMeta(
        epsg_code=28992,
        model_name="test-name",
        grid_settings=GridSettings(
            use_2d=True,
            use_1d_flow=True,
            use_2d_flow=True,
            use_0d_inflow=True,
            grid_space=20.0,
            dist_calc_points=25.0,
            kmax=4,
        ),
        tables_settings=TablesSettings(
            table_step_size=0.05,
            frict_coef=0.03,
            frict_coef_type=9,
        ),
    )
    quadtree_stats = QuadtreeStats(
        **{
            "lgrmin": 5,
            "kmax": 2,
            "mmax": np.array([2, 4], dtype=np.int32),
            "nmax": np.array([3, 5], dtype=np.int32),
            "dx": np.array([1.0, 2.0], dtype=np.float64),
            "dxp": 0.5,
            "x0p": 10.0,
            "y0p": 10.0,
        }
    )
    cross_sections = CrossSections(
        id=[0, 1, 2],
        content_pk=[3, 4, 6],
        width_1d=[0.2, 1.5, 3.1],
        count=[4, 2, -9999],
        offset=[0, 4, -9999],
    )
    cross_sections.tables = np.random.random((6, 2))
    nodes_embedded = Nodes(
        id=[0, 1],
        embedded_in=[1, 2],
        dmax=[2.3, 0.2],
    )
    obstacles = Obstacles(
        id=[0, 1],
        the_geom=[
            shapely.linestrings([[1, 1], [2, 2], [4, 4]]),
            shapely.linestrings([[1, 1], [2, 2], [3, 3]]),
        ],
    )
    breaches = PotentialBreaches(
        id=[0, 1],
        content_pk=[1, -9999],
        line_id=[3, 4],
        the_geom=shapely.points([[0, 2], [1, 1]]),
        maximum_breach_depth=[1.2, np.nan],
        levee_material=[1, -9999],
        code=["a", None],
        display_name=["aa", None],
    )
    surfaces = Surfaces(
        id=[0, 1],
        code=[b"1", b"2"],
        display_name=[b"d1", b"d2"],
        function=[b"f1", b"f2"],
        area=[1.0, 2.0],
        centroid_x=[1.0, 2.0],
        centroid_y=[1.0, 2.0],
        dry_weather_flow=[1.2, 1.2],
        nr_of_inhabitants=[1000.0, 2000.0],
        infiltration_flag=[True, False],
        outflow_delay=[1.1, 1.2],
        storage_limit=[10.0, 20.0],
        fb=[1.0, 0.0],
        fe=[1.0, 0.0],
        ka=[1.0, 0.0],
        kh=[1.0, 0.0],
        surface_class=None,
        surface_inclination=None,
        surface_sub_class=None,
    )

    surface_maps = SurfaceMaps(
        id=[1, 2, 3],
        imp=[1, 2, 1],
        fac=[1.0, 0.0, 0.5],
        nxc=[1.1, 2.1, 1.3],
        nyc=[1.3, 2.3, 1.5],
        pk=[1, 2, 3],
        cci=[1, 2, 2],
    )

    return Grid(
        nodes,
        lines,
        pumps,
        cross_sections,
        surfaces,
        surface_maps,
        nodes_embedded,
        obstacles,
        breaches,
        meta,
        quadtree_stats,
    )


@pytest.fixture
def crs_wkt_28992():
    """A current CRS matching EPSG:28992"""
    return CRS(CRS.from_epsg(28992).to_wkt())


@pytest.fixture
def crs_wkt_28992_legacy():
    """A non-current CRS matching the test DEM file"""
    return CRS.from_wkt(
        'PROJCS["Amersfoort / RD New",GEOGCS["Amersfoort",DATUM["Amersfoort",'
        'SPHEROID["Bessel 1841",6377397.155,299.1528128,AUTHORITY["EPSG","7004"]],'
        "TOWGS84[565.2369,50.0087,465.658,-0.406857,0.350733,-1.87035,4.0812],"
        'AUTHORITY["EPSG","6289"]],PRIMEM["Greenwich",0,AUTHORITY["EPSG","8901"]],'
        'UNIT["degree",0.0174532925199433,AUTHORITY["EPSG","9122"]],AUTHORITY["EPSG","4289"]],'
        'PROJECTION["Oblique_Stereographic"],PARAMETER["latitude_of_origin",52.1561605555556],'
        'PARAMETER["central_meridian",5.38763888888889],PARAMETER["scale_factor",0.9999079],'
        'PARAMETER["false_easting",155000],PARAMETER["false_northing",463000],'
        'UNIT["metre",1,AUTHORITY["EPSG","9001"]],AXIS["Easting",EAST],AXIS["Northing",NORTH],'
        'AUTHORITY["EPSG","28992"]]'
    )
