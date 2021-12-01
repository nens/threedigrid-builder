from threedigrid_builder.base import Breaches
from threedigrid_builder.base import GridSettings
from threedigrid_builder.base import Levees
from threedigrid_builder.base import Lines
from threedigrid_builder.base import Nodes
from threedigrid_builder.base import Pumps
from threedigrid_builder.base import TablesSettings
from threedigrid_builder.constants import NodeType
from threedigrid_builder.grid import CrossSections
from threedigrid_builder.grid import Grid
from threedigrid_builder.grid import GridMeta
from threedigrid_builder.grid import QuadtreeStats
from threedigrid_builder.interface.db import SQLite

import numpy as np
import os
import pathlib
import pygeos
import pytest


data_path = pathlib.Path(__file__).resolve().parents[0] / "data"


@pytest.fixture(scope="session")
def db():
    """Yields a threedigrid_builder.interface.db.SQLite object with access
    to the test v2_bergermeer.sqlite."""
    sqlite_path = data_path / "v2_bergermeer.sqlite"
    if not os.path.isfile(sqlite_path):
        pytest.skip("sample sqlite is not available", allow_module_level=True)
    return SQLite(sqlite_path)


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
        coordinates=[(1, 1), (2, 2), (3, 3)],
        bounds=[(0, 0, 1, 1), (0, 0, 0, 0), (0, 0, 0, 0)],
        node_type=[NodeType.NODE_2D_OPEN_WATER] + 2 * [NodeType.NODE_1D_NO_STORAGE],
    )
    lines = Lines(
        id=[0, 1, 2, 3, 4],
        dpumax=[1.2, 2.2, 3.3, 4.2, 5.1],
        line=[[0, 1], [1, 2], [2, 0], [0, 2], [2, 1]],
        line_geometries=[
            pygeos.linestrings([[1, 1], [2, 2]]),
            pygeos.linestrings([[1, 1], [2, 2], [3, 3]]),
            None,
            None,
            None,
        ],
        cross1=[0, -9999, 1, 1, 2],
        cross2=[-9999, -9999, -9999, -9999, 3],
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
    levees = Levees(
        id=[0, 1],
        the_geom=[
            pygeos.linestrings([[1, 1], [2, 2]]),
            pygeos.linestrings([[1, 1], [2, 2], [3, 3]]),
        ],
    )
    breaches = Breaches(
        id=[0, 1],
        coordinates=[[0, 0], [1, 1]],
        levl=[4, 3],
        levee_id=[1, 0],
    )
    return Grid(
        nodes,
        lines,
        pumps,
        cross_sections,
        nodes_embedded,
        levees,
        breaches,
        meta,
        quadtree_stats,
    )
