from threedigrid_builder.base import Breaches
from threedigrid_builder.base import Levees
from threedigrid_builder.base import Lines
from threedigrid_builder.constants import ContentType
from threedigrid_builder.constants import LineType
from threedigrid_builder.constants import Material
from threedigrid_builder.grid import ConnectedPoints

import numpy as np
import pygeos
import pytest


@pytest.fixture
def levees():
    return Levees(
        id=[1, 2],
        the_geom=[
            pygeos.linestrings([[0, 0], [0, 10], [10, 10]]),
            pygeos.linestrings([[5, 5], [8, 5]]),
        ],
        max_breach_depth=[4.0, np.nan],
        material=[Material.CLAY, -9999],
    )


@pytest.fixture
def connected_points():
    return ConnectedPoints(
        id=[0, 1, 2],
        levee_id=[1, 1, 2],
        exchange_level=[2.0, np.nan, np.nan],
    )


@pytest.fixture
def lines():
    return Lines(
        id=[0, 1, 2, 3, 4],
        the_geom=[
            pygeos.linestrings([[0, 0], [0, 10], [10, 10]]),
            pygeos.linestrings([[5, 5], [8, 5]]),
        ],
        max_breach_depth=[4.0, np.nan],
        material=[Material.CLAY, -9999],
    )


@pytest.fixture
def lines():
    return Lines(
        id=[0, 1, 2, 3, 4],
        the_geom=[
            pygeos.linestrings([[0, 0], [0, 10], [10, 10]]),
            pygeos.linestrings([[5, 5], [8, 5]]),
        ],
        max_breach_depth=[4.0, np.nan],
        material=[Material.CLAY, -9999],
    )
