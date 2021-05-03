from threedigrid_builder.constants import CalculationType
from threedigrid_builder.constants import CrossSectionShape
from threedigrid_builder.constants import FrictionType
from threedigrid_builder.constants import ManholeIndicator
from threedigrid_builder.grid import Channels
from threedigrid_builder.grid import ConnectionNodes
from threedigrid_builder.grid import CrossSectionDefinitions
from threedigrid_builder.grid import CrossSectionLocations
from threedigrid_builder.interface import SQLite
from unittest import mock

import numpy as np
import pygeos


def test_init():
    path = "/some/path"

    with mock.patch("threedigrid_builder.interface.db.ThreediDatabase") as db:
        sqlite = SQLite(path)

    db.assert_called_with(
        connection_settings={"db_path": path, "db_file": path}, db_type="spatialite"
    )

    assert sqlite.db is db.return_value


def test_get_channels(db):
    channels = db.get_channels()
    assert isinstance(channels, Channels)

    # some test samples
    assert len(channels.id) == 1175
    assert (
        pygeos.to_wkt(channels.the_geom[0])
        == "LINESTRING (109798 518896, 109815 518823, 109818 518812, 109822 518800, 109834 518750)"
    )
    assert channels.id[11] == 12
    assert channels.code[42] == "151"
    assert channels.calculation_type[494] == CalculationType.CONNECTED
    assert channels.dist_calc_points[580] == 100.0
    assert channels.connection_node_start_id[536] == 1377
    assert channels.connection_node_end_id[1163] == 1056
    assert channels.id[1174] == 666668899


def test_get_connection_nodes(db):
    connection_nodes = db.get_connection_nodes()
    assert isinstance(connection_nodes, ConnectionNodes)

    # some test samples
    assert len(connection_nodes.id) == 1360
    assert pygeos.to_wkt(connection_nodes.the_geom[0]) == "POINT (110404 517792)"
    assert connection_nodes.id[11] == 12
    assert connection_nodes.storage_area[39] == 0.64
    assert np.isnan(connection_nodes.storage_area[49])
    assert connection_nodes.code[494] == ""
    # manhole fields
    assert connection_nodes.manhole_id[10] == 11
    assert connection_nodes.manhole_id[100] == -9999
    assert connection_nodes.calculation_type[1] == CalculationType.CONNECTED
    assert connection_nodes.manhole_indicator[6] == ManholeIndicator.OUTLET
    assert connection_nodes.bottom_level[9] == -3.51
    assert connection_nodes.drain_level[1] == -0.82
    assert connection_nodes.surface_level[35] == -0.54
    assert connection_nodes.manhole_shape[40] == "00"
    assert connection_nodes.manhole_width[32] == 0.8


def test_get_cross_section_definitions(db):
    definitions = db.get_cross_section_definitions()
    assert isinstance(definitions, CrossSectionDefinitions)

    # some test samples
    assert len(definitions.id) == 11
    assert definitions.id[8] == 97
    assert definitions.shape[7] == CrossSectionShape.RECTANGLE
    assert definitions.height[10] == 0.4
    assert definitions.width[2] == 0.315


def test_get_cross_section_locations(db):
    locations = db.get_cross_section_locations()
    assert isinstance(locations, CrossSectionLocations)

    # some test samples
    assert len(locations.id) == 1175
    assert pygeos.to_wkt(locations.the_geom[96]) == "POINT (111104 521655)"
    assert locations.id[11] == 12
    assert locations.definition_id[365] == 98
    assert locations.channel_id[448] == 452
    assert locations.reference_level[691] == -3.0
    assert locations.bank_level[995] == -1.7
    assert locations.friction_type[1103] == FrictionType.MANNING
    assert locations.friction_value[1103] == 0.03


def test_get_grid_refinements(db):
    grid_refinements = db.get_grid_refinements()

    # some test samples
    assert len(grid_refinements.id) == 6
    assert (
        pygeos.to_wkt(grid_refinements.the_geom[3])
        == "LINESTRING (110173 517604, 110327 517527, 110461 517809, 110249 517909, 110147 517700, 110304 517616, 110368 517765, 110280 517816, 110242 517726, 110300 517701, 110316 517749)"
    )
    assert (
        pygeos.to_wkt(grid_refinements.the_geom[4])
        == "POLYGON ((108334 517481, 108701 517460, 108686 517171, 108324 517197, 108334 517481))"
    )
    assert grid_refinements.refinement_level[1] == 1
    assert grid_refinements.display_name[3] == "riolering"
    assert grid_refinements.display_name[5] == "test_polygon2"
    assert grid_refinements.code[5] == "2"
