from unittest import mock

import numpy as np
import pytest
import shapely
from shapely.testing import assert_geometries_equal
from threedi_schema import constants

from threedigrid_builder.base import GridSettings, Pumps, TablesSettings
from threedigrid_builder.constants import (
    BoundaryType,
    CalculationType,
    CrossSectionShape,
    FrictionType,
    InitializationType,
)
from threedigrid_builder.exceptions import SchematisationError
from threedigrid_builder.grid import (
    BoundaryConditions1D,
    BoundaryConditions2D,
    Channels,
    ConnectionNodes,
    CrossSectionDefinitions,
    CrossSectionLocations,
    Culverts,
    Obstacles,
    Orifices,
    Pipes,
    Weirs,
)
from threedigrid_builder.interface import SQLite
from threedigrid_builder.interface.db import map_cross_section_definition


def test_init(tmp_path):
    path = tmp_path / "some.sqlite"
    path.touch()

    with mock.patch(
        "threedigrid_builder.interface.db.ThreediDatabase"
    ) as db, mock.patch.object(SQLite, "get_version") as get_version:
        get_version.return_value = 229
        sqlite = SQLite(path)

    db.assert_called_with(path)

    assert sqlite.db is db.return_value


def test_init_no_file(tmp_path):
    path = tmp_path / "some.sqlite"

    with pytest.raises(FileNotFoundError):
        SQLite(path)


def test_init_bad_version(tmp_path):
    path = tmp_path / "some.sqlite"
    path.touch()

    with mock.patch(
        "threedigrid_builder.interface.db.ThreediDatabase"
    ), mock.patch.object(SQLite, "get_version") as get_version:
        with pytest.raises(SchematisationError):
            get_version.return_value = 171
            SQLite(path)


def test_get_version(db):
    assert db.get_version() == 300


def test_get_boundary_conditions_1d(db):
    boundary_conditions_1d = db.get_boundary_conditions_1d()
    assert isinstance(boundary_conditions_1d, BoundaryConditions1D)

    # some test samples
    assert len(boundary_conditions_1d) == 4
    assert boundary_conditions_1d.id[1] == 2
    assert boundary_conditions_1d.type[2] == BoundaryType.DISCHARGE
    assert boundary_conditions_1d.connection_node_id[3] == 59


def test_get_boundary_conditions_2d(db):
    boundary_conditions_2d = db.get_boundary_conditions_2d()
    assert isinstance(boundary_conditions_2d, BoundaryConditions2D)

    assert len(boundary_conditions_2d) == 0


def test_get_channels(db):
    channels = db.get_channels()
    assert isinstance(channels, Channels)

    # some test samples
    assert len(channels.id) == 1175
    assert_geometries_equal(
        channels.the_geom[0],
        shapely.from_wkt(
            "LINESTRING (109798 518896, 109815 518823, 109818 518812, 109822 518800, 109834 518750)"
        ),
        tolerance=1,
    )
    assert channels.id[11] == 12
    assert channels.code[42] == "151"
    assert channels.calculation_type[494] == CalculationType.CONNECTED
    assert channels.dist_calc_points[580] == 100.0
    assert channels.connection_node_start_id[536] == 1377
    assert channels.connection_node_end_id[1163] == 1056
    assert channels.id[1174] == 666668899
    assert channels.display_name[33] == "801"


def test_get_connection_nodes(db):
    connection_nodes = db.get_connection_nodes()
    assert isinstance(connection_nodes, ConnectionNodes)

    # some test samples
    assert len(connection_nodes.id) == 1360
    assert_geometries_equal(
        connection_nodes.the_geom[0],
        shapely.from_wkt("POINT (110404 517792)"),
        tolerance=1,
    )
    assert connection_nodes.id[11] == 12
    assert connection_nodes.storage_area[39] == 0.64
    assert np.isnan(connection_nodes.storage_area[49])
    assert connection_nodes.code[494] == ""
    assert connection_nodes.initial_waterlevel[169] == -0.4
    # manhole fields
    assert connection_nodes.calculation_type[1] == CalculationType.CONNECTED
    assert connection_nodes.manhole_indicator[6] == 1
    assert connection_nodes.bottom_level[9] == -3.51
    assert np.isnan(connection_nodes.bottom_level[100])
    assert connection_nodes.drain_level[1] == -0.82
    assert connection_nodes.surface_level[35] == -0.54
    assert connection_nodes.display_name[33] == "71512"


def test_get_cross_section_definitions(db):
    definitions = db.get_cross_section_definitions()
    assert isinstance(definitions, CrossSectionDefinitions)

    # some test samples
    assert len(definitions.id) == 1365
    assert definitions.id[8] == 8
    assert definitions.shape[7] == CrossSectionShape.RECTANGLE
    assert definitions.height[10] == 0.4
    assert definitions.width[2] == 3
    mask_tabulated_yz = (
        definitions.shape == constants.CrossSectionShape.TABULATED_YZ.value
    )
    idx_tabulated_yz = definitions.id_to_index(definitions.id[mask_tabulated_yz])[0]
    assert (
        definitions.cross_section_table[idx_tabulated_yz]
        == "0,2\n1,1\n2,0.5\n3,0\n4,1\n5,1.5\n6,2"
    )
    assert (
        definitions.cross_section_vegetation_table[idx_tabulated_yz]
        == "10,0.2,1,7\n5,0.1,0.5,7\n5,0.1,0.5,7\n5,0.3,0.5,7\n15,0.4,1,7\n20,0.1,1.5,7"
    )


def test_get_cross_section_locations(db):
    locations = db.get_cross_section_locations()
    assert isinstance(locations, CrossSectionLocations)

    # some test samples
    assert len(locations.id) == 1175
    assert_geometries_equal(
        locations.the_geom[96], shapely.from_wkt("POINT (111104 521655)"), tolerance=1
    )
    assert locations.id[11] == 12
    assert locations.channel_id[448] == 452
    assert locations.reference_level[691] == -3.0
    assert locations.bank_level[995] == -1.7
    assert locations.friction_type[1103] == FrictionType.MANNING
    assert locations.friction_value[1103] == 0.03


def test_get_grid_refinements(db):
    grid_refinements = db.get_grid_refinements()

    # some test samples
    assert len(grid_refinements.id) == 6
    assert_geometries_equal(
        grid_refinements.the_geom[3],
        shapely.from_wkt(
            "LINESTRING (110173 517604, 110327 517527, 110461 517809, 110249 517909, 110147 517700, 110304 517616, 110368 517765, 110280 517816, 110242 517726, 110300 517701, 110316 517749)"
        ),
        tolerance=1,
    )
    assert_geometries_equal(
        grid_refinements.the_geom[4],
        shapely.from_wkt(
            "POLYGON ((108334 517481, 108701 517460, 108686 517171, 108324 517197, 108334 517481))"
        ),
        tolerance=1,
    )
    assert grid_refinements.refinement_level[1] == 1
    assert grid_refinements.display_name[3] == "riolering"
    assert grid_refinements.display_name[5] == "test_polygon2"
    assert grid_refinements.code[5] == "2"


def test_get_obstacles(db):
    obstacles = db.get_obstacles()
    assert isinstance(obstacles, Obstacles)

    assert len(obstacles) == 3
    assert obstacles.id[1] == 2
    assert_geometries_equal(
        shapely.get_point(obstacles.the_geom[0], 0),
        shapely.from_wkt("POINT (110241 519070)"),
        tolerance=1,
    )
    assert obstacles.crest_level[1] == 0.0


def test_get_pipes(db):
    pipes = db.get_pipes()
    assert isinstance(pipes, Pipes)

    # some test samples
    assert len(pipes.id) == 42
    assert pipes.id[1] == 2
    assert pipes.code[2] == "71022_71023"
    assert pipes.calculation_type[3] == CalculationType.ISOLATED
    assert pipes.connection_node_start_id[4] == 37
    assert pipes.connection_node_end_id[5] == 35
    assert pipes.invert_level_start_point[16] == -3.91
    assert pipes.invert_level_end_point[19] == -3.62
    assert pipes.sewerage_type[24] == 2
    assert pipes.friction_type[28] == FrictionType.MANNING
    assert pipes.friction_value[36] == 0.0145
    assert pipes.display_name[33] == "71518_71517"
    assert_geometries_equal(
        pipes.the_geom[0],
        shapely.from_wkt("LINESTRING (110267.3 517868.8, 110264.3 517863.5)"),
        tolerance=1,
    )


def test_get_settings(db):
    result = db.get_settings()
    assert result["model_name"] == "simple_infil_no_grndwtr"

    g = result["grid_settings"]
    assert isinstance(g, GridSettings)
    # some settings copied over from SQLite:
    assert g.use_1d_flow is True
    assert g.use_2d_flow is True
    assert g.grid_space == 20.0
    assert g.dist_calc_points == 15.0
    assert g.kmax == 4
    assert g.embedded_cutoff_threshold == 0.05
    assert g.max_angle_1d_advection == 0.4 * np.pi
    # use_2d is based on the presence of dem_file:
    assert g.use_2d is True

    s = result["tables_settings"]
    assert isinstance(s, TablesSettings)
    assert s.table_step_size == 0.05
    assert s.frict_type == 2
    assert s.frict_coef == 0.03
    assert s.frict_coef_type == InitializationType.GLOBAL
    assert s.interception_global == 100.0
    assert s.interception_type == InitializationType.NO_AGG
    assert s.table_step_size_1d == 0.05
    assert s.maximum_table_step_size == 5.0
    assert s.manhole_storage_area is None

    # groundwater settings
    assert s.groundwater_impervious_layer_level == -5.0
    assert s.groundwater_impervious_layer_level_type == InitializationType.GLOBAL
    assert s.phreatic_storage_capacity == 0.25
    assert s.phreatic_storage_capacity_type == InitializationType.GLOBAL
    assert s.equilibrium_infiltration_rate == 100.0
    assert s.equilibrium_infiltration_rate_type == InitializationType.GLOBAL
    assert s.initial_infiltration_rate == 300.0
    assert s.initial_infiltration_rate_type == InitializationType.GLOBAL
    assert s.infiltration_decay_period == 0.1
    assert s.infiltration_decay_period_type == InitializationType.GLOBAL
    assert s.groundwater_hydro_connectivity == 1.0
    assert s.groundwater_hydro_connectivity_type == InitializationType.GLOBAL

    # there are interflow settings, but most are unset
    assert s.interflow_type == 0  # means: no interflow
    assert s.porosity is None
    assert s.porosity_type is None
    assert s.porosity_layer_thickness is None
    assert s.impervious_layer_elevation is None
    assert s.hydraulic_conductivity is None
    assert s.hydraulic_conductivity_type is None

    # there are no infiltration settings
    assert s.infiltration_rate is None
    assert s.infiltration_rate_type is None
    assert s.infiltration_surface_option == 0
    assert s.max_infiltration_capacity is None
    assert s.max_infiltration_capacity_type is None

    # there are no vegetation_drag settings
    assert s.vegetation_height is None
    assert s.vegetation_height_type is None
    assert s.vegetation_stem_count is None
    assert s.vegetation_stem_count_type is None
    assert s.vegetation_stem_diameter is None
    assert s.vegetation_stem_diameter_type is None
    assert s.vegetation_drag_coefficient is None
    assert s.vegetation_drag_coefficient is None


def test_get_pumps(db):
    pumps = db.get_pumps()
    assert isinstance(pumps, Pumps)
    # some test samples
    assert len(pumps) == 19
    assert pumps.id[11] == 13
    assert pumps.code[0] == "Rioolgemaal"
    assert pumps.capacity[12] == 0.288
    assert pumps.connection_node_start_id[13] == 1006
    assert pumps.connection_node_end_id[0] == -9999  # NULL handling
    assert pumps.connection_node_end_id[2] == 218
    assert pumps.type_[5] == 1
    assert pumps.start_level[0] == -4.0
    assert pumps.lower_stop_level[18] == -1.9
    assert np.isnan(pumps.upper_stop_level[15])
    assert pumps.display_name[10] == "KGM-JL-18"


def test_get_culverts(db):
    culverts = db.get_culverts()
    assert isinstance(culverts, Culverts)

    # some test samples
    assert len(culverts) == 92

    assert culverts.id[89] == 2000112
    assert culverts.code[35] == "500"
    assert_geometries_equal(
        culverts.the_geom[36],
        shapely.from_wkt("LINESTRING (108351 516428, 108357 516430)"),
        tolerance=1,
    )
    assert culverts.dist_calc_points[52] == 20.0
    assert culverts.connection_node_start_id[45] == 1220
    assert culverts.connection_node_end_id[61] == 1620
    assert culverts.calculation_type[4] == CalculationType.ISOLATED
    assert culverts.invert_level_start_point[21] == -2.39
    assert culverts.invert_level_end_point[83] == -3.28
    assert culverts.discharge_coefficient_negative[0] == 0.8
    assert culverts.discharge_coefficient_positive[0] == 0.8
    assert culverts.friction_type[28] == FrictionType.MANNING
    assert culverts.friction_value[36] == 0.03
    assert culverts.display_name[32] == "KDU-Q-18482"


def test_get_orifices(db):
    orifices = db.get_orifices()
    assert isinstance(orifices, Orifices)

    # no orifices in test dataset
    assert len(orifices) == 0


def test_get_weirs(db):
    weirs = db.get_weirs()
    assert isinstance(weirs, Weirs)

    # some test samples
    assert len(weirs) == 56

    assert weirs.id[13] == 17
    assert weirs.code[11] == "317"
    assert weirs.connection_node_start_id[40] == 1643
    assert weirs.connection_node_end_id[7] == 394
    assert weirs.crest_level[26] == -2.508
    assert weirs.crest_type[0] == CalculationType.SHORT_CRESTED
    assert weirs.discharge_coefficient_negative[0] == 0.0
    assert weirs.discharge_coefficient_positive[0] == 0.8
    assert weirs.friction_type[28] == FrictionType.MANNING
    assert weirs.friction_value[36] == 0.03
    assert weirs.display_name[33] == "KST-JL-76"
    assert weirs.sewerage[0] == 1
    assert_geometries_equal(
        weirs.the_geom[0],
        shapely.from_wkt("LINESTRING (110278.3 517669.1, 110276.3 517669.8)"),
        tolerance=1,
    )


def test_get_dem_average(db):
    dem_avg_areas = db.get_dem_average_areas()
    # no dem_average_areas in test dataset
    assert len(dem_avg_areas) == 0


def test_get_surfaces(db):
    surfaces = db.get_surfaces()
    # No surfaces in test dataset
    assert len(surfaces) == 0


def test_get_windshieldings(db):
    windshieldings = db.get_windshieldings()
    # No windshielding in test dataset
    assert len(windshieldings) == 0


def test_get_potential_breaches(db):
    potential_breaches = db.get_potential_breaches()
    assert len(potential_breaches) == 2


def test_get_exchange_lines(db):
    exchange_lines = db.get_exchange_lines()
    # No exchange lines in test dataset
    assert len(exchange_lines) == 0


def test_map_cross_section_definition():
    mapping = {"pipe": {0: 0, 1: 1}, "weir": {1: 0}, "cross_section_location": {0: 10}}
    pipes = Pipes(id=[0, 1], cross_section_definition_id=[10, 20])
    weirs = Weirs(id=[0, 1], cross_section_definition_id=[10, 20])
    cross_section_locations = CrossSectionLocations(id=[0], definition_id=[100])
    objects = [pipes, weirs, cross_section_locations]
    map_cross_section_definition(objects, mapping)
    np.testing.assert_array_equal(pipes.cross_section_definition_id, [0, 1])
    np.testing.assert_array_equal(weirs.cross_section_definition_id, [10, 0])
    np.testing.assert_array_equal(cross_section_locations.definition_id, [10])
