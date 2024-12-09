import numpy as np
import shapely

from threedigrid_builder.constants import ContentType
from threedigrid_builder.grid.zero_d import (
    SURFACE_CLASS_MAP,
    SURFACE_INCLINATION_MAP,
    SURFACE_TYPE_PROPERTIES,
    SurfaceParams,
    Surfaces,
)


def test_surface_params():
    """
    Check surface params conversions
    """
    surface_class = []
    surface_inclination = []

    expected = {}
    for surface_cls in SURFACE_CLASS_MAP:
        for surface_incl in SURFACE_INCLINATION_MAP:
            surface_class.append(surface_cls)
            surface_inclination.append(surface_incl)
            values_dict = SURFACE_TYPE_PROPERTIES[
                f"{SURFACE_CLASS_MAP[surface_cls]}_{SURFACE_INCLINATION_MAP[surface_incl]}"
            ]
            for k, v in values_dict.items():
                if k not in expected:
                    expected[k] = []
                expected[k].append(v)

    params = SurfaceParams.from_surface_class_and_inclination(
        np.array(surface_class), np.array(surface_inclination)
    )

    for k, v in expected.items():
        value_to_check = getattr(params, k)
        assert np.all(value_to_check == np.array(v, dtype=value_to_check.dtype))


def test_surfaces(grid_all):
    surface = Surfaces(
        id=[1, 2, 3],
        # function=[b"1", b"2", b"3"],
        outflow_delay=[1.0, 2.0, 2.0],
        surface_layer_thickness=[1.0, 2.0, 2.0],
        infiltration=[0.0, 1.0, 1.0],
        max_infiltration_capacity=[1.0, 2.0, 2.0],
        min_infiltration_capacity=[0.0, 1.0, 1.0],
        infiltration_decay_constant=[1.1, 1.2, 1.2],
        infiltration_recovery_constant=[0.9, 1.2, 1.2],
        surface_id=[1, 2, 2],
        code=[b"1", b"2", b"2"],
        display_name=[b"d1", b"d2", b"d2"],
        # nr_of_inhabitants=[1000.0, 2000.0, 2000.0],
        area=[100.0, 200.0, 200.0],
        # dry_weather_flow=[1.0, 2.0, 2.0],
        the_geom=[shapely.points([1.0, 2.0]), None, None],
        connection_node_id=[1, 2, 1],
        connection_node_the_geom=[
            shapely.points([1.0, 2.0]),
            shapely.points([2.0, 3.0]),
            shapely.points([1.0, 2.0]),
        ],
        percentage=[100.0, 50.0, 50.0],
    )

    grid_all.nodes.content_type[1:] = ContentType.TYPE_V2_CONNECTION_NODES.value
    grid_surfaces = surface.as_grid_surfaces()

    expected = (
        ("id", (2,), "int32"),
        ("code", (2,), "object"),
        ("display_name", (2,), "object"),
        ("function", (2,), "object"),
        ("area", (2,), "float64"),
        ("centroid_x", (2,), "float64"),
        ("centroid_y", (2,), "float64"),
        ("dry_weather_flow", (2,), "float64"),
        ("nr_of_inhabitants", (2,), "float64"),
        ("infiltration_flag", (2,), "bool"),
        ("outflow_delay", (2,), "float64"),
        ("storage_limit", (2,), "float64"),
        ("fb", (2,), "float64"),
        ("fe", (2,), "float64"),
        ("ka", (2,), "float64"),
        ("kh", (2,), "float64"),
        ("surface_class", (2,), "object"),
        ("surface_inclination", (2,), "object"),
        ("surface_sub_class", (2,), "object"),
    )

    for name, shape, dtype in expected:
        value = getattr(grid_surfaces, name)
        assert value.shape == shape
        assert value.dtype == np.dtype(dtype)

    # Second surface centroid is based on geometries of both coupled connection_node geoms
    assert np.all(grid_surfaces.centroid_x == np.array([1.0, 1.5]))
    assert np.all(grid_surfaces.centroid_y == np.array([2.0, 2.5]))

    # Check surface maps
    expected = (
        ("fac", (3,), "float64"),
        ("imp", (3,), "int32"),
        ("nxc", (3,), "float64"),
        ("nyc", (3,), "float64"),
        ("pk", (3,), "int32"),
        ("cci", (3,), "int32"),
    )

    grid_surface_maps = surface.as_surface_maps(grid_all.nodes, grid_all.nodes_embedded)

    for name, shape, dtype in expected:
        value = getattr(grid_surface_maps, name)
        assert value.shape == shape
        assert value.dtype == np.dtype(dtype)
