from threedigrid_builder.base import Lines
from threedigrid_builder.base import Nodes
from threedigrid_builder.interface import GridAdminOut

import h5py
import numpy as np
import pygeos
import pytest
import tempfile


@pytest.fixture(scope="session")
def h5_out():
    nodes = Nodes(
        id=[1, 2, 3],
    )
    lines = Lines(
        id=[1, 2, 3, 4, 5],
        line=[[1, 2], [2, 3], [3, 1], [1, 3], [3, 2]],
        line_geometries=[
            pygeos.linestrings([[1, 1], [2, 2]]),
            pygeos.linestrings([[1, 1], [2, 2], [3, 3]]),
            None,
            None,
            None,
        ],
    )
    quadtree_stats = {
        "lgrmin": 5,
        "kmax": 2,
        "mmax": np.array([2, 4], dtype=np.int32),
        "nmax": np.array([3, 5], dtype=np.int32),
        "dx": np.array([1.0, 2.0], dtype=np.float64),
        "dxp": 0.5,
        "x0p": 10.0,
        "y0p": 10.0,
    }

    with tempfile.NamedTemporaryFile(suffix=".h5") as tmpfile:
        path = tmpfile.name
        with GridAdminOut(path) as out:
            out.write_grid_characteristics(nodes, lines, 28992)
            out.write_grid_counts(nodes, lines)
            out.write_quadtree(quadtree_stats)
            out.write_nodes(nodes)
            out.write_lines(lines)

        with h5py.File(path, "r") as f:
            yield f


# obtained from bergermeer gridadmin.h5, edited:
# - 13280 nodes to 3 nodes
# - int64 to int32
@pytest.mark.parametrize(
    "dataset,shape,dtype",
    [
        ("bottom_level", (4,), "float64"),
        ("calculation_type", (4,), "int32"),
        ("cell_coords", (4, 4), "float64"),
        ("content_pk", (4,), "int32"),
        ("coordinates", (2, 4), "float64"),
        ("display_name", (4,), "|S64"),
        ("drain_level", (4,), "float64"),
        ("id", (4,), "int32"),
        ("initial_waterlevel", (4,), "float64"),
        ("is_manhole", (4,), "int32"),
        ("manhole_indicator", (4,), "int32"),
        ("node_type", (4,), "int32"),
        ("pixel_coords", (4, 4), "int32"),
        ("pixel_width", (4,), "int32"),
        ("shape", (4,), "|S4"),
        ("storage_area", (4,), "float64"),
        ("sumax", (4,), "float64"),
        ("surface_level", (4,), "float64"),
        ("width", (4,), "float64"),
        ("x_coordinate", (4,), "float64"),
        ("y_coordinate", (4,), "float64"),
        ("z_coordinate", (4,), "float64"),
        ("zoom_category", (4,), "int32"),
        ("code", (4,), "|S32"),  # added
    ],
)
def test_write_nodes(h5_out, dataset, shape, dtype):
    assert h5_out["nodes"][dataset].shape == shape
    assert h5_out["nodes"][dataset].dtype == np.dtype(dtype)


# obtained from bergermeer gridadmin.h5, edited:
# - 5 lines to 5 lines
# - int64 to int32
@pytest.mark.parametrize(
    "dataset,shape,dtype",
    [
        ("calculation_type", (6,), "int32"),
        ("code", (6,), "|S32"),
        ("connection_node_end_pk", (6,), "int32"),
        ("connection_node_start_pk", (6,), "int32"),
        ("content_pk", (6,), "int32"),
        ("content_type", (6,), "|S10"),
        ("crest_level", (6,), "float64"),
        ("crest_type", (6,), "int32"),
        ("cross_section_height", (6,), "int32"),
        ("cross_section_shape", (6,), "int32"),
        ("cross_section_width", (6,), "float64"),
        ("discharge_coefficient", (6,), "float64"),
        ("discharge_coefficient_negative", (6,), "float64"),
        ("discharge_coefficient_positive", (6,), "float64"),
        ("display_name", (6,), "|S64"),
        ("dist_calc_points", (6,), "float64"),
        ("friction_type", (6,), "int32"),
        ("friction_value", (6,), "float64"),
        ("id", (6,), "int32"),
        ("invert_level_end_point", (6,), "float64"),
        ("invert_level_start_point", (6,), "float64"),
        ("kcu", (6,), "int32"),
        ("lik", (6,), "int32"),
        ("line", (2, 6), "int32"),
        ("cross_pix_coords", (4, 6), "int32"),
        ("line_coords", (4, 6), "float64"),
        ("line_geometries", (6,), "object"),
        ("material", (6,), "int32"),
        ("sewerage", (6,), "int32"),
        ("sewerage_type", (6,), "int32"),
        ("zoom_category", (6,), "int32"),
        ("ds1d", (6,), "float64"),  # added
        ("dpumax", (6,), "float64"),  # added
        ("flod", (6,), "float64"),  # added
        ("flou", (6,), "float64"),  # added
        ("cross1", (6,), "int32"),  # added
        ("cross2", (6,), "int32"),  # added
        ("cross_weight", (6,), "float64"),  # added
    ],
)
def test_write_lines(h5_out, dataset, shape, dtype):
    assert h5_out["lines"][dataset].shape == shape
    assert h5_out["lines"][dataset].dtype == np.dtype(dtype)


@pytest.mark.parametrize(
    "dataset,shape,dtype",
    [
        ("lgrmin", (), "int32"),
        ("kmax", (), "int32"),
        ("mmax", (2,), "int32"),
        ("nmax", (2,), "int32"),
        ("dx", (2,), "float64"),
        ("dxp", (), "float64"),
        ("x0p", (), "float64"),
        ("y0p", (), "float64"),
    ],
)
def test_write_quadtree(h5_out, dataset, shape, dtype):
    assert h5_out["grid_coordinate_attributes"][dataset].shape == shape
    assert h5_out["grid_coordinate_attributes"][dataset].dtype == np.dtype(dtype)


@pytest.mark.parametrize(
    "dataset,shape,dtype",
    [
        ("infl1d", (), "i4"),
        ("ingrw1d", (), "i4"),
        ("jap1d", (), "i4"),
        ("l1dtot", (), "i4"),
        ("lgutot", (), "i4"),
        ("lgvtot", (), "i4"),
        ("liutot", (), "i4"),
        ("livtot", (), "i4"),
        ("n1dobc", (), "i4"),
        ("n1dtot", (), "i4"),
        ("n2dobc", (), "i4"),
        ("n2dtot", (), "i4"),
        ("ngr2bc", (), "i4"),
    ],
)
def test_write_grid_counts(h5_out, dataset, shape, dtype):
    assert h5_out["meta"][dataset].shape == shape
    assert h5_out["meta"][dataset].dtype == np.dtype(dtype)


@pytest.mark.parametrize(
    "attr,dtype",
    [
        ("epsg_code", "i4"),
        ("has_1d", "i4"),
        ("has_2d", "i4"),
        ("extent_1d", "float64"),
        ("extent_2d", "float64"),
        ("has_interception", "i4"),
        ("has_pumpstations", "i4"),
        ("has_simple_infiltration", "i4"),
        # ('model_name', "S"),  # For later concern.
        # ('model_slug', "S"),  # For later concern.
        # ('revision_hash', "S"),  # For later concern.
        ("revision_nr", "i4"),
        ("threedigrid_builder_version", "i4"),
    ],
)
def test_write_grid_characteristics(h5_out, attr, dtype):
    assert h5_out.attrs[attr].dtype == np.dtype(dtype)
