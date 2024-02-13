from copy import copy

import h5py
import numpy as np
import pytest
from numpy.testing import assert_equal

from threedigrid_builder.interface import GridAdminOut


@pytest.fixture(scope="session")
def h5_out(tmpdir_factory, grid_all):
    path = tmpdir_factory.mktemp("h5") / "gridadmin.h5"

    grid_all.meta.has_0d = True

    with GridAdminOut(path) as out:
        out.write(grid_all)

    with h5py.File(path, "r") as f:
        yield f


@pytest.fixture(scope="session")
def h5_out_1d(tmpdir_factory, grid_all):
    path = tmpdir_factory.mktemp("h5") / "gridadmin_1d.h5"

    grid_1d = copy(grid_all)
    grid_1d.quadtree_stats = None
    with GridAdminOut(path) as out:
        out.write(grid_1d)

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
        ("dmax", (4,), "float64"),  # added
        ("s1d", (4,), "float64"),  # added
        ("embedded_in", (4,), "int32"),  # added
        ("boundary_type", (4,), "int32"),  # added
        ("has_dem_averaged", (4,), "int32"),  # added
    ],
)
def test_write_nodes(h5_out, dataset, shape, dtype):
    assert h5_out["nodes"][dataset].shape == shape
    assert h5_out["nodes"][dataset].dtype == np.dtype(dtype)


@pytest.mark.parametrize(
    "dataset,shape,dtype",
    [
        ("bottom_level", (3,), "float64"),
        ("calculation_type", (3,), "int32"),
        ("cell_coords", (4, 3), "float64"),
        ("content_pk", (3,), "int32"),
        ("coordinates", (2, 3), "float64"),
        ("display_name", (3,), "|S64"),
        ("drain_level", (3,), "float64"),
        ("id", (3,), "int32"),
        ("initial_waterlevel", (3,), "float64"),
        ("is_manhole", (3,), "int32"),
        ("manhole_indicator", (3,), "int32"),
        ("node_type", (3,), "int32"),
        ("pixel_coords", (4, 3), "int32"),
        ("pixel_width", (3,), "int32"),
        ("shape", (3,), "|S4"),
        ("storage_area", (3,), "float64"),
        ("sumax", (3,), "float64"),
        ("surface_level", (3,), "float64"),
        ("width", (3,), "float64"),
        ("x_coordinate", (3,), "float64"),
        ("y_coordinate", (3,), "float64"),
        ("z_coordinate", (3,), "float64"),
        ("zoom_category", (3,), "int32"),
        ("code", (3,), "|S32"),  # added
        ("dmax", (3,), "float64"),  # added
        ("s1d", (3,), "float64"),  # added
        ("embedded_in", (3,), "int32"),  # added
        ("boundary_id", (3,), "int32"),  # added
        ("boundary_type", (3,), "int32"),  # added
        ("has_dem_averaged", (3,), "int32"),  # added
    ],
)
def test_write_nodes_embedded(h5_out, dataset, shape, dtype):
    assert h5_out["nodes_embedded"][dataset].shape == shape
    assert h5_out["nodes_embedded"][dataset].dtype == np.dtype(dtype)


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
        ("cross_section_height", (6,), "float64"),
        ("cross_section_shape", (6,), "int32"),
        ("cross_section_width", (6,), "float64"),
        ("discharge_coefficient", (6,), "float64"),
        ("discharge_coefficient_negative", (6,), "float64"),
        ("discharge_coefficient_positive", (6,), "float64"),
        ("display_name", (6,), "|S64"),
        ("dist_calc_points", (6,), "float64"),
        ("friction_type", (6,), "int32"),
        ("friction_value", (6,), "float64"),
        ("veg_coef1", (6,), "float64"),
        ("veg_coef2", (6,), "float64"),
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
        ("s1d", (6,), "float64"),  # added
        ("ds1d", (6,), "float64"),  # added
        ("dpumax", (6,), "float64"),  # added
        ("flod", (6,), "float64"),  # added
        ("flou", (6,), "float64"),  # added
        ("cross1", (6,), "int32"),  # added
        ("cross2", (6,), "int32"),  # added
        ("cross_weight", (6,), "float64"),  # added
        ("windshieldings", (8, 6), "float64"),
    ],
)
def test_write_lines(h5_out, dataset, shape, dtype):
    assert h5_out["lines"][dataset].shape == shape
    assert h5_out["lines"][dataset].dtype == np.dtype(dtype)


def test_line_geometries(h5_out):
    # line geometries are stored as a variable-length array [x, x, ..., y, y, ...]
    data = h5_out["lines"]["line_geometries"][:]
    assert np.isnan(data[0]).all()
    assert data[1].tolist() == [1, 2, 1, 2]
    assert data[2].tolist() == [1, 2, 3, 1, 2, 3]


def test_line_cross_mapping(h5_out):
    # cross ids should be mapped to (1 based) cross indexes
    assert_equal(h5_out["lines"]["cross1"][1:], [2, -9999, 1, 1, 2])
    assert_equal(h5_out["lines"]["cross2"][1:], [-9999, -9999, -9999, -9999, 3])


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
        ("lgrmin", (), "int32"),
        ("kmax", (), "int32"),
        ("mmax", (1,), "int32"),
        ("nmax", (1,), "int32"),
        ("dx", (1,), "float64"),
        ("dxp", (), "float64"),
        ("x0p", (), "float64"),
        ("y0p", (), "float64"),
    ],
)
def test_write_quadtree_pure_1d(h5_out_1d, dataset, shape, dtype):
    assert h5_out_1d["grid_coordinate_attributes"][dataset].shape == shape
    assert h5_out_1d["grid_coordinate_attributes"][dataset].dtype == np.dtype(dtype)


@pytest.mark.parametrize(
    "dataset,shape,dtype",
    [
        # ('ijmax', (), 'i4'),  deprecated
        # ('imax', (), 'i4'),  deprecated
        ("infl1d", (), "i4"),
        ("ingrw1d", (), "i4"),
        # ('jap1d', (), 'i4'),  deprecated
        # ('jmax', (), 'i4'),  deprecated
        ("l1dtot", (), "i4"),
        ("l2dtot", (), "i4"),
        # ('levnms', (), 'i4'),  deprecated
        # ('lgrmin', (), 'i4'),  deprecated
        ("lgrtot", (), "i4"),
        ("lgutot", (), "i4"),
        ("lgvtot", (), "i4"),
        # ('linall', (), 'i4'),  deprecated (can be derived from others)
        # ('lintot', (), 'i4'),  deprecated (can be derived from others)
        ("liutot", (), "i4"),
        ("livtot", (), "i4"),
        ("n1dobc", (), "i4"),
        ("n1dtot", (), "i4"),
        # ('n2dall', (), 'i4'),  deprecated (can be derived from others)
        ("n2dobc", (), "i4"),
        ("n2dtot", (), "i4"),
        ("ngr2bc", (), "i4"),
        ("ngrtot", (), "i4"),
        ("nob2dg", (), "i4"),
        ("nob2ds", (), "i4"),
        # ('nodall', (), 'i4'),  deprecated (can be derived from others)
        # ('nodobc', (), 'i4'),  deprecated (can be derived from others)
        # ('nodtot', (), 'i4'),  deprecated (can be derived from others)
    ],
)
def test_write_meta(h5_out, dataset, shape, dtype):
    assert h5_out["meta"][dataset].shape == shape
    assert h5_out["meta"][dataset].dtype == np.dtype(dtype)


@pytest.mark.parametrize(
    "attr,shape,dtype",
    [
        ("epsg_code", (), "S"),
        ("crs_wkt", (), "S"),  # added
        ("has_1d", (), "bool"),  # changed to bool
        ("has_2d", (), "bool"),  # changed to bool
        ("has_0d", (), "bool"),
        ("extent_1d", (4,), "float64"),
        ("extent_2d", (4,), "float64"),
        ("zero_dim_extent", (4,), "float64"),
        ("has_breaches", (), "bool"),  # changed to bool
        ("has_groundwater", (), "bool"),  # changed to bool
        ("has_groundwater_flow", (), "bool"),  # changed to bool
        ("has_interception", (), "bool"),  # changed to bool
        ("has_pumpstations", (), "bool"),  # changed to bool
        ("has_simple_infiltration", (), "bool"),  # changed to bool
        ("has_interflow", (), "bool"),  # added
        ("has_initial_waterlevels", (), "bool"),  # added
        ("model_name", (), "S"),
        ("model_slug", (), "S"),
        ("revision_hash", (), "S"),
        ("revision_nr", (), "int32"),  # changed to int
        ("threedi_version", (), "S"),
        ("threedicore_version", (), "S"),
    ],
)
def test_write_attrs(h5_out, attr, shape, dtype):
    assert h5_out.attrs[attr].shape == shape
    actual = h5_out.attrs[attr].dtype
    if dtype == "S":
        actual = actual.char
    assert np.dtype(actual) == np.dtype(dtype)


# obtained from bergermeer gridadmin.h5, edited:
# - 20 pumps to 4
# - int64 to int32
@pytest.mark.parametrize(
    "dataset,shape,dtype",
    [
        ("bottom_level", (4,), "float64"),
        ("capacity", (4,), "float64"),
        ("connection_node_end_pk", (4,), "int32"),
        ("connection_node_start_pk", (4,), "int32"),
        ("content_pk", (4,), "int32"),
        ("coordinates", (2, 4), "float64"),
        ("display_name", (4,), "|S64"),  # increased size from 24 to 64
        ("id", (4,), "int32"),
        ("lower_stop_level", (4,), "float64"),
        ("node1_id", (4,), "int32"),
        ("node2_id", (4,), "int32"),
        ("node_coordinates", (4, 4), "float64"),
        # ("nodp1d", (2, 3), "int32"), removed
        # ("p1dtyp", (4,), "int32"), removed
        ("start_level", (4,), "float64"),
        ("type", (4,), "int32"),
        ("zoom_category", (4,), "int32"),
        ("code", (4,), "|S32"),  # added
        ("upper_stop_level", (4,), "float64"),  # added
    ],
)
def test_write_pumps(h5_out, dataset, shape, dtype):
    assert h5_out["pumps"][dataset].shape == shape
    assert h5_out["pumps"][dataset].dtype == np.dtype(dtype)


@pytest.mark.parametrize(
    "dataset,shape,dtype",
    [
        ("id", (4,), "int32"),
        ("code", (4,), "|S32"),  # added
        ("shape", (4,), "int32"),
        ("content_pk", (4,), "int32"),
        ("width_1d", (4,), "float64"),
        ("offset", (4,), "int32"),
        ("count", (4,), "int32"),
        ("tables", (2, 6), "float64"),
    ],
)
def test_write_cross_sections(h5_out, dataset, shape, dtype):
    assert h5_out["cross_sections"][dataset].shape == shape
    assert h5_out["cross_sections"][dataset].dtype == np.dtype(dtype)


# obtained from bergermeer gridadmin.h5
# note that there is no dummy element in this group (just 2 levees)
@pytest.mark.parametrize(
    "dataset,shape,dtype",
    [
        ("id", (2,), "int32"),
        ("crest_level", (2,), "float64"),
        ("max_breach_depth", (2,), "float64"),
        ("coords", (2,), "object"),
    ],
)
def test_write_levees(h5_out, dataset, shape, dtype):
    assert h5_out["levees"][dataset].shape == shape
    assert h5_out["levees"][dataset].dtype == np.dtype(dtype)


def test_write_levees_coords(h5_out):
    # coords are stored as a variable-length array [x, x, ..., y, y, ...]
    data = h5_out["levees"]["coords"][:]
    assert data[0].tolist() == [1, 2, 4, 1, 2, 4]
    assert data[1].tolist() == [1, 2, 3, 1, 2, 3]


# obtained from bergermeer gridadmin.h5, edited:
# - 1 breach to 2
# - int64 to int32
@pytest.mark.parametrize(
    "dataset,shape,dtype",
    [
        ("id", (3,), "int32"),
        ("levl", (3,), "int32"),
        ("levbr", (3,), "float64"),
        ("levmat", (3,), "int32"),
        ("content_pk", (3,), "int32"),
        ("coordinates", (2, 3), "float64"),
        ("code", (3,), "|S32"),  # added
        ("display_name", (3,), "|S64"),  # added
        # ("llev", (1931,), "int32"),  dropped
        # ("kcu", (2,), "int32"),  dropped
        # ("seq_ids", (2,), "int32"),  dropped
    ],
)
def test_write_breaches(h5_out, dataset, shape, dtype):
    assert h5_out["breaches"][dataset].shape == shape
    assert h5_out["breaches"][dataset].dtype == np.dtype(dtype)


@pytest.mark.parametrize(
    "group,dataset,expected",
    [
        ("nodes", "id", 1),
        ("nodes", "dmax", 1.2),
        ("nodes_embedded", "embedded_in", 2),
        ("nodes_embedded", "dmax", 2.3),
        ("lines", "id", 1),
        ("lines", "dpumax", 1.2),
        ("pumps", "id", 1),
        ("pumps", "capacity", 0.1),
        ("cross_sections", "id", 1),
        ("cross_sections", "width_1d", 0.2),
        ("cross_sections", "offset", 0),  # reference to tables dataset, not increased
        ("lines", "line", [1, 2]),  # reference to node
        ("lines", "cross1", 2),  # reference to cross section
        ("lines", "cross2", -9999),  # reference to cross section
        ("pumps", "node1_id", 1),  # reference to node
        ("pumps", "node2_id", 2),  # reference to node
        ("breaches", "id", 1),
        ("breaches", "levl", 4),  # reference to line
    ],
)
def test_not_off_by_one(h5_out, group, dataset, expected):
    # gridadmin contains a dummy element at index 0 (so index 1 is the first)
    # references should also be increased by one
    assert h5_out[group][dataset][..., 1].tolist() == expected


@pytest.mark.parametrize(
    "dataset,shape,dtype",
    [
        ("area", (3,), "float64"),
        ("cci", (4,), "int32"),
        ("centroid_x", (3,), "float64"),
        ("centroid_y", (3,), "float64"),
        ("code", (3,), "|S100"),
        ("display_name", (3,), "|S250"),
        ("dry_weather_flow", (3,), "float64"),
        ("fac", (4,), "float64"),
        ("fb", (3,), "float64"),
        ("fe", (3,), "float64"),
        ("function", (3,), "|S64"),
        ("id", (3,), "int32"),
        ("imp", (4,), "int32"),
        ("infiltration_flag", (3,), "bool"),
        ("ka", (3,), "float64"),
        ("kh", (3,), "float64"),
        ("nr_of_inhabitants", (3,), "float64"),
        ("nxc", (4,), "float64"),
        ("nyc", (4,), "float64"),
        ("outflow_delay", (3,), "float64"),
        ("pk", (4,), "int32"),
        ("storage_limit", (3,), "float64"),
    ],
)
def test_write_surface(h5_out, dataset, shape, dtype):
    assert h5_out["surface"][dataset].shape == shape
    assert h5_out["surface"][dataset].dtype == np.dtype(dtype)


def test_encoding(h5_out):
    assert h5_out["nodes"]["display_name"][1].decode() == "foo"
    assert h5_out["nodes"]["display_name"][2].decode() == "金蟾"
    assert h5_out["nodes"]["display_name"][3].decode() == ""
