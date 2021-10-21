from threedigrid_builder.base import GridSettings
from threedigrid_builder.base import Lines
from threedigrid_builder.base import Nodes
from threedigrid_builder.base import Pumps
from threedigrid_builder.base import TablesSettings
from threedigrid_builder.grid import CrossSections
from threedigrid_builder.grid import GridMeta
from threedigrid_builder.grid import QuadtreeStats
from threedigrid_builder.interface import GridAdminOut

import h5py
import numpy as np
import pygeos
import pytest


@pytest.fixture(scope="session")
def h5_out(tmpdir_factory):
    nodes = Nodes(
        id=[0, 1, 2],
        dmax=[1.2, 2.2, 3.3],
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
        epsg_code=12432634,
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

    path = tmpdir_factory.mktemp("h5") / "gridadmin.h5"
    with GridAdminOut(path) as out:
        out.write_meta(meta)
        out.write_grid_counts(nodes, lines)
        out.write_quadtree(quadtree_stats)
        out.write_nodes(nodes)
        out.write_nodes_embedded(nodes_embedded)
        out.write_lines(lines)
        out.write_pumps(pumps)
        out.write_cross_sections(cross_sections)

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
        ("s1d", (6,), "float64"),  # added
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


def test_line_geometries(h5_out):
    # line geometries are stored as a variable-length array [x, x, ..., y, y, ...]
    data = h5_out["lines"]["line_geometries"][:]
    assert data[0].tolist() == [-9999, -9999]
    assert data[1].tolist() == [1, 2, 1, 2]
    assert data[2].tolist() == [1, 2, 3, 1, 2, 3]


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
def test_write_meta(h5_out, dataset, shape, dtype):
    assert h5_out["meta"][dataset].shape == shape
    assert h5_out["meta"][dataset].dtype == np.dtype(dtype)


@pytest.mark.parametrize(
    "attr,shape,dtype",
    [
        ("epsg_code", (), "int32"),  # changed to int
        ("has_1d", (), "bool"),  # changed to bool
        ("has_2d", (), "bool"),  # changed to bool
        ("extent_1d", (4,), "float64"),
        ("extent_2d", (4,), "float64"),
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


@pytest.mark.xfail
@pytest.mark.parametrize(
    "group,attr",
    [
        ("breaches", "prepared"),
        ("levees", "prepared"),
        ("lines", "channels_prepared"),
        ("lines", "culverts_prepared"),
        ("lines", "lines_prepared"),
        ("lines", "orifices_prepared"),
        ("lines", "pipes_prepared"),
        ("lines", "weirs_prepared"),
        ("nodes", "connectionnodes_prepared"),
        ("nodes", "manholes_prepared"),
        ("nodes", "prepared"),
        ("pumps", "prepared"),
    ],
)  # [(x, y) for x in list(f) for y in list(f[x].attrs)]
def test_write_sub_attrs(h5_out, group, attr):
    assert h5_out[group].attrs[attr] == 1


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
        ("lines", "cross1", 1),  # reference to cross section
        ("lines", "cross2", -9999),  # reference to cross section
        ("pumps", "node1_id", 1),  # reference to node
        ("pumps", "node2_id", 2),  # reference to node
    ],
)
def test_not_off_by_one(h5_out, group, dataset, expected):
    # gridadmin contains a dummy element at index 0 (so index 1 is the first)
    # references should also be increased by one
    assert h5_out[group][dataset][..., 1].tolist() == expected
