import pytest

from threedigrid_builder.interface import DictOut


@pytest.fixture(scope="session")
def dict_out(grid_all):
    with DictOut(path=None) as out:
        result = out.write(grid_all, geometry_format="wkt")

    return result


def test_nodes(dict_out):
    assert dict_out["nodes"]["id"].shape == (3,)
    assert dict_out["nodes"]["geometry"][0] == "POINT (1 1)"
    assert dict_out["nodes"]["node_type"][0] == "NODE_2D_OPEN_WATER"


def test_lines(dict_out):
    assert dict_out["lines"]["id"].shape == (5,)
    assert dict_out["lines"]["geometry"][0] == "LINESTRING (1 1, 2 2)"


def test_cells(dict_out):
    assert dict_out["cells"]["id"].shape == (1,)
    assert dict_out["cells"]["geometry"][0] == "POLYGON ((1 0, 1 1, 0 1, 0 0, 1 0))"


def test_breaches(dict_out):
    assert dict_out["breaches"]["id"].shape == (2,)
    assert dict_out["breaches"]["geometry"][0] == "POINT (0 2)"


def test_nodes_embedded(dict_out):
    assert dict_out["nodes_embedded"]["id"].shape == (2,)


def test_meta(dict_out):
    assert dict_out["meta"]["epsg_code"] == 28992


@pytest.mark.parametrize(
    "group,dataset,expected",
    [
        ("nodes", "id", 1),
        ("nodes", "dmax", 1.2),
        ("nodes_embedded", "embedded_in", 2),
        ("lines", "id", 1),
        ("lines", "dpumax", 1.2),
        ("lines", "node_1", 1),  # reference to node
        ("lines", "node_2", 2),  # reference to node
        ("lines", "cross_id1", 4),  # reference to cross section def id, not increased
        ("breaches", "id", 1),
        ("breaches", "line_id", 4),  # reference to line
    ],
)
def test_not_off_by_one(dict_out, group, dataset, expected):
    # dict does not contain a dummy element, but ids are increased by 1 for consistency
    # with the gridadmin
    assert dict_out[group][dataset][0] == expected
