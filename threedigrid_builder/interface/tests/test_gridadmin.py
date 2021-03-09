import h5py
import pytest
import numpy as np
import tempfile
import pygeos

from threedigrid_builder.base import Nodes
from threedigrid_builder.base import Lines
from threedigrid_builder.interface import GridAdminOut


@pytest.fixture(scope="session")
def h5_out():
    nodes = Nodes(id=[1, 2, 3])
    lines = Lines(
        id=[1, 2, 3, 4, 5],
        line=[[1, 2], [2, 3], [3, 1], [1, 3], [3, 2]],
        line_geometries=[
            pygeos.linestrings([[1, 1], [2, 2]]),
            pygeos.linestrings([[1, 1], [2, 2], [3, 3]]),
            None,
            None,
            None
        ]
    )

    with tempfile.NamedTemporaryFile(suffix=".h5") as tmpfile:
        path = tmpfile.name
        with GridAdminOut(path) as out:
            out.write_nodes(nodes, pixel_size=0.5)
            out.write_lines(lines)

        with h5py.File(path, "r") as f:
            yield f


# obtained from bergermeer gridadmin.h5, edited:
# - 13280 nodes to 3 nodes
# - int64 to int32
@pytest.mark.parametrize("dataset,shape,dtype", [
    ('bottom_level', (3,), 'float64'),
    ('calculation_type', (3,), 'int32'),
    ('cell_coords', (4, 3), 'float64'),
    ('content_pk', (3,), 'int32'),
    ('coordinates', (2, 3), 'float64'),
    ('display_name', (3,), '|S64'),
    ('drain_level', (3,), 'float64'),
    ('id', (3,), 'int32'),
    ('initial_waterlevel', (3,), 'float64'),
    ('is_manhole', (3,), 'int32'),
    ('manhole_indicator', (3,), 'int32'),
    ('node_type', (3,), 'int32'),
    ('pixel_coords', (4, 3), 'int32'),
    ('pixel_width', (3,), 'int32'),
    ('seq_id', (3,), 'int32'),
    ('shape', (3,), '|S4'),
    ('storage_area', (3,), '|S32'),
    ('sumax', (3,), 'float64'),
    ('surface_level', (3,), 'float64'),
    ('width', (3,), 'float64'),
    ('x_coordinate', (3,), 'float64'),
    ('y_coordinate', (3,), 'float64'),
    ('z_coordinate', (3,), 'float64'),
    ('zoom_category', (3,), 'int32'),
    ('code', (3,), '|S32'),  # added
])
def test_write_nodes(h5_out, dataset, shape, dtype):
    assert h5_out["nodes"][dataset].shape == shape
    assert h5_out["nodes"][dataset].dtype == np.dtype(dtype)


# obtained from bergermeer gridadmin.h5, edited:
# - 5 lines to 5 lines
# - int64 to int32
@pytest.mark.parametrize("dataset,shape,dtype", [  
    ('calculation_type', (5,), 'int32'),
    ('code', (5,), '|S32'),
    ('connection_node_end_pk', (5,), 'int32'),
    ('connection_node_start_pk', (5,), 'int32'),
    ('content_pk', (5,), 'int32'),
    ('content_type', (5,), '|S10'),
    ('crest_level', (5,), 'float64'),
    ('crest_type', (5,), 'int32'),
    ('cross_section_height', (5,), 'int32'),
    ('cross_section_shape', (5,), 'int32'),
    ('cross_section_width', (5,), 'float64'),
    ('discharge_coefficient', (5,), 'float64'),
    ('discharge_coefficient_negative', (5,), 'float64'),
    ('discharge_coefficient_positive', (5,), 'float64'),
    ('display_name', (5,), '|S64'),
    ('dist_calc_points', (5,), 'float64'),
    ('friction_type', (5,), 'int32'),
    ('friction_value', (5,), 'float64'),
    ('id', (5,), 'int32'),
    ('invert_level_end_point', (5,), 'float64'),
    ('invert_level_start_point', (5,), 'float64'),
    ('kcu', (5,), 'int32'),
    ('lik', (5,), 'int32'),
    ('line', (2, 5), 'int32'),
    ('line_coords', (4, 5), 'float64'),
    ('line_geometries', (5,), 'object'),
    ('material', (5,), 'int32'),
    ('sewerage', (5,), 'int32'),
    ('sewerage_type', (5,), 'int32'),
    ('zoom_category', (5,), 'int32'),
    ('ds1d', (5,), 'float64'),  # added
    ('dpumax', (5,), 'float64'),  # added
    ('flod', (5,), 'float64'),  # added
    ('flou', (5,), 'float64'),  # added
    ('cross1', (5,), 'int32'),  # added
    ('cross2', (5,), 'int32'),  # added
    ('cross_weight', (5,), 'float64'),  # added
])
def test_write_lines(h5_out, dataset, shape, dtype):
    assert h5_out["lines"][dataset].shape == shape
    assert h5_out["lines"][dataset].dtype == np.dtype(dtype)
