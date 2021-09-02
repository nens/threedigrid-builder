from numpy.testing import assert_array_equal
from threedigrid_builder.application import make_grid
from threedigrid_builder.constants import CrossSectionShape
from threedigrid_builder.constants import LineType
from threedigrid_builder.constants import NodeType
from unittest.mock import Mock

import h5py
import numpy as np
import pathlib


data_path = pathlib.Path(__file__).resolve().parent / "data"


def count_unique(arr):
    """Return a dict of counts per unique value in the input array"""
    return dict(zip(*np.unique(arr, return_counts=True)))


def test_integration(tmp_path):
    progress_callback = Mock()
    make_grid(
        data_path / "v2_bergermeer.sqlite",
        data_path / "dem_test_5m.tif",
        tmp_path / "gridadmin.h5",
        model_area_path=None,  # untested
        meta={
            "model_slug": "slug-123abc",
            "revision_hash": "123abc",
            "revision_nr": 24,
            "threedi_version": "1.2.3.dev",
        },
        progress_callback=progress_callback,
    )
    with h5py.File(tmp_path / "gridadmin.h5", "r") as f:
        ## NODES
        assert f["nodes"]["id"].shape == (7906,)
        assert_array_equal(f["nodes"]["id"][:], np.arange(f["nodes"]["id"].shape[0]))
        assert count_unique(f["nodes"]["node_type"]) == {
            -9999: 1,
            NodeType.NODE_2D_OPEN_WATER: 5374,
            NodeType.NODE_1D_NO_STORAGE: 2485,
            NodeType.NODE_1D_STORAGE: 42,
            NodeType.NODE_1D_BOUNDARIES: 4,
        }

        ## LINES
        assert f["lines"]["id"].shape == (15498,)
        assert_array_equal(f["lines"]["id"][:], np.arange(f["lines"]["id"].shape[0]))
        assert count_unique(f["lines"]["kcu"]) == {
            -9999: 1,
            LineType.LINE_1D_ISOLATED: 716,
            LineType.LINE_1D_CONNECTED: 1545,
            LineType.LINE_1D_SHORT_CRESTED: 56,
            LineType.LINE_1D_DOUBLE_CONNECTED: 219,
            LineType.LINE_1D2D_SINGLE_CONNECTED_CLOSED: 23,
            LineType.LINE_1D2D_SINGLE_CONNECTED_OPEN_WATER: 1512,
            LineType.LINE_1D2D_DOUBLE_CONNECTED_OPEN_WATER: 396,
            LineType.LINE_2D: 9537,
            LineType.LINE_2D_OBSTACLE: 1493,
        }
        assert count_unique(f["lines"]["content_type"]) == {
            b"": 12962,
            b"v2_channel": 2346,
            b"v2_culvert": 92,
            b"v2_pipe": 42,
            b"v2_weir": 56,
        }

        ## PUMPS
        assert f["pumps"]["id"].shape == (20,)
        assert_array_equal(f["pumps"]["id"][:], np.arange(f["pumps"]["id"].shape[0]))

        ## CROSS SECTIONS
        assert_array_equal(
            f["cross_sections"]["id"][:], np.arange(f["cross_sections"]["id"].shape[0])
        )
        assert count_unique(f["cross_sections"]["shape"]) == {
            -9999: 1,
            CrossSectionShape.CIRCLE: 8,
            CrossSectionShape.RECTANGLE: 3,
        }

        ## EMBEDDED NODES
        assert_array_equal(f["nodes_embedded"]["id"][:], [0])

    # progress increases
    args = [x[0] for x in progress_callback.call_args_list]
    assert all([b[0] > a[0] for (a, b) in zip(args[:-1], args[1:])])
