import pathlib
import shutil
from unittest.mock import Mock

import h5py
import numpy as np
import pytest
from numpy.testing import assert_almost_equal, assert_array_equal

import threedigrid_builder
from threedigrid_builder.application import make_grid
from threedigrid_builder.constants import CrossSectionShape, LineType, NodeType

project_root = pathlib.Path(__file__).resolve().parent.parent
unittests_data_path = project_root / "threedigrid_builder/tests/data"


def count_unique(arr):
    """Return a dict of counts per unique value in the input array"""
    return dict(zip(*np.unique(arr, return_counts=True)))


@pytest.mark.parametrize(
    "filename",
    ["v2_bergermeer.sqlite"],
)
def test_integration(tmp_path, filename):
    shutil.copyfile(unittests_data_path / filename, tmp_path / filename)
    progress_callback = Mock()
    make_grid(
        tmp_path / filename,
        unittests_data_path / "dem_test_5m.tif",
        tmp_path / "gridadmin.h5",
        meta={
            "model_slug": "slug-123abc",
            "revision_hash": "123abc",
            "revision_nr": 24,
            "threedi_version": "1.2.3.dev",
        },
        progress_callback=progress_callback,
        upgrade=True,
    )
    with h5py.File(tmp_path / "gridadmin.h5", "r") as f:
        ## REPROJECTION
        # spot-check a coordinates as reprojection WGS84 -> RD might mismatch
        # WGS84: 4.728282895,52.645792838
        idx = np.where(f["nodes"]["content_pk"][:] == 1)[0][0]
        assert_almost_equal(
            f["nodes"]["coordinates"][:, idx],
            [110404.2, 517792.3],
            decimal=3,
        )

        ## NODES
        assert f["nodes"]["id"].shape == (13280,)  # Inpy: (13280, )
        assert_array_equal(f["nodes"]["id"][:], np.arange(f["nodes"]["id"].shape[0]))
        assert count_unique(f["nodes"]["node_type"]) == {
            -9999: 1,
            NodeType.NODE_2D_OPEN_WATER: 5374,
            NodeType.NODE_2D_GROUNDWATER: 5374,
            NodeType.NODE_1D_NO_STORAGE: 2485,  # Inpy: 2527
            NodeType.NODE_1D_STORAGE: 42,  # Inpy: 0
            NodeType.NODE_1D_BOUNDARIES: 4,
        }
        assert np.count_nonzero(f["nodes"]["is_manhole"][:] == 1) == 42
        assert np.count_nonzero(f["nodes"]["content_pk"][:] > 0) == 1360
        assert (
            np.count_nonzero(np.isfinite(f["nodes"]["initial_waterlevel"][:])) == 2531
        )

        ## LINES
        assert f["lines"]["id"].shape == (31916,)  # Inpy: (31916, )
        assert_array_equal(f["lines"]["id"][:], np.arange(f["lines"]["id"].shape[0]))
        assert count_unique(f["lines"]["kcu"]) == {
            -9999: 1,  # Inpy: 0: 1
            LineType.LINE_2D_GROUNDWATER: 11037,
            LineType.LINE_1D_ISOLATED: 716,
            LineType.LINE_1D_CONNECTED: 1545,
            LineType.LINE_1D_SHORT_CRESTED: 56,
            LineType.LINE_1D_DOUBLE_CONNECTED: 219,
            LineType.LINE_1D2D_SINGLE_CONNECTED_CLOSED: 23,
            LineType.LINE_1D2D_SINGLE_CONNECTED_OPEN_WATER: 1512,  # Inpy: 1907
            LineType.LINE_1D2D_DOUBLE_CONNECTED_OPEN_WATER: 396,  # Inpy: 0
            # Inpy: LineType.LINE_1D2D_POSSIBLE_BREACH: 1,
            LineType.LINE_2D: 9544,  # Inpy: 9553
            LineType.LINE_2D_OBSTACLE: 1493,  # Inpy: 1484
            LineType.LINE_2D_VERTICAL: 5374,
        }
        assert count_unique(f["lines"]["content_type"]) == {
            b"": 29378,  # Inpy: 29380
            b"v2_channel": 2346,
            b"v2_culvert": 92,
            b"v2_pipe": 42,
            b"v2_weir": 56,
            b"v2_breach": 2,
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
            CrossSectionShape.TABULATED_TRAPEZIUM: 1,
            CrossSectionShape.TABULATED_YZ: 1,
        }

        ## EMBEDDED NODES
        assert_array_equal(f["nodes_embedded"]["id"][:], [0])

        ## LEVEES
        assert f["levees"]["id"][:].tolist() == [1, 2, 3]

        ## BREACHES
        assert f["breaches"]["id"][:].tolist() == [0, 1]
        assert_almost_equal(
            f["breaches"]["coordinates"][:, 1], [108988.84, 517280.88], decimal=2
        )

        ## COUNTS
        assert {ds: f["meta"][ds][()] for ds in f["meta"]} == {
            "infl1d": 1931,
            "ingrw1d": 0,
            "l1dtot": 2532,
            "l2dtot": 11037,
            "lgrtot": 5374,
            "lgutot": 5476,
            "lgvtot": 5561,
            "liutot": 5476,
            "livtot": 5561,
            "n1dobc": 4,
            "n1dtot": 2527,
            "n2dobc": 0,
            "n2dtot": 5374,
            "ngr2bc": 0,
            "ngrtot": 5374,
            "nob2dg": 0,
            "nob2ds": 0,
        }

        ## ATTRIBUTES
        attrs = {
            k: (v.decode() if isinstance(v, bytes) else v) for (k, v) in f.attrs.items()
        }
        assert_almost_equal(
            attrs.pop("extent_1d"), [105427.6, 511632.4, 115887.0, 523549.6], decimal=1
        )
        assert_almost_equal(
            attrs.pop("extent_2d"), [106314, 514912, 111114, 519872], decimal=1
        )

        assert np.all(np.isnan(attrs.pop("zero_dim_extent")))
        crs_wkt = attrs.pop("crs_wkt")
        assert "28992" in crs_wkt

        assert attrs == {
            "epsg_code": "28992",
            "has_1d": True,
            "has_2d": True,
            "has_breaches": True,
            "has_embedded": False,
            "has_groundwater": True,
            "has_groundwater_flow": True,
            "has_initial_waterlevels": True,
            "has_interception": True,
            "has_interflow": False,
            "has_max_infiltration_capacity": False,
            "has_pumpstations": True,
            "has_simple_infiltration": False,
            "has_vegetation": False,
            "model_name": "simple_infil_no_grndwtr",
            "model_slug": "slug-123abc",
            "revision_hash": "123abc",
            "revision_nr": 24,
            "threedi_tables_version": "",
            "threedi_version": "1.2.3.dev",
            "threedicore_version": "",
            "threedigrid_builder_version": threedigrid_builder.__version__,
            "has_0d": False,
        }

    # progress increases
    args = [x[0] for x in progress_callback.call_args_list]
    assert all([b[0] > a[0] for (a, b) in zip(args[:-1], args[1:])])
