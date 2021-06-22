from threedigrid_builder.application import make_grid
from unittest.mock import Mock

import h5py
import pathlib


data_path = pathlib.Path(__file__).resolve().parent / "data"


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
        assert f["nodes"]["id"].shape == (7906,)
        assert f["lines"]["id"].shape == (13574,)
        assert f["pumps"]["id"].shape == (20,)

    # progress increases
    args = [x[0] for x in progress_callback.call_args_list]
    assert all([b[0] > a[0] for (a, b) in zip(args[:-1], args[1:])])
