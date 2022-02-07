from threedigrid_builder.interface import GDALInterface

import numpy as np
import pytest


def test_does_not_exist(tmp_path):
    with pytest.raises(Exception):
        with GDALInterface(tmp_path / "i-do-not-exist.tif"):
            pass


def test_crs(dem_path, crs_wkt_28992_legacy):
    with GDALInterface(dem_path) as dem:
        assert dem.crs == crs_wkt_28992_legacy


def test_pixel_size(dem_path):
    with GDALInterface(dem_path) as dem:
        assert dem.pixel_size == 0.5


def test_read(dem_path):
    with GDALInterface(dem_path) as dem:
        result = dem.read()
        assert result["area_mask"].shape == (9517, 9726)
        assert result["area_mask"].dtype == np.int16
        assert np.count_nonzero(result["area_mask"] == 1) == 47180799
        assert result["pixel_size"] == 0.5
        assert result["width"] == 9517
        assert result["height"] == 9726
        assert result["bbox"] == (106314.0, 514912.0, 111072.5, 519775.0)
