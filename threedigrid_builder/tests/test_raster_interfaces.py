import numpy as np
import pytest
from osgeo import gdal

from threedigrid_builder.interface import GDALInterface


def test_does_not_exist(tmp_path):
    with pytest.raises(Exception):
        with GDALInterface(tmp_path / "i-do-not-exist.tif"):
            pass


@pytest.mark.skipif(
    int(gdal.VersionInfo()[:3]) >= 304,
    reason="this test is known to fail with GDAL versions below 3.4",
)
def test_crs_legacy(dem_path, crs_wkt_28992_legacy):
    with GDALInterface(dem_path) as dem:
        assert dem.crs == crs_wkt_28992_legacy


@pytest.mark.skipif(
    int(gdal.VersionInfo()[:3]) < 304,
    reason="this test is known to fail with GDAL 3.4 and newer",
)
def test_crs(dem_path, crs_wkt_28992):
    with GDALInterface(dem_path) as dem:
        assert dem.crs == crs_wkt_28992


def test_pixel_size(dem_path):
    with GDALInterface(dem_path) as dem:
        assert dem.pixel_size == 0.5


def test_read(dem_path):
    with GDALInterface(dem_path) as dem:
        result = dem.read()
        assert result["area_mask"].shape == (9517, 9726)
        assert result["area_mask"].dtype == np.int8
        assert np.count_nonzero(result["area_mask"] == 1) == 47180799
        assert result["pixel_size"] == 0.5
        assert result["width"] == 9517
        assert result["height"] == 9726
        assert result["bbox"] == (106314.0, 514912.0, 111072.5, 519775.0)
