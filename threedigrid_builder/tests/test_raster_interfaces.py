from threedigrid_builder.interface import GDALInterface
from threedigrid_builder.interface import RasterioInterface

import numpy as np
import pytest


interfaces = [pytest.param(RasterioInterface)]


@pytest.fixture(params=[RasterioInterface, GDALInterface])
def interface(request):
    try:
        request.param(None, None)
    except ImportError:
        pytest.skip(f"Interface {request.param} not available")
    else:
        yield request.param


def test_does_not_exist(interface, tmp_path):
    with pytest.raises(Exception):
        with interface(tmp_path / "i-do-not-exist.tif"):
            pass


def test_epsg_code(interface, dem_path):
    with interface(dem_path) as dem:
        assert dem.epsg_code == 28992


def test_pixel_size(interface, dem_path):
    with interface(dem_path) as dem:
        assert dem.pixel_size == 0.5


def test_read(interface, dem_path):
    with interface(dem_path) as dem:
        result = dem.read()
        assert result["area_mask"].shape == (9517, 9726)
        assert result["area_mask"].dtype == np.int16
        assert np.count_nonzero(result["area_mask"] == 1) == 47180799
        assert result["pixel_size"] == 0.5
        assert result["width"] == 9517
        assert result["height"] == 9726
        assert result["bbox"] == (106314.0, 514912.0, 111072.5, 519775.0)
