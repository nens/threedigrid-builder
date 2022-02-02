from pyproj.crs import CRS
from threedigrid_builder.interface import GDALInterface

import numpy as np
import pytest


def test_does_not_exist(tmp_path):
    with pytest.raises(Exception):
        with GDALInterface(tmp_path / "i-do-not-exist.tif"):
            pass


def test_epsg_code(dem_path):
    with GDALInterface(dem_path) as dem:
        assert dem.crs == CRS(
            'PROJCS["Amersfoort / RD New",GEOGCS["Amersfoort",DATUM["Amersfoort",SPHEROID["Bessel 1841",6377397.155,299.1528128,AUTHORITY["EPSG","7004"]],TOWGS84[565.2369,50.0087,465.658,-0.406857,0.350733,-1.87035,4.0812],AUTHORITY["EPSG","6289"]],PRIMEM["Greenwich",0,AUTHORITY["EPSG","8901"]],UNIT["degree",0.0174532925199433,AUTHORITY["EPSG","9122"]],AUTHORITY["EPSG","4289"]],PROJECTION["Oblique_Stereographic"],PARAMETER["latitude_of_origin",52.1561605555556],PARAMETER["central_meridian",5.38763888888889],PARAMETER["scale_factor",0.9999079],PARAMETER["false_easting",155000],PARAMETER["false_northing",463000],UNIT["metre",1,AUTHORITY["EPSG","9001"]],AXIS["Easting",EAST],AXIS["Northing",NORTH],AUTHORITY["EPSG","28992"]]'
        )


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
