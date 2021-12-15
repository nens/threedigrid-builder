from threedigrid_builder.base import RasterInterface
from threedigrid_builder.exceptions import SchematisationError

import numpy as np


try:
    from osgeo import gdal

    gdal.UseExceptions()
except ImportError:
    gdal = None

__all__ = ["GDALInterface"]

GT_TOLERANCE = 7


def get_epsg_code(sr):
    """
    Return epsg code from a osr.SpatialReference object
    """
    if sr.IsGeographic():
        raise SchematisationError(
            f"The supplied DEM file has geographic projection '{sr.GetName()}'"
        )
    if sr.GetAuthorityName("PROJCS") != "EPSG":
        raise SchematisationError(
            f"The supplied DEM file has a none-EPSG projection '{sr.GetName()}'"
        )
    return int(sr.GetAuthorityCode("PROJCS"))


class GDALInterface(RasterInterface):
    def __init__(self, *args, **kwargs):
        if gdal is None:
            raise ImportError("Cannot use GDALInterface if GDAL is not available.")
        super().__init__(*args, **kwargs)
        if self.model_area_path is not None:
            raise NotImplementedError("Model areas are not available with GDAL")

    def __enter__(self):
        self._dataset = gdal.Open(self.path.as_posix(), gdal.GA_ReadOnly)
        c, a, b, f, d, e = self._dataset.GetGeoTransform()
        self.set_transform((a, b, c, d, e, f))
        self.set_epsg_code(get_epsg_code(self._dataset.GetSpatialRef()))
        return self

    def __exit__(self, *args, **kwargs):
        del self._dataset

    def read(self):
        width, height, bbox, data = self._create_area_arr_from_dem()

        return {
            "pixel_size": self.pixel_size,
            "width": width,
            "height": height,
            "bbox": bbox,
            "area_mask": np.flipud(data).T.astype(
                dtype=np.int16, copy=False, order="F"
            ),
        }

    def _create_area_arr_from_dem(self):
        xpixel, _, xmin, _, ypixel, ymax = self.transform
        width, height = self._dataset.RasterXSize, self._dataset.RasterYSize
        xmax = xmin + width * xpixel
        ymin = ymax + height * ypixel

        if xmax < xmin:
            xmin, xmax = xmax, xmin
        if ymax < ymin:
            ymin, ymax = ymax, ymin

        data = self._dataset.GetRasterBand(1).GetMaskBand().ReadAsArray()
        data[data > 0] = 1
        return width, height, (xmin, ymin, xmax, ymax), data
