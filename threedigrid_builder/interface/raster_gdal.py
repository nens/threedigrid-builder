import numpy as np

from threedigrid_builder.base import RasterInterface

try:
    from osgeo import gdal

    gdal.UseExceptions()
except ImportError:
    gdal = None

__all__ = ["GDALInterface"]

GT_TOLERANCE = 7


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
        self.set_crs(self._dataset.GetProjection())
        return self

    def __exit__(self, *args, **kwargs):
        del self._dataset

    def read(self):
        width, height, bbox, mask = self._create_area_arr_from_dem()

        return {
            "pixel_size": self.pixel_size,
            "width": width,
            "height": height,
            "bbox": bbox,
            "area_mask": np.flipud(mask).T.astype(dtype=np.int8, copy=True, order="F"),
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

        band = self._dataset.GetRasterBand(1)
        size_j, size_i = band.GetBlockSize()
        nodata = band.GetNoDataValue()
        mask = np.zeros((height, width), dtype=np.int8)
        n_blocks_j = ((width - 1) // size_j) + 1
        n_blocks_i = ((height - 1) // size_i) + 1

        for i in range(n_blocks_i):
            for j in range(n_blocks_j):
                i1 = i * size_i
                j1 = j * size_j
                i2 = min(i1 + size_i, height)
                j2 = min(j1 + size_j, width)
                data = band.ReadAsArray(
                    xoff=j1, yoff=i1, win_xsize=j2 - j1, win_ysize=i2 - i1
                )
                _mask = np.isfinite(data)
                if nodata is not None and np.isfinite(nodata):
                    _mask &= data != nodata
                mask[i1:i2, j1:j2] = _mask
        return width, height, (xmin, ymin, xmax, ymax), mask
