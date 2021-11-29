from threedigrid_builder.base import RasterInterface
from threedigrid_builder.exceptions import SchematisationError

import json
import numpy as np


try:
    import rasterio
except ImportError:
    rasterio = None

__all__ = ["RasterioInterface"]


def get_epsg_code(crs):
    """
    Return epsg code from a osr.SpatialReference object
    """
    if crs.is_geographic:
        raise SchematisationError(
            f"The supplied DEM file has geographic projection '{str(crs)}'"
        )
    if not crs.is_epsg_code:
        raise SchematisationError(
            f"The supplied DEM file has a non-EPSG projection '{str(crs)}'"
        )
    return crs.to_epsg()


class RasterioInterface(RasterInterface):
    def __init__(self, *args, **kwargs):
        if rasterio is None:
            raise ImportError(
                "Cannot use RasterioInterface if rasterio is not available."
            )
        super().__init__(*args, **kwargs)

    def __enter__(self):
        self._env = rasterio.Env().__enter__()
        self._raster = rasterio.open(self.path, "r").__enter__()
        profile = self._raster.profile
        self.set_transform(profile["transform"][:6])
        self.set_epsg_code(get_epsg_code(profile["crs"]))
        return self

    def __exit__(self, *args, **kwargs):
        self._raster.__exit__(*args, **kwargs)
        self._env.__exit__(*args, **kwargs)

    def read(self):
        if self.model_area_path is not None:
            area_geometry = self._load_geometry(self.model_area_path)
            width, height, bbox, data = self._create_area_arr_from_geometry(
                area_geometry
            )
        else:
            width, height, bbox, data = self._create_area_arr_from_dem()

        return {
            "pixel_size": self.pixel_size,
            "width": width,
            "height": height,
            "bbox": bbox,
            "area_mask": np.flipud(data).T.astype(
                dtype=np.int32, copy=False, order="F"
            ),
        }

    def _create_area_arr_from_dem(self):
        data = self._raster.read_masks(1)
        data[data > 0] = 1
        return self._raster.width, self._raster.height, self._raster.bounds, data

    @staticmethod
    def _load_geometry(model_area_json):
        with model_area_json.open("r") as f:
            area_geometry = json.load(f)["features"][0]
        return area_geometry

    def _create_area_arr_from_geometry(self, area_geometry):
        import rasterio.features as rasterio_features

        dem = self._raster

        window = rasterio_features.geometry_window(
            dem, (area_geometry,), pixel_precision=3
        )
        left = dem.bounds.left + window.col_off * self.pixel_size
        top = dem.bounds.top + window.row_off * -self.pixel_size
        bbox = (
            left,
            top + window.height * -self.pixel_size,
            left + window.width * self.pixel_size,
            top,
        )
        self.transform = rasterio.transform.from_origin(
            left, top, self.pixel_size, self.pixel_size
        )

        data = rasterio.features.rasterize(
            (area_geometry["geometry"], 1),
            out_shape=(window.height, window.width),
            transform=self.transform,
            dtype=np.int32,
        )
        return window.width, window.height, bbox, data
