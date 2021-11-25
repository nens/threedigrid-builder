from threedigrid_builder.base import RasterInterface

import json
import numpy as np


try:
    from rasterio import features as rasterio_features

    import rasterio
except ImportError:
    rasterio = rasterio_features = None

__all__ = ["RasterioInterface"]


class RasterioInterface(RasterInterface):
    def __init__(self, *args, **kwargs):
        if rasterio is None:
            raise ImportError("Module 'rasterio' is not available")
        super().__init__(*args, **kwargs)

    def __enter__(self):
        if rasterio is None:
            raise ImportError("Module 'rasterio' is not available")
        self._env = rasterio.Env().__enter__()
        self._raster = rasterio.open(self.path, "r").__enter__()
        self.pixel_size = self._raster.profile["transform"][0]
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
        # self._i_off = window.col_off
        # self._j_off = dem.height - window.row_off - window.height
        # offset_model_area_dem = (self._i_off, self._j_off)
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

    # def create_model_area_tiff(self, out_filepath):

    #     with rasterio.Env():
    #         with rasterio.open(self.dem_filepath, "r") as dem:
    #             profile = dem.profile
    #             profile.update(
    #                 width=self.width,
    #                 height=self.height,
    #                 out_shape=(self.height, self.width),
    #                 transform=self.transform,
    #                 dtype=np.int32,
    #                 compress="lzw",
    #             )
    #             with rasterio.open(out_filepath, "w", **profile) as area_tiff:
    #                 area_tiff.write(self._area_raster.astype(np.int32), 1)
