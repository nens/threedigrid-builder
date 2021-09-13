import json
import numpy as np


try:
    from rasterio import features as rasterio_features

    import rasterio
except ImportError:
    rasterio = rasterio_features = None


class Subgrid:
    def __init__(self, filepath, model_area=None):

        self.filepath = filepath
        with rasterio.Env():
            with rasterio.open(filepath, "r") as dem:
                self.transform = dem.profile["transform"]
                self._dxp = self.transform[0]
                if model_area is not None:
                    area_geometry = self._load_geometry(model_area)
                    self._area_raster = self._create_area_arr_from_geometry(
                        dem, area_geometry
                    )
                else:
                    self._area_raster = self._create_area_arr_from_dem(dem)

    def __repr__(self):
        return f"<SubgridMeta object of raster file {self.filepath}>"

    @property
    def mask(self):
        return np.flipud(self._area_raster).T.astype(
            dtype=np.int32, copy=False, order="F"
        )

    @property
    def bbox(self):
        return self._bbox

    @property
    def pixel_size(self):
        return self._dxp

    @property
    def width(self):
        return self._width

    @property
    def height(self):
        return self._height

    @property
    def offset_model_area_dem(self):
        return (self._i_off, self._j_off)

    def _load_geometry(self, model_area_json):
        with open(model_area_json, "r") as f:
            area_geometry = json.load(f)["features"][0]
        return area_geometry

    def _create_area_arr_from_geometry(self, dem, area_geometry):

        window = rasterio_features.geometry_window(
            dem, (area_geometry,), pixel_precision=3
        )
        self._height = window.height
        self._width = window.width
        left = dem.bounds.left + window.col_off * self.pixel_size
        top = dem.bounds.top + window.row_off * -self.pixel_size
        self._bbox = (
            left,
            top + window.height * -self.pixel_size,
            left + window.width * self.pixel_size,
            top,
        )
        self._i_off = window.col_off
        self._j_off = dem.height - window.row_off - window.height
        self.transform = rasterio.transform.from_origin(
            left, top, self.pixel_size, self.pixel_size
        )

        data = rasterio.features.rasterize(
            (area_geometry["geometry"], 1),
            out_shape=(window.height, window.width),
            transform=self.transform,
            dtype=np.int32,
        )
        return data

    def _create_area_arr_from_dem(self, dem):
        self._bbox = dem.bounds
        self._height = dem.height
        self._width = dem.width
        self._i_off = self._j_off = 0
        data = dem.read_masks(1)
        data[data > 0] = 1

        return data

    def get_meta(self):
        """Return meta information of subgrid from raster."""
        return {
            "pixel_size": self.pixel_size,
            "width": self.width,
            "height": self.height,
            "bbox": self.bbox,
            "area_mask": self.mask,
        }

    def create_model_area_tiff(self, out_filepath):

        with rasterio.Env():
            with rasterio.open(self.dem_filepath, "r") as dem:
                profile = dem.profile
                profile.update(
                    width=self.width,
                    height=self.height,
                    out_shape=(self.height, self.width),
                    transform=self.transform,
                    dtype=np.int32,
                    compress="lzw",
                )
                with rasterio.open(out_filepath, "w", **profile) as area_tiff:
                    area_tiff.write(self._area_raster.astype(np.int32), 1)
