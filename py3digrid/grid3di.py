"""

"""
from gridwrapper.m_quadtree import QuadTreeFortran
import numpy as np
import logging
import rasterio
from rasterio import features as rasterio_features
import json

from enum import Enum
from enum import unique

logger = logging.getLogger(__name__)


@unique
class RefineType(Enum):
    LINESTRING = 1
    POLYGON = 2


class QuadTree:

    def __init__(self, subgrid, max_ref_lvl, min_grid_size):

        min_pix_cell = min_grid_size / subgrid.pixel_size

        self._quadtree = QuadTreeFortran(
            subgrid.origin[0],
            subgrid.origin[1],
            subgrid.pixel_size,
            subgrid.width,
            subgrid.height,
            min_pix_cell,
            max_ref_lvl
        )

    @property
    def origin(self):
        return self._quadtree.get_origin()

    @property
    def active_node_count(self):
        return self._quadtree.active_node_count()

    @property
    def quadtree(self):
        return self._quadtree

    def set_refinement(self, refinement):
        return self.quadtree.set_refinement(
            refinement['id'],
            refinement['geometry'].T,
            refinement['refinement_level'],
            RefineType[refinement['geometry_type'].upper()].value
        )

    def create_active_grid(self, model_area_pix):
        self.quadtree.make_quadtree()
        self.quadtree.set_active_2d_comp_cells(model_area_pix)


class SubgridArea:

    def __init__(self, dem_filepath, model_area=None):

        self.dem_filepath = dem_filepath
        with rasterio.Env():
            with rasterio.open(dem_filepath, 'r') as dem:
                self.transform = dem.profile['transform']
                self._dxp = self.transform[0]
                if model_area is not None:
                    area_geometry = self._load_geometry(model_area)
                    self._area_raster = self._create_area_arr_from_geometry(
                        dem, area_geometry
                    )
                else:
                    self._area_raster = self._create_area_arr_from_dem(dem)

    @property
    def area_pix(self):
        return np.flipud(self._area_raster).T

    @property
    def origin(self):
        return (self._left, self._bottom)

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
        with open(model_area_json, 'r') as f:
            area_geometry = json.load(f)['features'][0]
        return area_geometry

    def _create_area_arr_from_geometry(self, dem, area_geometry):

        window = rasterio_features.geometry_window(
                dem, (area_geometry, ), pixel_precision=3
            )
        self._height = window.height
        self._width = window.width
        self._left = dem.bounds.left + window.col_off * self.pixel_size
        top = dem.bounds.top + window.row_off * -self.pixel_size
        self._bottom = top + window.height * -self.pixel_size
        self._i_off = window.col_off
        self._j_off = dem.height - window.row_off - window.height
        self.transform = rasterio.transform.from_origin(
            self._left, top, self.pixel_size, self.pixel_size)

        data = rasterio.features.rasterize(
            (area_geometry['geometry'], 1),
            out_shape=(window.height, window.width),
            transform=self.transform,
            dtype=np.int32
        )
        return data

    def _create_area_arr_from_dem(self, dem):
        self._left = dem.bounds.left
        self._bottom = dem.bounds.bottom
        self._height = dem.height
        self._width = dem.width
        self._i_off = self._j_off = 0
        data = dem.read_masks(1)
        data[data > 0] = 1

        return data

    def create_model_area_tiff(self, out_filepath):

        with rasterio.Env():
            with rasterio.open(self.dem_filepath, 'r') as dem:
                profile = dem.profile
                profile.update(
                    width=self.width,
                    height=self.height,
                    out_shape=(self.height, self.width),
                    transform=self.transform,
                    dtype=np.int32,
                    compress='lzw'
                )
                with rasterio.open(out_filepath, 'w', **profile) as area_tiff:
                    area_tiff.write(self._area_raster.astype(np.int32), 1)
