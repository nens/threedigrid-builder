"""

"""
from gridwrapper.m_quadtree import QuadTreeFortran
import logging

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
