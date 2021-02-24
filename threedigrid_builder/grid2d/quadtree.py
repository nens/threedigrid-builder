import numpy as np
import pygeos
from ..lib.quadtree import set_refinement, create_quadtree


__all__ = ["QuadTree"]


class QuadTree:

    def __init__(
        self,
        origin,
        num_refine_levels,
        min_gridsize,
        subgrid_width,
        subgrid_height,
        pixel_size,
        refinements
    ):

        self.lgrmin = min_gridsize / pixel_size
        self.kmax = num_refine_levels
        self.mmax = np.empty((self.kmax,), dtype=np.int32, order='F')
        self.nmax = np.empty((self.kmax,), dtype=np.int32, order='F')
        self.dx = np.empty((self.kmax,), dtype=np.float64, order='F')

        max_grid_x_pix = self._determine_max_quadtree_pixels(
            subgrid_width
        )
        max_grid_y_pix = self._determine_max_quadtree_pixels(
            subgrid_height
        )
        self.mmax[:] = max_grid_x_pix / (self.lgrmin * 2 ** np.arange(0, self.kmax))
        self.nmax[:] = max_grid_y_pix / (self.lgrmin * 2 ** np.arange(0, self.kmax))
        self.dx[:] = self.lgrmin * 2 ** np.arange(0, self.kmax) * pixel_size
        ur_corner = (
            origin[0] + max_grid_x_pix * pixel_size,
            origin[1] + max_grid_y_pix * pixel_size
        )
        self.bbox = np.array(tuple(self.origin) + ur_corner)
        self.lg = np.full(
            (self.mmax[0], self.nmax[0]),
            self.kmax,
            dtype=np.int32,
            order='F'
        )
    
        self.set_refinements(refinements)
        self.create_quadtree()

    def _determine_max_quadtree_pixels(self, side):
        max_pix_largest_cell = self.lgrmin * 2 ** (self.kmax - 1)
        return (max_pix_largest_cell * (side / (max_pix_largest_cell))) \
            + max_pix_largest_cell

    def set_refinements(self, refinements):
        for i in range(len(refinements['id'])):
            geom = np.asfortranarray(
                pygeos.get_coordinates(refinements['the_geom'][i])
            )
            set_refinement(
                id=refinements['id'][i],
                geom=geom,
                level=refinements['refinement_level'][i],
                type=pygeos.get_type_id(refinements['the_geom'][i]),
                origin=np.array(self.origin),
                bbox=self.bbox,
                mmax=self.mmax,
                nmax=self.nmax,
                dx=self.dx,
                lg=self.lg
            )

    def create_quadtree(self):
        return create_quadtree(self.kmax, self.mmax, self.nmax, self.lg)
