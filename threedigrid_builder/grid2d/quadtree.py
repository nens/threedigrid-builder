from ..lib.quadtree import create_quadtree
from ..lib.quadtree import set_refinement

import numpy as np
import pygeos


__all__ = ["QuadTree"]


class QuadTree:
    """Defines active cell levels for computational grid.
    """
    def __init__(
        self, subgrid_meta, num_refine_levels, min_gridsize, refinements
    ):

        self._lgrmin = min_gridsize / subgrid_meta.pixel_size
        self._kmax = num_refine_levels
        self._mmax = np.empty((self.kmax,), dtype=np.int32, order='F')
        self._nmax = np.empty((self.kmax,), dtype=np.int32, order='F')
        self._dx = np.empty((self.kmax,), dtype=np.float64, order='F')

        max_grid_x_pix = self._determine_max_quadtree_pixels(
            subgrid_meta.width
        )
        max_grid_y_pix = self._determine_max_quadtree_pixels(
            subgrid_meta.height
        )

        pix_grid_levels = self._lgrmin * 2 ** np.arange(0, self.kmax)
        self._mmax[:] = max_grid_x_pix / pix_grid_levels
        self._nmax[:] = max_grid_y_pix / pix_grid_levels
        self._dx[:] = pix_grid_levels * subgrid_meta.pixel_size
        self.bbox = np.array(subgrid_meta.bbox)

        # Array with dimensions of smallest active grid level and contains
        # map of active grid level for each quadtree cell.
        self.lg = np.full(
            (self.mmax[0], self.nmax[0]),
            self.kmax,
            dtype=np.int32,
            order='F'
        )

        self.apply_refinements(refinements)
        self.active_cells = create_quadtree(
            subgrid_meta.area_pix, self.kmax, self.mmax, self.nmax, self.lg
        )

    @property
    def min_cell_pixels(self):
        """Returns minimum number of pixels in smallles computational_cell.
        """
        return self._lgrmin

    @property
    def kmax(self):
        """Returns maximum number of active grid levels in quadtree.
        """
        return self._kmax

    @property
    def mmax(self):
        """Returns array with column dimensions of every
        active grid level [0:kmax].
        """
        return self._mmax

    @property
    def nmax(self):
        """Returns array with row dimensions of every
        active grid level [0:kmax].
        """
        return self._nmax

    @property
    def dx(self):
        """Returns array with cell widths at every
        active grid level [0:kmax].
        """
        return self._dx

    def _determine_max_quadtree_pixels(self, side):
        """Compute pixels dimension of active quadtree. (extend raster
        dimension to quadtree cell bound.)

        Args:
          number of pixels of raster dimension.

        Returns:
          Value of pixel dimension of quadtree dimension.
        """
        max_pix_largest_cell = self.min_cell_pixels * 2 ** (self.kmax - 1)
        return (max_pix_largest_cell * (side / (max_pix_largest_cell))) \
            + max_pix_largest_cell

    def apply_refinements(self, refinements):
        """Set active grid levels for based on refinement dict and 
        filling lg variable for refinement locations.

        Args:
          dict of refinemnets containing at least
          - geometry
          - refinement_level
          - id

        """
        for i in range(len(refinements['id'])):
            geom = np.asfortranarray(
                pygeos.get_coordinates(refinements['the_geom'][i])
            )
            set_refinement(
                id=refinements['id'][i],
                geom=geom,
                level=refinements['refinement_level'][i],
                type=pygeos.get_type_id(refinements['the_geom'][i]),
                bbox=self.bbox,
                mmax=self.mmax,
                nmax=self.nmax,
                dx=self.dx,
                lg=self.lg
            )