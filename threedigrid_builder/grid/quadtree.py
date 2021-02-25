from threedigrid_builder.base import Lines
from threedigrid_builder.base import Nodes
from threedigrid_builder.constants import NodeType
from threedigrid_builder.grid.fwrapper import create_quadtree
from threedigrid_builder.grid.fwrapper import set_refinement
from threedigrid_builder.grid.fwrapper import set_2d_computational_nodes
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

        self._apply_refinements(refinements)
        self.active_cells = create_quadtree(
            self.kmax,
            self.mmax,
            self.nmax,
            self.min_cell_pixels,
            subgrid_meta.area_pix,
            self.lg
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

    def _apply_refinements(self, refinements):
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

    def get_nodes(self, subgrid_meta):
        """Compute 2D openwater Nodes based computed Quadtree

        Args:
          SubgridMeta object for passing active model_area to node computation.

        Return:
          Nodes object
        """

        # Create all arrays for filling in external Fortran routine.
        id = np.arange(self.active_cells)
        nodk = np.empty(len(id), dtype=np.int32, order='F')
        nodm = np.empty(len(id), dtype=np.int32, order='F')
        nodn = np.empty(len(id), dtype=np.int32, order='F')
        quad_nod = np.empty(self.lg.shape, dtype=np.int32, order='F')
        bounds = np.empty((len(id), 4), dtype=np.float64, order='F')
        coords = np.empty((len(id), 2), dtype=np.float64, order='F')
        
        # Node type is always openwater at first init
        node_type = np.full(
            (len(id)), NodeType.NODE_2D_OPEN_WATER, dtype='O', order='F'
        )

        set_2d_computational_nodes(
            np.array([self.bbox[0], self.bbox[1]]),
            self.kmax,
            self.mmax,
            self.nmax,
            self.dx,
            self.lg,
            nodk,
            nodm,
            nodn,
            quad_nod,
            bounds,
            coords,
        )

        return Nodes(
            id=id,
            node_type=node_type,
            nodk=nodk,
            nodm=nodm,
            nodn=nodn,
            bounds=bounds,
            coordinates=coords
        )
