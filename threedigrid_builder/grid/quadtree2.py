import itertools
import logging
import math

import numpy as np

from threedigrid_builder.base import Lines, Nodes, Quarters
from threedigrid_builder.constants import LineType, NodeType
from threedigrid_builder.exceptions import SchematisationError
from threedigrid_builder.interface import (
    DictOut,
    GDALInterface,
    GeopackageOut,
    GridAdminOut,
    SQLite,
)

from .fgrid import m_cells, m_quadtree

logger = logging.getLogger(__name__)


__all__ = ["QuadTree2"]


def reduce_refinement_levels(refinements, num_refine_levels):
    """Determine the lowest refinement level that is actually in use.

    Args:
      refinements (GridRefinements): all gridrefinements polygon and linestrings
        with refinement levels. refinement_level may be changed inplace.
      num_refine_levels (int): the number of refinement levels supplied in the settings

    Returns:
      num_refine_levels (int), possibly lower than the input one
    """
    if refinements is not None and len(refinements) > 0:
        min_level = refinements.refinement_level.min().item()
    else:
        min_level = num_refine_levels

    min_level -= 1
    if refinements is not None:
        refinements.refinement_level -= min_level
    return num_refine_levels - min_level


class QuadTree2:
    """Defines active cell levels for computational grid."""

    def __init__(
        self,
        dem_path,
        num_refine_levels,
        min_gridsize,
        use_2d_flow,
        refinements,
        max_mem=128 * 1024**2,
    ):
        self.dem_path = dem_path
        self.num_refine_levels = num_refine_levels
        self.min_gridsize = min_gridsize
        self.use_2d_flow = use_2d_flow
        self.refinements = refinements
        self.max_mem = max_mem

        with GDALInterface(self.dem_path) as dem_interface:
            self.pixel_size = dem_interface.pixel_size
            self.width, self.height, bbox = dem_interface.get_width_height_bbox()
            self.origin = (bbox[0], bbox[1])
            self.dataset_chunk = dem_interface._dataset.GetRasterBand(1).GetBlockSize()

        min_num_pix = self.min_gridsize / self.pixel_size
        if min_num_pix % 2 == 0:
            self.lgrmin = int(min_num_pix)
        else:
            raise SchematisationError(
                f"Smallest 2D grid cell does not contain an even number of pixels. "
                f"Smallest 2D grid cell size: {min_gridsize}m. Pixel size: {self.pixel_size}m."
            )

        # Maximum number of active grid levels in quadtree.
        self.kmax = reduce_refinement_levels(refinements, num_refine_levels)
        self.lgrmin *= 2 ** (num_refine_levels - self.kmax)

        # Array with cell widths at every active grid level [0:kmax]
        self.mmax = np.empty((self.kmax,), dtype=np.int32, order="F")
        # Array with row dimensions of every active grid level [0:kmax].
        self.nmax = np.empty((self.kmax,), dtype=np.int32, order="F")
        # Array with cell widths at every active grid level [0:kmax].
        self.dx = np.empty((self.kmax,), dtype=np.float64, order="F")

        lvl_multiplr = 2 ** np.arange(self.kmax - 1, -1, -1, dtype=np.int32)

        # Determine number of largest cells that fit over subgrid extent.
        self.lgrmax = self.lgrmin * lvl_multiplr[0]
        max_large_cells_col = int(math.ceil(self.width / (self.lgrmax)))
        max_large_cells_row = int(math.ceil(self.height / (self.lgrmax)))

        # Calculate grid level dimensions.
        self.mmax[:] = max_large_cells_col * lvl_multiplr
        self.nmax[:] = max_large_cells_row * lvl_multiplr
        self.dx[:] = self.lgrmin * lvl_multiplr[::-1] * self.pixel_size

        # Array with dimensions of smallest active grid level and contains
        # map of active grid level for each quadtree cell.
        self.lg = self._apply_refinements(refinements)

        # Array with dimensions of smallest active grid level and contains
        # idx of active grid level for each quadtree cell.
        self.quad_idx = np.empty(
            (self.mmax[0], self.nmax[0]), dtype=np.int32, order="F"
        )
        self.n_cells = np.array(0, dtype=np.int32, order="F")
        self.n_lines_u = np.array(0, dtype=np.int32, order="F")
        self.n_lines_v = np.array(0, dtype=np.int32, order="F")

    def __repr__(self):
        return (
            f"<Quadtree object with {self.kmax} refinement levels "
            f"and {self.n_cells} active computational cells>"
        )

    def _apply_refinements(self, refinements):
        """Set active grid levels for based on refinements and
        setting lg variable for refinement locations.

        Args:
          refinements (GridRefinements): all gridrefinements polygon and linestrings
            with refinement levels.

        Returns:
            lg (ndarray): Array with spatial refinement locations.
        """

        if refinements is None:
            return np.full(
                (self.mmax[0], self.nmax[0]), self.kmax, dtype=np.int32, order="F"
            )

        lg = refinements.rasterize(
            origin=self.origin,
            height=self.nmax[0],
            width=self.mmax[0],
            cell_size=self.dx[0],
            no_data_value=self.kmax,
        )
        return np.asfortranarray(lg.T)

    def make_quadtree(self):
        m_quadtree.balance_quadtree(self.kmax, self.mmax, self.nmax, self.lg)
        tile = self.determine_tile(self.max_mem, byte_count=1)
        for (x0, y0, x1, y1), mask in self.iter_by_tile(*tile):
            # Check DEM orientation, seems like top row is now flipped to last column in fortran
            import ipdb

            ipdb.set_trace()
            print(x0, y0, x1, y1)

        return []

    def determine_tile(
        self, max_mem: int = 128 * 1024**2, byte_count: int = 1
    ) -> tuple[int, int]:
        """Determine tile size for DEM reading based on maximum memory and precision.

        Args:
          max_mem (int): maximum memory to use for each tile in bytes.
          byte_count (int): number of bytes per data point (1 for int8, 4 for float32).
        """
        # compute how many pixels with size of byte_count fit in max_mem
        max_cells = int(max_mem / byte_count)

        # compute the shape of a native chunk (round up to whole grid cells)
        dataset_chunk = np.asarray(self.dataset_chunk)
        chunk = np.ceil(dataset_chunk / self.lgrmax) * self.lgrmax

        # compute the shape of the dataset (round up to whole grid cells)
        dataset_shape = np.array([self.width, self.height])
        shape = np.ceil(dataset_shape / self.lgrmax) * self.lgrmax

        n = max_cells / chunk.prod()  # number of chunks that fit in memory
        if n >= 2:
            # first multiply the dimension that is closest to full image shape
            smallest, largest = np.argsort(chunk / shape)
            chunk[largest] = min(int(n) * chunk[largest], shape[largest])
            # then multiply the other dimension if necessary
            n = max_cells / chunk.prod()
            if n >= 2:
                chunk[smallest] = min(int(n) * chunk[smallest], shape[smallest])

        width, height = chunk.astype(int).tolist()
        return width, height

    def iter_by_tile(self, width: int, height: int):
        """Iterate over DEM by tiles of given width and height.

        Args:
          width (int): tile width in pixels.
          height (int): tile height in pixels.
        """
        # determine the width of the largest cell
        if width % self.lgrmax != 0 or height % self.lgrmax != 0:
            raise ValueError(
                "width and height should be a multiple of {}".format(self.lgrmax)
            )

        i1 = 0
        j1 = 0
        i2 = int(np.ceil(float(self.width) / width))
        j2 = int(np.ceil(float(self.height) / height))

        # yield the tiles
        with GDALInterface(self.dem_path) as dem_interface:
            for i, j in itertools.product(range(i1, i2), range(j1, j2)):
                x1, y1, x2, y2 = (
                    i * width,
                    j * height,
                    (i + 1) * width,
                    (j + 1) * height,
                )
                mask = dem_interface.bounded_dem_mask_read((x1, y1, x2, y2))
                yield (x1, y1, x2, y2), mask
