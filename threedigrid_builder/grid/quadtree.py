import itertools
import logging
import math

import numpy as np

from threedigrid_builder.base import Lines, Nodes, Quarters
from threedigrid_builder.constants import LineType, NodeType
from threedigrid_builder.exceptions import SchematisationError

from .fgrid import m_cells, m_quadtree

logger = logging.getLogger(__name__)


__all__ = ["QuadTree"]


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


class QuadTree:
    """Defines active cell levels for computational grid."""

    lgrmin: int
    pixel_size: float
    origin: float
    kmax: int
    mmax: int
    nmax: int
    dx: float
    lg: int
    quad_idx: int
    transform: float

    def __init__(
        self, subgrid_meta, num_refine_levels, min_gridsize, use_2d_flow, refinements
    ):
        pixel_size = subgrid_meta["pixel_size"]
        min_num_pix = min_gridsize / pixel_size
        if min_num_pix % 2 == 0:
            self.lgrmin = int(min_num_pix)
        else:
            raise SchematisationError(
                f"Smallest 2D grid cell does not contain an even number of pixels. "
                f"Smallest 2D grid cell size: {min_gridsize}m. Pixel size: {pixel_size}m."
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
        self.origin = (subgrid_meta["bbox"][0], subgrid_meta["bbox"][1])
        self.pixel_size = subgrid_meta["pixel_size"]

        lvl_multiplr = 2 ** np.arange(self.kmax - 1, -1, -1, dtype=np.int32)

        # Determine number of largest cells that fit over subgrid extent.
        max_pix_largest_cell = self.min_cell_pixels * lvl_multiplr[0]
        max_large_cells_col = int(
            math.ceil(subgrid_meta["width"] / (max_pix_largest_cell))
        )
        max_large_cells_row = int(
            math.ceil(subgrid_meta["height"] / (max_pix_largest_cell))
        )

        # Calculate grid level dimensions.
        self.mmax[:] = max_large_cells_col * lvl_multiplr
        self.nmax[:] = max_large_cells_row * lvl_multiplr
        self.dx[:] = self.lgrmin * lvl_multiplr[::-1] * subgrid_meta["pixel_size"]

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

        m_quadtree.make_quadtree(
            self.kmax,
            self.mmax,
            self.nmax,
            self.min_cell_pixels,
            use_2d_flow,
            subgrid_meta["area_mask"],
            self.lg,
            self.quad_idx,
            self.n_cells,
            self.n_lines_u,
            self.n_lines_v,
        )

        if self.n_cells.sum() == 0:
            raise SchematisationError(
                "There are no 2D cells because the raster contains no active pixels"
            )

    def __repr__(self):
        return f"<Quadtree object with {self.kmax} refinement levels and {self.n_cells} active computational cells>"  # NOQA

    @property
    def min_cell_pixels(self):
        """Returns minimum number of pixels in smallles computational_cell."""
        return self.lgrmin

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

    def get_nodes_lines(self, area_mask, node_id_counter, line_id_counter):
        """Compute 2D openwater Nodes based on computed Quadtree.

        Return:
          Nodes object
        """

        # Create all node arrays for filling in external Fortran routine.

        # We need a counter for all nodes.
        id_n = itertools.islice(node_id_counter, self.n_cells)
        nodk = np.empty((self.n_cells,), dtype=np.int32, order="F")
        nodm = np.empty((self.n_cells,), dtype=np.int32, order="F")
        nodn = np.empty((self.n_cells,), dtype=np.int32, order="F")
        bounds = np.empty((self.n_cells, 4), dtype=np.float64, order="F")
        coords = np.empty((self.n_cells, 2), dtype=np.float64, order="F")
        pixel_coords = np.empty((self.n_cells, 4), dtype=np.int32, order="F")

        # Node type is always openwater at first init
        node_type = np.full(
            (self.n_cells,), NodeType.NODE_2D_OPEN_WATER, dtype="O", order="F"
        )

        # Create all line array for filling in external Fortran routine
        total_lines = self.n_lines_u + self.n_lines_v

        # Line type is always openwater at first init
        kcu = np.full((total_lines,), LineType.LINE_2D_U, dtype="i4", order="F")
        kcu[self.n_lines_u : total_lines] = LineType.LINE_2D_V

        # Node connection array
        line = np.empty((total_lines, 2), dtype=np.int32, order="F")
        cross_pix_coords = np.full((total_lines, 4), -9999, dtype=np.int32, order="F")

        m_cells.set_2d_computational_nodes_lines(
            np.array([self.origin[0], self.origin[1]], dtype=np.float64),
            self.min_cell_pixels,
            self.kmax,
            self.mmax,
            self.nmax,
            self.dx,
            self.lg,
            nodk,
            nodm,
            nodn,
            self.quad_idx,
            bounds,
            coords,
            pixel_coords,
            area_mask,
            line,
            cross_pix_coords,
            self.n_lines_u,
            self.n_lines_v,
        )

        idx = line[np.arange(total_lines), np.argmin(nodk[line[:, :]], axis=1)]
        lik = nodk[idx]
        lim = nodm[idx]
        lin = nodn[idx]

        nodes = Nodes(
            id=id_n,
            node_type=node_type,
            nodk=nodk,
            nodm=nodm,
            nodn=nodn,
            bounds=bounds,
            coordinates=coords,
            pixel_coords=pixel_coords,
            has_dem_averaged=0,
        )

        if np.any(line[:, 0] == line[:, 1]):
            raise RuntimeError("Quadtree produced self-connected lines")

        lines = Lines(
            id=itertools.islice(line_id_counter, len(line)),
            kcu=kcu,
            line=line,
            lik=lik,
            lim=lim,
            lin=lin,
            cross_pix_coords=cross_pix_coords,
        )

        return nodes, lines

    def get_quarters_admin(self, nodes, lines):
        """Compute 2D openwater quarter administration based on computed Quadtree.
        Determine the connecting flowlines for a quarter in u and v direction.
        Determine the adjacent calculation nodes for a quarter in u and v direction.

        Return:
          Quarters object
        """

        n_2d_nodes = np.count_nonzero(
            np.isin(nodes.node_type, NodeType.NODE_2D_OPEN_WATER)
        )
        n_2d_bnd_nodes = np.count_nonzero(
            np.isin(nodes.node_type, NodeType.NODE_2D_BOUNDARIES)
        )

        quarter_line = np.full((4 * n_2d_nodes, 2), -9999, dtype=np.int32, order="F")
        neighbour_node = np.full((4 * n_2d_nodes, 2), -9999, dtype=np.int32, order="F")

        m_cells.set_quarter_admin(
            nodes.nodk,
            nodes.nodm,
            nodes.nodn,
            lines.line,
            lines.kcu,
            quarter_line,
            neighbour_node,
            self.n_lines_u,
            self.n_lines_v,
            n_2d_bnd_nodes,
        )

        return Quarters(
            id=np.arange(0, 4 * n_2d_nodes),
            line=quarter_line,
            neighbour_node=neighbour_node,
        )
