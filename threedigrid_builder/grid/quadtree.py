from threedigrid_builder.base import Lines
from threedigrid_builder.base import Nodes
from threedigrid_builder.constants import LineType
from threedigrid_builder.constants import NodeType
from threedigrid_builder.grid.fwrapper import create_quadtree
from threedigrid_builder.grid.fwrapper import set_2d_computational_nodes_lines
from threedigrid_builder.grid.fwrapper import set_refinement

import itertools
import math
import numpy as np
import pygeos


__all__ = ["QuadTree"]


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

    def __init__(self, subgrid_meta, num_refine_levels, min_gridsize, refinements):

        self.lgrmin = int(min_gridsize / subgrid_meta["pixel_size"])
        # Maximum number of active grid levels in quadtree.
        self.kmax = num_refine_levels
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
        self.lg = np.full(
            (self.mmax[0], self.nmax[0]), self.kmax, dtype=np.int32, order="F"
        )

        if refinements is not None:
            self._apply_refinements(np.array(subgrid_meta["bbox"]), refinements)
        # Array with dimensions of smallest active grid level and contains
        # idx of active grid level for each quadtree cell.
        self.quad_idx = np.empty(
            (self.mmax[0], self.nmax[0]), dtype=np.int32, order="F"
        )

        self.n_cells, self.n_lines = create_quadtree(
            self.kmax,
            self.mmax,
            self.nmax,
            self.min_cell_pixels,
            subgrid_meta["area_mask"],
            self.quad_idx,
            self.lg,
        )

    def __repr__(self):
        return f"<Quadtree object with {self.kmax} refinement levels and {self.n_cells} active computational cells>"  # NOQA

    @property
    def min_cell_pixels(self):
        """Returns minimum number of pixels in smallles computational_cell."""
        return self.lgrmin

    def _apply_refinements(self, subgrid_bbox, refinements):
        """Set active grid levels for based on refinement dict and
        filling lg variable for refinement locations.

        Args:
          dict of refinemnets containing at least
          - geometry
          - refinement_level
          - id

        """
        for i in range(len(refinements.id)):
            geom = np.asfortranarray(pygeos.get_coordinates(refinements.the_geom[i]))
            set_refinement(
                id=refinements.id[i],
                geom=geom,
                level=refinements.refinement_level[i],
                type=pygeos.get_type_id(refinements.the_geom[i]),
                bbox=subgrid_bbox,
                mmax=self.mmax,
                nmax=self.nmax,
                dx=self.dx,
                lg=self.lg,
            )

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
        total_lines = self.n_lines[0] + self.n_lines[1]

        # Line type is always openwater at first init
        kcu = np.full((total_lines,), LineType.LINE_2D_U, dtype="i4", order="F")
        kcu[self.n_lines[0] : self.n_lines[0] + self.n_lines[1]] = LineType.LINE_2D_V

        # Node connection array
        line = np.empty((total_lines, 2), dtype=np.int32, order="F")
        cross_pix_coords = np.full((total_lines, 4), -9999, dtype=np.int32, order="F")

        set_2d_computational_nodes_lines(
            np.array([self.origin[0], self.origin[1]]),
            self.min_cell_pixels,
            self.kmax,
            self.mmax,
            self.nmax,
            self.dx,
            self.lg,
            self.quad_idx,
            area_mask,
            nodk,
            nodm,
            nodn,
            bounds,
            coords,
            pixel_coords,
            line,
            cross_pix_coords,
            self.n_lines[0],
            self.n_lines[1],
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
