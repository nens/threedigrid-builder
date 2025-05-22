import itertools

import numpy as np

from threedigrid_builder.base import Lines, Nodes
from threedigrid_builder.constants import LineType, NodeType

from .fgrid import m_clone

__all__ = ["Clone"]


class Clone:
    """Defines active clone cells for computational grid."""

    def __init__(
        self,
        clone_array,
        clone_mask,
        quadtree,
        line,
        cross_pix_coords,
        nodk,
        nodm,
        nodn,
        bounds,
        coords,
        pixel_coords,
        area_mask,
    ):
        n_clones = clone_array.max() + 1
        if n_clones > 0:
            self.n_cells = np.copy(quadtree.n_cells)
            self.cell_numbering = np.zeros(
                (self.n_cells), dtype=np.int32, order="F"
            )  # For the new numbering of the quadtree cells
            self.clone_numbering = np.zeros(
                (n_clones), dtype=np.int32, order="F"
            )  # For the numbering of the clone cells
            self.clones_in_cell = np.zeros(
                (self.n_cells), dtype=np.int32, order="F"
            )  # Total number of clones in each quadtree cell
            self.n_newlines_u = np.copy(
                quadtree.n_lines_u
            )  # Total number of new 2D lines in u direction
            self.n_newlines_v = np.copy(
                quadtree.n_lines_v
            )  # Total number of new 2D lines in v direction
            self.n_interclone_lines = np.array(
                0, dtype=np.int32, order="F"
            )  # Total number of the lines linking 2 clones with the same quadtree cell
            lgrmin = quadtree.lgrmin

            ## To find the active clone cells and renumbering the whole cells (including clones)
            # TODO: Use -9999 for non-existing clone numbering.
            m_clone.find_active_clone_cells(
                self.n_cells,
                clone_array,
                self.cell_numbering,
                self.clone_numbering,
                self.clones_in_cell,
            )

            ## Count the new lines
            m_clone.make_clones(
                lgrmin,
                area_mask,
                clone_mask,
                self.clones_in_cell,
                self.cell_numbering,
                self.clone_numbering,
                line,
                nodk,
                nodm,
                nodn,
                self.n_newlines_u,
                self.n_newlines_v,
                self.n_interclone_lines,
            )

            # this includes all lines with new numbering for the cells
            total_line_number = (
                self.n_newlines_u + self.n_newlines_v + self.n_interclone_lines
            )
            self.line_new = np.empty(
                (total_line_number, 2),
                dtype=np.int32,
                order="F",
            )

            self.cross_pix_coords_new = np.full(
                (total_line_number, 4), -9999, dtype=np.int32, order="F"
            )

            ## Creat the new line administration
            m_clone.set_lines_with_clones(
                quadtree.n_lines_u,
                quadtree.n_lines_v,
                self.n_interclone_lines,
                lgrmin,
                area_mask,
                clone_mask,
                clone_array,
                self.clones_in_cell,
                self.cell_numbering,
                self.clone_numbering,
                line,
                nodk,
                nodm,
                nodn,
                self.line_new,
                self.cross_pix_coords_new,
                cross_pix_coords,
            )

            self.nodk_new = np.empty((self.n_cells), dtype=np.int32, order="F")
            self.nodm_new = np.empty((self.n_cells), dtype=np.int32, order="F")
            self.nodn_new = np.empty((self.n_cells), dtype=np.int32, order="F")
            self.bounds = np.empty((self.n_cells, 4), dtype=np.float64, order="F")
            self.coords = np.empty((self.n_cells, 2), dtype=np.float64, order="F")
            self.pixel_coords = np.empty((self.n_cells, 4), dtype=np.int32, order="F")

            m_clone.reset_nod_parameters(
                clone_array,
                nodk,
                nodm,
                nodn,
                self.nodk_new,
                self.nodm_new,
                self.nodn_new,
                bounds,
                coords,
                pixel_coords,
                self.bounds,
                self.coords,
                self.pixel_coords,
            )

        else:
            print("No clone cells are created")

    def update(self, quadtree, clone, nodes, lines, node_id_counter, line_id_counter):
        quadtree.n_cells = clone.n_cells
        quadtree.n_lines_u = clone.n_newlines_u
        quadtree.n_lines_v = clone.n_newlines_v
        total_lines = clone.n_newlines_u + clone.n_newlines_v + clone.n_interclone_lines
        # lines.line = clone.line_new
        lines.kcu = np.full(
            (total_lines,),
            LineType.LINE_2D_U,
            dtype="i4",
            order="F",
        )
        lines.kcu[
            quadtree.n_lines_u : quadtree.n_lines_u + quadtree.n_lines_v
        ] = LineType.LINE_2D_V
        lines.kcu[
            quadtree.n_lines_u + quadtree.n_lines_v : total_lines
        ] = LineType.LINE_INTERCLONE
        # nodes.nodk = clone.nodk_new
        # nodes.nodm = clone.nodm_new
        # nodes.nodn = clone.nodn_new

        node_type = np.full(
            (clone.n_cells,), NodeType.NODE_2D_OPEN_WATER, dtype="O", order="F"
        )

        idx = clone.line_new[
            np.arange(total_lines),
            np.argmin(clone.nodk_new[clone.line_new[:, :]], axis=1),
        ]
        # lines.lik = nodes.nodk[idx]
        # lines.lim = nodes.nodm[idx]
        # lines.lin = nodes.nodn[idx]

        nodes = Nodes(
            id=itertools.islice(node_id_counter, clone.n_cells),
            node_type=node_type,
            nodk=clone.nodk_new,
            nodm=clone.nodm_new,
            nodn=clone.nodn_new,
            bounds=clone.bounds,
            coordinates=clone.coords,
            pixel_coords=clone.pixel_coords,
            has_dem_averaged=0,
        )

        lines = Lines(
            id=itertools.islice(line_id_counter, len(clone.line_new)),
            kcu=lines.kcu,
            line=clone.line_new,
            lik=nodes.nodk[idx],
            lim=nodes.nodm[idx],
            lin=nodes.nodn[idx],
            cross_pix_coords=clone.cross_pix_coords_new,
        )
        return nodes, lines

    # if self.n_clone_cells.sum() == 0
    #     raise a message to show no clones were created
