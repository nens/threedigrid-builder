from threedigrid_builder.base import array_of
from threedigrid_builder.base import Lines
from threedigrid_builder.base import Nodes
from threedigrid_builder.constants import BoundaryType
from threedigrid_builder.constants import CalculationType
from threedigrid_builder.constants import ContentType
from threedigrid_builder.constants import LineType
from threedigrid_builder.constants import NodeType
from threedigrid_builder.exceptions import SchematisationError

import itertools
import numpy as np
import pygeos


__all__ = ["BoundaryConditions1D", "BoundaryConditions2D"]


class BoundaryCondition1D:
    id: int
    boundary_type: BoundaryType
    connection_node_id: int


@array_of(BoundaryCondition1D)
class BoundaryConditions1D:
    def apply(self, grid):
        """Apply 1D boundary conditions to a Grid

        Fields on nodes/lines that have a boundary condition be adjusted:
        - nodes.calculation_type: set (overridden) to BOUNDARY_NODE
        - nodes.boundary_id: from BoundaryConditions1D.id
        - nodes.boundary_type: from BoundaryConditions1D.boundary_type
        - nodes.node_type: set to NODE_1D_BOUNDARIES
        - lines.boundary_id: from BoundaryConditions1D.id

        Note:
        - it is assumed that nodes.content_pk is sorted for connection nodes
        - this method should go before 1D-2D connections and node embedding
        """
        # lookup corresponding nodes
        cn_idx = np.where(
            grid.nodes.content_type == ContentType.TYPE_V2_CONNECTION_NODES
        )[0]
        idx = cn_idx[
            np.searchsorted(grid.nodes.content_pk[cn_idx], self.connection_node_id)
        ]
        if not np.all(grid.nodes.content_pk[idx] == self.connection_node_id):
            missing = np.where(grid.nodes.content_pk[idx] != self.connection_node_id)[0]
            raise SchematisationError(
                f"1D boundary conditions {self.id[missing].tolist()} refer to missing "
                f"connection nodes {self.connection_node_id[missing].tolist()}."
            )

        # lookup the lines connecting to these nodes
        ids = grid.nodes.index_to_id(idx)
        line_idx, start_or_end = np.where(np.isin(grid.lines.line, ids))
        # lookup the boundary condition node and other node corresponding to these lines
        bc_node_id = grid.lines.line[line_idx, start_or_end]
        # this trick transforms 0 to 1 and 1 to 0:
        start_or_end_inv = (~(start_or_end.astype(bool))).astype(int)
        other_node_id = grid.lines.line[line_idx, start_or_end_inv]

        if not np.isin(ids, bc_node_id).all():
            missing = ~np.isin(ids, bc_node_id)
            raise SchematisationError(
                f"1D boundary conditions {self.id[missing].tolist()} refer to connection nodes "
                f"that have no attached objects ({self.connection_node_id[missing].tolist()})."
            )
        if len(ids) != len(bc_node_id):
            sorted_cn_id, counts = np.unique(bc_node_id, return_counts=True)
            duplicate_cn_id = sorted_cn_id[counts > 1]
            raise SchematisationError(
                f"Only one attached object is allowed for connection nodes "
                f"{duplicate_cn_id.tolist()} with a 1D boundary condition."
            )
        if np.isin(ids, other_node_id).any():
            both_sides = np.isin(ids, other_node_id)
            raise SchematisationError(
                f"1D boundary conditions {self.id[both_sides].tolist()} are too close to other "
                f"boundary conditions."
            )

        # set node attributes
        grid.nodes.calculation_type[idx] = CalculationType.BOUNDARY_NODE
        grid.nodes.boundary_id[idx] = self.id
        grid.nodes.boundary_type[idx] = self.boundary_type
        grid.nodes.node_type[idx] = NodeType.NODE_1D_BOUNDARIES
        # set line attributes
        grid.lines.is_1d_boundary[line_idx] = 1


class BoundaryCondition2D:
    id: int
    boundary_type: BoundaryType
    the_geom: pygeos.Geometry


@array_of(BoundaryCondition2D)
class BoundaryConditions2D:
    def get_nodes_and_lines(
        self, nodes, cell_tree, quadtree, node_id_counter, line_id_counter
    ):
        # we pad quad_idx to prevent out of bound errors later
        quad_idx = np.pad(quadtree.quad_idx, 1, mode="constant", constant_values=0)
        cells_done = np.empty((0,), dtype=int)

        boundary_nodes = Nodes(id=[])
        boundary_lines = Lines(id=[])
        for bc_idx in range(len(self)):
            bc_geom = self.the_geom[bc_idx]

            x1, y1, x2, y2 = pygeos.bounds(bc_geom)
            is_horizontal = (x2 - x1) > (y2 - y1)

            node_idx = np.sort(cell_tree.query(bc_geom))
            if node_idx.size == 0:
                # Pragmatic fix; because of reprojection issues, the boundary condition
                # can be right outside the cells.

                # Start by getting the closest point inside the model area (all cells)
                nearest_cell = cell_tree.geometries[cell_tree.nearest(bc_geom)[1][0]]
                shortest_line = pygeos.shortest_line(nearest_cell, bc_geom)
                x1, _, x2, _ = pygeos.bounds(nearest_cell)
                if pygeos.length(shortest_line) < (x2 - x1):
                    cell_coords, _ = pygeos.get_coordinates(shortest_line)

                    # Adjust the x / y coordinate to the x / y coordinate of the point
                    bc_coords = pygeos.get_coordinates(bc_geom)
                    bc_coords[:, int(is_horizontal)] = cell_coords[int(is_horizontal)]
                    bc_geom = pygeos.set_coordinates(bc_geom, bc_coords)

                    # Query again for intersecting cells
                    node_idx = np.sort(cell_tree.query(bc_geom))

            nodk = nodes.nodk[node_idx]
            sz = 2 ** (nodk - 1)
            # nodm and nodn are 1-based indices into quad_idx. So we subtract 1 from it.
            nodm = nodes.nodm[node_idx] - 1
            nodn = nodes.nodn[node_idx] - 1
            # Resized nodm, nodn (called x, y) so that they index into quad_idx.
            # The +1 accounts for the padding of quad_idx.
            x = nodm * sz + 1
            y = nodn * sz + 1
            # just to be sure, can be removed later:
            assert np.all(quad_idx[x, y] == (node_idx + 1))

            if is_horizontal:  # horizontal BC, looking at vertical neigbours
                # A cell has no bottom (upstream) neighbor if:
                # - the nodgrid entry before (under) at the left & right sides is 0
                before = (quad_idx[x, y - 1] == 0) & (quad_idx[x + sz - 1, y - 1] == 0)
                # A cell has no top (downstream) neighbor if:
                # - the nodgrid entry after (above) at the left & right sides is 0
                after = (quad_idx[x, y + sz] == 0) & (quad_idx[x + sz - 1, y + sz] == 0)
            else:  # vertical BC, looking at horizontal neigbours
                # A cell has no left (upstream) neighbor if:
                # - the nodgrid entry before (to the left) at bottom & top sides are 0
                before = (quad_idx[x - 1, y] == 0) & (quad_idx[x - 1, y + sz - 1] == 0)
                # A cell has no right (downstream) neighbor if:
                # - the nodgrid entry after (to the right) at bottom & top sides is 0
                after = (quad_idx[x + sz, y] == 0) & (quad_idx[x + sz, y + sz - 1] == 0)

            # The edge of the boundary is whichever side has the most edges
            n_edges_before = np.count_nonzero(before)
            n_edges_after = np.count_nonzero(after)
            if n_edges_after == 0 and n_edges_before == 0:
                raise SchematisationError(
                    f"2D boundary condition {self.id[bc_idx]} does not touch any edge "
                    f"cell."
                )
            elif n_edges_after == n_edges_before:
                s = "top and bottom" if is_horizontal else "right and left"
                raise SchematisationError(
                    f"2D boundary condition {self.id[bc_idx]} touches cells that have "
                    f"equal numbers of {s} edges."
                )
            is_before = n_edges_before > n_edges_after  # is_upstream

            # Filter the nodes to those that have edges
            if is_horizontal and is_before:
                kcu = LineType.LINE_2D_BOUNDARY_SOUTH
                node_idx = node_idx[before]
                edge_coord = nodes.bounds[node_idx, 1]
                cross_pix_coords_cols = (0, 1, 2, 1)
            elif is_horizontal:  # and not is_before
                kcu = LineType.LINE_2D_BOUNDARY_NORTH
                node_idx = node_idx[after]
                edge_coord = nodes.bounds[node_idx, 3]
                cross_pix_coords_cols = (0, 3, 2, 3)
            elif is_before:  # and not is_horizontal
                kcu = LineType.LINE_2D_BOUNDARY_WEST
                node_idx = node_idx[before]
                edge_coord = nodes.bounds[node_idx, 0]
                cross_pix_coords_cols = (0, 1, 0, 3)
            else:  # not is_horizontal and not is_before
                kcu = LineType.LINE_2D_BOUNDARY_EAST
                node_idx = node_idx[after]
                edge_coord = nodes.bounds[node_idx, 2]
                cross_pix_coords_cols = (2, 1, 2, 3)

            if edge_coord.size > 1 and np.any(edge_coord[1:] != edge_coord[:-1]):
                coords = np.unique(edge_coord)
                axis = "y" if is_horizontal else "x"
                raise SchematisationError(
                    f"2D boundary condition {self.id[bc_idx]} touches cells with "
                    f"different edge coordinates ({axis}={coords.tolist()})."
                )

            # Copy the nodes and change the attributes into a boundary node
            node_idx = node_idx[~np.in1d(node_idx, cells_done)]
            if len(node_idx) == 0:
                continue

            new_nodes = nodes[node_idx]
            new_nodes.id[:] = list(itertools.islice(node_id_counter, len(new_nodes)))
            new_nodes.boundary_id[:] = self.id[bc_idx]
            new_nodes.boundary_type[:] = self.boundary_type[bc_idx]
            new_nodes.node_type[:] = NodeType.NODE_2D_BOUNDARIES

            # Get the cross_pix_coords before resetting the pixel_coords
            cross_pix_coords = np.array(
                [new_nodes.pixel_coords[:, x] for x in cross_pix_coords_cols]
            ).T

            # a boundary cell has no real place on the pixels:
            new_nodes.pixel_coords[:] = -9999

            # the bounds and coordinates are shifted:
            size = new_nodes.bounds[:, 2] - new_nodes.bounds[:, 0]
            axes = (1, 3) if is_horizontal else (0, 2)
            add_or_subtract = -1 if is_before else 1
            new_nodes.bounds[:, axes[0]] += size * add_or_subtract
            new_nodes.bounds[:, axes[1]] += size * add_or_subtract
            new_nodes.coordinates[:, axes[0]] += size * add_or_subtract
            if is_horizontal:
                new_nodes.nodn[:] += add_or_subtract
            else:
                new_nodes.nodm[:] += add_or_subtract

            # Create new lines
            if is_before:
                line = np.array([new_nodes.id, nodes.index_to_id(node_idx)]).T
            else:
                line = np.array([nodes.index_to_id(node_idx), new_nodes.id]).T

            new_lines = Lines(
                id=itertools.islice(line_id_counter, len(new_nodes)),
                kcu=kcu,
                line=line,
                cross_pix_coords=cross_pix_coords,
                lik=new_nodes.nodk,
            )

            cells_done = np.append(cells_done, node_idx)
            boundary_nodes += new_nodes
            boundary_lines += new_lines

        return boundary_nodes, boundary_lines
