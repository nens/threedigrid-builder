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
        - lines.kcu: set to LINE_1D_BOUNDARY

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
        grid.lines.kcu[line_idx] = LineType.LINE_1D_BOUNDARY


class BoundaryCondition2D:
    id: int
    boundary_type: BoundaryType
    the_geom: pygeos.Geometry


@array_of(BoundaryCondition2D)
class BoundaryConditions2D:
    def get_nodes(self, grid, node_id_counter):
        min_length = grid.quadtree.dx[0] / 2

        nodes = Nodes(id=[])
        for bc_idx in range(len(self)):
            x1, y1, x2, y2 = pygeos.bounds(self.the_geom[bc_idx])
            is_horizontal = (x2 - x1) > (y2 - y1)
            length = x2 - x1 if is_horizontal else y2 - y1
            if length < min_length:
                raise SchematisationError(
                    f"Boundary condition {self.id[bc_idx]} is too small."
                )

            node_idx = np.sort(grid.cell_tree.query(self.the_geom[bc_idx]))
            nodk = grid.nodes.nodk[node_idx]
            sz = 2 ** (nodk - 1)
            nodm = grid.nodes.nodm[node_idx] - 1  # nodm indexes 1-based into quad_idx
            nodn = grid.nodes.nodn[node_idx] - 1  # nodn indexes 1-based into quad_idx
            quad_idx = grid.quadtree.quad_idx  # indexes in quad_idx are 1-based
            assert np.all(quad_idx[nodm * sz, nodn * sz] == (node_idx + 1))

            if is_horizontal:
                b = nodn > 0
                b[b] = quad_idx[nodm[b] * sz[b], nodn[b] * sz[b] - 1] != 0
                b[b] = quad_idx[(nodm[b] + 1) * sz[b] - 1, nodn[b] * sz[b] - 1] != 0
                a = nodn + sz < quad_idx.shape[1]
                a[a] = quad_idx[nodm[a] * sz[a], (nodn[a] + 1) * sz[a]] != 0
                a[a] = quad_idx[(nodm[a] + 1) * sz[a] - 1, (nodn[a] + 1) * sz[a]] != 0
            else:
                b = nodm > 0
                b[b] = quad_idx[nodm[b] * sz[b] - 1, nodn[b] * sz[b]] != 0
                b[b] = quad_idx[nodm[b] * sz[b] - 1, (nodn[b] + 1) * sz[b] - 1] != 0
                a = (nodm + 1) * sz < quad_idx.shape[0]
                a[a] = quad_idx[(nodm[a] + 1) * sz[a], nodn[a] * sz[a]] != 0
                a[a] = quad_idx[(nodm[a] + 1) * sz[a], (nodn[a] + 1) * sz[a] - 1] != 0

            # the edge of the boundary is whichever side has the most edges
            if np.count_nonzero(a) > np.count_nonzero(b):
                node_idx = node_idx[a]
                kcu = (
                    LineType.LINE_2D_BOUNDARY_SOUTH
                    if is_horizontal
                    else LineType.LINE_2D_BOUNDARY_WEST
                )
                is_before = False
            else:
                node_idx = node_idx[b]
                kcu = (
                    LineType.LINE_2D_BOUNDARY_NORTH
                    if is_horizontal
                    else LineType.LINE_2D_BOUNDARY_EAST
                )
                is_before = True

            # check the edge location (coordinate should be the same for all cells)
            bounds = grid.nodes.bounds[node_idx]
            if kcu == LineType.LINE_2D_BOUNDARY_WEST:
                edge_coord = bounds[:, 0]
            elif kcu == LineType.LINE_2D_BOUNDARY_EAST:
                edge_coord = bounds[:, 2]
            elif kcu == LineType.LINE_2D_BOUNDARY_NORTH:
                edge_coord = bounds[:, 3]
            elif kcu == LineType.LINE_2D_BOUNDARY_SOUTH:
                edge_coord = bounds[:, 1]

            if edge_coord.size > 1 and np.any(edge_coord[1:] != edge_coord[:-1]):
                coords = np.unique(edge_coord)
                axis = "y" if is_horizontal else "x"
                raise SchematisationError(
                    f"2D boundary condition {self.id[bc_idx]} touches cells with "
                    f"different edge coordinates ({axis}={coords.tolist()})."
                )

            # now create the new nodes
            new_nodes = grid.nodes[node_idx]
            new_nodes.id[:] = list(itertools.islice(node_id_counter, len(new_nodes)))
            new_nodes.boundary_id[:] = self.id[bc_idx]
            new_nodes.boundary_type[:] = self.boundary_type[bc_idx]
            new_nodes.node_type[:] = NodeType.NODE_2D_BOUNDARIES
            new_nodes.nodk[:] = -9999
            new_nodes.nodm[:] = -9999
            new_nodes.nodn[:] = -9999
            size = new_nodes.bounds[:, 2] - new_nodes.bounds[:, 0]
            size_px = new_nodes.pixel_coords[:, 0] - new_nodes.pixel_coords[:, 2]
            if is_horizontal:
                axes = (1, 3)
            else:
                axes = (0, 2)
            if is_before:
                size *= -1
                size_px *= -1
            new_nodes.bounds[:, axes[0]] -= size
            new_nodes.bounds[:, axes[1]] -= size
            new_nodes.coordinates[:, axes[0]] -= size
            new_nodes.pixel_coords[:, axes[0]] -= size_px
            new_nodes.pixel_coords[:, axes[1]] -= size_px

            nodes += new_nodes

        return nodes
