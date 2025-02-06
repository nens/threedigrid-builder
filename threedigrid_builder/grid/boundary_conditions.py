import itertools
from typing import Iterator, Tuple

import numpy as np
import shapely

from threedigrid_builder.base import Array, Lines, Nodes
from threedigrid_builder.constants import (
    BoundaryType,
    CalculationType,
    ContentType,
    LineType,
    NodeType,
)
from threedigrid_builder.exceptions import SchematisationError

__all__ = ["BoundaryConditions1D", "BoundaryConditions2D"]


GROUNDWATER_BOUNDARY_TYPES = frozenset(
    {
        BoundaryType.GROUNDWATERLEVEL,
        BoundaryType.GROUNDWATERDISCHARGE,
        BoundaryType.GROUNDWATER_TOTAL_DISCHARGE_2D,
    }
)


class BoundaryCondition1D:
    id: int
    type: BoundaryType
    connection_node_id: int


class BoundaryConditions1D(Array[BoundaryCondition1D]):
    def apply(self, grid):
        """Apply 1D boundary conditions to a Grid

        Fields on nodes/lines that have a boundary condition be adjusted:
        - nodes.calculation_type: set (overridden) to BOUNDARY_NODE
        - nodes.boundary_id: from BoundaryConditions1D.id
        - nodes.boundary_type: from BoundaryConditions1D.type
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
        grid.nodes.boundary_type[idx] = self.type
        grid.nodes.node_type[idx] = NodeType.NODE_1D_BOUNDARIES
        # set line attributes
        grid.lines.is_1d_boundary[line_idx] = 1


class BoundaryCondition2D:
    id: int
    type: BoundaryType
    the_geom: shapely.Geometry


class BoundaryConditions2D(Array[BoundaryCondition2D]):
    def get_intersecting_node_idx(
        self, idx: int, cell_tree: shapely.STRtree
    ) -> np.ndarray:
        bc_geom = self.the_geom[idx]

        x1, y1, x2, y2 = shapely.bounds(bc_geom)
        is_horizontal = (x2 - x1) > (y2 - y1)

        node_idx = np.sort(cell_tree.query(bc_geom))
        if node_idx.size == 0:
            # Pragmatic fix; because of reprojection issues, the boundary condition
            # can be right outside the cells.

            # Start by getting the closest point inside the model area (all cells)
            nearest_cell = cell_tree.geometries[cell_tree.nearest(bc_geom)]
            shortest_line = shapely.shortest_line(nearest_cell, bc_geom)
            x1, _, x2, _ = shapely.bounds(nearest_cell)
            if shapely.length(shortest_line) < (x2 - x1):
                cell_coords, _ = shapely.get_coordinates(shortest_line)

                # Adjust the x / y coordinate to the x / y coordinate of the point
                bc_coords = shapely.get_coordinates(bc_geom)
                bc_coords[:, int(is_horizontal)] = cell_coords[int(is_horizontal)]
                bc_geom = shapely.set_coordinates(bc_geom, bc_coords)

                # Query again for intersecting cells
                node_idx = np.sort(cell_tree.query(bc_geom))

        return node_idx, is_horizontal

    def get_neighoring_node_idx(
        self,
        idx: int,
        nodes: Nodes,
        node_idx: np.ndarray,
        quad_idx: np.ndarray,
        is_horizontal: bool,
    ) -> np.ndarray:
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
                f"2D boundary condition {self.id[idx]} does not touch any edge "
                f"cell."
            )
        elif n_edges_after == n_edges_before:
            s = "top and bottom" if is_horizontal else "right and left"
            raise SchematisationError(
                f"2D boundary condition {self.id[idx]} touches cells that have "
                f"equal numbers of {s} edges."
            )

        is_before = n_edges_before > n_edges_after  # is_upstream

        return node_idx[before if is_before else after], is_before

    def adapt_node_idx_groundwater(self, idx: int, nodes: Nodes, node_idx):
        is_groundwater = self.boundary_type_is_groundwater(self.type[idx])
        if is_groundwater:
            n_groundwater_cells = nodes.n_groundwater_cells
            if n_groundwater_cells == 0:
                raise SchematisationError(
                    f"2D boundary condition {self.id[idx]} has groundwater type while "
                    f"there are no groundwater cells."
                )

            node_idx += n_groundwater_cells
        return is_groundwater

    @staticmethod
    def get_kcu(is_horizontal: bool, is_before: bool, is_groundwater: bool) -> LineType:
        return {
            (True, True, False): LineType.LINE_2D_BOUNDARY_SOUTH,
            (True, False, False): LineType.LINE_2D_BOUNDARY_NORTH,
            (False, True, False): LineType.LINE_2D_BOUNDARY_WEST,
            (False, False, False): LineType.LINE_2D_BOUNDARY_EAST,
            (True, True, True): LineType.LINE_2D_GROUNDWATER_BOUNDARY_SOUTH,
            (True, False, True): LineType.LINE_2D_GROUNDWATER_BOUNDARY_NORTH,
            (False, True, True): LineType.LINE_2D_GROUNDWATER_BOUNDARY_WEST,
            (False, False, True): LineType.LINE_2D_GROUNDWATER_BOUNDARY_EAST,
        }[(is_horizontal, is_before, is_groundwater)]

    @staticmethod
    def get_edge_coord_col(kcu: LineType) -> int:
        return {
            LineType.LINE_2D_BOUNDARY_SOUTH: 1,
            LineType.LINE_2D_BOUNDARY_NORTH: 3,
            LineType.LINE_2D_BOUNDARY_WEST: 0,
            LineType.LINE_2D_BOUNDARY_EAST: 2,
            LineType.LINE_2D_GROUNDWATER_BOUNDARY_SOUTH: 1,
            LineType.LINE_2D_GROUNDWATER_BOUNDARY_NORTH: 3,
            LineType.LINE_2D_GROUNDWATER_BOUNDARY_WEST: 0,
            LineType.LINE_2D_GROUNDWATER_BOUNDARY_EAST: 2,
        }[kcu]

    @staticmethod
    def get_cross_pix_coords_cols(kcu: LineType) -> Tuple[int, int, int, int]:
        return {
            LineType.LINE_2D_BOUNDARY_SOUTH: (0, 1, 2, 1),
            LineType.LINE_2D_BOUNDARY_NORTH: (0, 3, 2, 3),
            LineType.LINE_2D_BOUNDARY_WEST: (0, 1, 0, 3),
            LineType.LINE_2D_BOUNDARY_EAST: (2, 1, 2, 3),
            LineType.LINE_2D_GROUNDWATER_BOUNDARY_SOUTH: (0, 1, 2, 1),
            LineType.LINE_2D_GROUNDWATER_BOUNDARY_NORTH: (0, 3, 2, 3),
            LineType.LINE_2D_GROUNDWATER_BOUNDARY_WEST: (0, 1, 0, 3),
            LineType.LINE_2D_GROUNDWATER_BOUNDARY_EAST: (2, 1, 2, 3),
        }[kcu]

    @staticmethod
    def is_horizontal(kcu: LineType) -> bool:
        return kcu in {
            LineType.LINE_2D_BOUNDARY_SOUTH,
            LineType.LINE_2D_BOUNDARY_NORTH,
            LineType.LINE_2D_GROUNDWATER_BOUNDARY_SOUTH,
            LineType.LINE_2D_GROUNDWATER_BOUNDARY_NORTH,
        }

    @staticmethod
    def is_before(kcu: LineType) -> bool:
        return kcu in {
            LineType.LINE_2D_BOUNDARY_SOUTH,
            LineType.LINE_2D_BOUNDARY_WEST,
            LineType.LINE_2D_GROUNDWATER_BOUNDARY_SOUTH,
            LineType.LINE_2D_GROUNDWATER_BOUNDARY_WEST,
        }

    @staticmethod
    def is_groundwater(kcu: LineType) -> bool:
        return kcu in {
            LineType.LINE_2D_GROUNDWATER_BOUNDARY_SOUTH,
            LineType.LINE_2D_GROUNDWATER_BOUNDARY_NORTH,
            LineType.LINE_2D_GROUNDWATER_BOUNDARY_WEST,
            LineType.LINE_2D_GROUNDWATER_BOUNDARY_EAST,
        }

    @staticmethod
    def boundary_type_is_groundwater(type: BoundaryType) -> bool:
        return type in {
            BoundaryType.GROUNDWATERLEVEL,
            BoundaryType.GROUNDWATERDISCHARGE,
            BoundaryType.GROUNDWATER_TOTAL_DISCHARGE_2D,
        }

    def check_edge_coord(
        self, edge_coord: np.ndarray, idx: int, is_horizontal: bool
    ) -> None:
        if edge_coord.size > 1 and np.any(edge_coord[1:] != edge_coord[:-1]):
            coords = np.unique(edge_coord)
            axis = "y" if is_horizontal else "x"
            raise SchematisationError(
                f"2D boundary condition {self.id[idx]} touches cells with "
                f"different edge coordinates ({axis}={coords.tolist()})."
            )

    def create_boundary_cells(
        self,
        idx: int,
        kcu: LineType,
        nodes: Nodes,
        node_idx: np.ndarray,
        node_id_counter: Iterator[int],
    ) -> Nodes:
        boundary_cells = nodes[node_idx]
        boundary_cells.id[:] = list(
            itertools.islice(node_id_counter, len(boundary_cells))
        )
        boundary_cells.boundary_id[:] = self.id[idx]
        boundary_cells.boundary_type[:] = self.type[idx]
        boundary_cells.node_type[:] = (
            NodeType.NODE_2D_GROUNDWATER_BOUNDARIES
            if self.is_groundwater(kcu)
            else NodeType.NODE_2D_BOUNDARIES
        )
        boundary_cells.pixel_coords[:] = -9999

        # the bounds and coordinates are shifted:
        size = boundary_cells.bounds[:, 2] - boundary_cells.bounds[:, 0]
        axes = (1, 3) if self.is_horizontal(kcu) else (0, 2)
        add_or_subtract = -1 if self.is_before(kcu) else 1
        boundary_cells.bounds[:, axes[0]] += size * add_or_subtract
        boundary_cells.bounds[:, axes[1]] += size * add_or_subtract
        boundary_cells.coordinates[:, axes[0]] += size * add_or_subtract
        if self.is_horizontal(kcu):
            boundary_cells.nodn[:] += add_or_subtract
        else:
            boundary_cells.nodm[:] += add_or_subtract
        return boundary_cells

    def create_boundary_lines(
        self,
        kcu: LineType,
        nodes: Nodes,
        node_idx: np.ndarray,
        boundary_cells: Nodes,
        line_id_counter: Iterator[int],
    ) -> Lines:
        if self.is_before(kcu):
            line = np.array([boundary_cells.id, nodes.index_to_id(node_idx)]).T
        else:
            line = np.array([nodes.index_to_id(node_idx), boundary_cells.id]).T

        cross_pix_coords = np.array(
            [
                nodes.pixel_coords[node_idx, x]
                for x in self.get_cross_pix_coords_cols(kcu)
            ]
        ).T
        return Lines(
            id=itertools.islice(line_id_counter, len(boundary_cells)),
            kcu=kcu,
            line=line,
            cross_pix_coords=cross_pix_coords,
            lik=boundary_cells.nodk,
        )

    def get_nodes_and_lines(
        self,
        nodes: Nodes,
        cell_tree: shapely.STRtree,
        quadtree,
        node_id_counter: Iterator[int],
        line_id_counter: Iterator[int],
    ) -> Tuple[Nodes, Lines]:
        # we pad quad_idx to prevent out of bound errors later
        quad_idx = np.pad(quadtree.quad_idx, 1, mode="constant", constant_values=0)

        cells_done = np.empty((0,), dtype=int)

        boundary_nodes = Nodes(id=[])
        boundary_lines = Lines(id=[])
        for bc_idx in range(len(self)):
            node_idx, is_horizontal = self.get_intersecting_node_idx(bc_idx, cell_tree)
            node_idx, is_before = self.get_neighoring_node_idx(
                bc_idx, nodes, node_idx, quad_idx, is_horizontal
            )
            is_groundwater = self.adapt_node_idx_groundwater(bc_idx, nodes, node_idx)

            kcu = self.get_kcu(is_horizontal, is_before, is_groundwater)
            edge_coord = nodes.bounds[node_idx, self.get_edge_coord_col(kcu)]
            self.check_edge_coord(edge_coord, bc_idx, is_horizontal)

            node_idx = node_idx[~np.in1d(node_idx, cells_done)]
            if len(node_idx) == 0:
                continue
            cells_done = np.append(cells_done, node_idx)

            cells_ = self.create_boundary_cells(
                bc_idx, kcu, nodes, node_idx, node_id_counter
            )
            lines_ = self.create_boundary_lines(
                kcu, nodes, node_idx, cells_, line_id_counter
            )

            boundary_nodes += cells_
            boundary_lines += lines_

        return boundary_nodes, boundary_lines
