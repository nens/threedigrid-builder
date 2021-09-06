from threedigrid_builder.base import array_of
from threedigrid_builder.constants import BoundaryType
from threedigrid_builder.constants import CalculationType
from threedigrid_builder.constants import ContentType
from threedigrid_builder.constants import LineType
from threedigrid_builder.constants import NodeType
from threedigrid_builder.exceptions import SchematisationError

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

        Calculation nodes that have a boundary condition will be moved to the end of the
        node list. The following fields are adjusted:
        - id: the nodes are moved to the end, the ids are renumbered accordingly
        - calculation_type: set (overridden) to ISOLATED
        - boundary_type: from BoundaryConditions1D
        - node_type: set to NODE_1D_BOUNDARIES

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
            missing = np.where(grid.nodes.content_pk[idx] == self.connection_node_id)[0]
            raise SchematisationError(
                f"1D boundary conditions {self.id[missing]} refer to non-existing "
                f"connection nodes {self.connection_node_id[missing]}."
            )

        # check the lines connecting to these nodes
        ids = grid.nodes.index_to_id(idx)
        line_idx, start_or_end = np.where(np.isin(grid.lines.line, ids))
        bc_node_id = grid.lines.line[line_idx, start_or_end]
        other_node_id = grid.lines.line[
            line_idx, (~start_or_end.astype(bool)).astype(int)
        ]

        if not np.isin(ids, bc_node_id).all():
            missing = ~np.isin(ids, bc_node_id)
            raise SchematisationError(
                f"1D boundary conditions {self.id[missing]} refer to connection nodes "
                f"that have no attached objects ({self.connection_node_id[missing]})."
            )
        if len(ids) != len(bc_node_id):
            sorted_cn_id, counts = np.unique(bc_node_id, return_counts=True)
            duplicate_cn_id = sorted_cn_id[counts > 1]
            raise SchematisationError(
                f"Only one attached object is allowed for connection nodes "
                f"{duplicate_cn_id} with a 1d boundary condition."
            )
        if np.isin(ids, other_node_id).any():
            both_sides = np.isin(ids, other_node_id)
            raise SchematisationError(
                f"1D boundary conditions {self.id[both_sides]} are too close to other "
                f"boundary conditions."
            )

        ## set node attributes
        grid.nodes.calculation_type[idx] = CalculationType.ISOLATED
        grid.nodes.content_pk[idx] = self.id
        grid.nodes.content_type[idx] = ContentType.TYPE_V2_1D_BOUNDARY_CONDITIONS
        grid.nodes.boundary_type[idx] = self.boundary_type
        grid.nodes.node_type[idx] = NodeType.NODE_1D_BOUNDARIES
        ## move the nodes to the end
        x = np.arange(len(grid.nodes))
        grid.nodes.reorder(np.concatenate([np.delete(x, idx), idx]))
        ## set line attributes
        # note: the LINE_1D_BOUNDARY will be mapped back to LINE_1D_ISOLATED later
        grid.lines.kcu[line_idx] = LineType.LINE_1D_BOUNDARY
        ## map lines.line according to node reordering
        # create a mapping with new node idx on the position of the old node indices
        is_bc = np.isin(np.arange(len(grid.nodes)), idx)
        new_idx = np.full_like(grid.nodes.id, fill_value=-9999)
        new_idx[~is_bc] = np.arange(len(grid.nodes) - len(self))
        new_idx[idx] = np.arange(len(self)) + len(grid.nodes) - len(self)
        grid.lines.line[:] = np.take(new_idx, grid.lines.line)


class BoundaryCondition2D:
    id: int
    boundary_type: BoundaryType
    the_geom: pygeos.Geometry


@array_of(BoundaryCondition2D)
class BoundaryConditions2D:
    pass
