from threedigrid_builder.base import Nodes
from threedigrid_builder.constants import CalculationType
from threedigrid_builder.constants import ContentType
from threedigrid_builder.exceptions import SchematisationError

import numpy as np
import pygeos


__all__ = ["embed_nodes"]


def take(arr, idx):
    """Take values at idx from arr, filling NaN where idx is -9999"""
    result = np.take(arr, idx, mode="clip")
    result[idx == -9999] = np.nan
    return result


def embed_nodes(grid):
    """Integrate embedded connection nodes into the 2D cells.

    This changes the grid (nodes as well as lines) inplace.

    Returned are the removed (embedded) nodes.
    """
    is_embedded = (grid.nodes.calculation_type == CalculationType.EMBEDDED) & (
        grid.nodes.content_type == ContentType.TYPE_V2_CONNECTION_NODES
    )
    embedded_idx = np.where(is_embedded)[0]
    n_embedded = len(embedded_idx)
    if n_embedded == 0:
        return Nodes(id=[])

    # The query_bulk returns 2 1D arrays: one with indices into the embedded nodes
    # and one with indices into the cells.
    idx = grid.cell_tree.query_bulk(pygeos.points(grid.nodes.coordinates[embedded_idx]))
    # Address cases of no or multiple cells per embedded node
    unique_embedded_idx, first_unique_index = np.unique(idx[0], return_index=True)
    # multiple cells: just take the first one
    idx = idx[:, first_unique_index]
    combs_embedded_idx = embedded_idx[idx[0]]
    combs_cell_idx = idx[1]
    if len(unique_embedded_idx) != n_embedded:
        missing_idx = embedded_idx[~np.isin(embedded_idx, combs_embedded_idx)]
        missing_cn_id = grid.nodes.content_pk[missing_idx].tolist()
        raise SchematisationError(
            f"Embedded connection nodes {missing_cn_id} are outside the 2D cells"
        )

    # create a mapping with new node ids on the position of the old node indices
    new_ids = np.full_like(grid.nodes.id, fill_value=-9999)
    new_ids[~is_embedded] = np.arange(len(grid.nodes) - n_embedded)
    new_ids[combs_embedded_idx] = grid.nodes.id[combs_cell_idx]

    # cut out the embedded nodes
    embedded_nodes = grid.nodes[combs_embedded_idx]
    embedded_nodes.id[:] = np.arange(len(embedded_nodes))
    embedded_nodes.embedded_in[:] = new_ids[combs_cell_idx]
    grid.nodes = grid.nodes[~is_embedded]

    # map node indices (fixes node replacement and node renumbering in one go)
    grid.nodes.id[:] = np.arange(len(grid.nodes))
    grid.lines.line[:] = np.take(new_ids, grid.lines.line)

    # check if there are no cells connecting to itself
    is_self_connected = grid.lines.line[:, 0] == grid.lines.line[:, 1]
    if np.any(is_self_connected):
        mask = np.isin(
            embedded_nodes.embedded_in, grid.lines.line[is_self_connected, 0]
        )
        raise SchematisationError(
            f"Embedded connection nodes {embedded_nodes.content_pk[mask].tolist()} connect to "
            f"one another within the same 2D cell."
        )

    return embedded_nodes
