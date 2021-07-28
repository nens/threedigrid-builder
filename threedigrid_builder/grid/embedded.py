import numpy as np
import pygeos

from threedigrid_builder.constants import CalculationType
from threedigrid_builder.exceptions import SchematisationError

__all__ = ["embed_nodes"]


def embed_nodes(nodes, lines, cell_tree):
    is_embedded = nodes.calculation_type == CalculationType.EMBEDDED
    idx = cell_tree.query_bulk(pygeos.points(nodes.coordinates[is_embedded]))
    cell_ids = nodes.id[idx[0]]
    embedded_ids = nodes.id[is_embedded][idx[1]]

    # go through an N x M array where N is the number of cells with embedded nodes and
    # M the max number of neighbors for 1 node
    unique_vals, counts = np.unique(cell_ids, return_counts=True)
    n = unique_vals.shape[0]
    m = counts.max()

    embedded_ids_arr = np.full((n, m), -9999, dtype=embedded_ids.dtype)
    embedded_ids_arr

    # some checks
    # if len(unique_vals) != np.count_nonzero(is_embedded):
    #     ids_err = set(nodes.id[is_embedded]) - set(embedded_ids)
    #     content_pks = nodes.content_pk[nodes.id_to_index(ids_err)]
    #     raise SchematisationError(
    #         f"Connection nodes {content_pks.tolist()} are embedded outside of the 2D domain."
    #     )

   


    pass
