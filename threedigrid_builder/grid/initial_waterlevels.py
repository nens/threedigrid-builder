import numpy as np
import shapely

from threedigrid_builder.base import search
from threedigrid_builder.constants import ContentType


def _compute_for_interpolated_nodes(nodes, cn_idx, objects):
    """Compute the initial waterlevels by interpolating between connection nodes

    This function is to be used for interpolated 1D nodes. It should be called
    after setting bottom levels (dmax) on connection nodes.

    Args:
        nodes (Nodes)
        cn_idx (ndarray of int): masks nodes so that you have only connection nodes
        objects (LinearObjects): channels, pipes, or culverts

    Returns:
        tuple of:
        - mask into nodes
        - initial waterlevel for each interpolated node in the mask
    """
    is_type = nodes.content_type == objects.content_type
    if not is_type.any():
        return is_type, None
    if shapely.is_missing(objects.the_geom).any():
        raise ValueError(
            f"{objects.__class__.__name__} found without a geometry. Call "
            f"set_geometries first."
        )
    lengths = shapely.length(objects.the_geom)

    idx = objects.id_to_index(nodes.content_pk[is_type])
    weights = nodes.s1d[is_type] / lengths[idx]
    if np.any(weights < 0.0) or np.any(weights > 1.0):
        raise ValueError("Encountered nodes outside of the linear object bounds")

    # convert connection_node_start_id to node indexes
    left_cn_idx = search(
        nodes.content_pk,
        objects.connection_node_start_id[idx],
        mask=cn_idx,
        assume_ordered=True,
    )
    right_cn_idx = search(
        nodes.content_pk,
        objects.connection_node_end_id[idx],
        mask=cn_idx,
        assume_ordered=True,
    )

    left = nodes.initial_waterlevel[left_cn_idx]
    right = nodes.initial_waterlevel[right_cn_idx]

    # Some logic to handle 'dry' initial waterlevels:
    # - (both sides nan): interpolated nodes become dry as well
    # - (one side nan): interpolate from initial_waterlevel to dmax
    # - (both value): interpolate between two initial_waterlevels
    has_value_left = np.isfinite(left)
    has_value_right = np.isfinite(right)
    left = np.where(
        ~has_value_left & has_value_right,
        nodes.dmax[left_cn_idx],
        left,
    )
    right = np.where(
        ~has_value_right & has_value_left,
        nodes.dmax[right_cn_idx],
        right,
    )

    return is_type, weights * right + (1 - weights) * left


def compute_initial_waterlevels(nodes, connection_nodes, channels, pipes, culverts):
    """Apply initial waterlevels (global or per connection nodes) to all 1D nodes.

    Bottom levels (dmax) should be set already.

    Args:
        connection_nodes (ConnectionNodes): used to map ids to indices
        channels (Channels)
        pipes (Pipes)
        culverts (Culverts)

    """
    if not np.any(np.isfinite(connection_nodes.initial_waterlevel)):
        return
    is_1d = nodes.content_type == ContentType.TYPE_V2_CONNECTION_NODES
    cn_idx = np.where(is_1d)[0]

    # First set the connection node initial waterlevels
    nodes.initial_waterlevel[cn_idx] = connection_nodes.initial_waterlevel[
        connection_nodes.id_to_index(nodes.content_pk[cn_idx])
    ]

    # Then interpolate
    for objects in (channels, pipes, culverts):
        is_type, initial_waterlevel = _compute_for_interpolated_nodes(
            nodes, cn_idx, objects
        )
        if initial_waterlevel is not None:
            nodes.initial_waterlevel[is_type] = initial_waterlevel
            is_1d |= is_type

    # Replace dry values with bottom levels (dmax)
    is_1d = np.where(is_1d)[0]
    is_dry = ~np.isfinite(nodes.initial_waterlevel[is_1d]) | (
        nodes.initial_waterlevel[is_1d] < nodes.dmax[is_1d]
    )
    nodes.initial_waterlevel[is_1d[is_dry]] = nodes.dmax[is_1d[is_dry]]
