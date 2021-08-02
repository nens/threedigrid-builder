from .linear import counts_to_ranges
from threedigrid_builder.base import Lines
from threedigrid_builder.base import Nodes
from threedigrid_builder.constants import CalculationType
from threedigrid_builder.constants import ContentType
from threedigrid_builder.exceptions import SchematisationError

import itertools
import numpy as np
import pygeos


__all__ = ["embed_nodes", "embed_channels"]


def take(arr, idx):
    """Take values at idx from arr, filling NaN where idx is -9999"""
    result = np.take(arr, idx, mode="clip")
    result[idx == -9999] = np.nan
    return result


def embed_nodes(grid, embedded_node_id_counter):
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
    embedded_nodes.id[:] = list(
        itertools.islice(embedded_node_id_counter, len(embedded_nodes))
    )
    embedded_nodes.node_type[:] = -9999  # clear the node_type; these are not calc nodes
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


def embed_channels(
    channels,
    connection_nodes,
    cell_tree,
    embedded_node_id_counter,
    line_id_counter,
    connection_node_offset=0,
):
    """Create embedded nodes for channels

    All channels are expected to be of EMBEDDED calculation type.
    """
    # variable shorthands:
    # - ch_ channels
    # - ech_ embedded channels
    if len(channels) == 0:
        return Nodes(id=[]), Lines(id=[])

    if cell_tree is None or len(cell_tree) == 0:
        raise SchematisationError(
            f"Channels {channels.id} have an embedded calculation type "
            f"while there is no 2D domain."
        )

    # The query_bulk returns 2 1D arrays: one with indices into the embedded nodes
    # and one with indices into the cells.
    idx = cell_tree.query_bulk(channels.the_geom, "intersects")
    # Get the velocity points (unique per channel) where the channel crosses a cell edge
    _points = pygeos.intersection(
        channels.the_geom[idx[0]],
        pygeos.get_exterior_ring(cell_tree.geometries[idx[1]]),
    )
    line_vp, line_ch_idx = pygeos.get_parts(
        pygeos.extract_unique_points(
            pygeos.geometrycollections(_points, indices=idx[0])
        ),
        return_index=True,
    )

    # Measure the location of the velocity points along the channels and sort by it
    line_s = pygeos.line_locate_point(channels.the_geom[line_ch_idx], line_vp)
    _sorter = np.lexsort((line_s, line_ch_idx))
    line_s = line_s[_sorter]
    line_ch_idx = line_ch_idx[_sorter]

    # EDGE CASES
    # 1) the start/end node of a channel is exactly at a cell edge
    mask = (line_s > 0.0) & (line_s < pygeos.length(channels.the_geom)[line_ch_idx])
    line_s = line_s[mask]
    line_ch_idx = line_ch_idx[mask]
    # 2) the velocity point 'touches' a cell edge but does not cross it
    pnt_a = pygeos.line_interpolate_point(channels.the_geom[line_ch_idx], line_s - 1e-7)
    pnt_b = pygeos.line_interpolate_point(channels.the_geom[line_ch_idx], line_s + 1e-7)
    _line_idx_a, line_cell_idx_a = cell_tree.query_bulk(pnt_a)
    _line_idx_b, line_cell_idx_b = cell_tree.query_bulk(pnt_b)
    if len(line_cell_idx_a) != len(line_s) or len(line_cell_idx_b) != len(line_s):
        # 2.2) A side of a velocity point is parallel to the cell edge.
        # For example, the following situation:
        #  line_s = [5., 11.]
        #  _line_idx_a, line_cell_idx_a = [0, 1, 1], [6, 6, 8]
        #  _line_idx_b, line_cell_idx_b = [0, 0, 1], [6, 8, 8]
        # Means that line 0 goes from cell 6 to ? and line 1 goes from cell ? to 8
        #
        # Such an undetermined crossings should be merged so that
        #  line_s = [8.]  # halfway
        #  _line_idx_a, line_cell_idx_a = [0], [6]
        #  _line_idx_b, line_cell_idx_b = [0], [8]
        #
        # In general we can have lines that are
        # - n -> ? (type X), ? -> n (type Y), ? -> ? (type Z)
        # Theoretically you would expect first X, then 0-x times Z, then Y.
        #
        # We are going to do the following:
        # 1. Fix all lines of type X by finding the next Y and merging with that
        # 2. Remove all lines of type Z en Y.
        type_yz = np.where(np.bincount(_line_idx_a) > 1)[0]
        type_xz = np.where(np.bincount(_line_idx_b) > 1)[0]
        type_y = np.setdiff1d(type_yz, type_xz)
        type_x = np.setdiff1d(type_xz, type_yz)
        for i in type_x:
            to_fix = np.where(_line_idx_b == i)[0]
            for j in itertools.count(i + 1):
                if j >= len(line_s) or line_ch_idx[i] != line_ch_idx[j]:
                    # This is a type x crossing without a type Y afterwards. This could
                    # happen when the channel end is on the cell edge. Ignore this
                    # crossing.
                    type_yz = np.append(type_yz, i)
                    break
                if j in type_y:
                    line_cell_idx_b[to_fix[0]] = line_cell_idx_b[_line_idx_b == j]
                    _line_idx_b[to_fix[1:]] = j  # so that this is deleted later
                    line_s[i] = 0.5 * (line_s[i] + line_s[j])
                    break
        line_s = np.delete(line_s, type_yz)
        line_ch_idx = np.delete(line_ch_idx, type_yz)
        line_cell_idx_a = line_cell_idx_a[~np.isin(_line_idx_a, type_yz)]
        line_cell_idx_b = line_cell_idx_b[~np.isin(_line_idx_b, type_yz)]
    line_touches = np.where(line_cell_idx_a == line_cell_idx_b)[0]
    if len(line_touches) > 0:
        line_s = np.delete(line_s, line_touches)
        line_ch_idx = np.delete(line_ch_idx, line_touches)
        line_cell_idx_a = np.delete(line_cell_idx_a, line_touches)
        line_cell_idx_b = np.delete(line_cell_idx_b, line_touches)
    # 3) there is no velocity point (embedded channel completely in 1 cell)
    # --> add the velocity point halfway the channel (it is not at a cell edge)
    ch_n_segments = np.bincount(line_ch_idx, minlength=len(channels))
    ch_no_segments = np.where(ch_n_segments == 0)[0]
    if len(ch_no_segments) > 0:
        insert_where = np.searchsorted(line_ch_idx, ch_no_segments)
        line_s = np.insert(
            line_s,
            insert_where,
            0.5 * pygeos.length(channels.the_geom[ch_no_segments]),
        )
        line_ch_idx = np.insert(line_ch_idx, insert_where, ch_no_segments)
        line_cell_idx_a = np.insert(line_cell_idx_a, insert_where, -9999)
        line_cell_idx_b = np.insert(line_cell_idx_b, insert_where, -9999)
        ch_n_segments[ch_no_segments] = 1

    # The virtual nodes are halfway
    ch_start, ch_end = counts_to_ranges(ch_n_segments)
    node_s = (np.delete(line_s, ch_start) + np.delete(line_s, ch_end - 1)) / 2
    node_ch_idx = np.delete(line_ch_idx, ch_start)
    node_point = pygeos.line_interpolate_point(channels.the_geom[node_ch_idx], node_s)
    node_cell_idx = np.delete(line_cell_idx_a, ch_start)

    # Now we create the virtual nodes
    embedded_nodes = Nodes(
        id=itertools.islice(embedded_node_id_counter, len(node_ch_idx)),
        coordinates=pygeos.get_coordinates(node_point),
        content_type=ContentType.TYPE_V2_CHANNEL,
        content_pk=channels.id[node_ch_idx],
        calculation_type=CalculationType.EMBEDDED,
        s1d=node_s,
        embedded_in=node_cell_idx,  # cell_tree is constructed such that idx == ids
    )

    lines = channels.get_lines(
        connection_nodes,
        None,
        embedded_nodes,
        line_id_counter,
        connection_node_offset=connection_node_offset,
        line_id_attr="embedded_in",  # take .embedded_in for filling lines.line
    )
    # override the velocity point locations (the defaulted to the line midpoint, while
    # for embedded channels we force them to the cell edges)
    lines.s1d[:] = line_s

    return embedded_nodes, lines
