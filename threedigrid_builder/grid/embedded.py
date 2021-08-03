from .linear import counts_to_ranges
from threedigrid_builder.base import Lines
from threedigrid_builder.base import Nodes
from threedigrid_builder.constants import CalculationType
from threedigrid_builder.constants import ContentType
from threedigrid_builder.exceptions import SchematisationError

import itertools
import numpy as np
import pygeos


# This is for judging from what cell to what cell a channel-cell edge intersection goes:
VPOINT_EPSILON = 1e-7


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


class EmbeddedChannels:
    """Embedded channels

    Attributes:
        channels (Channels): the channels (input)
        vpoint_s (ndarray of float): positions of velocity points along channels
        vpoint_ch_idx (ndarray of int): index into channels corresponding to vpoint_s
        vpoint_line (length-2 iterable of ndarrays of int): cell ids belonging
          to the velocity points. Each velocity point has two cell ids because the
          point is on a cell edge by definition.
    """

    def __init__(self, channels, vpoint_s, vpoint_ch_idx, vpoint_line=None):
        self.channels = channels
        self.vpoint_s = vpoint_s
        self.vpoint_ch_idx = vpoint_ch_idx
        self.vpoint_line = vpoint_line

    @property
    def counts(self):
        """The number of lines per channel"""
        return np.bincount(self.vpoint_ch_idx, minlength=len(self.channels))

    def __getitem__(self, index):
        """Get velocity points at ``index``"""
        return self.__class__(
            self.channels,
            self.vpoint_s[index],
            self.vpoint_ch_idx[index],
            self.vpoint_line[:, index] if self.vpoint_line is not None else None,
        )

    def delete(self, index):
        """Delete velocity points at ``index``"""
        if self.vpoint_line is not None:
            vpoint_line = np.delete(self.vpoint_line, index, axis=1)
        else:
            vpoint_line = None
        return self.__class__(
            self.channels,
            np.delete(self.vpoint_s, index),
            np.delete(self.vpoint_ch_idx, index),
            vpoint_line,
        )

    def insert(self, where, vpoint_s, vpoint_ch_idx, vpoint_line=-9999):
        """Insert velocity points at ``where``"""
        if self.vpoint_line is not None:
            vpoint_line = np.insert(self.vpoint_line, where, vpoint_line, axis=1)
        else:
            vpoint_line = None
        return self.__class__(
            self.channels,
            np.insert(self.vpoint_s, where, vpoint_s),
            np.insert(self.vpoint_ch_idx, where, vpoint_ch_idx),
            vpoint_line,
        )

    @classmethod
    def from_cell_tree(cls, channels, cell_tree):
        """Construct the velocity points in the embedded channels from a cell tree."""
        # The query_bulk returns 2 1D arrays: one with indices into the channels
        # and one with indices into the cells.
        idx = cell_tree.query_bulk(channels.the_geom, "intersects")
        # Get the intersections between channels and cell edges: these are the vpoints
        _points = pygeos.intersection(
            channels.the_geom[idx[0]],
            pygeos.get_exterior_ring(cell_tree.geometries[idx[1]]),
        )
        # Get unique points per channel (there will be duplicates).
        # This also converts linestring intersections; every vertex becomes a point.
        multipoints = pygeos.extract_unique_points(
            pygeos.geometrycollections(_points, indices=idx[0])
        )
        # Flatten into a 'ragged array' structure (points + indices into channels)
        vpoints, vpoint_ch_idx = pygeos.get_parts(multipoints, return_index=True)
        # Measure the location along the channels
        vpoint_s = pygeos.line_locate_point(channels.the_geom[vpoint_ch_idx], vpoints)
        return cls(channels, vpoint_s, vpoint_ch_idx)

    def get_nodes(self, node_id_counter):
        """Get the (virtual) Nodes corresponding to the velocity points"""
        ch_start, ch_end = counts_to_ranges(self.counts)
        node_s = (
            np.delete(self.vpoint_s, ch_start) + np.delete(self.vpoint_s, ch_end - 1)
        ) / 2
        node_ch_idx = np.delete(self.vpoint_ch_idx, ch_start)
        node_point = pygeos.line_interpolate_point(
            self.channels.the_geom[node_ch_idx], node_s
        )
        node_cell_id = np.delete(self.vpoint_line[0], ch_start)
        return Nodes(
            id=itertools.islice(node_id_counter, len(node_ch_idx)),
            coordinates=pygeos.get_coordinates(node_point),
            content_type=ContentType.TYPE_V2_CHANNEL,
            content_pk=self.channels.id[node_ch_idx],
            calculation_type=CalculationType.EMBEDDED,
            s1d=node_s,
            embedded_in=node_cell_id,
        )

    def sort(self):
        """Sort the velocity points by channel and then by position on the channel"""
        return self[np.lexsort((self.vpoint_s, self.vpoint_ch_idx))]

    def clean_vpoints_channel_ending(self):
        """Solve an edge case where a channel ending is at a cell edge.

        The EmbeddedChannels construction erroneously puts velocity points at channel
        endings if the ending is precisely at the cell edge. These are filtered out.
        """
        ch_lengths = pygeos.length(self.channels.the_geom)
        return self[
            (self.vpoint_s > 0.0) & (self.vpoint_s < ch_lengths[self.vpoint_ch_idx])
        ]

    def set_cell_ids(self, cell_tree):
        """Compute 2 cell ids (from and to) for each velocity point (line), inplace.

        This also deals with edge cases in which a line lies tangent to a cell edge. If
        a line goes exactly over a cell edge, this cell is not included. In other words,
        an embedded channel is only embedded in a cell if the channel intersects with
        the cell interior.
        """
        # For each vpoint, see from which to which cell it goes. This is done by
        # going a bit before and after the vpoint and using the cell_tree.
        pnt_a = pygeos.line_interpolate_point(
            self.channels.the_geom[self.vpoint_ch_idx], self.vpoint_s - VPOINT_EPSILON
        )
        pnt_b = pygeos.line_interpolate_point(
            self.channels.the_geom[self.vpoint_ch_idx], self.vpoint_s + VPOINT_EPSILON
        )
        vp_idx_a, cell_idx_a = cell_tree.query_bulk(pnt_a)
        vp_idx_b, cell_idx_b = cell_tree.query_bulk(pnt_b)
        cell_idx_a, cell_idx_b = self._fix_tangent_lines(
            vp_idx_a, cell_idx_a, vp_idx_b, cell_idx_b
        )
        self.vpoint_line = np.array([cell_idx_a, cell_idx_b], dtype=int)

    def clean_vpoints_not_crossing(self):
        """Remove velocity points that only touch (but not cross) a cell edge"""
        line_touches = self.vpoint_line[0] == self.vpoint_line[1]
        if np.any(line_touches):
            return self[~line_touches]
        else:
            return self

    def merge_close_points(self, threshold):
        """Merge velocity points that are closer than ``threshold`` together"""
        line_ds = np.diff(self.vpoint_s)
        line_too_short = np.where(
            (line_ds < threshold) & (np.diff(self.vpoint_ch_idx) == 0)
        )[0]
        if len(line_too_short) == 0:
            return self
        # Get consecutive ranges of too short segments
        line_too_short += 1
        _is_start = (line_too_short - np.roll(line_too_short, 1)) > 1
        _is_start[0] = True
        _first = line_too_short[_is_start] - 1
        _last = line_too_short[np.roll(_is_start, -1)]
        # Adjust the first of each consecutive range (change vpoint_s and vpoint_line)
        copy = self[:]
        copy.vpoint_s[_first] = (copy.vpoint_s[_first] + copy.vpoint_s[_last]) / 2
        copy.vpoint_line[1][_first] = copy.vpoint_line[0][_last]
        # Delete the rest
        return copy.delete(line_too_short)

    def fix_zero_points(self):
        """Add a velocity point for channels without any"""
        ch_no_lines = np.where(self.counts == 0)[0]
        if len(ch_no_lines) == 0:
            return self

        insert_where = np.searchsorted(self.vpoint_ch_idx, ch_no_lines)
        return self.insert(
            insert_where,
            vpoint_s=0.5 * pygeos.length(self.channels.the_geom[ch_no_lines]),
            vpoint_ch_idx=ch_no_lines,
        )

    def _fix_tangent_lines(self, vpoint_idx_a, cell_idx_a, vpoint_idx_b, cell_idx_b):
        """Fix velocity points whose line around it is on a cell edge, inplace.

        Args:
            vpoint_idx_a (ndarray of int): velocity point index
            cell_idx_a (ndarray of int): index of the cell before the velocity point
            vpoint_idx_b (ndarray of int): velocity point index
            cell_idx_b (ndarray of int): index of the cell after the velocity point

        Returns:
            tuple of cell_idx_a, cell_idx_b with precisly 2 cell ids per velocity point

        Example:
            There is a channel that goes from cell 6 to the edge between 6 and 8 and
            then to cell 8. This channel will give 2 velocity points:
            Vpoint 0 goes from cell 6 to 6/8 and vpoint 1 goes from cell 6/8 to 8

            >>> vpoint_s = [5., 11.]
            >>> vpoint_idx_a, vpoint_line_a = [0, 1, 1], [6, 6, 8]
            >>> vpoint_idx_b, vpoint_line_b = [0, 0, 1], [6, 8, 8]

            These vpoints will be merged into:

            >>> vpoint_s = [8.]  # halfway
            >>> vpoint_idx_a, cell_idx_a = [0], [6]
            >>> vpoint_idx_b, cell_idx_b = [0], [8]

        Notes:
            In general we can have lines that are
            - n -> ? (type X), ? -> n (type Y), ? -> ? (type Z)
            Theoretically you would expect first X, then 0-x times Z, then Y.

            This method does the following:
            1. Fix all lines of type X by finding the next Y and merging with that
            2. Remove all lines of type Z en Y.
        """
        if len(cell_idx_a) == len(self.vpoint_s) and len(cell_idx_b) == len(
            self.vpoint_s
        ):
            return cell_idx_a, cell_idx_b

        type_yz = np.where(np.bincount(vpoint_idx_a) > 1)[0]
        type_xz = np.where(np.bincount(vpoint_idx_b) > 1)[0]
        type_y = np.setdiff1d(type_yz, type_xz)
        type_x = np.setdiff1d(type_xz, type_yz)
        for i in type_x:
            to_fix = np.where(vpoint_idx_b == i)[0]
            for j in itertools.count(i + 1):
                if (
                    j >= len(self.vpoint_s)
                    or self.vpoint_ch_idx[i] != self.vpoint_ch_idx[j]
                ):
                    # This is a type x crossing without a type Y afterwards. This could
                    # happen when the channel end is on the cell edge. Ignore this
                    # crossing.
                    type_yz = np.append(type_yz, i)
                    break
                if j in type_y:
                    cell_idx_b[to_fix[0]] = cell_idx_b[vpoint_idx_b == j]
                    vpoint_idx_b[to_fix[1:]] = j  # so that this is deleted later
                    self.vpoint_s[i] = 0.5 * (self.vpoint_s[i] + self.vpoint_s[j])
                    break

        self.vpoint_s = np.delete(self.vpoint_s, type_yz)
        self.vpoint_ch_idx = np.delete(self.vpoint_ch_idx, type_yz)
        cell_idx_a = cell_idx_a[~np.isin(vpoint_idx_a, type_yz)]
        cell_idx_b = cell_idx_b[~np.isin(vpoint_idx_b, type_yz)]
        return cell_idx_a, cell_idx_b


def embed_channels(
    channels,
    connection_nodes,
    cell_tree,
    embedded_node_id_counter,
    line_id_counter,
    connection_node_offset=0,
    embedded_cutoff_threshold=None,
):
    """Create embedded nodes for channels

    All channels are expected to be of EMBEDDED calculation type.
    """
    if len(channels) == 0:
        return Nodes(id=[]), Lines(id=[])

    if cell_tree is None or len(cell_tree) == 0:
        raise SchematisationError(
            f"Channels {channels.id} have an embedded calculation type "
            f"while there is no 2D domain."
        )

    embedded_channels = EmbeddedChannels.from_cell_tree(channels, cell_tree)
    embedded_channels = embedded_channels.clean_vpoints_channel_ending()
    embedded_channels = embedded_channels.sort()
    embedded_channels.set_cell_ids(cell_tree)
    embedded_channels = embedded_channels.clean_vpoints_not_crossing()
    if embedded_cutoff_threshold:
        embedded_channels = embedded_channels.merge_close_points(
            embedded_cutoff_threshold
        )
    embedded_channels = embedded_channels.fix_zero_points()

    # Now we create the virtual nodes
    embedded_nodes = embedded_channels.get_nodes(embedded_node_id_counter)
    # And construct the lines (also connecting to connection nodes)
    lines = channels.get_lines(
        connection_nodes,
        None,
        embedded_nodes,
        line_id_counter,
        connection_node_offset=connection_node_offset,
        line_id_attr="embedded_in",  # take .embedded_in for filling lines.line
    )
    # Override the velocity point locations (the defaulted to the line midpoint, while
    # for embedded channels we force them to the cell edges)
    lines.s1d[:] = embedded_channels.vpoint_s

    return embedded_nodes, lines
