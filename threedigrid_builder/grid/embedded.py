import itertools
from collections import defaultdict, deque

import numpy as np
import shapely

from threedigrid_builder.base import Nodes
from threedigrid_builder.constants import CalculationType, ContentType
from threedigrid_builder.exceptions import SchematisationError

from .linear import counts_to_ranges

# This is for judging from what cell to what cell a channel-cell edge intersection goes:
VPOINT_EPSILON = 1e-7


__all__ = ["embed_nodes", "embed_linear_objects"]


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

    # The query returns 2 1D arrays: one with indices into the embedded nodes
    # and one with indices into the cells.
    idx = grid.cell_tree.query(shapely.points(grid.nodes.coordinates[embedded_idx]))
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
    if grid.pumps is not None:
        mask = grid.pumps.line != -9999
        grid.pumps.line[mask] = np.take(new_ids, grid.pumps.line[mask])
    if grid.surface_maps is not None:
        grid.surface_maps.cci[:] = np.take(new_ids, grid.surface_maps.cci)

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


class EmbeddedObjects:
    """Embedded linear objects

    Attributes:
        objects (LinearObjects): the linear objects (input)
        vpoint_s (ndarray of float): positions of velocity points along objects
        vpoint_ch_idx (ndarray of int): index into objects corresponding to vpoint_s
        vpoint_line (length-2 iterable of ndarrays of int): cell ids belonging
          to the velocity points. Each velocity point has two cell ids because the
          point is on a cell edge by definition.
    """

    def __init__(
        self, objects, vpoint_s, vpoint_ch_idx, vpoint_line=None, vpoint_line_s=None
    ):
        self.objects = objects
        self.vpoint_s = vpoint_s
        self.vpoint_ch_idx = vpoint_ch_idx
        self.vpoint_line = vpoint_line
        self.vpoint_line_s = vpoint_line_s

    @property
    def counts(self):
        """The number of lines per channel"""
        return np.bincount(self.vpoint_ch_idx, minlength=len(self.objects))

    def __repr__(self):
        return (
            f"<{self.__class__.__name__} object with "
            f"{len(self.objects)} {self.objects.__class__.__name__} and "
            f"{len(self.vpoint_s)} velocity points>"
        )

    def __getitem__(self, index):
        """Get velocity points at ``index``"""
        return self.__class__(
            self.objects,
            self.vpoint_s[index],
            self.vpoint_ch_idx[index],
            self.vpoint_line[:, index] if self.vpoint_line is not None else None,
            self.vpoint_line_s[index] if self.vpoint_line_s is not None else None,
        )

    def delete(self, index):
        """Delete velocity points at ``index``"""
        if self.vpoint_line is not None:
            vpoint_line = np.delete(self.vpoint_line, index, axis=1)
        else:
            vpoint_line = None
        if self.vpoint_line_s is not None:
            vpoint_line_s = np.delete(self.vpoint_line_s, index)
        else:
            vpoint_line_s = None
        return self.__class__(
            self.objects,
            np.delete(self.vpoint_s, index),
            np.delete(self.vpoint_ch_idx, index),
            vpoint_line,
            vpoint_line_s,
        )

    def insert(
        self, where, vpoint_s, vpoint_ch_idx, vpoint_line=-9999, vpoint_line_s=np.nan
    ):
        """Insert velocity points at ``where``"""
        if self.vpoint_line is not None:
            vpoint_line = np.insert(self.vpoint_line, where, vpoint_line, axis=1)
        else:
            vpoint_line = None
        if self.vpoint_line_s is not None:
            vpoint_line_s = np.insert(self.vpoint_line_s, where, vpoint_line_s, axis=1)
        else:
            vpoint_line_s = None
        return self.__class__(
            self.objects,
            np.insert(self.vpoint_s, where, vpoint_s),
            np.insert(self.vpoint_ch_idx, where, vpoint_ch_idx),
            vpoint_line,
            vpoint_line_s,
        )

    @classmethod
    def from_cell_tree(cls, objects, cell_tree):
        """Construct the velocity points in the embedded objects from a cell tree."""
        # The query returns 2 1D arrays: one with indices into the objects
        # and one with indices into the cells.
        idx = cell_tree.query(objects.the_geom, "intersects")
        # Get the intersections between objects and cell edges: these are the vpoints
        _points = shapely.intersection(
            objects.the_geom[idx[0]],
            shapely.get_exterior_ring(cell_tree.geometries[idx[1]]),
        )
        # Get unique points per channel (there will be duplicates).
        # This also converts linestring intersections; every vertex becomes a point.
        multipoints = shapely.extract_unique_points(
            shapely.geometrycollections(_points, indices=idx[0])
        )
        # Flatten into a 'ragged array' structure (points + indices into objects)
        vpoints, vpoint_ch_idx = shapely.get_parts(multipoints, return_index=True)
        # Measure the location along the objects
        vpoint_s = shapely.line_locate_point(objects.the_geom[vpoint_ch_idx], vpoints)
        return cls(objects, vpoint_s, vpoint_ch_idx)

    def get_nodes(self, node_id_counter):
        """Get the (virtual) Nodes corresponding to the velocity points"""
        ch_start, ch_end = counts_to_ranges(self.counts)
        node_s = (
            np.delete(self.vpoint_s, ch_start) + np.delete(self.vpoint_s, ch_end - 1)
        ) / 2
        self.vpoint_line_s = self.vpoint_s.copy()
        int_nodes = np.delete(np.arange(len(self.vpoint_line_s)), ch_start)
        self.vpoint_line_s[int_nodes] -= node_s
        node_ch_idx = np.delete(self.vpoint_ch_idx, ch_start)
        node_point = shapely.line_interpolate_point(
            self.objects.the_geom[node_ch_idx], node_s
        )
        node_cell_id = np.delete(self.vpoint_line[0], ch_start)
        return Nodes(
            id=itertools.islice(node_id_counter, len(node_ch_idx)),
            coordinates=shapely.get_coordinates(node_point),
            content_type=self.objects.content_type,
            content_pk=self.objects.id[node_ch_idx],
            calculation_type=CalculationType.EMBEDDED,
            s1d=node_s,
            embedded_in=node_cell_id,
        )

    def sort(self):
        """Sort the velocity points by channel and then by position on the channel"""
        return self[np.lexsort((self.vpoint_s, self.vpoint_ch_idx))]

    def clean_vpoints_channel_ending(self):
        """Solve an edge case where a channel ending is at a cell edge.

        The EmbeddedObjects construction erroneously puts velocity points at channel
        endings if the ending is precisely at the cell edge. These are filtered out.
        """
        ch_lengths = shapely.length(self.objects.the_geom)
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
        pnt_a = shapely.line_interpolate_point(
            self.objects.the_geom[self.vpoint_ch_idx], self.vpoint_s - VPOINT_EPSILON
        )
        pnt_b = shapely.line_interpolate_point(
            self.objects.the_geom[self.vpoint_ch_idx], self.vpoint_s + VPOINT_EPSILON
        )
        vp_idx_a, cell_idx_a = cell_tree.query(pnt_a)
        vp_idx_b, cell_idx_b = cell_tree.query(pnt_b)

        # List velocity points that are outside of the model.
        n_vpoints = len(self.vpoint_s)
        no_intersct_a = np.bincount(vp_idx_a, minlength=n_vpoints) == 0
        no_intersct_b = np.bincount(vp_idx_b, minlength=n_vpoints) == 0
        if np.any(no_intersct_a | no_intersct_b):
            bad_objects_ids = self.vpoint_ch_idx[no_intersct_a | no_intersct_b].tolist()
            raise SchematisationError(
                f"{self.objects.__class__.__name__} {sorted(set(bad_objects_ids))} are "
                f"not completely inside the 2D cell."
            )

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
        copy.vpoint_line[1, _first] = copy.vpoint_line[1, _last]
        # Delete the rest
        return copy.delete(line_too_short)

    def fix_zero_points(self):
        """Add a velocity point for objects without any"""
        ch_no_lines = np.where(self.counts == 0)[0]
        if len(ch_no_lines) == 0:
            return self

        insert_where = np.searchsorted(self.vpoint_ch_idx, ch_no_lines)
        return self.insert(
            insert_where,
            vpoint_s=0.5 * shapely.length(self.objects.the_geom[ch_no_lines]),
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
            tuple of cell_idx_a, cell_idx_b, each with same length as self.vpoint_s

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
            To solve this connectivity problem, we recast the data structure into a
            directed graph containing both the cells and the velocity points. A standard
            shortest path algorithm is then used to get the path that crosses the
            least amount of cells.
        """
        n_vpoints = len(self.vpoint_s)
        if len(cell_idx_a) == n_vpoints and len(cell_idx_b) == n_vpoints:
            return cell_idx_a, cell_idx_b

        # Check which vpoints could have multiple cell transitions
        vpoint_has_mult_a = np.bincount(vpoint_idx_a, minlength=n_vpoints) > 1
        vpoint_has_mult_b = np.bincount(vpoint_idx_b, minlength=n_vpoints) > 1
        mult_vpoint_idx = np.where(vpoint_has_mult_a | vpoint_has_mult_b)[0]
        mult_ch_idx = self.vpoint_ch_idx[mult_vpoint_idx]

        # List consecutive ranges of vpoints with multiple cell transitions
        is_first = ((mult_vpoint_idx - np.roll(mult_vpoint_idx, 1)) > 1) | (
            mult_ch_idx != np.roll(mult_ch_idx, 1)
        )
        is_first[0] = True
        is_last = np.roll(is_first, -1)

        # Iterate over consecutive ranges (and keep a global mask of vpoints to keep)
        vpoints_to_keep = np.ones(n_vpoints, dtype=bool)
        for _first, _last in zip(mult_vpoint_idx[is_first], mult_vpoint_idx[is_last]):
            vpoints_to_keep[_first : _last + 1] = False

            # Construct the graph (e.g. {"c2": ["v1", "v2"], "v2": ["c1"], ...}
            graph = defaultdict(list)
            mask_a = np.where((vpoint_idx_a >= _first) & (vpoint_idx_a <= _last))[0]
            for vpoint_idx, cell_idx in zip(vpoint_idx_a[mask_a], cell_idx_a[mask_a]):
                graph[f"c{cell_idx}"].append(f"v{vpoint_idx}")
            mask_b = np.where((vpoint_idx_b >= _first) & (vpoint_idx_b <= _last))[0]
            for vpoint_idx, cell_idx in zip(vpoint_idx_b[mask_b], cell_idx_b[mask_b]):
                graph[f"v{vpoint_idx}"].append(f"c{cell_idx}")

            # Get the shortest path (e.g. ['c2', 'v2', 'c1'])
            # Note: There could be multiple cell options to start/end with, loop over
            # combinations and take the shortest outcome.
            path = None
            for a, b in itertools.product(
                cell_idx_a[vpoint_idx_a == _first], cell_idx_b[vpoint_idx_b == _last]
            ):
                _path = shortest_path(graph, f"c{a}", f"c{b}")
                if path is None or len(_path) < len(path):
                    path = _path
            n_links = int(len(path) / 2)  # len(path) is always odd, but int() floors

            # Reduce the size of cell_idx_a/b, vpoint_idx_a/b
            vpoint_idx_a = np.delete(vpoint_idx_a, mask_a[n_links:])
            cell_idx_a = np.delete(cell_idx_a, mask_a[n_links:])
            vpoint_idx_b = np.delete(vpoint_idx_b, mask_b[n_links:])
            cell_idx_b = np.delete(cell_idx_b, mask_b[n_links:])

            # Overwrite the cell_idx_a/b, vpoint_idx_a/b to connect the right cells
            for i in range(n_links):
                vpoint = int(path[2 * i + 1][1:])
                vpoints_to_keep[vpoint] = True
                vpoint_idx_a[mask_a[i]] = vpoint
                cell_idx_a[mask_a[i]] = int(path[2 * i][1:])
                vpoint_idx_b[mask_b[i]] = vpoint
                cell_idx_b[mask_b[i]] = int(path[2 * i + 2][1:])

        self.vpoint_s = self.vpoint_s[vpoints_to_keep]
        self.vpoint_ch_idx = self.vpoint_ch_idx[vpoints_to_keep]
        return cell_idx_a, cell_idx_b


def embed_linear_objects(
    objects,
    cell_tree,
    embedded_cutoff_threshold,
    embedded_node_id_counter,
):
    """Create embedded nodes for linear objects

    All linear objects are expected to be of EMBEDDED calculation type.
    """
    if len(objects) == 0:
        return Nodes(id=[]), np.empty((0,)), np.empty((0,))

    if cell_tree is None or len(cell_tree) == 0:
        raise SchematisationError(
            f"{objects.__class__.__name__} {objects.id.tolist()} have an embedded "
            f"calculation type while there is no 2D domain."
        )

    embedded_objects = EmbeddedObjects.from_cell_tree(objects, cell_tree)
    embedded_objects = embedded_objects.clean_vpoints_channel_ending()
    embedded_objects = embedded_objects.sort()
    embedded_objects.set_cell_ids(cell_tree)
    embedded_objects = embedded_objects.clean_vpoints_not_crossing()
    if embedded_cutoff_threshold:
        embedded_objects = embedded_objects.merge_close_points(
            embedded_cutoff_threshold
        )
    embedded_objects = embedded_objects.fix_zero_points()

    # Now we create the virtual nodes
    embedded_nodes = embedded_objects.get_nodes(embedded_node_id_counter)

    return embedded_nodes, embedded_objects.vpoint_s, embedded_objects.vpoint_line_s


def _bfs(graph, start):
    """See shortest_path"""
    queue, enqueued = deque([(None, start)]), set([start])
    while queue:
        parent, n = queue.popleft()
        yield parent, n
        new = set(graph[n]) - enqueued
        enqueued |= new
        queue.extend([(n, child) for child in new])


def shortest_path(graph, start, end):
    """Find the shortest path in a directed graph using Breath First Search (BFS)

    Source:
       https://code.activestate.com/recipes/576675/, August 2021, MIT Licensed
    """
    if start == end:
        return [start]
    parents = {}
    for parent, child in _bfs(graph, start):
        parents[child] = parent
        if child == end:
            revpath = [end]
            while True:
                parent = parents[child]
                revpath.append(parent)
                if parent == start:
                    break
                child = parent
            return list(reversed(revpath))
