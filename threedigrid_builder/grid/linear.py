from threedigrid_builder.base import Lines
from threedigrid_builder.base import Nodes
from threedigrid_builder.constants import CalculationType
from threedigrid_builder.constants import NodeType

import itertools
import numpy as np
import pygeos


COORD_EQUAL_ATOL = 1e-8  # the distance below which coordinates are considered equal


def sanitize_linestrings(linestrings, points_1, points_2):
    """Make sure that the linestrings go from points_1 to points_2.

    The linestring is reversed if the distance from the linestring start to point_1
    and v.v. is lowered.

    Then, point_1 and/or point_2 are added to the linestring if the distance is larger
    than the tolerance of 1e-8.

    This function acts inplace on linestrings.
    """
    has_geom = np.where(pygeos.is_geometry(linestrings))[0]
    if len(has_geom) == 0:
        return

    # reverse the geometries where necessary
    node_start = points_1[has_geom]
    node_end = points_2[has_geom]
    line_start = pygeos.get_point(linestrings[has_geom], 0)
    line_end = pygeos.get_point(linestrings[has_geom], -1)

    dist_start_start = pygeos.distance(node_start, line_start)
    dist_end_end = pygeos.distance(node_end, line_end)
    dist_start_end = pygeos.distance(node_start, line_end)
    dist_end_start = pygeos.distance(node_end, line_start)
    needs_reversion = has_geom[
        (dist_start_start > dist_start_end) & (dist_end_end > dist_end_start)
    ]
    if len(needs_reversion) > 0:
        linestrings[needs_reversion] = pygeos.reverse(linestrings[needs_reversion])
        line_start[needs_reversion] = pygeos.get_point(linestrings[needs_reversion], 0)
        line_end[needs_reversion] = pygeos.get_point(linestrings[needs_reversion], -1)
        dist_start_start[needs_reversion] = dist_start_end[needs_reversion]
        dist_end_end[needs_reversion] = dist_end_start[needs_reversion]

    # for the remaining geometries; add point if necessary
    add_first = has_geom[dist_start_start > COORD_EQUAL_ATOL]
    if len(add_first) > 0:
        linestrings[add_first] = prepend_point(
            linestrings[add_first], points_1[add_first]
        )
    add_end = has_geom[dist_end_end > COORD_EQUAL_ATOL]
    if len(add_end) > 0:
        linestrings[add_end] = append_point(linestrings[add_end], points_2[add_end])


class BaseLinear:
    content_type = None  # to be defined by subclasses

    def set_geometries(self, connection_nodes):
        """Set the_geom from connection nodes where necessary.

        This also reverses geometries if necessary.

        Args:
            connection_nodes (ConnectionNodes): to take the coordinates from
        """
        # get the start and end points
        points_1 = connection_nodes.the_geom[
            connection_nodes.id_to_index(self.connection_node_start_id)
        ]
        points_2 = connection_nodes.the_geom[
            connection_nodes.id_to_index(self.connection_node_end_id)
        ]
        sanitize_linestrings(self.the_geom, points_1, points_2)

        # construct the geometries where necessary
        has_no_geom = pygeos.is_missing(self.the_geom)
        coordinates = np.empty((np.count_nonzero(has_no_geom), 2, 2))
        coordinates[:, 0, 0] = pygeos.get_x(points_1[has_no_geom])
        coordinates[:, 0, 1] = pygeos.get_y(points_1[has_no_geom])
        coordinates[:, 1, 0] = pygeos.get_x(points_2[has_no_geom])
        coordinates[:, 1, 1] = pygeos.get_y(points_2[has_no_geom])
        self.the_geom[has_no_geom] = pygeos.linestrings(coordinates)

    def interpolate_nodes(self, node_id_counter, global_dist_calc_points):
        """Compute nodes on each linear object with constant intervals

        The following fields are expected to be filled on self:
        - id
        - the_geom
        - dist_calc_points
        - calculation_type

        Args:
            node_id_counter (iterable): an iterable yielding integers
            global_dist_calc_points (float): Default node interdistance.

        Returns:
            tuple of nodes (Nodes)
            nodes has data in the following columns:
            - id: counter generated from node_id_counter
            - coordinates: (x, y) coordinates of the node
            - content_pk: the id of the linear object from which the node originates
            - node_type: NodeType.NODE_1D_NO_STORAGE
            - calculation_type: the calculation type copied from the linear object
            - s1d: distance (along the linestring) to the start of the linestring
            The nodes are ordered by content_pk and then by position on the linestring.
        """
        # Do not add nodes for embedded objects
        not_embedded = self.calculation_type != CalculationType.EMBEDDED
        if np.all(not_embedded):
            not_embedded = slice(None)  # faster than bool mask

        if pygeos.is_missing(self.the_geom).any():
            raise ValueError(
                f"{self.__class__.__name__} encountered without a geometry."
            )

        # normalize global dist_calc_points
        if (
            global_dist_calc_points is None
            or np.isnan(global_dist_calc_points)
            or global_dist_calc_points <= 0.0
        ):
            global_dist_calc_points = np.inf  # means: no interpolation

        # insert default dist_calc_points where necessary
        dists = self.dist_calc_points[not_embedded].copy()
        dists[np.isnan(dists)] = global_dist_calc_points
        dists[dists <= 0] = global_dist_calc_points

        # interpolate the node geometries
        points, index, dist_to_start = segmentize(self.the_geom[not_embedded], dists)

        # construct the nodes with available attributes
        nodes = Nodes(
            id=itertools.islice(node_id_counter, len(points)),
            coordinates=pygeos.get_coordinates(points),
            content_type=self.content_type,
            content_pk=self.id[not_embedded][index],
            node_type=NodeType.NODE_1D_NO_STORAGE,
            calculation_type=self.calculation_type[not_embedded][index],
            s1d=dist_to_start,
        )
        return nodes

    def get_embedded(
        self, cell_tree, embedded_cutoff_threshold, embedded_node_id_counter
    ):
        """ """
        from threedigrid_builder.grid import embed_linear_objects

        return embed_linear_objects(
            self[self.calculation_type == CalculationType.EMBEDDED],
            cell_tree,
            embedded_cutoff_threshold,
            embedded_node_id_counter,
        )

    def get_lines(
        self,
        connection_nodes,
        nodes,
        line_id_counter,
        connection_node_offset=0,
        embedded_mode=False,
    ):
        """Compute the grid lines for the linear objects.

        The following fields are expected to be filled on self:
        - id
        - the_geom
        - connection_node_start_id
        - connection_node_end_id
        - calculation_type
        - cross_section_definition_id (pipes and culverts only)
        - display_name

        Args:
            connection_nodes (ConnectionNodes): used to map ids to indices
            nodes (Nodes): interpolated nodes (see interpolate_nodes)
            line_id_counter (iterable): an iterable yielding integers
            connection_node_offset (int): offset to give connection node
              indices in the returned lines.line. Default 0.
            embedded_mode (bool): whether to get lines for embedded objects
              if embedded_mode is true, the line.line come from node.embedded_in

        Returns:
            Lines with data in the following columns:
            - id: counter generated from line_id_counter
            - line: 2 node ids per line (or a different attribute, see line_id_attr)
            - content_pk: the id of the linear from which this line originates
            - s1d: the positon of the line center along the linear object
            - ds1d: the arclength of the line
            - kcu: the calculation_type of the linear object
            - line_geometries: the linestrings (segments of self.the_geom)
            These attributes are filled for pipes/culverts only:
            - cross_id1 & cross_id2: the id of the cross section definition
            - cross_weight: 1.0 (which means that cross_id2 should be ignored)
            - frict_type1 & frict_type2: the friction type (both are equal)
            - frict_value1 & frict_value2: the friction value (both are equal)
            - invert_level_start_point: copied from pipe/culvert
            - invert_level_end_point: copied from pipe/culvert
            - dpumax: largest of the two invert levels
            - discharge_coefficient_positive & _positive: culverts only
            The lines are ordered by content_pk and then by position on the linestring.
        """
        if embedded_mode:
            objs = self[self.calculation_type == CalculationType.EMBEDDED]
        else:
            objs = self[self.calculation_type != CalculationType.EMBEDDED]

        # count the number of segments per object
        node_line_idx = objs.id_to_index(nodes.content_pk)
        segment_counts = np.bincount(node_line_idx, minlength=len(objs)) + 1

        # cut the channel geometries into segment geometries
        start_s, end_s, segment_idx = segment_start_end(
            objs.the_geom, segment_counts, nodes.s1d
        )
        segments = line_substring(objs.the_geom, start_s, end_s, segment_idx)

        # Copy properties to segments
        display_name = np.take(self.display_name, segment_idx)
        zoom_category = np.take(self.zoom_category, segment_idx)
        connection_node_start_id = np.take(self.connection_node_start_id, segment_idx)
        connection_node_end_id = np.take(self.connection_node_end_id, segment_idx)
        dist_calc_points = np.take(self.dist_calc_points, segment_idx)

        # set the right node indices for each segment
        first_idx, last_idx = counts_to_ranges(segment_counts)
        last_idx -= 1  # convert slice end into last index
        line = np.full((len(segments), 2), -9999, dtype=np.int32)

        # convert connection_node_start_id to index and put it at first segments' start
        line[first_idx, 0] = (
            connection_nodes.id_to_index(objs.connection_node_start_id)
            + connection_node_offset
        )
        # convert connection_node_end_id to index and put it at last segments' end
        line[last_idx, 1] = (
            connection_nodes.id_to_index(objs.connection_node_end_id)
            + connection_node_offset
        )
        # set node indices to line start where segment start is not a conn. node
        mask = np.ones(len(segments), dtype=bool)
        mask[first_idx] = False
        line[mask, 0] = nodes.id if not embedded_mode else nodes.embedded_in
        # set node indices to line end where segment end is not a conn. node
        mask = np.ones(len(segments), dtype=bool)
        mask[last_idx] = False
        line[mask, 1] = nodes.id if not embedded_mode else nodes.embedded_in

        # map to obj ids
        content_pk = objs.id[segment_idx]

        # conditionally add the cross section definition (for pipes and culverts only)
        try:
            cross_id = objs.cross_section_definition_id[segment_idx]
            cross_weight = 1.0
        except AttributeError:
            cross_id = -9999
            cross_weight = np.nan

        # conditionally add the invert levels (for pipes and culverts only)
        try:
            invert_start = self.compute_bottom_level(content_pk, start_s)
            invert_end = self.compute_bottom_level(content_pk, end_s)
            dpumax = np.maximum(invert_start, invert_end)
        except AttributeError:
            invert_end = invert_start = dpumax = np.nan

        # conditionally add friction type and value (for pipes and culverts only)
        try:
            frict_type = objs.friction_type[segment_idx]
            frict_value = objs.friction_value[segment_idx]
        except AttributeError:
            frict_type = -9999
            frict_value = np.nan

        # Conditionally add discharge coefficients, (for culverts only). If one culvert has
        # multiple segments positive coefficient goes onto the first segment and negative
        # coefficient goes onto last segment (otherwise we have to much energy losses.)
        try:
            dc_positive = np.full((len(segments)), 1.0, dtype=np.float64)
            dc_negative = np.full((len(segments)), 1.0, dtype=np.float64)
            dc_positive[first_idx] = objs.discharge_coefficient_positive
            dc_negative[last_idx] = objs.discharge_coefficient_negative
        except AttributeError:
            dc_positive = 1.0
            dc_negative = 1.0

        # Conditionally add material (pipes only)
        try:
            material = objs.material[segment_idx]
        except AttributeError:
            material = -9999

        # Conditionally add sewerage_type (pipes only)
        try:
            sewerage_type = objs.sewerage_type[segment_idx]
        except AttributeError:
            sewerage_type = -9999

        # construct the result
        return Lines(
            id=itertools.islice(line_id_counter, len(segments)),
            line_geometries=segments,
            line=line,
            content_type=objs.content_type,
            content_pk=content_pk,
            s1d=(start_s + end_s) / 2,
            ds1d_half=(end_s - start_s) / 2,
            ds1d=end_s - start_s,
            kcu=objs.calculation_type[segment_idx],
            cross_id1=cross_id,
            cross_id2=cross_id,
            cross_weight=cross_weight,
            frict_type1=frict_type,
            frict_value1=frict_value,
            frict_type2=frict_type,
            frict_value2=frict_value,
            invert_level_start_point=invert_start,
            invert_level_end_point=invert_end,
            dpumax=dpumax,
            discharge_coefficient_positive=dc_positive,
            discharge_coefficient_negative=dc_negative,
            display_name=display_name,
            zoom_category=zoom_category,
            connection_node_start_id=connection_node_start_id,
            connection_node_end_id=connection_node_end_id,
            dist_calc_points=dist_calc_points,
            material=material,
            sewerage_type=sewerage_type,
        )

    def compute_bottom_level(self, ids, s):
        """Compute the bottom level by interpolating between invert levels

        This function is to be used on nodes and lines between on pipes and culverts.

        Args:
            ids (ndarray of int): the id (content_pk) of the object
            s (ndarray of float): the position of the node/line along the object

        Returns:
            an array of the same shape as ids and ds containing the interpolated values
        """
        if pygeos.is_missing(self.the_geom).any():
            raise ValueError(
                f"{self.__class__.__name__} found without a geometry. Call "
                f"set_geometries first."
            )
        lengths = pygeos.length(self.the_geom)
        idx = self.id_to_index(ids)
        weights = s / lengths[idx]
        if np.any(weights < 0.0) or np.any(weights > 1.0):
            raise ValueError("Encountered nodes outside of the linear object bounds")

        left = self.invert_level_start_point[idx]
        right = self.invert_level_end_point[idx]

        return weights * right + (1 - weights) * left

    def compute_drain_level(self, ids, s, connection_nodes):
        """Compute the drain level by interpolating between Manhole drain levels

        This function is to be used for interpolated nodes on pipes and culverts.

        Args:
            ids (ndarray of int): the id (content_pk) of the object
            s (ndarray of float): the position of the node measured along the object
            connection_nodes (ConnectionNodes): to obtain drain levels

        Returns:
            an array of the same shape as ids and ds containing the interpolated values
        """
        if pygeos.is_missing(self.the_geom).any():
            raise ValueError(
                f"{self.__class__.__name__} found without a geometry. Call "
                f"set_geometries first."
            )
        lengths = pygeos.length(self.the_geom)
        idx = self.id_to_index(ids)
        weights = s / lengths[idx]
        if np.any(weights < 0.0) or np.any(weights > 1.0):
            raise ValueError("Encountered nodes outside of the linear object bounds")

        # convert connection_node_start_id to index to get the drain level
        left_cn_idx = connection_nodes.id_to_index(self.connection_node_start_id)[idx]
        right_cn_idx = connection_nodes.id_to_index(self.connection_node_end_id)[idx]

        left = connection_nodes.drain_level[left_cn_idx]
        right = connection_nodes.drain_level[right_cn_idx]

        result = weights * right + (1 - weights) * left

        # fix cases with nan on one side
        result[np.isnan(result)] = left[np.isnan(result)]
        result[np.isnan(result)] = right[np.isnan(result)]
        return result


def counts_to_ranges(counts):
    """Convert an array of list-of-lists counts to ranges.

    ``counts`` define lengths of lists that are concatenated into a 1D array.
    The output ranges index into that array. For example, ``arr[start[i]:stop[i]]`` will
    give all values belonging to list `i`.

    Args:
        counts (ndarray of int): list lengths

    Returns:
        tuple of start (ndarray of int), stop (ndarray of int)

    Example:
    >>> counts_to_ranges([3, 2, 0, 1, 0])
    (array([0, 3, 5, 5, 6]), array([3, 5, 5, 6, 6]))
    """
    if len(counts) == 0:
        return np.empty((0,), dtype=int), np.empty((0,), dtype=int)
    stop = np.cumsum(counts)
    start = np.roll(stop, 1)
    start[0] = 0
    return start, stop


def counts_to_row_index(counts):
    """Convert an array of list-of-lists counts into row indices into a 2D array.

    ``counts`` define lengths of lists that are concatenated into a 1D array. The output
    of this function assigns a row index to every element in the 1D array. The row index
    is counting over the outer list.

    Args:
        counts (ndarray of int): list lengths

    Returns:
        row index (ndarray of int)

    Example:
    >>> counts_to_row_index([3, 2, 0, 1, 0])
    array([0, 0, 0, 1, 1, 3])
    """
    (n,) = counts.shape
    return np.repeat(np.arange(n), counts)


def counts_to_column_index(counts):
    """Convert an array of list-of-lists counts into row indices into a 2D array.

    ``counts`` define lengths of lists that are concatenated into a 1D array. The output
    of this function assigns a column index to every element in the 1D array. The column
    index is counting over elements in a single list.

    Args:
        counts (ndarray of int): list lengths

    Returns:
        column index (ndarray of int)

    Example:
    >>> counts_to_indices([3, 2, 0, 1, 0])
    array([0, 1, 2, 0, 1, 0])
    """
    if len(counts) == 0:
        return np.empty((0,), dtype=int)
    start, stop = counts_to_ranges(counts)
    return np.arange(stop[-1]) - np.repeat(start, counts)


def segmentize(linestrings, desired_segment_size):
    """Return points that divide linestrings into segments of equal length.

    Args:
        linestrings (ndarray of pygeos.Geometry): linestrings to segmentize
        desired_segment_size (ndarray of float): the desired size of the segments; the
           actual size will depend on the linestring length and is computed by rounding
           ``line length / size`` to the nearest integer.
           Inf inputs will lead to no segmentation for that line (n_segments=1)

    Returns:
        nodes: the points where segments connect.
          this excludes the start and end of the input linestrings.
        line_idx: indices mapping nodes to input linestrings
        dist_to_start: the location of the node measured along the linestring
    """
    # compute number of nodes to add per channel
    length = pygeos.length(linestrings)
    n_segments = np.maximum(np.round(length / desired_segment_size).astype(int), 1)
    segment_size = length / n_segments
    n_nodes = n_segments - 1

    # get the distance to the start of each channel
    i = counts_to_row_index(n_nodes)  # e.g. [0, 0, 0, 1, 1, 3]
    j = counts_to_column_index(n_nodes)  # e.g. [0, 1, 2, 0, 1, 0]
    dist_to_start = (j + 1) * segment_size[i]

    nodes = pygeos.line_interpolate_point(
        linestrings[i],
        dist_to_start,  # note: this only copies geometry pointers
    )
    return nodes, i, dist_to_start


def line_substring(linestrings, start, end, index=None):
    """Divide linestrings into segments given [start, end] measured along the line.

    Note that there is no bound check done on start and end. Expect bogus results or
    errors if start is negative or if end is larger than the linestring length.

    Args:
        linestrings (ndarray of pygeos.Geometry): Linestrings to segmentize
        start (ndarray of float): The start of the segment, measured along the line.
        end (ndarray of float): The end of the segment, measured along the line.
        index (ndarray of int, optional): An optional linestring index per start/end.
           If not given, then all input arrays should be equal sized.

    Returns:
        segments: the segments (sublinestrings of the input linestrings)
    """
    (n_lines,) = linestrings.shape
    if index is None:
        n_segments = n_lines
        index = np.arange(n_lines)
    else:
        (n_segments,) = index.shape

    if n_segments == 0:
        return np.empty((0,), dtype=object)

    coords, coord_line_idx = pygeos.get_coordinates(linestrings, return_index=True)
    line_n_coords = pygeos.get_num_coordinates(linestrings)

    ## compute the length of each vertex-vertex distance
    coord_ds = np.sqrt(np.sum((coords - np.roll(coords, 1, axis=0)) ** 2, axis=1))
    coord_s = np.cumsum(coord_ds)
    # compute corresponding length of each given start and end
    line_first_coord_idx, line_end_coord_idx = counts_to_ranges(line_n_coords)
    start_s = start + coord_s[line_first_coord_idx][index]
    end_s = end + coord_s[line_first_coord_idx][index]

    ## interpolate the start & end points along the linestrings
    start_insert_at = np.searchsorted(coord_s, start_s - COORD_EQUAL_ATOL)
    end_insert_at = np.searchsorted(coord_s, end_s - COORD_EQUAL_ATOL)
    # find points that would be located exactly at an existing coordinate
    start_eq_coord = np.abs(coord_s[start_insert_at] - start_s) < COORD_EQUAL_ATOL
    end_eq_coord = np.abs(coord_s[end_insert_at] - end_s) < COORD_EQUAL_ATOL
    # insert the existing coords and interpolate the others
    start_coords = np.empty((n_segments, 2))
    end_coords = np.empty((n_segments, 2))
    start_coords[start_eq_coord] = coords[start_insert_at[start_eq_coord]]
    end_coords[end_eq_coord] = coords[end_insert_at[end_eq_coord]]
    for i in range(2):
        start_coords[~start_eq_coord, i] = np.interp(
            start_s[~start_eq_coord], coord_s, coords[:, i]
        )
        end_coords[~end_eq_coord, i] = np.interp(
            end_s[~end_eq_coord], coord_s, coords[:, i]
        )

    ## fill an array of coordinates for the segmetns
    # the number of coordinates per segment are:
    # + 2 (line start and end)
    # + coordinates in between (end_insert_at - start_insert_at)
    # + correction where start_eq_coord (then there is 1 less "in between")
    # + no correction where end_eq_coord
    segment_n_coords = 2 + end_insert_at - start_insert_at - start_eq_coord
    segment_first_coord_idx, segment_end_coord_idx = counts_to_ranges(segment_n_coords)

    # fill in the start & end coords
    segment_coords = np.empty((segment_n_coords.sum(), 2))
    segment_coords[segment_first_coord_idx] = start_coords
    segment_coords[segment_end_coord_idx - 1] = end_coords

    # fill the remaining segment coordinates with original coordinates
    i = counts_to_row_index(segment_n_coords - 2)  # e.g. [0, 0, 0, 1, 1, 3]
    j = counts_to_column_index(segment_n_coords - 2)  # e.g. [0, 1, 2, 0, 1, 0]
    segment_coords_to_fill = segment_first_coord_idx[i] + j + 1
    coords_to_fill_with = start_insert_at[i] + start_eq_coord[i] + j
    segment_coords[segment_coords_to_fill] = coords[coords_to_fill_with]

    # construct the segments
    segments = pygeos.linestrings(
        segment_coords,
        indices=counts_to_row_index(segment_n_coords),
    )

    return segments


def segment_start_end(linestrings, segment_counts, dist_to_start):
    """Utility function to use the output of segmentize as input of line_substrings.

    Args:
        linestrings (ndarray of pygeos.Geometry): linestrings to segmentize
        segment_counts (ndarray of int): the number of segments per linestring
        dist_to_start (ndarray of float): the location of added nodes, measured along
          the linestring. The length of this array is such that:
          ``len(linestrings) + len(dist_to_start) = segment_counts.sum()``

    Returns:
        tuple of:
        - segment_start (ndarray of float): The segment start, measured along the line.
        - segment_end (ndarray of float): The segment end, measured along the line.
        - segment_index (ndarray of int): Indices mapping segments to input linestrings.
    """
    start_idx, last_idx = counts_to_ranges(segment_counts)
    last_idx -= 1

    start_s = np.full(segment_counts.sum(), np.nan)
    end_s = np.full(segment_counts.sum(), np.nan)
    # every first start_s is 0.0, the rest are added nodes with an s1d
    start_s[start_idx] = 0.0
    start_s[np.isnan(start_s)] = dist_to_start
    # every last end_s equals to the length, the rest are added nodes with an s1d
    end_s[last_idx] = pygeos.length(linestrings)
    end_s[np.isnan(end_s)] = dist_to_start
    return start_s, end_s, counts_to_row_index(segment_counts)


def _add_point(linestrings, points, at=0):
    """Append or prepend points to linestrings"""
    counts = pygeos.get_num_coordinates(linestrings)
    coords, index = pygeos.get_coordinates(linestrings, return_index=True)
    start, end = counts_to_ranges(counts)
    if at == 0:
        insert_at = start
    elif at == -1:
        insert_at = end
    else:
        raise ValueError(f"Cannot insert at {at}")
    new_coords = np.insert(coords, insert_at, pygeos.get_coordinates(points), axis=0)
    new_index = np.insert(index, insert_at, np.arange(len(insert_at)))
    return pygeos.linestrings(new_coords, indices=new_index)


def prepend_point(linestrings, points):
    return _add_point(linestrings, points, at=0)


def append_point(linestrings, points):
    return _add_point(linestrings, points, at=-1)
