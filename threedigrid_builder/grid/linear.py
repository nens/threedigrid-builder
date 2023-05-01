import itertools
from typing import Optional

import numpy as np
import shapely

from threedigrid_builder.base import Lines, LineStrings, Nodes, PointsOnLine
from threedigrid_builder.base.linestrings import counts_to_ranges
from threedigrid_builder.constants import CalculationType, NodeType

COORD_EQUAL_ATOL = 1e-8  # the distance below which coordinates are considered equal


class BaseLinear:
    content_type = None  # to be defined by subclasses

    @property
    def linestrings(self):
        return LineStrings(objects=self)

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

        self.linestrings.sanitize(points_1, points_2)

        # construct the geometries where necessary
        has_no_geom = shapely.is_missing(self.the_geom)
        coordinates = np.empty((np.count_nonzero(has_no_geom), 2, 2))
        coordinates[:, 0, 0] = shapely.get_x(points_1[has_no_geom])
        coordinates[:, 0, 1] = shapely.get_y(points_1[has_no_geom])
        coordinates[:, 1, 0] = shapely.get_x(points_2[has_no_geom])
        coordinates[:, 1, 1] = shapely.get_y(points_2[has_no_geom])
        self.the_geom[has_no_geom] = shapely.linestrings(coordinates)

    def interpolate_nodes(
        self,
        node_id_counter,
        global_dist_calc_points,
        fixed_nodes: Optional[PointsOnLine] = None,
    ):
        """Compute nodes on each linear object with constant intervals

        The following fields are expected to be filled on self:
        - id
        - the_geom
        - dist_calc_points
        - calculation_type

        Args:
            node_id_counter (iterable): an iterable yielding integers
            global_dist_calc_points (float): Default node interdistance.
            fixed_nodes: optional fixed points that will become a node

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
        if shapely.is_missing(self.the_geom).any():
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

        if fixed_nodes is None:
            fixed_nodes = PointsOnLine.empty(self.linestrings)
        else:
            fixed_nodes = fixed_nodes[~(fixed_nodes.at_start | fixed_nodes.at_end)]

        # insert default dist_calc_points where necessary
        dists = self.dist_calc_points.copy()
        dists[np.isnan(dists)] = global_dist_calc_points
        dists[dists <= 0] = global_dist_calc_points
        # Do not add nodes for embedded objects:
        dists[self.calculation_type == CalculationType.EMBEDDED] = np.inf

        # interpolate the node geometries
        sublinestrings = self.linestrings.segmentize(fixed_nodes)
        points = sublinestrings.interpolate_points(dists[sublinestrings.linestring_idx])

        # fixed_nodes also become nodes
        points = points.merge_with(fixed_nodes)

        # construct the nodes with available attributes
        return Nodes(
            id=itertools.islice(node_id_counter, len(points)),
            coordinates=shapely.get_coordinates(points.the_geom),
            content_type=self.content_type,
            content_pk=self.id[points.linestring_idx],
            node_type=NodeType.NODE_1D_NO_STORAGE,
            calculation_type=self.calculation_type[points.linestring_idx],
            s1d=points.s1d,
            breach_ids=np.array([points.content_pk, points.secondary_content_pk]).T,
        )

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

        # Cut the linestrings up into segments (LinesOnLine)
        node_points = PointsOnLine.from_s1d(
            objs.linestrings,
            nodes.s1d,
            objs.id_to_index(nodes.content_pk),
            content_pk=nodes.content_pk,
        )
        segments = objs.linestrings.segmentize(node_points)
        segment_idx = segments.linestring_idx

        # Copy properties to segments
        display_name = np.take(objs.display_name, segment_idx)
        zoom_category = np.take(objs.zoom_category, segment_idx)
        connection_node_start_id = np.take(objs.connection_node_start_id, segment_idx)
        connection_node_end_id = np.take(objs.connection_node_end_id, segment_idx)
        dist_calc_points = np.take(objs.dist_calc_points, segment_idx)

        # set the right node indices for each segment
        first_idx, last_idx = counts_to_ranges(np.bincount(segments.linestring_idx))
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
            invert_start = self.compute_bottom_level(content_pk, segments.s1d_start)
            invert_end = self.compute_bottom_level(content_pk, segments.s1d_end)
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

        # conditionall add hydraulic conductivity (for channels and pipes only)
        try:
            exchange_thickness = objs.exchange_thickness
            hydraulic_resistance_in = (
                objs.hydraulic_conductivity_in[segment_idx]
                / exchange_thickness[segment_idx]
            )
            hydraulic_resistance_out = (
                objs.hydraulic_conductivity_out[segment_idx]
                / exchange_thickness[segment_idx]
            )
        except AttributeError:
            hydraulic_resistance_in = np.nan
            hydraulic_resistance_out = np.nan

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
            line_geometries=segments.the_geom,
            line=line,
            content_type=objs.content_type,
            content_pk=content_pk,
            s1d=segments.s1d,
            ds1d_half=segments.ds1d / 2,
            ds1d=segments.ds1d,
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
            hydraulic_resistance_in=hydraulic_resistance_in,
            hydraulic_resistance_out=hydraulic_resistance_out,
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
        if shapely.is_missing(self.the_geom).any():
            raise ValueError(
                f"{self.__class__.__name__} found without a geometry. Call "
                f"set_geometries first."
            )
        lengths = shapely.length(self.the_geom)
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
        if shapely.is_missing(self.the_geom).any():
            raise ValueError(
                f"{self.__class__.__name__} found without a geometry. Call "
                f"set_geometries first."
            )
        lengths = shapely.length(self.the_geom)
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

    @property
    def has_groundwater_exchange(self):
        return np.full(len(self), False, dtype=bool)

    def apply_has_groundwater_exchange(self, nodes: Nodes, lines: Lines):
        content_pk = self.id[self.has_groundwater_exchange]

        if len(content_pk) == 0:
            return

        node_ids = np.unique(
            lines.line[
                (lines.content_type == self.content_type)
                & np.isin(lines.content_pk, content_pk)
            ]
        )

        nodes.has_groundwater_exchange[
            np.isin(nodes.id, node_ids)
            & (nodes.calculation_type != CalculationType.BOUNDARY_NODE)
        ] = True
