import itertools
import logging
from typing import Iterator, List

import numpy as np
import shapely

from threedigrid_builder.base import LineHalfs, Lines, Nodes, replace, search
from threedigrid_builder.constants import CalculationType, ContentType, LineType
from threedigrid_builder.exceptions import SchematisationError

from .connection_nodes import ConnectionNodes
from .exchange_lines import ExchangeLines
from .obstacles import ObstacleAffectsType, Obstacles
from .potential_breaches import PotentialBreaches

logger = logging.getLogger(__name__)

__all__ = ["Lines1D2D"]


class Lines1D2D(Lines):
    @classmethod
    def create(cls, nodes: Nodes, line_id_counter: Iterator[int]) -> "Lines1D2D":
        """Create the 1D-2D lines

        Sets: id, line[:, 1], content_type, content_pk
        """
        is_single = nodes.calculation_type == CalculationType.CONNECTED
        is_double = nodes.calculation_type == CalculationType.DOUBLE_CONNECTED
        node_idx = np.concatenate(
            [np.where(is_single | is_double)[0], np.where(is_double)[0]]
        )
        node_idx.sort()

        # Merge the arrays
        line = np.full((len(node_idx), 2), -9999, dtype=np.int32)
        line[:, 1] = nodes.index_to_id(node_idx)
        return cls(
            id=itertools.islice(line_id_counter, len(node_idx)),
            line=line,
            content_type=nodes.content_type[node_idx],
            content_pk=nodes.content_pk[node_idx],
        )

    @property
    def side_2d(self):
        return shapely.points(self.line_coords[:, :2])

    @property
    def side_1d(self):
        return shapely.points(self.line_coords[:, 2:])

    def assign_connection_nodes_to_channels_from_breaches(
        self, nodes: Nodes, potential_breaches: PotentialBreaches
    ) -> None:
        """Replace the content_pk of connection nodes using the nodes.breach_ids.

        This is required for exchange line assignment.

        Requires: content_type, line[:, 1]
        Sets: content_pk, content_type
        """
        is_connection_node = np.where(
            self.content_type == ContentType.TYPE_V2_CONNECTION_NODES
        )[0]
        node_ids = self.line[is_connection_node, 1]
        breach_id = nodes.breach_ids[nodes.id_to_index(node_ids), 0]
        has_breach = breach_id != -9999
        channel_id = potential_breaches.channel_id[
            potential_breaches.id_to_index(breach_id[has_breach])
        ]
        self.content_pk[is_connection_node[has_breach]] = channel_id
        self.content_type[is_connection_node[has_breach]] = ContentType.TYPE_V2_CHANNEL

    def assign_connection_nodes_to_channels_from_lines(
        self, nodes: Nodes, lines: Lines
    ):
        """Replace the content_pk of connection nodes on connected objects.

        The connection node will be assigned to the first 'connected' channel
        that is attached to it. Double connected channels get priority.

        This is required for exchange line assignment.

        Requires: content_type, line[:, 1]
        Sets: content_pk, content_type
        """
        PRIORITY = {LineType.LINE_1D_DOUBLE_CONNECTED: 0, LineType.LINE_1D_CONNECTED: 1}

        is_connection_node = np.where(
            self.content_type == ContentType.TYPE_V2_CONNECTION_NODES
        )[0]
        node_ids = self.line[is_connection_node, 1]

        line_halfs = LineHalfs.for_connection_nodes(
            nodes,
            lines,
            line_mask=(
                (lines.content_type == ContentType.TYPE_V2_CHANNEL)
                & (np.isin(lines.kcu, list(PRIORITY.keys())))
            ),
            node_mask=np.isin(nodes.id, node_ids),
        )
        # Order the line_halfs by (channel_id, calc type) and pick the first
        line_halfs.reorder(
            np.lexsort(
                [
                    line_halfs.content_pk,
                    replace(line_halfs.kcu, PRIORITY),
                    line_halfs.node_id,
                ]
            )
        )
        channel_id_per_node = line_halfs.first(line_halfs.content_pk)
        # Use the node_id -> channel_id map to set channel_id on self
        idx = search(
            line_halfs.get_reduce_id(),
            node_ids,
            check_exists=False,
            assume_ordered=True,
        )
        mask = idx != -9999

        self.content_pk[is_connection_node[mask]] = channel_id_per_node[idx[mask]]
        self.content_type[is_connection_node[mask]] = ContentType.TYPE_V2_CHANNEL

    def assign_exchange_lines(self, exchange_lines: ExchangeLines) -> None:
        """Assign exchange lines to the 1D-2D lines.

        Resets all content_type / content_pk.

        Requires: content_pk, content_type
        Sets: content_pk, content_type
        """
        line_idx_1, line_idx_2 = self.split_in_two(self.line[:, 1])
        self._assign_exchange_lines(line_idx_1, exchange_lines, is_primary=True)
        self._assign_exchange_lines(line_idx_2, exchange_lines, is_primary=False)

        # reset content_pk / content_type where not set to exchange line
        has_no_exc = self.content_type != ContentType.TYPE_V2_EXCHANGE_LINE
        self.content_type[has_no_exc] = -9999
        self.content_pk[has_no_exc] = -9999

    def _assign_exchange_lines(
        self, idx, exchange_lines: ExchangeLines, is_primary: bool
    ) -> None:
        if len(idx) == 0:
            return
        is_channel = self.content_type[idx] == ContentType.TYPE_V2_CHANNEL
        channel_id = self.content_pk[idx[is_channel]]
        content_pk = exchange_lines.get_for_channel_id(
            channel_id, is_primary=is_primary
        )
        self.content_pk[idx[is_channel]] = content_pk
        self.content_type[idx[is_channel]] = np.where(
            content_pk != -9999, ContentType.TYPE_V2_EXCHANGE_LINE, -9999
        )

    def get_1d_node_idx(self, nodes: Nodes):
        """Get the 1D node index based on line[:, 1]"""
        return nodes.id_to_index(self.line[:, 1])

    def assign_line_coords(self, nodes: Nodes):
        """Copy the 1D node location to the 2D side for non-exchange line lines.

        Returns an array of shapely Point geometries
        Requires: line[:, 1], content_pk, content_type
        Sets: line_coords (side_1d and side_2d)
        """
        self.line_coords[:, 2:] = nodes.coordinates[self.get_1d_node_idx(nodes)]
        self.line_coords[:, :2] = self.line_coords[:, 2:]

    def assign_2d_side_from_exchange_lines(self, exchange_lines: ExchangeLines):
        """Compute the closest Point on exchange line for each 1D2D line

        Requires: line_coords[:, 2:] (side_1d), content_pk, content_type
        Sets: line_coords[:, :2] (side_2d)
        """
        # Create an array of node indices & corresponding exchange line indices
        has_exc = self.content_type == ContentType.TYPE_V2_EXCHANGE_LINE

        # Find the closest points on the exchange lines
        exc_geoms = exchange_lines.the_geom[
            exchange_lines.id_to_index(self.content_pk[has_exc])
        ]
        shapely.prepare(exc_geoms)
        line_geom = shapely.shortest_line(exc_geoms, self.side_1d[has_exc])

        self.line_coords[has_exc, :2] = shapely.get_coordinates(
            shapely.get_point(line_geom, 0)
        )

    def assign_breaches(self, nodes: Nodes, potential_breaches: PotentialBreaches):
        """Assign breaches to the 1D-2D lines.

        It is assumed that the breach ids match the calculation type of the nodes.
        This is done in PotentialBreachPoints.match_breach_ids_with_calculation_types.

        If a node is double connected, the 1 or 2 breaches that are present need to be
        matched against existing lines. The matching minimizes the total distance of
        the 2D ends of the potential breach to the 2D end of the 1D-2D line (which
        is possibly derived from an exchange line).

        Requires: line_coords[:, :2]  (from the exchange line or 1D node location)
        Sets: line_coords[:, :2], line_geometries, content_pk, content_type
        """
        node_idx = self.get_1d_node_idx(nodes)
        breach_counts = np.count_nonzero(nodes.breach_ids[node_idx, :] != -9999, axis=1)
        line_idx_with_breaches = np.where(breach_counts > 0)[0]

        if len(line_idx_with_breaches) == 0:
            return

        # define a function that computes the distance between breach and line 2D sides
        breach_2d_side = potential_breaches.side_2d
        line_2d_side = self.side_2d

        def dist(breach_idx, line_idx):
            return shapely.distance(breach_2d_side[breach_idx], line_2d_side[line_idx])

        # a mapping with a breach_id for each line
        line_to_breach = {}

        # loop over lines with breaches (i is a line index)
        for i in np.where(breach_counts > 0)[0]:
            if i in line_to_breach:
                continue  # happens for double connected nodes

            # find all lines that are connected to the node
            line_idx = np.where(self.line[:, 1] == self.line[i, 1])[0]
            if len(line_idx) == 1:
                l1, l2 = line_idx[0], -9999
            else:
                l1, l2 = line_idx
            b1, b2 = potential_breaches.id_to_index(nodes.breach_ids[node_idx[i], :])

            if l2 == -9999 and b2 == -9999:
                # 1 line, 1 breach
                line_to_breach[l1] = b1
            elif l2 == -9999 and b2 != -9999:
                raise ValueError(
                    f"Node {nodes.id[node_idx[i]]} has two breaches assigned while it is not double connected"
                )
            elif b2 == -9999:
                # 2 lines, 1 breach; find closest line
                if dist(b1, l1) <= dist(b1, l2):
                    line_to_breach[l1] = b1
                else:
                    line_to_breach[l2] = b1
            elif b2 != -9999:
                # 2 lines, 2 breach; find the combination that minimizes total distance
                if (dist(b1, l1) + dist(b2, l2)) <= (dist(b1, l2) + dist(b2, l1)):
                    line_to_breach[l1], line_to_breach[l2] = b1, b2
                else:
                    line_to_breach[l1], line_to_breach[l2] = b2, b1

        line_idx, breach_idx = np.array(list(line_to_breach.items()), dtype=np.int32).T

        self.line_coords[line_idx, :2] = shapely.get_coordinates(
            breach_2d_side[breach_idx]
        )
        self.line_geometries[line_idx] = potential_breaches.the_geom[breach_idx]
        self.content_type[line_idx] = ContentType.TYPE_V2_BREACH
        self.content_pk[line_idx] = potential_breaches.id[breach_idx]

    def assign_2d_node(self, cell_tree: shapely.STRtree) -> None:
        """Assigns the 2D node id based on the line_coords

        Requires: line_coords[:, :2] (the 2D side)
        Sets: line[:, 0], line_coords[:] (clears)
        """
        # The query returns 2 1D arrays: one with indices into the supplied node
        # geometries and one with indices into the tree of cells.
        idx = cell_tree.query(self.side_2d)
        # Address edge cases of multiple 1D-2D lines per node: just take the one
        _, unique_matches = np.unique(idx[0], return_index=True)
        line_idx, cell_idx = idx[:, unique_matches]
        self.line[line_idx, 0] = cell_idx
        self.line_coords[:] = np.nan

    def check_unassigned(self, nodes, is_groundwater: bool = False) -> None:
        """Checks 1D-2D lines where any of the required nodes is set to null, represented as -9999
        This is the case when the nodes are outside the 2D domain.
        """
        invalid_rows = self.line[:, 0] == -9999
        if invalid_rows.any():
            invalid_node_ids = self.line[invalid_rows, 1]
            invalid_nodes = nodes.id_to_index(invalid_node_ids)
            invalid_nodes_formatted = nodes.format_message(invalid_nodes)

            if is_groundwater:
                raise SchematisationError(
                    "The following objects have groundwater exchange properties but "
                    f"are (partially) outside of the 2D model domain: {invalid_nodes_formatted}."
                )

            raise SchematisationError(
                "The following objects are connected but are (partially) outside of "
                f"the 2D model domain: {invalid_nodes_formatted}."
            )

    def transfer_2d_node_to_groundwater(self, offset: int):
        """Transfers the 1D-2D line to a groundwater node

        Sets: line[:, 0]
        """
        has_2d_node = np.where(self.line[:, 0] != -9999)[0]
        if offset == 0:
            # no groundwater cells; reset the node id (will raise error later)
            self.line[has_2d_node, 0] = -9999
        else:
            self.line[has_2d_node, 0] += offset

    def assign_kcu(self, mask, is_closed) -> None:
        """Set kcu where it is not set already"""
        node_id = self.line[mask, 1]
        node_id_unique, counts = np.unique(node_id, return_counts=True)
        line_is_double = np.isin(self.line[mask, 1], node_id_unique[counts == 2])
        self.kcu[mask] = np.choose(
            line_is_double * 2 + is_closed,
            choices=[
                LineType.LINE_1D2D_SINGLE_CONNECTED_OPEN_WATER,
                LineType.LINE_1D2D_SINGLE_CONNECTED_CLOSED,
                LineType.LINE_1D2D_DOUBLE_CONNECTED_OPEN_WATER,
                LineType.LINE_1D2D_DOUBLE_CONNECTED_CLOSED,
            ],
        )

    def assign_dpumax_from_breaches(self, potential_breaches: PotentialBreaches):
        """Set the dpumax based exchange_lines.exchange_level"""
        has_breach = self.content_type == ContentType.TYPE_V2_BREACH
        self.assign_dpumax(
            has_breach,
            potential_breaches.exchange_level[
                potential_breaches.id_to_index(self.content_pk[has_breach])
            ],
        )

    def assign_dpumax_from_exchange_lines(self, exchange_lines: ExchangeLines):
        """Set the dpumax based exchange_lines.exchange_level"""
        has_exc = self.content_type == ContentType.TYPE_V2_EXCHANGE_LINE
        self.assign_dpumax(
            has_exc,
            exchange_lines.exchange_level[
                exchange_lines.id_to_index(self.content_pk[has_exc])
            ],
        )

    def assign_dpumax_from_obstacles_open(self, obstacles: Obstacles):
        self._assign_dpumax_from_obstacles(
            obstacles,
            affects_type=ObstacleAffectsType.AFFECTS_1D2D_OPEN_WATER,
            line_types=[
                LineType.LINE_1D2D_SINGLE_CONNECTED_OPEN_WATER,
                LineType.LINE_1D2D_DOUBLE_CONNECTED_OPEN_WATER,
            ],
        )

    def assign_dpumax_from_obstacles_closed(self, obstacles: Obstacles):
        self._assign_dpumax_from_obstacles(
            obstacles,
            affects_type=ObstacleAffectsType.AFFECTS_1D2D_CLOSED,
            line_types=[
                LineType.LINE_1D2D_SINGLE_CONNECTED_CLOSED,
                LineType.LINE_1D2D_DOUBLE_CONNECTED_CLOSED,
            ],
        )

    def _assign_dpumax_from_obstacles(
        self,
        obstacles: Obstacles,
        affects_type: ObstacleAffectsType,
        line_types: List[LineType],
    ):
        """Set the dpumax based on intersected obstacles

        Only for open connections. Requires line_geometries.
        """
        is_included = np.isin(self.kcu, line_types)
        is_included_idx = np.where(is_included)[0]
        obstacle_crest_levels, obstacle_idx = obstacles.compute_dpumax(
            self, where=is_included_idx, affects_type=affects_type
        )
        self.assign_dpumax(is_included, obstacle_crest_levels)
        mask = obstacle_idx != -9999

        self.assign_ds1d_half_from_obstacles(
            obstacles, is_included_idx[mask], obstacle_idx[mask]
        )

    def assign_ds1d_half_from_obstacles(
        self, obstacles: Obstacles, line_idx, obstacle_idx
    ):
        if len(line_idx) == 0:
            return
        line_geoms = self.line_geometries[line_idx]
        # compute the intersections (use shortest_line and not intersects to
        # account for the case that the obstacle and the line intersect multiple times)
        points = shapely.get_point(
            shapely.shortest_line(obstacles.the_geom[obstacle_idx], line_geoms), 0
        )
        # compute the distance to the start of the line
        self.ds1d_half[line_idx] = shapely.line_locate_point(line_geoms, points)

    def assign_dpumax(self, mask, dpumax) -> None:
        """Assign dpumax only where it is not set already"""
        is_nan = np.isnan(self.dpumax)
        self.dpumax[mask & is_nan] = dpumax[is_nan[mask]]

    def assign_ds1d(self, nodes: Nodes) -> None:
        """Sets the length (ds1d) based on the 2D cell with, where not yet set.

        Requires: line[:, 0]
        Sets: ds1d
        """
        has_no_ds1d = np.isnan(self.ds1d)
        cell_idx = nodes.id_to_index(self.line[has_no_ds1d, 0])
        self.ds1d[has_no_ds1d] = nodes.bounds[cell_idx, 2] - nodes.bounds[cell_idx, 0]

    def assign_ds1d_half(self) -> None:
        """Sets the velocity point location on the line (ds1d_half), where not yet set.

        Requires: ds1d
        Sets: ds1d_half
        """
        has_no_ds1d_half = np.isnan(self.ds1d_half)
        self.ds1d_half[has_no_ds1d_half] = (
            shapely.length(self.line_geometries[has_no_ds1d_half]) / 2
        )

    def output_breaches(self, breaches: PotentialBreaches) -> PotentialBreaches:
        """Create the actual breaches based on the 1D-2D lines

        This results in 'named' breaches (which are present in the input potential_breaches).
        Other 1D-2D lines may be breached as well, but they have no metadata.
        """
        line_idx = np.where(self.content_type == ContentType.TYPE_V2_BREACH)[0]
        breach_idx = breaches.id_to_index(self.content_pk[line_idx])

        # Order by breach id
        sorter = np.argsort(breach_idx)
        line_idx = line_idx[sorter]
        breach_idx = breach_idx[sorter]

        # Filter by breaches that miss metadata
        mask = (breaches.levee_material[breach_idx] != -9999) & np.isfinite(
            breaches.maximum_breach_depth[breach_idx]
        )
        line_idx = line_idx[mask]
        breach_idx = breach_idx[mask]

        return PotentialBreaches(
            id=range(len(breach_idx)),
            line_id=self.id[line_idx],
            content_pk=breaches.id[breach_idx],
            maximum_breach_depth=breaches.maximum_breach_depth[breach_idx],
            levee_material=breaches.levee_material[breach_idx],
            the_geom=self.get_velocity_points(line_idx),
            code=breaches.code[breach_idx],
            display_name=breaches.display_name[breach_idx],
        )

    @classmethod
    def create_groundwater(
        cls, nodes: Nodes, line_id_counter: Iterator[int]
    ) -> "Lines1D2D":
        """Create a 1D-2D groundwater lines for each 1D node.

        Sets: id, line[:, 1], kcu
        """
        node_idx = np.where(nodes.has_groundwater_exchange)[0]
        line = np.full((len(node_idx), 2), -9999, dtype=np.int32)
        line[:, 1] = nodes.index_to_id(node_idx)
        return cls(
            id=itertools.islice(line_id_counter, len(node_idx)),
            line=line,
            kcu=LineType.LINE_1D2D_GROUNDWATER,
        )

    def assign_groundwater_exchange(self, nodes: Nodes, cn: ConnectionNodes):
        """Compute the hydraulic_conductivity values for groundwater exchange manholes
        when present. Set all dpumax values for 1d2d groundwater exchange lines.

        Sets:
        - hydraulic_conductivity_in for manholes
        - hydraulic_conductivity_out for manholes
        - dpumax (from 1D nodes dmax)
        """

        node_idx = self.get_1d_node_idx(nodes)
        idx = np.where(
            nodes.content_type[node_idx] == ContentType.TYPE_V2_CONNECTION_NODES
        )[0]
        cn_idx = cn.id_to_index(nodes.content_pk[node_idx[idx]])

        mask = cn.has_groundwater_exchange[cn_idx]
        idx = idx[mask]
        cn_idx = cn_idx[mask]

        self.hydraulic_resistance_in[idx] = (
            cn.hydraulic_conductivity_in[cn_idx] / cn.exchange_thickness[cn_idx]
        )
        self.hydraulic_resistance_out[idx] = (
            cn.hydraulic_conductivity_out[cn_idx] / cn.exchange_thickness[cn_idx]
        )

        self.dpumax = nodes[self.get_1d_node_idx(nodes)].dmax
