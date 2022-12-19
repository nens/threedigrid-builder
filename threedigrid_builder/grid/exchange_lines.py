import itertools
from typing import Iterator

import numpy as np
import pygeos

from threedigrid_builder.base import Array, Lines, Nodes, search
from threedigrid_builder.constants import CalculationType, ContentType, LineType

from .levees import PotentialBreaches
from .obstacles import Obstacles

__all__ = ["ExchangeLines", "Lines1D2D"]


class ExchangeLine:
    id: int
    the_geom: pygeos.Geometry
    channel_id: int
    exchange_level: float


class ExchangeLines(Array[ExchangeLine]):
    pass


class Lines1D2D(Lines):
    @classmethod
    def create(cls, nodes: Nodes, line_id_counter: Iterator[int]) -> "Lines1D2D":
        """Create the 1D-2D lines, having only the 1D node index"""
        is_single = nodes.calculation_type == CalculationType.CONNECTED
        is_double = nodes.calculation_type == CalculationType.DOUBLE_CONNECTED
        node_idx = np.concatenate(
            [np.where(is_single | is_double)[0], np.where(is_double)[0]]
        )
        node_idx.sort()

        # Merge the arrays
        line = np.full((len(node_idx), 2), -9999, dtype=np.int32)
        line[:, 1] = nodes.index_to_id(node_idx)
        return Lines1D2D(id=itertools.islice(line_id_counter, len(node_idx)), line=line)

    @property
    def side_2d(self):
        return pygeos.points(self.line_coords[:, :2])

    @property
    def side_1d(self):
        return pygeos.points(self.line_coords[:, 2:])

    def assign_exchange_lines(
        self, nodes: Nodes, exchange_lines: ExchangeLines
    ) -> None:
        """Assign exchange lines to the 1D-2D lines.

        Requires: line[:, 1]
        Sets: content_pk, content_type
        """
        for (line_idx, exchange_lines_idx) in zip(
            self.split_in_two(self.line[:, 1]),
            exchange_lines.split_in_two(exchange_lines.channel_id),
        ):
            if len(line_idx) == 0 or len(exchange_lines_idx) == 0:
                continue

            node_idx = nodes.id_to_index(self.line[line_idx, 1])

            is_channel = nodes.content_type[node_idx] == ContentType.TYPE_V2_CHANNEL
            idx = search(
                exchange_lines.channel_id,
                nodes.content_pk[node_idx[is_channel]],
                assume_ordered=True,  # split_in_two orders
                check_exists=False,
                mask=exchange_lines_idx,
            )
            self.content_pk[line_idx[is_channel]] = exchange_lines.index_to_id(idx)
            self.content_type[line_idx[is_channel]] = np.where(
                idx != -9999, ContentType.TYPE_V2_EXCHANGE_LINE, -9999
            )

    def get_1d_node_idx(self, nodes: Nodes):
        """Get the 1D node index based on line[:, 1]"""
        return nodes.id_to_index(self.line[:, 1])

    def assign_2d_side(self, nodes: Nodes, exchange_lines: ExchangeLines):
        """Compute a Point on the (user-requested) 2D side for each line

        If there is no exchange line this is the 1D node location. If there is
        an exchange line, it is the closest point on the exchange line.

        Returns an array of pygeos Point geometries
        Requires: line[:, 1], content_pk, content_type
        Sets: line_coords[:, :2] (the 2D side)
        """
        # Create an array of node indices & corresponding exchange line indices
        has_exc = self.content_type == ContentType.TYPE_V2_EXCHANGE_LINE

        # Collect the 1D sides of the 1D2D line
        coords_1d = nodes.coordinates[self.get_1d_node_idx(nodes)]

        # Find the closest points on the exchange lines
        exc_geoms = exchange_lines.the_geom[
            exchange_lines.id_to_index(self.content_pk[has_exc])
        ]
        pygeos.prepare(exc_geoms)
        line_geom = pygeos.shortest_line(exc_geoms, pygeos.points(coords_1d[has_exc]))

        self.line_coords[~has_exc, :2] = coords_1d[~has_exc]
        self.line_coords[has_exc, :2] = pygeos.get_coordinates(
            pygeos.get_point(line_geom, 0)
        )

    def assign_breaches(self, nodes: Nodes, potential_breaches: PotentialBreaches):
        """Assign breaches to the 1D-2D lines.

        Nodes may have 0, 1, or 2 'breach_ids' assigned, but this assignment did not
        take into account the calculation type (isolated/connected/double connected).
        Here, the nodes.breach_ids are matched against the actually present 1D-2D lines.

        The matching minimizes the total distance of the 2D ends of the potential breach
        to the 2D end of the 1D-2D line (which is possibly derived from an exchange line).

        Requires: line_coords[:, :2]  (from the exchange line or 1D node location)
        Sets: line_coords[:, :2], content_pk, content_type  (updated from the potential breach)
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
            return pygeos.distance(breach_2d_side[breach_idx], line_2d_side[line_idx])

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
                # 1 line, 2 breaches; find closest breach
                if dist(b1, l1) <= dist(b2, l1):
                    line_to_breach[l1] = b1
                else:
                    line_to_breach[l1] = b2
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

        self.line_coords[line_idx, :2] = pygeos.get_coordinates(
            breach_2d_side[breach_idx]
        )
        self.content_type[line_idx] = ContentType.TYPE_V2_BREACH
        self.content_pk[line_idx] = potential_breaches.id[breach_idx]

    def assign_2d_node(self, cell_tree: pygeos.STRtree) -> None:
        """Assigns the 2D node id based on the line_coords

        Requires: line_coords[:, :2] (the 2D side)
        Sets: line[:, 0], line_coords[:, :2] (clears)
        """
        # The query_bulk returns 2 1D arrays: one with indices into the supplied node
        # geometries and one with indices into the tree of cells.
        idx = cell_tree.query_bulk(self.side_2d)
        # Address edge cases of multiple 1D-2D lines per node: just take the one
        _, unique_matches = np.unique(idx[0], return_index=True)
        line_idx, cell_idx = idx[:, unique_matches]
        self.line[line_idx, 0] = cell_idx
        self.line_coords[:, :2] = np.nan

    def assign_kcu(self, mask, is_closed) -> None:
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

    def assign_dpumax_from_obstacles(self, obstacles: Obstacles):
        """Set the dpumax based on intersected obstacles

        Only for (open) connections that are computed via an exchange line.
        """
        has_exc = self.content_type == ContentType.TYPE_V2_EXCHANGE_LINE
        self.assign_dpumax(
            has_exc, obstacles.compute_dpumax(self, where=np.where(has_exc)[0])[0]
        )

    def assign_dpumax(self, mask, dpumax) -> None:
        """Assign dpumax only where it is not set already"""
        is_nan = np.isnan(self.dpumax)
        self.dpumax[mask & is_nan] = dpumax[is_nan[mask]]

    def assign_ds1d(self, nodes: Nodes) -> None:
        """Sets the length (ds1d) based on the 2D cell with

        Requires: line[:, 0]
        Sets: ds1d, ds1d_half
        """
        cell_idx = nodes.id_to_index(self.line[:, 0])
        self.ds1d[:] = nodes.bounds[cell_idx, 2] - nodes.bounds[cell_idx, 0]
        self.ds1d_half[:] = 0.5 * self.ds1d
