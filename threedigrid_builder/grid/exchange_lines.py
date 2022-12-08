import itertools
from typing import Iterator

import numpy as np
import pygeos

from threedigrid_builder.base import Array, Lines, Nodes, search
from threedigrid_builder.constants import CalculationType, ContentType, LineType

__all__ = ["ExchangeLines", "Lines1D2D"]


class ExchangeLine:
    id: int
    the_geom: pygeos.Geometry
    channel_id: int


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

    def compute_2d_side(self, nodes: Nodes, exchange_lines: ExchangeLines):
        """Compute a Point on the (user-requested) 2D side for each line

        If there is no exchange line this is the 1D node location. If there is
        an exchange line, it is the closest point on the exchange line.

        Returns an array of pygeos Point geometries
        Requires: line[:, 1], content_pk, content_type
        Sets: nothing
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

        result = np.empty(len(self), dtype=object)
        result[~has_exc] = pygeos.points(coords_1d[~has_exc])
        result[has_exc] = pygeos.get_point(line_geom, 0)
        return result

    def assign_2d_node(self, points_2d_side, cell_tree: pygeos.STRtree) -> None:
        """Assigns the 2D node id based on the line_coords

        Sets: line[:, 0]
        """
        # The query_bulk returns 2 1D arrays: one with indices into the supplied node
        # geometries and one with indices into the tree of cells.
        idx = cell_tree.query_bulk(points_2d_side)
        # Address edge cases of multiple 1D-2D lines per node: just take the one
        _, unique_matches = np.unique(idx[0], return_index=True)
        line_idx, cell_idx = idx[:, unique_matches]
        self.line[line_idx, 0] = cell_idx

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

    def assign_dpumax(self, mask, dpumax) -> None:
        self.dpumax[mask] = dpumax

    def assign_ds1d(self, nodes: Nodes) -> None:
        """Sets the length (ds1d) based on the 2D cell with

        Requires: line[:, 0]
        Sets: ds1d, ds1d_half
        """
        cell_idx = nodes.id_to_index(self.line[:, 0])
        self.ds1d[:] = nodes.bounds[cell_idx, 2] - nodes.bounds[cell_idx, 0]
        self.ds1d_half[:] = 0.5 * self.ds1d
