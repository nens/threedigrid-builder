import itertools
from typing import Iterator

import numpy as np
import shapely

from threedigrid_builder.base import Lines, Nodes
from threedigrid_builder.constants import LineType, NodeType
from threedigrid_builder.exceptions import SchematisationError


def get_nodes(nodes, node_id_counter):
    """Add groundwater nodes from open water nodes."""
    open_water_nodes = nodes[np.isin(nodes.node_type, NodeType.NODE_2D_OPEN_WATER)]
    if len(open_water_nodes) == 0:
        raise SchematisationError(
            "Unable to create groundwater nodes, no open water nodes found."
        )
    id_n = itertools.islice(node_id_counter, len(open_water_nodes))
    groundwater_nodes = Nodes(
        id=id_n,
        node_type=NodeType.NODE_2D_GROUNDWATER,
        bounds=open_water_nodes.bounds,
        pixel_coords=open_water_nodes.pixel_coords,
        nodk=open_water_nodes.nodk,
        nodm=open_water_nodes.nodm,
        nodn=open_water_nodes.nodn,
        coordinates=open_water_nodes.coordinates,
    )
    node_id_offset = groundwater_nodes.id[0] - open_water_nodes.id[0]
    return groundwater_nodes, node_id_offset


def get_vertical_lines(nodes, node_id_offset, line_id_counter):
    """Add vertical lines between open water - and groundwater nodes."""
    open_water_nodes = nodes[np.isin(nodes.node_type, NodeType.NODE_2D_OPEN_WATER)]

    id_n = itertools.islice(line_id_counter, len(open_water_nodes))
    return Lines(
        id=id_n,
        kcu=LineType.LINE_2D_VERTICAL,
        line=np.stack(
            (open_water_nodes.id, open_water_nodes.id + node_id_offset), axis=1
        ),
        lik=open_water_nodes.nodk,
        lim=open_water_nodes.nodm,
        lin=open_water_nodes.nodn,
    )


def get_lines(lines, node_id_offset, line_id_counter):
    """Add groundwater lines from open water - and obstacle lines."""
    open_water_lines = lines[
        np.isin(
            lines.kcu,
            (
                LineType.LINE_2D_U,
                LineType.LINE_2D_V,
                LineType.LINE_2D_OBSTACLE_U,
                LineType.LINE_2D_OBSTACLE_V,
            ),
        )
    ]
    id_n = itertools.islice(line_id_counter, len(open_water_lines))
    return Lines(
        id=id_n,
        kcu=LineType.LINE_2D_GROUNDWATER,
        line=open_water_lines.line + node_id_offset,
        lik=open_water_lines.lik,
        lim=open_water_lines.lim,
        lin=open_water_lines.lin,
        line_coords=open_water_lines.line_coords,
        cross_pix_coords=open_water_lines.cross_pix_coords,
    )


class Lines1D2DGroundwater(Lines):
    @classmethod
    def create(
        cls, nodes: Nodes, line_id_counter: Iterator[int]
    ) -> "Lines1D2DGroundwater":
        """Create the 1D-2D groundwater lines

        Sets: id, line[:, 1], content_type, content_pk
        """
        node_idx = np.where(np.isfinite(nodes.groundwater_exchange[:, 0]))[0]
        line = np.full((len(node_idx), 2), -9999, dtype=np.int32)
        line[:, 1] = nodes.index_to_id(node_idx)
        return cls(
            id=itertools.islice(line_id_counter, len(node_idx)),
            line=line,
            content_type=nodes.content_type[node_idx],
            content_pk=nodes.content_pk[node_idx],
            groundwater_exchange=nodes.groundwater_exchange[node_idx],
        )

    def assign_2d_node(self, nodes: Nodes, cell_tree: shapely.STRtree) -> None:
        """Assigns the 2D node id based on the node coordinate.

        Requires: line[:, 1]
        Sets: line[:, 0]
        """
        node_idx = nodes.id_to_index(self.line[:, 1])
        side_1d = shapely.points(nodes.coordinates[node_idx])

        # The query returns 2 1D arrays: one with indices into the supplied node
        # geometries and one with indices into the tree of cells.
        idx = cell_tree.query(side_1d)
        # Address edge cases of multiple 1D-2D lines per node: just take the one
        _, unique_matches = np.unique(idx[0], return_index=True)
        line_idx, cell_idx = idx[:, unique_matches]
        self.line[line_idx, 0] = cell_idx + nodes.n_groundwater_cells
