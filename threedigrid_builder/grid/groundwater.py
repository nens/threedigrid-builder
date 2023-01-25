import itertools

import numpy as np

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
