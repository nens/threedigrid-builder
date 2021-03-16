from threedigrid_builder.base import Lines
from threedigrid_builder.base import Nodes
from threedigrid_builder.constants import CalculationType
from threedigrid_builder.constants import ContentType
from threedigrid_builder.constants import LineType
from threedigrid_builder.grid import cross_sections

import numpy as np


__all__ = ["Grid"]

# For assigning calculation types to connection nodes:
LINE_TYPE_PRECEDENCE = (
    LineType.LINE_1D_ISOLATED,
    LineType.LINE_1D_CONNECTED,
    LineType.LINE_1D_EMBEDDED,
)
LINE_TYPE_MAPPING = {
    LineType.LINE_1D_ISOLATED: CalculationType.ISOLATED,
    LineType.LINE_1D_CONNECTED: CalculationType.CONNECTED,
    LineType.LINE_1D_EMBEDDED: CalculationType.EMBEDDED,
}


class Grid:
    def __init__(self, nodes: Nodes, lines: Lines, quadtree_statistics=None):
        if not isinstance(nodes, Nodes):
            raise TypeError(f"Expected Nodes instance, got {type(nodes)}")
        if not isinstance(lines, Lines):
            raise TypeError(f"Expected Lines instance, got {type(lines)}")
        self.nodes = nodes
        self.lines = lines
        self.epsg_code = None  # Grid is aware of its epsg_code
        self.quadtree_statistics = quadtree_statistics

    def __add__(self, other):
        """Concatenate two grids without renumbering nodes."""
        if self.__class__ is not other.__class__:
            raise TypeError(
                "Cannot concatenate {self} with {other} as they are not of "
                "equal types."
            )
        return self.__class__(
            nodes=self.nodes + other.nodes, lines=self.lines + other.lines
        )

    def __repr__(self):
        return f"<Grid object with {len(self.nodes)} nodes and {len(self.lines)} lines>"

    @classmethod
    def from_quadtree(cls, quadtree, area_mask, node_id_counter, line_id_counter):
        """Construct the 2D grid based on the quadtree object.

        Args:
            quadtree (QuadTree): grid of refinement levels and active cell idx
              are used.
            area_mask (numpy.ndarray) : Array with raster of active pixels.
            node_id_counter (iterable): an iterable yielding integers
            line_id_counter (iterable): an iterable yielding integers

        Returns:
            Grid with data in the following columns:
            - nodes.id: ids generated based on counter
            - nodes.coordinates: node coordinates
            - nodes.bounds: node bounds
            - nodes.node_type: Node_type (initially NODE_2D_OPEN_WATER)
            - nodes.nodk: Grid refinement level of node
            - nodes.nodm: horizontal index of node at grid refinement level nodk
            - nodes.nodn: vertical index of node at grid refinement level nodk
            - lines.id: ids generated based on counter
            - lines.line: lines between connecting nodes.
            - lines.lik:  Grid refinement level of line (smallest at refinements.)
            - lines.lim:  horizontal index of line at grid refimenent level lik
            - lines.lin:  vertical index line at grid refimenent level lik
        """

        nodes, lines = quadtree.get_nodes_lines(
            area_mask, node_id_counter, line_id_counter
        )

        quadtree_statistics = {
            "lgrmin": quadtree.lgrmin,
            "kmax": quadtree.kmax,
            "mmax": quadtree.mmax,
            "nmax": quadtree.nmax,
            "dx": quadtree.dx,
            "bbox": quadtree.bbox,
        }

        return cls(nodes=nodes, lines=lines, quadtree_statistics=quadtree_statistics)

    @classmethod
    def from_connection_nodes(cls, connection_nodes, node_id_counter):
        """Construct a grid (only nodes) for the connection nodes

        Args:
            connection_nodes (ConnectionNodes): id and the_geom are used
            node_id_counter (iterable): an iterable yielding integers

        Returns:
            Grid with data in the following columns:
            - nodes.id: ids generated based on counter
            - nodes.coordinates: node coordinates (from self.the_geom)
            - nodes.content_type: TYPE_V2_CONNECTION_NODES
            - nodes.content_pk: the user-supplied id
            - nodes.node_type: NODE_1D_NO_STORAGE / NODE_1D_STORAGE
            - nodes.calculation_type: only if set on Manhole
        """
        return cls(connection_nodes.get_nodes(node_id_counter), Lines(id=[]))

    @classmethod
    def from_channels(
        cls,
        connection_nodes,
        channels,
        global_dist_calc_points,
        node_id_counter,
        line_id_counter,
        connection_node_offset=0,
    ):
        """Construct a grid for the channels

        Args:
            connection_nodes (ConnectionNodes): used to map ids to indices
            channels (Channels)
            global_dist_calc_points (float): Default node interdistance.
            node_id_counter (iterable): an iterable yielding integers
            line_id_counter (iterable): an iterable yielding integers
            connection_node_offset (int): offset to give connection node
              indices in the returned lines.line. Default 0.

        Returns:
            Grid with data in the following columns:
            - nodes.id: counter generated here starting from node_id_offset
            - nodes.coordinates
            - nodes.content_type: ContentType.TYPE_V2_CHANNEL
            - nodes.content_pk: the id of the Channel from which this node originates
            - nodes.node_type: NODE_1D_NO_STORAGE
            - nodes.calculation_type: from the channel
            - lines.id: 0-based counter generated here
            - lines.line: lines between connection nodes and added channel
              nodes. The indices are offset using the respective parameters.
            - lines.content_type: ContentType.TYPE_V2_CHANNEL
            - lines.content_pk: the id of the Channel from which this line originates
            - lines.line_type: from the channel's calculation_type
        """
        nodes, segment_size = channels.interpolate_nodes(
            node_id_counter, global_dist_calc_points
        )
        lines = channels.get_lines(
            connection_nodes,
            nodes,
            line_id_counter,
            segment_size=segment_size,
            connection_node_offset=connection_node_offset,
        )
        return cls(nodes, lines)

    def set_channel_weights(self, locations, channels):
        """Set cross section weights to channel lines.

        The attributes lines.cross1, lines.cross2, lines.cross_weight are
        changed in place for lines whose content_type equals TYPE_V2_CHANNEL.

        Args:
            locations (CrossSectionLocations)
            channels (Channels): Used to lookup the channel geometry
        """
        cross_sections.compute_weights(self.lines, locations, channels)

    def set_calculation_types(self):
        """Set the calculation types for connection nodes that do not yet have one.

        The calculation_type of a connection nodes is based on
          - connection_node.manhole.calculation_type (set in ConnectionNode.get_nodes)
          - if not present, the calculation_type of adjacent lines is taken. if
            these are different the precedence is: 1 > 2 > 0
          - if not present, then isolated
        """
        node_mask = (
            self.nodes.content_type == ContentType.TYPE_V2_CONNECTION_NODES
        ) & (self.nodes.calculation_type == -9999)
        line_mask = (self.lines.content_type == ContentType.TYPE_V2_CHANNEL) | (
            self.lines.content_type == ContentType.TYPE_V2_PIPE
        )
        lines_line = self.lines.line[line_mask]
        lines_line_type = self.lines.line_type[line_mask]

        for node_idx in np.where(node_mask)[0]:
            node_id = self.nodes.id[node_idx]
            # fine the lines that connect to this node
            line_idx = np.where(lines_line == node_id)[0]
            # fine the corresponding types and convert to a set
            line_types = set(lines_line_type[line_idx])
            # set the node's line type if there is any line with a type
            for line_type in LINE_TYPE_PRECEDENCE:
                if line_type in line_types:
                    self.nodes.calculation_type[node_idx] = LINE_TYPE_MAPPING[line_type]
                    break
            else:
                # if there was no line with a type; set ISOLATED
                self.nodes.calculation_type[node_idx] = CalculationType.ISOLATED

    def set_bottom_levels(self, connection_nodes, channels):
        """Set the bottom levels (dmax and dpumax) for 1D nodes and lines

        The types are based on:
        1. channels: crosssection locations
        2. pipes, connection nodes:
          - from the manhole.bottom_level
          - if not present (from neighboring points? ?)
        3. pipes, interpolated nodes:
          - interpolate between invert level start & end
        4. dpumax = greatest of 2 neighboring nodes dmax
           except for pipes with no interpolated nodes, then take pipe's own level
        """

    def finalize(self, epsg_code=None):
        """Finalize the Grid, computing and setting derived attributes"""
        self.lines.set_line_coords(self.nodes)
        self.lines.fix_line_geometries()
        self.lines.set_discharge_coefficients()
        self.epsg_code = epsg_code
