from threedigrid_builder.base import Lines
from threedigrid_builder.base import Nodes
from threedigrid_builder.constants import ContentType
from threedigrid_builder.grid import cross_sections

import itertools
import pygeos


__all__ = ["Grid"]


class Grid:
    def __init__(self, nodes: Nodes, lines: Lines):
        if not isinstance(nodes, Nodes):
            raise TypeError(f"Expected Nodes instance, got {type(nodes)}")
        if not isinstance(lines, Lines):
            raise TypeError(f"Expected Lines instance, got {type(lines)}")
        self.nodes = nodes
        self.lines = lines
        self.epsg_code = None  # a Grid is aware of its projection
        self.pixel_size = None  # a Grid is aware of the pixel size in the DEM

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
    def from_quadtree(cls, quadtree, subgrid_meta):
        """Construct the 2D grid based on the quadtree object.
        """

        nodes, lines = quadtree.get_nodes_lines(subgrid_meta["area_mask"])

        return cls(nodes=nodes, lines=lines)

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
            - nodes.content_type: ContentType.TYPE_V2_CONNECTION_NODES
            - nodes.content_pk: the user-supplied id
        """
        nodes = Nodes(
            id=itertools.islice(node_id_counter, len(connection_nodes)),
            coordinates=pygeos.get_coordinates(connection_nodes.the_geom),
            content_type=ContentType.TYPE_V2_CONNECTION_NODES,
            content_pk=connection_nodes.id,
        )
        return cls(nodes, Lines(id=[]))

    @classmethod
    def from_channels(
        cls, connection_nodes, channels, global_dist_calc_points, node_id_counter
    ):
        """Construct a grid for the channels

        Args:
            connection_nodes (ConnectionNodes): used to map ids to indices
            channels (Channels)
            global_dist_calc_points (float): Default node interdistance.
            node_id_counter (iterable): an iterable yielding integers

        Returns:
            Grid with data in the following columns:
            - nodes.id: counter generated here starting from node_id_offset
            - nodes.coordinates
            - nodes.content_type: ContentType.TYPE_V2_CHANNEL
            - nodes.content_pk: the id of the Channel from which this node originates
            - lines.id: 0-based counter generated here
            - lines.line: lines between connetion nodes and added channel
              nodes. The indices are offset using the respective parameters.
            - lines.content_type: ContentType.TYPE_V2_CHANNEL
            - lines.content_pk: the id of the Channel from which this line originates
        """
        nodes, segment_size = channels.interpolate_nodes(
            node_id_counter, global_dist_calc_points
        )
        lines = channels.get_lines(connection_nodes, nodes, segment_size=segment_size)
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

    def finalize(self, epsg_code=None, pixel_size=None):
        """Finalize the Grid, computing and setting derived attributes"""
        self.lines.set_line_coords(self.nodes)
        self.lines.fix_line_geometries()
        self.epsg_code = epsg_code
        self.pixel_size = pixel_size