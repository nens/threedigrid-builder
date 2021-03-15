from threedigrid_builder.base import Lines
from threedigrid_builder.base import Nodes
from threedigrid_builder.constants import NodeType
from threedigrid_builder.constants import LineType
from threedigrid_builder.constants import ContentType
from threedigrid_builder.grid import cross_sections

import itertools
import numpy as np
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
        self.meta = {}  # Grid has meta information based on nodes and lines.
        self.attrs = (
            {}
        )  # Grid has attribute information based on attributes from input.

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
            - lines.id: 0-based counter generated here
            - lines.line: lines between connection nodes and added channel
              nodes. The indices are offset using the respective parameters.
            - lines.content_type: ContentType.TYPE_V2_CHANNEL
            - lines.content_pk: the id of the Channel from which this line originates
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

    def set_meta(self):
        self.meta["n2dtot"] = np.count_nonzero(
            self.nodes.node_type == NodeType.NODE_2D_OPEN_WATER
        )
        self.meta["n1dtot"] = np.count_nonzero(
            self.nodes.node_type == NodeType.NODE_1D_NO_STORAGE
        ) + np.count_nonzero(self.nodes.node_type == NodeType.NODE_1D_STORAGE)

        self.meta["liutot"] = np.count_nonzero(
            self.lines.line_type == LineType.LINE_2D_U
        )
        self.meta["livtot"] = np.count_nonzero(
            self.lines.line_type == LineType.LINE_2D_V
        )

        line_types_1d = [
            LineType.LINE_1D_EMBEDDED,
            LineType.LINE_1D_ISOLATED,
            LineType.LINE_1D_CONNECTED,
            LineType.LINE_1D_LONG_CRESTED,
            LineType.LINE_1D_SHORT_CRESTED,
            LineType.LINE_1D_DOUBLE_CONNECTED,
        ]
        self.meta["l1dtot"] = sum(
            [np.count_nonzero(self.lines.line_type == x) for x in line_types_1d]
        )

        line_types_1d2d = [
            LineType.LINE_1D2D_SINGLE_CONNECTED_WITH_STORAGE,
            LineType.LINE_1D2D_SINGLE_CONNECTED_WITHOUT_STORAGE,
            LineType.LINE_1D2D_DOUBLE_CONNECTED_WITH_STORAGE,
            LineType.LINE_1D2D_DOUBLE_CONNECTED_WITHOUT_STORAGE,
            LineType.LINE_1D2D_POSSIBLE_BREACH,
            LineType.LINE_1D2D_ACTIVE_BREACH,
        ]
        self.meta["infl1d"] = sum(
            [np.count_nonzero(self.lines.line_type == x) for x in line_types_1d2d]
        )

        line_types_1d2d_gw = [
            LineType.LINE_1D2D_GROUNDWATER_57,
            LineType.LINE_1D2D_GROUNDWATER_58,
        ]
        self.meta["ingrw1d"] = sum(
            [np.count_nonzero(self.lines.line_type == x) for x in line_types_1d2d_gw]
        )

    def set_attributes(self, epsg_code=None):

        self.attrs["epsg_code"] = epsg_code

        is_1d = (self.nodes.node_type == NodeType.NODE_1D_NO_STORAGE) + (
            self.nodes.node_type == NodeType.NODE_1D_STORAGE
        )
        if is_1d.any():
            extent_1d = np.array(
                [
                    np.amin(self.nodes.coordinates[is_1d, 0]),
                    np.amin(self.nodes.coordinates[is_1d, 1]),
                    np.amax(self.nodes.coordinates[is_1d, 0]),
                    np.amax(self.nodes.coordinates[is_1d, 1]),
                ]
            )
            self.attrs["extent_1d"] = extent_1d
            self.attrs["has_1d"] = 1
        else:
            self.attrs["extent_1d"] = np.array([-9999.0, -9999.0, -9999.0, -9999.0])
            self.attrs["has_1d"] = 0

        is_2d = self.nodes.node_type == NodeType.NODE_2D_OPEN_WATER
        if is_2d.any:
            extent_2d = np.array(
                [
                    np.amin(self.nodes.coordinates[is_2d, 0]),
                    np.amin(self.nodes.coordinates[is_2d, 1]),
                    np.amax(self.nodes.coordinates[is_2d, 0]),
                    np.amax(self.nodes.coordinates[is_2d, 1]),
                ]
            )
            self.attrs["extent_2d"] = extent_2d
            self.attrs["has_2d"] = 1
        else:
            self.attrs["extent_2d"] = np.array([-9999.0, -9999.0, -9999.0, -9999.0])
            self.attrs["has_2d"] = 0
        self.attrs["has_interception"] = 0
        self.attrs["has_pumpstations"] = 0
        self.attrs["has_simple_infiltration"] = 0
        self.attrs["model_name"] = "..."
        self.attrs["model_slug"] = "..."
        self.attrs["revision_hash"] = "..."
        self.attrs["revision_nr"] = 0
        self.attrs["threedigrid_builder_version"] = 0

    def finalize(self, epsg_code=None):
        """Finalize the Grid, computing and setting derived attributes"""
        self.lines.set_line_coords(self.nodes)
        self.lines.fix_line_geometries()
        self.set_meta()
        self.set_attributes(epsg_code=epsg_code)
