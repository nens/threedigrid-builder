from . import connection_nodes
from . import cross_sections
from threedigrid_builder.base import Lines
from threedigrid_builder.base import Nodes
from threedigrid_builder.constants import ContentType


__all__ = ["Grid"]


class Grid:
    def __init__(self, nodes: Nodes, lines: Lines, epsg_code=None, quadtree_stats=None):
        if not isinstance(nodes, Nodes):
            raise TypeError(f"Expected Nodes instance, got {type(nodes)}")
        if not isinstance(lines, Lines):
            raise TypeError(f"Expected Lines instance, got {type(lines)}")
        self.nodes = nodes
        self.lines = lines
        self.epsg_code = epsg_code  # Grid is aware of its epsg_code
        self.quadtree_stats = quadtree_stats

    def __add__(self, other):
        """Concatenate two grids without renumbering nodes."""
        if self.__class__ is not other.__class__:
            raise TypeError(
                "Cannot concatenate {self} with {other} as they are not of "
                "equal types."
            )
        new_attrs = {}
        for k, v in other.__dict__.items():
            if isinstance(v, Nodes) or isinstance(v, Lines):
                new_attrs[k] = getattr(self, k) + v
            else:
                new_attrs[k] = v if v is not None else getattr(self, k)
        return self.__class__(**new_attrs)

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

        # Some general quadtree grid statistics we need in the .h5 later on.
        quadtree_stats = {
            "lgrmin": quadtree.lgrmin,
            "kmax": quadtree.kmax,
            "mmax": quadtree.mmax,
            "nmax": quadtree.nmax,
            "dx": quadtree.dx,
            "dxp": quadtree.pixel_size,
            "x0p": quadtree.origin[0],
            "y0p": quadtree.origin[1],
        }

        return cls(nodes=nodes, lines=lines, quadtree_stats=quadtree_stats)

    @classmethod
    def from_connection_nodes(
        cls, connection_nodes, node_id_counter, connection_node_offset=0
    ):
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
            - lines.kcu: from the channel's calculation_type
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
        # Mask the lines to only the Channel lines
        line_mask = self.lines.content_type == ContentType.TYPE_V2_CHANNEL
        cross1, cross2, cross_weight = cross_sections.compute_weights(
            self.lines.content_pk[line_mask],
            self.lines.ds1d[line_mask],
            locations,
            channels,
        )
        self.lines.cross1[line_mask] = cross1
        self.lines.cross2[line_mask] = cross2
        self.lines.cross_weight[line_mask] = cross_weight

    def set_calculation_types(self):
        """Set the calculation types for connection nodes that do not yet have one.

        ConnectionNodes, Channels and Pipes should be present in the grid.
        """
        connection_nodes.set_calculation_types(self.nodes, self.lines)

    def set_bottom_levels(self, locations, channels, pipes, weirs, culverts):
        """Set the bottom levels (dmax and dpumax) for 1D nodes and lines

        Note that there should not be 2D nodes & lines in the grid yet.

        The levels are based on:
        1. channel nodes: interpolate between crosssection locations
        2. connection nodes: see connection_nodes.set_bottom_levels
        3. (not implemented) pipes, interpolated nodes:
          - interpolate between invert level start & end
        4. lines: dpumax = greatest of the two neighboring nodes dmax
          - except for channels with no interpolated nodes: take reference level, but
            only if that is higher than the two neighboring nodes.
        """
        # Channels, interpolated nodes
        mask = self.nodes.content_type == ContentType.TYPE_V2_CHANNEL
        self.nodes.dmax[mask] = cross_sections.compute_bottom_level(
            self.nodes.content_pk[mask], self.nodes.ds1d[mask], locations, channels
        )

        # Connection nodes: complex logic based on the connected objects
        connection_nodes.set_bottom_levels(
            self.nodes, self.lines, locations, channels, pipes, weirs, culverts
        )

        # Lines: based on the nodes
        self.lines.set_bottom_levels(self.nodes, allow_nan=False)

        # Fix channel lines: set dpumax of channel lines that have no interpolated nodes
        cross_sections.fix_dpumax(self.lines, self.nodes, locations)

    def finalize(self, epsg_code=None):
        """Finalize the Grid, computing and setting derived attributes"""
        self.lines.set_line_coords(self.nodes)
        self.lines.fix_line_geometries()
        self.lines.set_discharge_coefficients()
        self.epsg_code = epsg_code
