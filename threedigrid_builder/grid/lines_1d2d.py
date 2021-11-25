from threedigrid_builder.base import array_of
from threedigrid_builder.base import Breaches
from threedigrid_builder.base import Lines
from threedigrid_builder.base import search
from threedigrid_builder.constants import CalculationType
from threedigrid_builder.constants import ContentType
from threedigrid_builder.constants import LineType
from threedigrid_builder.constants import NodeType
from threedigrid_builder.exceptions import SchematisationError

import itertools
import numpy as np
import pygeos


__all__ = ["ConnectedPoints"]


HAS_1D2D_CONTENT_TYPES = (
    ContentType.TYPE_V2_CONNECTION_NODES,
    ContentType.TYPE_V2_CHANNEL,
    ContentType.TYPE_V2_PIPE,
    ContentType.TYPE_V2_CULVERT,
)

HAS_1D2D_CALC_TYPES = (
    CalculationType.CONNECTED,
    CalculationType.DOUBLE_CONNECTED,
)


class ConnectedPoint:
    id: int
    the_geom: pygeos.Geometry
    exchange_level: float
    content_type: ContentType
    content_pk: int
    node_number: int
    calculation_point_id: int  # only for error messages
    levee_id: int


@array_of(ConnectedPoint)
class ConnectedPoints:
    def get_node_index(self, nodes, channels, pipes, culverts):
        """Find the node index to which connected points belong.

        This method uses the ConnectedPoint's content_pk, content_type and node_number
        to find the index into grid.nodes to which a connected point belongs.

        The supplied nodes should be ordered by content_pk (per content_type) and then
        by position on the linear object (channel / pipe / culvert). 1D boundary node
        ids must be sorted by boundary_id, but manholes do not have to be sorted by
        manhole_id.
        """
        node_idx = np.empty_like(self.id)

        # Find manholes
        current_selection = self.content_type == ContentType.TYPE_V2_MANHOLE
        if current_selection.any():
            try:
                node_idx[current_selection] = search(
                    nodes.manhole_id,
                    self.content_pk[current_selection],
                    mask=(nodes.manhole_id != -9999),
                    assume_ordered=False,
                    check_exists=True,
                )
            except KeyError as e:
                bad = self.calculation_point_id[
                    np.where(current_selection)[0][e.indices]
                ]
                raise SchematisationError(
                    f"Calculation points {np.unique(bad).tolist()} refer to "
                    f"non-existing manholes."
                )

        # Find 1D boundary conditions
        current_selection = (
            self.content_type == ContentType.TYPE_V2_1D_BOUNDARY_CONDITIONS
        )
        if current_selection.any():
            try:
                node_idx[current_selection] = search(
                    nodes.boundary_id,
                    self.content_pk[current_selection],
                    mask=(nodes.node_type == NodeType.NODE_1D_BOUNDARIES),
                    assume_ordered=True,
                    check_exists=True,
                )
            except KeyError as e:
                bad = self.calculation_point_id[
                    np.where(current_selection)[0][e.indices]
                ]
                raise SchematisationError(
                    f"Calculation points {np.unique(bad).tolist()} refer to "
                    f"non-existing 1D boundary conditions."
                )

        HAS_1D2D_NO_CONN_NODES = {
            ContentType.TYPE_V2_CHANNEL: channels,
            ContentType.TYPE_V2_PIPE: pipes,
            ContentType.TYPE_V2_CULVERT: culverts,
        }
        for content_type, objs in HAS_1D2D_NO_CONN_NODES.items():
            current_selection = self.content_type == content_type
            if not current_selection.any():
                continue

            # Handle node_number == 1 seperately; they may not have interpolated nodes
            is_first_node = current_selection & (self.node_number == 1)
            if is_first_node.any():
                try:
                    obj_idx = objs.id_to_index(
                        self.content_pk[is_first_node], check_exists=True
                    )
                except KeyError as e:
                    bad = self.calculation_point_id[
                        np.where(is_first_node)[0][e.indices]
                    ]
                    raise SchematisationError(
                        f"Calculation points {np.unique(bad).tolist()} refer to "
                        f"non-existing {objs.__class__.__name__.lower()}."
                    )
                node_idx[is_first_node] = search(
                    nodes.content_pk,
                    objs.connection_node_start_id[obj_idx],
                    mask=(nodes.content_type == ContentType.TYPE_V2_CONNECTION_NODES),
                    assume_ordered=True,
                    check_exists=False,
                )

            # Search the nodes for the ones with node_number > 1
            current_selection = current_selection & (self.node_number > 1)
            if not current_selection.any():
                continue
            _node_idx = search(
                nodes.content_pk,
                self.content_pk[current_selection],
                mask=(nodes.content_type == content_type),
                assume_ordered=True,
                check_exists=False,  # happens if there are no interpol. nodes
            )
            node_idx[current_selection] = _node_idx

        # Handle node_number > 1 cases (these are interpolated nodes or end nodes)
        is_linear = np.isin(self.content_type, list(HAS_1D2D_NO_CONN_NODES.keys()))
        bad_node_number = is_linear & (self.node_number < 1)
        if bad_node_number.any():
            bad = self.calculation_point_id[bad_node_number]
            raise SchematisationError(
                f"Calculation points {np.unique(bad).tolist()} have node "
                f"numbers below 1."
            )
        is_interpolated_node = is_linear & (self.node_number > 1)
        node_idx[is_interpolated_node] += self.node_number[is_interpolated_node] - 2

        # Identify end nodes by checking if the the node_idx + node_number - 2 leads
        # to a different (the next) channel/pipe/culvert
        is_end_node = is_interpolated_node & (
            (np.take(nodes.content_type, node_idx, mode="clip") != self.content_type)
            | (np.take(nodes.content_pk, node_idx, mode="clip") != self.content_pk)
            | (node_idx >= len(nodes))
        )
        # check if the previous one is indeed an interpolated node
        node_idx[is_end_node] -= 1
        bad_node_number = (
            is_linear
            & (self.node_number > 2)
            & (
                (
                    np.take(nodes.content_type, node_idx, mode="clip")
                    != self.content_type
                )
                | (np.take(nodes.content_pk, node_idx, mode="clip") != self.content_pk)
                | (node_idx >= len(nodes))
            )
        )
        if bad_node_number.any():
            bad = self.calculation_point_id[bad_node_number]
            raise SchematisationError(
                f"Calculation points {np.unique(bad).tolist()} have too "
                f"large node numbers."
            )

        # handle end nodes (find the connection node id)
        for content_type, objs in HAS_1D2D_NO_CONN_NODES.items():
            _is_end_node = is_end_node & (self.content_type == content_type)
            if not _is_end_node.any():
                continue
            try:
                obj_idx = objs.id_to_index(
                    self.content_pk[_is_end_node], check_exists=True
                )
            except KeyError as e:
                bad = self.calculation_point_id[np.where(_is_end_node)[0][e.indices]]
                raise SchematisationError(
                    f"Calculation points {np.unique(bad).tolist()} refer to "
                    f"non-existing {objs.__class__.__name__.lower()}."
                )
            node_idx[_is_end_node] = search(
                nodes.content_pk,
                objs.connection_node_end_id[obj_idx],
                mask=(nodes.content_type == ContentType.TYPE_V2_CONNECTION_NODES),
                assume_ordered=True,
                check_exists=False,
            )

        return node_idx

    def get_line_mappings(self, nodes, cp_node_idx):
        """Return a node index and a connected point index per 1D2D line.

        Returns 3 arrays of shape (n_lines, ):
        - an int array with a node index
        - an int array with a connected point index (-9999 if there is none)
        - a bool array with True where a line is double connected
        """
        # Initialize the arrays
        is_single = nodes.calculation_type == CalculationType.CONNECTED
        is_double = nodes.calculation_type == CalculationType.DOUBLE_CONNECTED
        n_lines = np.count_nonzero(is_single) + 2 * np.count_nonzero(is_double)
        if n_lines == 0:
            return None, None, None

        line1_node_idx = np.where(is_single | is_double)[0]
        line1_cp_idx = np.full_like(line1_node_idx, fill_value=-9999)
        line2_node_idx = np.where(is_double)[0]
        line2_cp_idx = np.full_like(line2_node_idx, fill_value=-9999)

        # Get the first occurences of connected points
        cp1_node_idx, cp1_idx, counts = np.unique(
            cp_node_idx, return_index=True, return_counts=True
        )
        if np.any(counts > 2):
            bad_ids = np.unique(self.calculation_point_id[cp1_node_idx[counts > 2]])
            raise SchematisationError(
                f"Calculation points {bad_ids.tolist()} have too many connected points."
            )
        # Get the second occurences of connected points by inverting the _1 array
        cp2_idx = np.delete(np.arange(len(cp_node_idx)), cp1_idx)
        cp2_node_idx = np.take(cp_node_idx, cp2_idx)

        # Lookup the first and second occurences seperately
        try:
            idx = search(
                line1_node_idx,
                cp1_node_idx,
                assume_ordered=True,
                check_exists=True,
            )
        except KeyError as e:
            bad = cp1_idx[e.indices]
            raise SchematisationError(
                f"Connected points {self.id[bad]} (calculation points "
                f"{self.calculation_point_id[bad]}) refer to objects that do not have a "
                f"'connected' calculation type."
            )
        line1_cp_idx[idx] = cp1_idx

        try:
            idx = search(
                line2_node_idx,
                cp2_node_idx,
                assume_ordered=True,
                check_exists=True,
            )
        except KeyError as e:
            bad = cp2_idx[e.indices]
            raise SchematisationError(
                f"Connected points {self.id[bad]} (calculation points "
                f"{self.calculation_point_id[bad]}) are a second reference to objects that "
                f"do not have a 'double connected' but only a 'connected' calculation "
                f"type."
            )
        line2_cp_idx[idx] = cp2_idx

        # Merge the arrays
        line_node_idx = np.concatenate([line1_node_idx, line2_node_idx])
        line_cp_idx = np.concatenate([line1_cp_idx, line2_cp_idx])

        # Sort them by node index
        sorter = np.argsort(line_node_idx)
        line_node_idx = line_node_idx[sorter]
        line_cp_idx = line_cp_idx[sorter]
        line_is_double = np.isin(line_node_idx, np.where(is_double)[0])
        return line_node_idx, line_cp_idx, line_is_double

    def get_lines(
        self,
        cell_tree,
        nodes,
        connection_nodes,
        channels,
        pipes,
        locations,
        culverts,
        line_id_counter,
    ):
        """Compute 1D-2D flowlines for (double) connected 1D nodes.

        The line will connect to the cell in which the ConnectedPoint geometry is
        located. If the ConnectedPoint is unavailable, the geometry of the 1D node
        itself is taken instead. For edge cases the line will connect to the cell with
        the lowest id.

        The following line attributes will be set:
        - kcu (line_type): according to the get_1d2d_properties method on the
          corresponding 1D objects (ConnectionNodes, Channels, Pipes, Culverts)
        - dpumax (bottom_level): to the exchange_level of the ConnectedPoint or, if
          there is no ConnectedPoint or its exchange_level is empty, according to
          the get_1d2d_properties method
        - content_type: TYPE_V2_ADDED_CALCULATION_POINT for lines that come from a
          ConnectedPoint
        - content_pk: the ConnectedPoint's id if available

        Args:
          cell_tree (pygeos.STRtree): An STRtree containing the cells. The indices into
            the STRtree must be equal to the indices into the nodes.
          nodes (Grid): Aall 1D and 2D nodes to compute 1D-2D lines for.
            Requires nodes from channels, pipes, culverts, manholes and 1d boundary
            conditions.
            Be sure that the nodes already have a dmax and calculation_type set.
          connection_nodes (ConnectionNodes): for looking up manhole_id and drain_level
          channels (Channels)
          pipes (Pipes)
          locations (CrossSectionLocations): for looking up bank_level
          culverts (Culverts)
          line_id_counter (iterable): An iterable yielding integers

        Returns:
          Lines with data in the following columns:
          - id: integer generated from line_id_counter
          - kcu: LINE_1D2D_* type (see above)
          - dpumax: based on drain_level or dmax (see above)
        """
        # Create an array of node indices that need a 1D2D connection

        line_node_idx, line_cp_idx, line_is_double = self.get_line_mappings(
            nodes, self.get_node_index(nodes, channels, pipes, culverts)
        )
        if line_node_idx is None:
            return Lines(id=[])

        n_lines = len(line_node_idx)
        line_has_cp = line_cp_idx != -9999

        # Collect the 2D sides of the 1D2D line (user supplied or 1D node coordinate)
        points = np.empty(n_lines, dtype=object)
        points[line_has_cp] = self.the_geom[line_cp_idx[line_has_cp]]
        points[~line_has_cp] = pygeos.points(
            nodes.coordinates[line_node_idx[~line_has_cp]]
        )

        # The query_bulk returns 2 1D arrays: one with indices into the supplied node
        # geometries and one with indices into the tree of cells.
        idx = cell_tree.query_bulk(points)
        # Address edge cases of multiple 1D-2D lines per node: just take the one
        _, first_unique_index = np.unique(idx[0], return_index=True)
        idx = idx[:, first_unique_index]
        n_lines = idx.shape[1]
        # Error if there is a node without a 1D-2D line
        if n_lines != line_node_idx.shape[0]:
            # The code in this if clause is only for pretty error formatting.
            out_of_bounds = np.delete(line_node_idx, idx[0])
            types = nodes.content_type[out_of_bounds]
            pks = nodes.content_pk[out_of_bounds]
            pretty_names = ("connection nodes", "channels", "pipes", "culverts")
            object_pk_list = [
                f"{pretty_name} {sorted(set(pks[types == content_type]))}"
                for content_type, pretty_name in zip(
                    HAS_1D2D_CONTENT_TYPES, pretty_names
                )
                if np.any(types == content_type)
            ]
            raise SchematisationError(
                f"The following object(s) have a connected calculation type but are "
                f"(partially) outside of the 2D calculation cells: "
                f"{', '.join(object_pk_list)}."
            )
        if n_lines == 0:
            return Lines(id=[])
        node_idx = line_node_idx[idx[0]]  # convert to node indexes
        node_id = nodes.index_to_id(node_idx)  # convert to node ids
        cell_id = nodes.index_to_id(idx[1])  # convert to cell ids

        # Identify different types of objects and dispatch to the associated functions
        is_ch = nodes.content_type[node_idx] == ContentType.TYPE_V2_CHANNEL
        is_cn = nodes.content_type[node_idx] == ContentType.TYPE_V2_CONNECTION_NODES
        is_pipe = nodes.content_type[node_idx] == ContentType.TYPE_V2_PIPE
        is_culvert = nodes.content_type[node_idx] == ContentType.TYPE_V2_CULVERT

        is_closed = np.zeros(n_lines, dtype=bool)
        dpumax = np.full(n_lines, fill_value=np.nan, dtype=np.float64)

        is_closed[is_cn], dpumax[is_cn] = connection_nodes.get_1d2d_properties(
            nodes,
            node_idx[is_cn],
            channels,
            locations,
        )
        is_closed[is_ch], dpumax[is_ch] = channels.get_1d2d_properties(
            nodes, node_idx[is_ch], locations
        )
        is_closed[is_pipe], dpumax[is_pipe] = pipes.get_1d2d_properties(
            nodes,
            node_idx[is_pipe],
            connection_nodes,
        )
        is_closed[is_culvert], dpumax[is_culvert] = culverts.get_1d2d_properties(
            nodes,
            node_idx[is_culvert],
            connection_nodes,
        )

        # Override the exchange_level if there is a connected point with one
        if line_has_cp.any():
            exchange_levels = np.take(self.exchange_level, line_cp_idx, mode="clip")
            use_exchange_level = line_has_cp & np.isfinite(exchange_levels)
            dpumax[use_exchange_level] = exchange_levels[use_exchange_level]

        # Make content_type and content_pk for tracing the ConnectedPoints
        content_type = np.where(
            line_has_cp,
            ContentType.TYPE_V2_ADDED_CALCULATION_POINT,
            -9999,
        )
        content_pk = np.full(n_lines, fill_value=-9999, dtype=np.int32)
        content_pk[line_has_cp] = self.index_to_id(line_cp_idx[line_has_cp])

        # map "is_closed" to "kcu" (including double/single connected properties)
        # map the two binary arrays on numbers 0, 1, 2, 3
        options = line_is_double * 2 + is_closed
        kcu = np.choose(
            options,
            choices=[
                LineType.LINE_1D2D_SINGLE_CONNECTED_OPEN_WATER,
                LineType.LINE_1D2D_SINGLE_CONNECTED_CLOSED,
                LineType.LINE_1D2D_DOUBLE_CONNECTED_OPEN_WATER,
                LineType.LINE_1D2D_DOUBLE_CONNECTED_CLOSED,
            ],
        )

        return Lines(
            id=itertools.islice(line_id_counter, n_lines),
            line=np.array([cell_id, node_id]).T,
            kcu=kcu,
            dpumax=dpumax,
            content_type=content_type,
            content_pk=content_pk,
        )

    def get_breaches(self, lines, levees):
        """Every ConnectedPoint that has a levee_id will get an associated Breach"""
        # get lines that have a connected point
        (line_idx,) = np.where(
            lines.content_type == ContentType.TYPE_V2_ADDED_CALCULATION_POINT
        )
        # refine to those connected points that have a levee
        conn_pnt_idx = self.id_to_index(lines.content_pk[line_idx])
        has_levee = self.levee_id[conn_pnt_idx] != -9999
        if not np.any(has_levee):
            return Breaches(id=[])
        line_idx = line_idx[has_levee]
        conn_pnt_idx = conn_pnt_idx[has_levee]
        levee_idx = levees.id_to_index(self.levee_id[conn_pnt_idx])

        # only consider levees that have the correct properties set
        mask = (np.isfinite(levees.max_breach_depth) & (levees.material != -9999))[
            levee_idx
        ]
        if not np.any(mask):
            return Breaches(id=[])
        line_idx = line_idx[mask]
        levee_idx = levee_idx[mask]

        # compute the intersections (use shortest_line and not intersects to
        # account for the possibility that the levee may not intersect the line)
        points = pygeos.get_point(
            pygeos.shortest_line(
                levees.the_geom[levee_idx], lines.line_geometries[line_idx]
            ),
            0,
        )

        return Breaches(
            id=range(len(line_idx)),
            levl=lines.index_to_id(line_idx),
            levee_id=levees.index_to_id(levee_idx),
            levmat=levees.material[levee_idx],
            levbr=levees.max_breach_depth[levee_idx],
            content_pk=lines.content_pk[line_idx],
            coordinates=np.array([pygeos.get_x(points), pygeos.get_y(points)]).T,
        )
