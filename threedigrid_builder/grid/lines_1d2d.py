from threedigrid_builder.base import array_of
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
    user_ref: str  # useful for debugging
    calc_pnt_id: int  # useful for debugging


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
                bad = np.where(current_selection)[0][e.indices]
                raise SchematisationError(
                    f"Connected points {self.id[bad]} (calculation points "
                    f"{self.calc_pnt_id[bad]}) refer to non-existing objects."
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
                bad = np.where(current_selection)[0][e.indices]
                raise SchematisationError(
                    f"Connected points {self.id[bad]} (calculation points "
                    f"{self.calc_pnt_id[bad]}) refer to non-existing objects."
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
                    bad = np.where(is_first_node)[0][e.indices]
                    raise SchematisationError(
                        f"Connected points {self.id[bad]} (calculation points "
                        f"{self.calc_pnt_id[bad]}) refer to non-existing objects."
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
            raise SchematisationError(
                f"Connected points {self.id[bad_node_number]} have a node number "
                f"that is out of bounds."
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
            raise SchematisationError(
                f"Connected points {self.id[bad_node_number]} have a node number "
                f"that is out of bounds."
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
                bad = np.where(_is_end_node)[0][e.indices]
                raise SchematisationError(
                    f"Connected points {self.id[bad]} (calculation points "
                    f"{self.calc_pnt_id[bad]}) refer to non-existing objects."
                )
            node_idx[_is_end_node] = search(
                nodes.content_pk,
                objs.connection_node_end_id[obj_idx],
                mask=(nodes.content_type == ContentType.TYPE_V2_CONNECTION_NODES),
                assume_ordered=True,
                check_exists=False,
            )

        return node_idx

    def get_line_mapping(self, nodes, conn_point_node_idx):
        """Return a node index and a connected point index per 1D2D line.

        Returns integer arrays and one boolean array:
        - an integer array with a node index for every 1D2D line
        - an integer array with a connected point index for every 1D2D line
        - a boolean array with whether a 1D2D is double connected

        If there is no ConnectedPoint for a 1D2D line, the connected point index is -9999.
        """
        node_is_single = np.where(nodes.calculation_type == CalculationType.CONNECTED)[
            0
        ]
        node_is_double = np.where(
            nodes.calculation_type == CalculationType.DOUBLE_CONNECTED
        )[0]
        n_lines = len(node_is_single) + 2 * len(node_is_double)
        if n_lines == 0:
            return None, None, None

        line_conn_point_idx = np.full((n_lines,), fill_value=-9999, dtype=int)

        # Merge the indices, doubling the double connected ones
        line_node_idx = np.concatenate([node_is_single, np.repeat(node_is_double, 2)])
        line_is_double = np.ones(len(line_node_idx), dtype=bool)
        line_is_double[: len(node_is_single)] = False
        line_is_first = np.ones(len(line_node_idx), dtype=bool)
        line_is_first[len(node_is_single) + 1 : len(line_node_idx) : 2] = False
        line_is_second = np.zeros(len(line_node_idx), dtype=bool)
        line_is_second[len(node_is_single) + 1 : len(line_node_idx) : 2] = True

        # Sort them in ascending node order
        sorter = np.argsort(line_node_idx)
        line_node_idx = line_node_idx[sorter]
        line_is_double = line_is_double[sorter]
        line_is_first = line_is_first[sorter]
        line_is_second = line_is_second[sorter]

        conn_point_node_idx_1, conn_point_idx_1, counts = np.unique(
            conn_point_node_idx, return_index=True, return_counts=True
        )
        if np.any(counts > 2):
            raise SchematisationError(
                "Some calculation nodes have too many connected points."
            )
        conn_point_idx_2 = np.delete(
            np.arange(len(conn_point_node_idx)), conn_point_idx_1
        )
        conn_point_node_idx_2 = conn_point_node_idx[conn_point_idx_2]
        try:
            idx = search(
                line_node_idx,
                conn_point_node_idx_1,
                mask=line_is_first,
                assume_ordered=True,
                check_exists=True,
            )
        except KeyError as e:
            bad = conn_point_idx_1[e.indices]
            raise SchematisationError(
                f"Connected points {self.id[bad]} (calculation points "
                f"{self.calc_pnt_id[bad]}) refer to objects that do not have a "
                f"'connected' calculation type."
            )
        line_conn_point_idx[idx] = conn_point_idx_1

        try:
            idx = search(
                line_node_idx,
                conn_point_node_idx_2,
                mask=line_is_second,
                assume_ordered=True,
                check_exists=True,
            )
        except KeyError as e:
            bad = conn_point_idx_1[e.indices]
            raise SchematisationError(
                f"Connected points {self.id[bad]} (calculation points "
                f"{self.calc_pnt_id[bad]}) are a second reference to objects that "
                f"do not have a 'double connected' but only a 'connected' calculation "
                f"type."
            )
        line_conn_point_idx[idx] = conn_point_idx_2

        return line_node_idx, line_conn_point_idx, line_is_double

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

        The line type and bottom level (kcu and dpumax) are set according to the
        "get_1d2d_properties" on the corresponding objects (ConnectionNodes, Channels,
        Pipes, Culverts).

        If the 2D bottom level will turn out higher, this will be corrected later by
        threedi-tables.

        The line will connect to the cell in which the 1D node is located. For edge cases,
        the line will connect to the cell with the lowest id. Users may want to influence
        which cells are connected to. This is currently unimplemented.

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

        line_node_idx, line_cp_idx, line_is_double = self.get_line_mapping(
            nodes, self.get_node_index(nodes, channels, pipes, culverts)
        )
        if line_node_idx is None:
            return Lines(id=[])

        # The query_bulk returns 2 1D arrays: one with indices into the supplied node
        # geometries and one with indices into the tree of cells.
        idx = cell_tree.query_bulk(pygeos.points(nodes.coordinates[line_node_idx]))
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
            nodes, node_idx[is_cn]
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
            line=np.array([node_id, cell_id]).T,
            kcu=kcu,
            dpumax=dpumax,
        )
