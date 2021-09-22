from threedigrid_builder.base import array_of
from threedigrid_builder.base import Lines
from threedigrid_builder.base import search
from threedigrid_builder.constants import CalculationType
from threedigrid_builder.constants import ContentType, NodeType
from threedigrid_builder.constants import LineType
from threedigrid_builder.exceptions import SchematisationError

import itertools
import numpy as np
import pygeos


__all__ = ["get_1d2d_lines", "ConnectedPoints"]


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


@array_of(ConnectedPoint)
class ConnectedPoints:
    def get_node_ids(self, grid):
        """Find the node ids to which connected points belong

        The supplied nodes should be ordered by content_pk (per content_type) and then
        by position on the linear object (channel / pipe / culvert). 1D boundary node
        ids must be sorted by boundary_id, but manholes do not have to be sorted by
        manhole_id. 
        """
        nodes = grid.nodes

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
                conn_pnt_ids = self.id[current_selection[e.indices]].tolist()
                raise SchematisationError(
                    f"Connected points {conn_pnt_ids} refer to non-existing objects."
                )

        # Find 1D boundary conditions
        current_selection = self.content_type == ContentType.TYPE_V2_1D_BOUNDARY_CONDITIONS
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
                conn_pnt_ids = self.id[current_selection[e.indices]].tolist()
                raise SchematisationError(
                    f"Connected points {conn_pnt_ids} refer to non-existing objects."
                )

        HAS_1D2D_NO_CONN_NODES = (
            ContentType.TYPE_V2_CHANNEL,
            ContentType.TYPE_V2_PIPE,
            ContentType.TYPE_V2_CULVERT,
        )
        for content_type in HAS_1D2D_NO_CONN_NODES:
            current_selection = self.content_type == content_type
            if not current_selection.any():
                continue

            # Search the nodes for the ones with the correct content_type & content_pk
            try:
                _node_idx = search(
                    nodes.content_pk,
                    self.content_pk[current_selection],
                    mask=(nodes.content_type == content_type),
                    assume_ordered=True,
                    check_exists=True,
                )
            except KeyError as e:
                conn_pnt_ids = self.id[current_selection[e.indices]].tolist()
                raise SchematisationError(
                    f"Connected points {conn_pnt_ids} refer to non-existing objects."
                )
            node_idx[current_selection] = _node_idx

        # handle node_number
        is_linear = np.isin(self.content_type, HAS_1D2D_NO_CONN_NODES)
        mask = is_linear & (self.node_number > 1)
        node_idx[mask] += self.node_number[mask] - 2

        # first and last nodes are actually connection nodes
        is_first_node = is_linear & (self.node_number == 1)

        # identify last nodes by checking if the the node_idx + node_number - 2 leads
        # to a different (the next) channel/pipe/culvert
        is_last_node = is_linear & (
            (np.take(nodes.content_type, node_idx, mode="wrap") != self.content_type)
            | (np.take(nodes.content_pk, node_idx, mode="wrap") != self.content_pk)
            | (node_idx >= len(nodes))
        )
        # check if the previous one is indeed an interpolated node
        node_idx[is_last_node] -= 1
        bad_node_number = is_linear & (
            (np.take(nodes.content_type, node_idx, mode="wrap") != self.content_type)
            | (np.take(nodes.content_pk, node_idx, mode="wrap") != self.content_pk)
            | (node_idx >= len(nodes))
            | (self.node_number < 1)
        )
        if bad_node_number.any():
            raise SchematisationError(
                f"Connected points {self.id[bad_node_number]} have a node number "
                f"that is out of bounds."
            )

        node_ids = nodes.index_to_id(node_idx)

        # handle first nodes (connect to a connection node)
        if is_first_node.any():
            line_idx = search(
                grid.lines.line[:, 1],
                node_ids[is_first_node],
                assume_ordered=False,
                check_exists=True,
            )
            node_ids[is_first_node] = grid.lines.line[line_idx, 0]

        # handle last nodes (connect to a connection node)
        if is_last_node.any():
            line_idx = search(
                grid.lines.line[:, 0],
                node_ids[is_last_node],
                assume_ordered=False,
                check_exists=True,
            )
            node_ids[is_last_node] = grid.lines.line[line_idx, 1]

        return node_ids


def get_1d2d_lines(
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
        nodes (Nodes): All 1D and 2D nodes to compute 1D-2D lines for. Be sure that
          the 1D nodes already have a dmax and calculation_type set.
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
    connected_idx = np.where(
        np.isin(nodes.content_type, HAS_1D2D_CONTENT_TYPES)
        & np.isin(nodes.calculation_type, HAS_1D2D_CALC_TYPES)
    )[0]
    n_connected_nodes = connected_idx.shape[0]
    if n_connected_nodes == 0:
        return Lines(id=[])

    # The query_bulk returns 2 1D arrays: one with indices into the supplied node
    # geometries and one with indices into the tree of cells.
    idx = cell_tree.query_bulk(pygeos.points(nodes.coordinates[connected_idx]))
    # Address edge cases of multiple 1D-2D lines per node: just take the one
    _, first_unique_index = np.unique(idx[0], return_index=True)
    idx = idx[:, first_unique_index]
    n_lines = idx.shape[1]
    # Error if there is a node without a 1D-2D line
    if n_lines != n_connected_nodes:
        # The code in this if clause is only for pretty error formatting.
        out_of_bounds = np.delete(connected_idx, idx[0])
        types = nodes.content_type[out_of_bounds]
        pks = nodes.content_pk[out_of_bounds]
        pretty_names = ("connection nodes", "channels", "pipes", "culverts")
        object_pk_list = [
            f"{pretty_name} {sorted(set(pks[types == content_type]))}"
            for content_type, pretty_name in zip(HAS_1D2D_CONTENT_TYPES, pretty_names)
            if np.any(types == content_type)
        ]
        raise SchematisationError(
            f"The following object(s) have a connected calculation type but are "
            f"(partially) outside of the 2D calculation cells: "
            f"{', '.join(object_pk_list)}."
        )
    if n_lines == 0:
        return Lines(id=[])
    node_idx = connected_idx[idx[0]]  # convert to node indexes
    node_id = nodes.index_to_id(node_idx)  # convert to node ids
    cell_id = nodes.index_to_id(idx[1])  # convert to cell ids

    # create a 'duplicator' array that duplicates lines from double connected nodes
    is_double = nodes.calculation_type[node_idx] == CalculationType.DOUBLE_CONNECTED
    duplicator = np.ones(n_lines, dtype=int)
    duplicator[is_double] = 2
    duplicator = np.repeat(np.arange(n_lines), duplicator)

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
    options = is_double * 2 + is_closed
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
        id=itertools.islice(line_id_counter, len(duplicator)),
        line=np.array([node_id, cell_id]).T[duplicator],
        kcu=kcu[duplicator],
        dpumax=dpumax[duplicator],
    )
