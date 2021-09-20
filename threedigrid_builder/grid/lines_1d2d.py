from threedigrid_builder.base import array_of
from threedigrid_builder.base import Lines
from threedigrid_builder.constants import CalculationType
from threedigrid_builder.constants import ContentType
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
        by position on the linear object (channel / pipe / culvert).
        """
        # TODO map manhole and boundary ids to connection node ids
        result = np.empty_like(self.id)

        for content_type in HAS_1D2D_CONTENT_TYPES - {
            ContentType.TYPE_V2_CONNECTION_NODES
        }:
            mask = self.content_type == content_type
            if not mask.any():
                continue

            node_idx = np.where(grid.nodes.content_type == content_type)[0]
            node_idx = node_idx[
                np.searchsorted(grid.nodes.content_pk[node_idx], self.content_pk[mask])
            ]
            missing = grid.nodes.content_pk[node_idx] != self.content_pk[mask]
            if missing.any():
                raise SchematisationError(
                    f"Connected points {self.id[missing]} refer to non-existing objects"
                )
            result[mask] = grid.nodes.index_to_id(node_idx)

        # print(result)
        return result


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
