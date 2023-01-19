import dataclasses

import numpy as np
import shapely

from threedigrid_builder.base import Lines, Nodes, OutputInterface
from threedigrid_builder.constants import (
    BoundaryType,
    CalculationType,
    ContentType,
    LineType,
    Material,
    NodeType,
)
from threedigrid_builder.grid import Grid, GridMeta, PotentialBreaches

__all__ = ["DictOut"]


NODE_FIELDS = (
    "id",
    "node_type",
    "calculation_type",
    "content_type",
    "content_pk",
    "dmax",  # bottom level
    "dimp",  # bottom level for groundwater
    "storage_area",
    "boundary_id",  # referring to the id of the boundary condition
    "boundary_type",
    "manhole_id",  # referring to the id of a manhole
    "initial_waterlevel",
)


CELL_FIELDS = (
    "id",
    "node_type",
    "boundary_id",  # referring to the id of the boundary condition
    "boundary_type",
    "has_dem_averaged",  # Boolean attribute to tell if dem is averaged for node.
)

LINE_FIELDS = (
    "id",
    "kcu",  # calculation type of the line
    "node_1",  # node id
    "node_2",  # node id
    "ds1d",  # arclength
    "ds1d_half",
    "content_type",
    "content_pk",
    "dpumax",  # bottom_level at the velocity point
    "invert_level_start_point",  # bottom level at line start
    "invert_level_end_point",  # bottom level at line end
    "flod",  # obstacle height
    "cross_id1",  # the id of the cross section definition
    "cross_id2",  # the id of the cross section definition
    "cross_weight",
    "frict_type1",
    "frict_type2",
    "frict_value1",
    "frict_value2",
    "discharge_coefficient_positive",
    "discharge_coefficient_negative",
)


EMBEDDED_NODE_FIELDS = (
    "id",
    "content_type",
    "content_pk",
    "embedded_in",
)


META_FIELDS = (
    "epsg_code",
    "threedigrid_builder_version",
    "has_1d",
    "has_2d",
    "has_embedded",
    "has_breaches",
    "has_groundwater",
    "has_groundwater_flow",
    "has_interception",
    "has_pumpstations",
    "has_simple_infiltration",
    "has_max_infiltration_capacity",
    "has_interflow",
    "has_initial_waterlevels",
)


def _enum_to_str(arr, enum_type):
    result = np.full_like(arr, "", dtype=object)
    result[arr != -9999] = [enum_type(x).name for x in arr[arr != -9999]]
    return result


def increase(arr):
    """Increase arr by one where arr is not -9999"""
    return np.add(arr, 1 * (arr != -9999), dtype=arr.dtype)


class DictOut(OutputInterface):
    def __init__(self, path):
        if path is not None:
            raise ValueError("The dict writer cannot write to a path")
        super().__init__(path)

    @staticmethod
    def available():
        return True

    def __enter__(self):
        return self

    def __exit__(self, *args, **kwargs):
        pass

    def write(self, grid: Grid, geometry_format="wkb"):
        """Write a grid to a dictionary of layers.

        A layer is a dictionary of 1D ndarrays. The following layers are returned:

        - nodes: all computational nodes as Point geometries.
        - cells: a subselection of the nodes as Polygon geometries
        - lines: all flowlines as LineString geometries. lines interconnect nodes.
        - breaches: potential breach locations (Point geometries)
        - nodes_embedded: (virtual) nodes that are embedded into computational nodes
        - meta: metadata (e.g. "epsg_code"). the fields are scalars, not 1D ndarrays

        Args:
            grid (Grid)
            geometry_format {"native", "wkb", "wkt"}: optionally serialize geometries

        Returns:
            dict of dicts of 1D ndarrays
        """
        if geometry_format == "native":
            geom_serializer = None
        elif geometry_format == "wkb":
            geom_serializer = shapely.to_wkb
        elif geometry_format == "wkt":
            geom_serializer = shapely.to_wkt
        else:
            raise ValueError(f"Unknown geometry format '{geometry_format}'")

        nodes, cells = self.get_nodes_cells(grid.nodes)
        nodes_emb = self.get_embedded(grid.nodes_embedded)
        lines = self.get_lines(grid.lines)
        breaches = self.get_breaches(grid.breaches)
        meta = self.get_meta(grid.meta)
        result = {
            "nodes": nodes,
            "cells": cells,
            "lines": lines,
            "breaches": breaches,
            "nodes_embedded": nodes_emb,
            "meta": meta,
        }
        if geom_serializer is not None:
            for key in result.keys():
                if result[key] is None or "geometry" not in result[key]:
                    continue
                result[key]["geometry"] = geom_serializer(result[key]["geometry"])
        return result

    def get_nodes_cells(self, nodes: Nodes):
        """Create "nodes" and "cells" dictionaries

        Args:
            nodes (Nodes)

        Returns:
            tuple 2 dicts of 1D ndarrays (nodes, cells)
        """
        node_data = nodes.to_dict()

        is_2d = np.isin(
            nodes.node_type, (NodeType.NODE_2D_OPEN_WATER, NodeType.NODE_2D_BOUNDARIES)
        )

        # construct points from nodes.coordinates
        node_geometries = np.empty(len(nodes), dtype=object)
        coordinates = node_data.pop("coordinates")
        has_coord = np.isfinite(coordinates).all(axis=1)
        node_geometries[has_coord] = shapely.points(coordinates[has_coord])

        # construct cells from nodes.bounds
        bounds = node_data.pop("bounds")[is_2d]
        has_cell = np.isfinite(bounds).all(axis=1)
        cell_geometries = np.empty(bounds.shape[0], dtype=object)
        cell_geometries[has_cell] = shapely.box(*bounds[has_cell].T)

        # convert enums to strings
        node_data["node_type"] = _enum_to_str(node_data["node_type"], NodeType)
        node_data["content_type"] = _enum_to_str(node_data["content_type"], ContentType)
        node_data["calculation_type"] = _enum_to_str(
            node_data["calculation_type"], CalculationType
        )
        node_data["boundary_type"] = _enum_to_str(
            node_data["boundary_type"], BoundaryType
        )

        # go from 0-based to 1-based indexing
        for field in ("id", "embedded_in"):
            node_data[field] = increase(node_data[field])

        node_data_filt = {field: node_data[field] for field in NODE_FIELDS}
        node_data_filt["geometry"] = node_geometries

        cell_data = {field: node_data[field][is_2d] for field in CELL_FIELDS}
        cell_data["geometry"] = cell_geometries

        return node_data_filt, cell_data

    def get_embedded(self, nodes_embedded: Nodes):
        """Creates the "nodes_embedded" dictionary

        Args:
            nodes_embedded (Nodes)

        Returns:
            dicts of 1D ndarrays
        """
        if nodes_embedded is None or len(nodes_embedded) == 0:
            return None
        node_data = nodes_embedded.to_dict()

        # construct points from nodes.coordinates
        node_geometries = np.empty(len(nodes_embedded), dtype=object)
        coordinates = node_data.pop("coordinates")
        has_coord = np.isfinite(coordinates).all(axis=1)
        node_geometries[has_coord] = shapely.points(coordinates[has_coord])

        # convert enums to strings
        node_data["node_type"] = _enum_to_str(node_data["node_type"], NodeType)
        node_data["content_type"] = _enum_to_str(node_data["content_type"], ContentType)
        node_data["calculation_type"] = _enum_to_str(
            node_data["calculation_type"], CalculationType
        )
        node_data["boundary_type"] = _enum_to_str(
            node_data["boundary_type"], BoundaryType
        )

        # go from 0-based to 1-based indexing
        for field in ("id", "embedded_in"):
            node_data[field] = increase(node_data[field])

        node_data_filt = {field: node_data[field] for field in EMBEDDED_NODE_FIELDS}
        node_data_filt["geometry"] = node_geometries

        return node_data_filt

    def get_lines(self, lines: Lines):
        """Get "lines" dictionary

        Args:
            lines (Lines)

        Returns:
            lines dict of 1D ndarrays
        """
        line_data = lines.to_dict()

        # convert enums to strings
        line_data["kcu"] = _enum_to_str(line_data["kcu"], LineType)
        line_data["content_type"] = _enum_to_str(line_data["content_type"], ContentType)

        # go from 0-based to 1-based indexing
        for field in ("id", "line"):
            line_data[field] = increase(line_data[field])

        # cast lines.line to 2 1D arrays
        line_data["node_1"], line_data["node_2"] = line_data.pop("line").T

        geometries = line_data["line_geometries"]
        line_data = {field: line_data[field] for field in LINE_FIELDS}
        line_data["geometry"] = geometries
        return line_data

    def get_breaches(self, breaches: PotentialBreaches):
        """Get "breaches" dictionary

        Args:
            breaches (Breaches)

        Returns:
            breaches dict of 1D ndarrays
        """
        if len(breaches) == 0:
            return

        # sort by id
        return {
            "id": increase(breaches.id),
            "line_id": increase(breaches.line_id),
            "content_pk": breaches.content_pk,
            "maximum_breach_depth": breaches.maximum_breach_depth,
            "levee_material": _enum_to_str(breaches.levee_material, Material),
            "geometry": breaches.the_geom,
            "code": breaches.code,
            "display_name": breaches.display_name,
        }

    def get_meta(self, meta: GridMeta):
        meta_data = dataclasses.asdict(meta)
        meta_data = {field: meta_data[field] for field in META_FIELDS}
        return meta_data
