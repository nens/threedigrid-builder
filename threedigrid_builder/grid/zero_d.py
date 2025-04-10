from dataclasses import dataclass
from typing import Dict, Optional

import numpy as np
import shapely

from threedigrid_builder.base import Array, surfaces
from threedigrid_builder.base.nodes import Nodes
from threedigrid_builder.constants import ContentType

__all__ = ["Surfaces"]


SURFACE_TYPE_PROPERTIES = {
    # Explanation:
    # outflow_delay: c (reaction factor) (/min)
    # surface_storage: berging (mm)
    # infiltration: boolean flag
    # fb: Max infiltration capacity (mm/h)
    # fe: Min infiltration capacity (mm/h)
    # ka: Time factor reduction (afname) of infiltration capacity (/h)
    # kh: Time factor recovery (herstel) of infiltration capacity (/h)
    # gesloten verharding, hellend
    "gvh_hel": {
        "outflow_delay": 0.5 / 60.0,
        "surface_storage": 0.0 / 1000.0,
        "infiltration": False,
        "fb": 0.0 / 60.0 / 60.0 / 1000.0,  # for consistency
        "fe": 0.0 / 60.0 / 60.0 / 1000.0,
        "ka": 0.0 / 60.0 / 60.0,
        "kh": 0.0 / 60.0 / 60.0,
    },
    # gesloten verharding, vlak
    "gvh_vla": {
        "outflow_delay": 0.2 / 60.0,
        "surface_storage": 0.5 / 1000.0,
        "infiltration": False,
        "fb": 0.0 / 60.0 / 60.0 / 1000.0,
        "fe": 0.0 / 60.0 / 60.0 / 1000.0,
        "ka": 0.0 / 60.0 / 60.0,
        "kh": 0.0 / 60.0 / 60.0,
    },
    # gesloten verharding, vlak uitgestrekt
    "gvh_vlu": {
        "outflow_delay": 0.1 / 60.0,
        "surface_storage": 1.0 / 1000.0,
        "infiltration": False,
        "fb": 0.0 / 60.0 / 60.0 / 1000.0,
        "fe": 0.0 / 60.0 / 60.0 / 1000.0,
        "ka": 0.0 / 60.0 / 60.0,
        "kh": 0.0 / 60.0 / 60.0,
    },
    # open verharding, hellend
    "ovh_hel": {
        "outflow_delay": 0.5 / 60.0,
        "surface_storage": 0.0 / 1000.0,
        "infiltration": True,
        "fb": 2.0 / 60.0 / 60.0 / 1000.0,
        "fe": 0.5 / 60.0 / 60.0 / 1000.0,
        "ka": 3.0 / 60.0 / 60.0,
        "kh": 0.1 / 60.0 / 60.0,
    },
    # open verharding, vlak
    "ovh_vla": {
        "outflow_delay": 0.2 / 60.0,
        "surface_storage": 0.5 / 1000.0,
        "infiltration": True,
        "fb": 2.0 / 60.0 / 60.0 / 1000.0,
        "fe": 0.5 / 60.0 / 60.0 / 1000.0,
        "ka": 3.0 / 60.0 / 60.0,
        "kh": 0.1 / 60.0 / 60.0,
    },
    # open verharding, vlak uitgestrekt
    "ovh_vlu": {
        "outflow_delay": 0.1 / 60.0,
        "surface_storage": 1.0 / 1000.0,
        "infiltration": True,
        "fb": 2.0 / 60.0 / 60.0 / 1000.0,
        "fe": 0.5 / 60.0 / 60.0 / 1000.0,
        "ka": 3.0 / 60.0 / 60.0,
        "kh": 0.1 / 60.0 / 60.0,
    },
    # dak, hellend
    "dak_hel": {
        "outflow_delay": 0.5 / 60.0,
        "surface_storage": 0.0 / 1000.0,
        "infiltration": False,
        "fb": 0.0 / 60.0 / 60.0 / 1000.0,
        "fe": 0.0 / 60.0 / 60.0 / 1000.0,
        "ka": 0.0 / 60.0 / 60.0,
        "kh": 0.0 / 60.0 / 60.0,
    },
    # dak, vlak
    "dak_vla": {
        "outflow_delay": 0.2 / 60.0,
        "surface_storage": 2.0 / 1000.0,
        "infiltration": False,
        "fb": 0.0 / 60.0 / 60.0 / 1000.0,
        "fe": 0.0 / 60.0 / 60.0 / 1000.0,
        "ka": 0.0 / 60.0 / 60.0,
        "kh": 0.0 / 60.0 / 60.0,
    },
    # dak, vlak uitgestrekt
    "dak_vlu": {
        "outflow_delay": 0.1 / 60.0,
        "surface_storage": 4.0 / 1000.0,
        "infiltration": False,
        "fb": 0.0 / 60.0 / 60.0 / 1000.0,
        "fe": 0.0 / 60.0 / 60.0 / 1000.0,
        "ka": 0.0 / 60.0 / 60.0,
        "kh": 0.0 / 60.0 / 60.0,
    },
    # onverhard, hellend
    "onv_hel": {
        "outflow_delay": 0.5 / 60.0,
        "surface_storage": 2.0 / 1000.0,
        "infiltration": True,
        "fb": 5.0 / 60.0 / 60.0 / 1000.0,
        "fe": 1.0 / 60.0 / 60.0 / 1000.0,
        "ka": 3.0 / 60.0 / 60.0,
        "kh": 0.1 / 60.0 / 60.0,
    },
    # onverhard, vlak
    "onv_vla": {
        "outflow_delay": 0.2 / 60.0,
        "surface_storage": 4.0 / 1000.0,
        "infiltration": True,
        "fb": 5.0 / 60.0 / 60.0 / 1000.0,
        "fe": 1.0 / 60.0 / 60.0 / 1000.0,
        "ka": 3.0 / 60.0 / 60.0,
        "kh": 0.1 / 60.0 / 60.0,
    },
    # onverhard, vlak uitgestrekt
    "onv_vlu": {
        "outflow_delay": 0.1 / 60.0,
        "surface_storage": 6.0 / 1000.0,
        "infiltration": True,
        "fb": 5.0 / 60.0 / 60.0 / 1000.0,
        "fe": 1.0 / 60.0 / 60.0 / 1000.0,
        "ka": 3.0 / 60.0 / 60.0,
        "kh": 0.1 / 60.0 / 60.0,
    },
}


SURFACE_CLASS_MAP = dict(
    [
        ("pand", "dak"),
        ("half verhard", "onv"),
        ("onverhard", "onv"),
        ("open verharding", "ovh"),
        ("gesloten verharding", "gvh"),
    ]
)

SURFACE_INCLINATION_MAP = dict(
    [("hellend", "hel"), ("vlak", "vla"), ("uitgestrekt", "vlu")]
)


@dataclass
class SurfaceParams:
    outflow_delay: np.ndarray
    surface_storage: np.ndarray
    infiltration: np.ndarray
    fb: np.ndarray
    fe: np.ndarray
    ka: np.ndarray
    kh: np.ndarray

    @classmethod
    def from_surface_class_and_inclination(
        cls, surface_class: np.ndarray, surface_inclination: np.ndarray
    ) -> "SurfaceParams":
        FIELDS = (
            ("outflow_delay", 0.0, np.float32),
            ("surface_storage", 0.0, np.float64),
            ("infiltration", False, bool),
            ("fb", 0.0, np.float64),
            ("fe", 0.0, np.float64),
            ("ka", 0.0, np.float64),
            ("kh", 0.0, np.float64),
        )

        NO_INFILTRATION_SKIP = ("fb", "fe", "ka", "kh")

        kwargs = {}
        for k, v, dtype in FIELDS:
            kwargs[k] = np.full(surface_class.shape, v, dtype)

        for sc in np.unique(surface_class):
            for si in np.unique(surface_inclination):
                sc_lookup = SURFACE_CLASS_MAP[sc]
                si_lookup = SURFACE_INCLINATION_MAP[si]
                type_props = SURFACE_TYPE_PROPERTIES[sc_lookup + "_" + si_lookup]
                mask = (surface_class == sc) & (surface_inclination == si)

                for k, _, _ in FIELDS:
                    if not type_props["infiltration"] and k in NO_INFILTRATION_SKIP:
                        continue
                    kwargs[k][mask] = type_props[k]

        return SurfaceParams(**kwargs)


class BaseSurface:
    surface_id: int
    code: str
    display_name: str
    area: float
    the_geom: shapely.Geometry

    connection_node_id: int
    connection_node_the_geom: shapely.Geometry
    percentage: float


class Surface(BaseSurface):
    """
    Datamodel based on "surface" table
    """

    id: int
    outflow_delay: float
    surface_layer_thickness: float
    infiltration: bool
    max_infiltration_capacity: float
    min_infiltration_capacity: float
    infiltration_decay_constant: float
    infiltration_recovery_constant: float


def fill_missing_centroids(
    centroids: np.ndarray, surface_ids: np.ndarray, connection_node_the_geom: np.ndarray
):
    """
    Fill any missing centroids with (centroids of the) connection node geometries
    """
    no_centroid_mask = shapely.is_missing(centroids)

    if np.any(no_centroid_mask):
        # Try to fill empty centroids
        empty_surface_ids = surface_ids[no_centroid_mask]
        s_unique = np.unique(empty_surface_ids)
        sort_idx = np.argsort(s_unique)
        lookup = sort_idx[np.searchsorted(s_unique, empty_surface_ids, sorter=sort_idx)]

        computed_centroids = shapely.centroid(
            shapely.multipoints(
                geometries=connection_node_the_geom[no_centroid_mask], indices=lookup
            )
        )

        centroids[no_centroid_mask] = computed_centroids[lookup]


class BaseSurfaces:
    @property
    def unique_surfaces_mask(self) -> np.ndarray:
        if not hasattr(self, "_unique_surfaces_mask"):
            _, self._unique_surfaces_mask = np.unique(
                self.surface_id, return_index=True
            )
        return self._unique_surfaces_mask

    def as_grid_surfaces(
        self, extra_fields: Optional[Dict] = None
    ) -> surfaces.Surfaces:
        if extra_fields is None:
            extra_fields = {}

        centroids = shapely.centroid(self.the_geom)
        no_centroid_mask = shapely.is_missing(centroids)

        if np.any(no_centroid_mask):
            fill_missing_centroids(
                centroids, self.surface_id, self.connection_node_the_geom
            )

        centroid_coords = shapely.get_coordinates(centroids)

        return surfaces.Surfaces(
            id=np.arange(1, len(self.unique_surfaces_mask) + 1, dtype=int),
            code=self.code[self.unique_surfaces_mask],
            display_name=self.display_name[self.unique_surfaces_mask],
            area=self.area[self.unique_surfaces_mask],
            centroid_x=centroid_coords[:, 0][self.unique_surfaces_mask],
            centroid_y=centroid_coords[:, 1][self.unique_surfaces_mask],
            **extra_fields
        )

    def as_surface_maps(
        self, nodes: Nodes, nodes_embedded: Nodes
    ) -> surfaces.SurfaceMaps:
        search_array = self.surface_id[self.unique_surfaces_mask]
        sort_idx = np.argsort(self.surface_id[self.unique_surfaces_mask])
        imp_lookup = sort_idx[
            np.searchsorted(search_array, self.surface_id, sorter=sort_idx)
        ]
        connection_node_coords = shapely.get_coordinates(self.connection_node_the_geom)

        # Find connection_node (calc) node id's both in nodes and embedded nodes
        connection_node_mask = (
            nodes.content_type == ContentType.TYPE_V2_CONNECTION_NODES.value
        )
        embedded_mask = (
            nodes_embedded.content_type == ContentType.TYPE_V2_CONNECTION_NODES.value
        )

        cc_node_ids = np.concatenate(
            (
                nodes.content_pk[connection_node_mask],
                nodes_embedded.content_pk[embedded_mask],
            )
        )
        node_ids = np.concatenate(
            (nodes.id[connection_node_mask], nodes_embedded.embedded_in[embedded_mask])
        )

        sort_idx = np.argsort(cc_node_ids)
        cc_lookup = sort_idx[
            np.searchsorted(cc_node_ids, self.connection_node_id, sorter=sort_idx)
        ]

        return surfaces.SurfaceMaps(
            id=np.arange(1, len(self.connection_node_id) + 1, dtype=int),
            pk=self.connection_node_id,
            nxc=connection_node_coords[:, 0],
            nyc=connection_node_coords[:, 1],
            fac=self.percentage / 100.0,
            imp=imp_lookup,
            cci=node_ids[cc_lookup],
        )


class Surfaces(Array[Surface], BaseSurfaces):
    def as_grid_surfaces(self) -> surfaces.Surfaces:
        _, unique_surfaces_mask = np.unique(self.surface_id, return_index=True)

        infiltration_flag = self.infiltration[unique_surfaces_mask].astype("bool")
        fb = self.max_infiltration_capacity[unique_surfaces_mask] / (
            60.0 * 60.0 * 1000.0
        )
        fe = self.min_infiltration_capacity[unique_surfaces_mask] / (
            60.0 * 60.0 * 1000.0
        )
        ka = self.infiltration_decay_constant[unique_surfaces_mask] / (60.0 * 60.0)
        kh = self.infiltration_recovery_constant[unique_surfaces_mask] / (60.0 * 60.0)

        mask = np.invert(infiltration_flag)

        fb[mask] = 0.0
        fe[mask] = 0.0
        ka[mask] = 0.0
        kh[mask] = 0.0

        extra_fields = dict(
            outflow_delay=self.outflow_delay[unique_surfaces_mask] / 60.0,
            storage_limit=self.surface_layer_thickness[unique_surfaces_mask] / 1000.0,
            infiltration_flag=self.infiltration[unique_surfaces_mask].astype("bool"),
            fb=fb,
            fe=fe,
            ka=ka,
            kh=kh,
        )

        return super().as_grid_surfaces(extra_fields)
