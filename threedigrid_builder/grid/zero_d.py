from dataclasses import dataclass
from threedi_modelchecker.threedi_model.models import ImperviousSurface
from threedigrid_builder.base import array_of
from threedigrid_builder.constants import InfiltrationSurfaceOption, SurfaceClass
from threedigrid_builder.base.surfaces import Surfaces as BaseSurfaces
from threedigrid_builder.constants import ContentType

import pygeos
import numpy as np

__all__ = ["Surfaces", "ImperviousSurfaces"]


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
    'gvh_hel': {
        'outflow_delay': 0.5 / 60.0,
        'surface_storage': 0.0 / 1000.0,
        'infiltration': False,
        'fb': 0.0 / 60.0 / 60.0 / 1000.0,   # for consistency
        'fe': 0.0 / 60.0 / 60.0 / 1000.0,
        'ka': 0.0 / 60.0 / 60.0,
        'kh': 0.0 / 60.0 / 60.0
    },
    # gesloten verharding, vlak
    'gvh_vla': {
        'outflow_delay': 0.2 / 60.0,
        'surface_storage': 0.5 / 1000.0,
        'infiltration': False,
        'fb': 0.0 / 60.0 / 60.0 / 1000.0,
        'fe': 0.0 / 60.0 / 60.0 / 1000.0,
        'ka': 0.0 / 60.0 / 60.0,
        'kh': 0.0 / 60.0 / 60.0
    },
    # gesloten verharding, vlak uitgestrekt
    'gvh_vlu': {
        'outflow_delay': 0.1 / 60.0,
        'surface_storage': 1.0 / 1000.0,
        'infiltration': False,
        'fb': 0.0 / 60.0 / 60.0 / 1000.0,
        'fe': 0.0 / 60.0 / 60.0 / 1000.0,
        'ka': 0.0 / 60.0 / 60.0,
        'kh': 0.0 / 60.0 / 60.0
    },
    # open verharding, hellend
    'ovh_hel': {
        'outflow_delay': 0.5 / 60.0,
        'surface_storage': 0.0 / 1000.0,
        'infiltration': True,
        'fb': 2.0 / 60.0 / 60.0 / 1000.0,
        'fe': 0.5 / 60.0 / 60.0 / 1000.0,
        'ka': 3.0 / 60.0 / 60.0,
        'kh': 0.1 / 60.0 / 60.0
    },
    # open verharding, vlak
    'ovh_vla': {
        'outflow_delay': 0.2 / 60.0,
        'surface_storage': 0.5 / 1000.0,
        'infiltration': True,
        'fb': 2.0 / 60.0 / 60.0 / 1000.0,
        'fe': 0.5 / 60.0 / 60.0 / 1000.0,
        'ka': 3.0 / 60.0 / 60.0,
        'kh': 0.1 / 60.0 / 60.0
    },
    # open verharding, vlak uitgestrekt
    'ovh_vlu': {
        'outflow_delay': 0.1 / 60.0,
        'surface_storage': 1.0 / 1000.0,
        'infiltration': True,
        'fb': 2.0 / 60.0 / 60.0 / 1000.0,
        'fe': 0.5 / 60.0 / 60.0 / 1000.0,
        'ka': 3.0 / 60.0 / 60.0,
        'kh': 0.1 / 60.0 / 60.0
    },
    # dak, hellend
    'dak_hel': {
        'outflow_delay': 0.5 / 60.0,
        'surface_storage': 0.0 / 1000.0,
        'infiltration': False,
        'fb': 0.0 / 60.0 / 60.0 / 1000.0,
        'fe': 0.0 / 60.0 / 60.0 / 1000.0,
        'ka': 0.0 / 60.0 / 60.0,
        'kh': 0.0 / 60.0 / 60.0
    },
    # dak, vlak
    'dak_vla': {
        'outflow_delay': 0.2 / 60.0,
        'surface_storage': 2.0 / 1000.0,
        'infiltration': False,
        'fb': 0.0 / 60.0 / 60.0 / 1000.0,
        'fe': 0.0 / 60.0 / 60.0 / 1000.0,
        'ka': 0.0 / 60.0 / 60.0,
        'kh': 0.0 / 60.0 / 60.0
    },
    # dak, vlak uitgestrekt
    'dak_vlu': {
        'outflow_delay': 0.1 / 60.0,
        'surface_storage': 4.0 / 1000.0,
        'infiltration': False,
        'fb': 0.0 / 60.0 / 60.0 / 1000.0,
        'fe': 0.0 / 60.0 / 60.0 / 1000.0,
        'ka': 0.0 / 60.0 / 60.0,
        'kh': 0.0 / 60.0 / 60.0
    },
    # onverhard, hellend
    'onv_hel': {
        'outflow_delay': 0.5 / 60.0,
        'surface_storage': 2.0 / 1000.0,
        'infiltration': True,
        'fb': 5.0 / 60.0 / 60.0 / 1000.0,
        'fe': 1.0 / 60.0 / 60.0 / 1000.0,
        'ka': 3.0 / 60.0 / 60.0,
        'kh': 0.1 / 60.0 / 60.0
    },
    # onverhard, vlak
    'onv_vla': {
        'outflow_delay': 0.2 / 60.0,
        'surface_storage': 4.0 / 1000.0,
        'infiltration': True,
        'fb': 5.0 / 60.0 / 60.0 / 1000.0,
        'fe': 1.0 / 60.0 / 60.0 / 1000.0,
        'ka': 3.0 / 60.0 / 60.0,
        'kh': 0.1 / 60.0 / 60.0
    },
    # onverhard, vlak uitgestrekt
    'onv_vlu': {
        'outflow_delay': 0.1 / 60.0,
        'surface_storage': 6.0 / 1000.0,
        'infiltration': True,
        'fb': 5.0 / 60.0 / 60.0 / 1000.0,
        'fe': 1.0 / 60.0 / 60.0 / 1000.0,
        'ka': 3.0 / 60.0 / 60.0,
        'kh': 0.1 / 60.0 / 60.0
    }
}


SURFACE_CLASS_MAP = dict([
    ("pand", "dak"),
    ("half_verhard", "onv"),
    ("onverhard", "onv"),
    ("open verharding", "ovh"),
    ("gesloten verharding", "gvh")   
])

SURFACE_INCLINATION_MAP = dict([
    ("hellend", "hel"),
    ("vlak", "vla"),
    ("uitgestrekt", "vlu")
])


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
    def from_surface_class_and_inclination(cls, surface_class: np.ndarray, surface_inclination: np.ndarray) -> "SurfaceParams":
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
        for k,v, dtype in FIELDS:
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


class Surface:
    id: int
    function: str
    code: str
    display_name: str
    nr_of_inhabitants: float
    area: float
    dry_weather_flow: float
    the_geom: pygeos.Geometry

    # Surface Parameters   
    outflow_delay: float
    surface_layer_thickness: float
    infiltration: float
    max_infiltration_capacity: float
    min_infiltration_capacity: float
    infiltration_decay_constant: float
    infiltration_recovery_constant: float


class ImperviousSurface:
    id: int
    code: str
    display_name: str
    surface_inclination: str
    surface_class: SurfaceClass
    surface_sub_class: str
    nr_of_inhabitants: float
    area: float
    dry_weather_flow: float
    the_geom: pygeos.Geometry
    connection_node_id: int
    connection_node_the_geom: pygeos.Geometry
    percentage: float


@array_of(Surface)
class Surfaces:
    def apply(self, grid):
        pass    


@array_of(ImperviousSurface)
class ImperviousSurfaces:
    def apply(self, grid):
        centroids = pygeos.get_coordinates(pygeos.centroid(self.the_geom))
        connection_node_coords = pygeos.get_coordinates(self.connection_node_the_geom)
         
        grid.surfaces.id = np.arange(1, len(self.id) + 1, dtype=int)
        grid.surfaces.centroid_x = centroids[:, 0]
        grid.surfaces.centroid_y = centroids[:, 1]
        grid.surfaces.code = self.code
        grid.surfaces.display_name = self.display_name
        grid.surfaces.surface_inclination = self.surface_inclination
        grid.surfaces.surface_class = self.surface_class
        grid.surfaces.surface_sub_class = self.surface_sub_class
        grid.surfaces.nr_of_inhabitants = self.nr_of_inhabitants
        grid.surfaces.area = self.area
        grid.surfaces.dry_weather_flow = self.dry_weather_flow

        params = SurfaceParams.from_surface_class_and_inclination(
            self.surface_class, self.surface_inclination)

        grid.surfaces.outflow_delay = params.outflow_delay
        grid.surfaces.storage_limit = params.surface_storage
        grid.surfaces.infiltration_flag = params.infiltration
        grid.surfaces.fb = params.fb
        grid.surfaces.fe = params.fe
        grid.surfaces.ka = params.ka
        grid.surfaces.kh = params.kh

        grid.surfaces.cid = self.connection_node_id
        grid.surfaces.pk = self.connection_node_id
        grid.surfaces.nxc = connection_node_coords[:, 0]
        grid.surfaces.nyc = connection_node_coords[:, 1]
        grid.surfaces.fac = self.percentage / 100.0
        grid.surfaces.imp = grid.surfaces.id

        grid.surfaces.function = np.full(self.id.shape, b"")

        # Find connection_node (calc) node id's
        connection_node_mask = grid.nodes.content_type == ContentType.TYPE_V2_CONNECTION_NODES.value
        cc_node_ids = grid.nodes.content_pk[connection_node_mask]
        node_ids = grid.nodes.id[connection_node_mask]

        sort_idx = np.argsort(cc_node_ids)
        lookup = sort_idx[np.searchsorted(
            cc_node_ids, self.connection_node_id, sorter=sort_idx)]

        grid.surfaces.cci = node_ids[lookup]
