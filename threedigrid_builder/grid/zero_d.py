from threedi_modelchecker.threedi_model.models import ImperviousSurface
from threedigrid_builder.base import array_of
from threedigrid_builder.constants import SurfaceClass
from threedigrid_builder.base.surfaces import Surfaces as BaseSurfaces

import pygeos
import numpy as np

__all__ = ["Surfaces", "ImperviousSurfaces"]


_NONE_TO_DEFAULT_BASE_FIELDS = {
    "function": "",
    "code": "",
    "display_name": "",
    "zoom_category": 0,
    "nr_of_inhabitants": 0.0,
    "area": 0.0,
    "dry_weather_flow": 0.0,

    # sewerage
    "surface_class": "",
    "surface_sub_class": "",
    "surface_inclination": "",
}


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


@array_of(Surface)
class Surfaces:
    def apply(self, grid):
        pass    


@array_of(ImperviousSurface)
class ImperviousSurfaces:
    def apply(self, grid):
        centroids = pygeos.get_coordinates(pygeos.centroid(self.the_geom))
        centroid_x, centroid_y = centroids[:, 0], centroids[:, 1]

        grid.surfaces.id = self.id
        grid.surfaces.centroid_x = centroid_x
        grid.surfaces.centroid_y = centroid_y
        grid.surfaces.code = self.code
        grid.surfaces.display_name = self.display_name
        grid.surfaces.surface_inclination = self.surface_inclination
        grid.surfaces.surface_class = self.surface_class
        grid.surfaces.surface_sub_class = self.surface_sub_class
        grid.surfaces.nr_of_inhabitants = self.nr_of_inhabitants
        grid.surfaces.area = self.area
        grid.surfaces.dry_weather_flow = self.dry_weather_flow


        _SEWERAGE_M_DEF = np.dtype(
            [
                # ("id", np.int32),
                # ("surface_class", np.dtype("S128")),
                # ("code", np.dtype("S100")),
                #("display_name", np.dtype("S255")),
                #("surface_sub_class", np.dtype("S128")),
                #("nr_of_inhabitants", np.float64),
                #("area", np.float64),
                #("dry_weather_flow", np.float32),
                #("surface_inclination", np.dtype("S64")),
                #("centroid_x", np.float64),
                #("centroid_y", np.float64),
                ("outflow_delay", np.float32),
                ("storage_limit", np.float64),
                ("infiltration_flag", np.bool_),
                ("fb", np.float64),
                ("fe", np.float64),
                ("ka", np.float64),
                ("kh", np.float64),
            ]
        )


        pass