from .array import array_of
import numpy as np


__all__ = ["Surfaces"]


class Surface:
    id: int
    code: str
    display_name: str
    function: str
    area: float
    centroid_x: float
    centroid_y: float
    dry_weather_flow: float
    nr_of_inhabitants: float
    infiltration_flag: bool
    outflow_delay: float
    storage_limit: float
    fac: float
    fb: float
    fe: float
    imp: int
    ka: float
    kh: float
    nxc: float
    nyc: float
    pk: int
    cci: int
    cid: int
    surface_class: str
    surface_inclination: str
    surface_sub_class: str


@array_of(Surface)
class Surfaces:
    """Zero-d surfaces."""

    def get_extent(self):
        if self.centroid_x is None or self.centroid_y is None:
            return None

        extent = (
            np.amin(self.centroid_x),
            np.amin(self.centroid_y),
            np.amax(self.centroid_x),
            np.amax(self.centroid_y),
        )
        if any(np.isnan(val) for val in extent):
            raise ValueError("Not all surface centroids have coordinates.")
        return extent
