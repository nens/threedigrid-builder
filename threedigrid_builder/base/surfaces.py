import numpy as np

from .array import Array

__all__ = ["Surfaces", "SurfaceMaps"]


class Surface:
    # Fields for every surface
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
    fb: float  # max_infiltration_capacity
    fe: float  # min_infiltration_capacity
    ka: float  # infiltration_decay_constant
    kh: float  # infiltration_recovery_constant

    # Only filled in for impervious surfaces
    surface_class: str
    surface_inclination: str
    surface_sub_class: str


class SurfaceMap:
    # Fields for mapping surface to connection node mapping
    # (surfaces can be connected to multiple connection_nodes)
    id: int
    fac: float  # Fraction of water on surface that is going to connection-node (0-1)
    imp: int  # 0-based index to surface array's (see Surface above)
    nxc: float  # Connection node x coordinate
    nyc: float  # Connection node y coordinate
    pk: int  # connection_node id
    cci: int  # node id


class Surfaces(Array[Surface]):
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


class SurfaceMaps(Array[SurfaceMap]):
    """Zero-d surfaces maps."""

    pass
