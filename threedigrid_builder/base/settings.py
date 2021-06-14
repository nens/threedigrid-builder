from dataclasses import dataclass
from dataclasses import fields
from threedigrid_builder.constants import InitializationType
from typing import Optional


__all__ = ["GridSettings", "TablesSettings"]


@dataclass
class GridSettings:
    """Settings necessary for threedigrid-builder."""

    ## from GlobalSettings
    use_2d: bool
    use_1d_flow: bool
    use_2d_flow: bool
    grid_space: float
    dist_calc_points: float
    kmax: int
    embedded_cutoff_threshold: float = 0.05
    max_angle_1d_advection: float = 90.0

    @classmethod
    def from_dict(cls, dct):
        """Construct skipping unknown fields and None values"""
        class_fields = {f.name for f in fields(cls)}
        return cls(
            **{k: v for k, v in dct.items() if k in class_fields and v is not None}
        )


@dataclass
class TablesSettings:
    """Settings necessary for threedi-tables."""

    ## from GlobalSettings
    table_step_size: float
    frict_coef: float
    frict_coef_type: InitializationType
    frict_type: int = 4
    interception_global: Optional[float] = None
    interception_type: Optional[InitializationType] = None
    table_step_size_1d: float = None  # actual default is set in __post_init__
    table_step_size_volume_2d: float = None  # actual default  is set in __post_init__

    # TODO --> https://github.com/nens/threedigrid-builder/issues/86
    manhole_storage_area: Optional[float] = None

    # obstacle detection could be a tool in QGis?
    # dem_obstacle_detection: bool
    # dem_obstacle_height: float

    ## from Groundwater
    groundwater_impervious_layer_level: Optional[float] = None
    groundwater_impervious_layer_level_type: Optional[InitializationType] = None
    phreatic_storage_capacity: Optional[float] = None
    phreatic_storage_capacity_type: Optional[InitializationType] = None
    equilibrium_infiltration_rate: Optional[float] = None
    equilibrium_infiltration_rate_type: Optional[InitializationType] = None
    initial_infiltration_rate: Optional[float] = None
    initial_infiltration_rate_type: Optional[InitializationType] = None
    infiltration_decay_period: Optional[float] = None
    infiltration_decay_period_type: Optional[InitializationType] = None
    groundwater_hydro_connectivity: Optional[float] = None
    groundwater_hydro_connectivity_type: Optional[InitializationType] = None

    ## from Interflow
    interflow_type: int = 0  # InterflowType (0 means: no interflow)
    porosity: Optional[float] = None
    porosity_type: Optional[InitializationType] = None
    porosity_layer_thickness: Optional[float] = None
    impervious_layer_elevation: Optional[float] = None
    hydraulic_conductivity: Optional[float] = None
    hydraulic_conductivity_type: Optional[InitializationType] = None

    ## from SimpleInfiltration
    infiltration_rate: Optional[float] = None
    infiltration_rate_type: Optional[InitializationType] = None
    infiltration_surface_option: Optional[int] = None  # InfiltrationSurfaceOption
    max_infiltration_capacity_type: Optional[InitializationType] = None

    def __post_init__(self):
        # set some defaults
        if self.table_step_size_1d is None:
            self.table_step_size_1d = self.table_step_size
        if self.table_step_size_volume_2d is None:
            self.table_step_size_volume_2d = self.table_step_size

    @classmethod
    def from_dict(cls, dct):
        """Construct skipping unknown fields and None values"""
        class_fields = {f.name for f in fields(cls)}
        return cls(
            **{k: v for k, v in dct.items() if k in class_fields and v is not None}
        )
