from dataclasses import dataclass
from dataclasses import fields
from threedigrid_builder.constants import InitializationType


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
    frict_type: int
    frict_coef: float
    frict_coef_type: InitializationType
    interception_global: float
    interception_type: InitializationType
    table_step_size_1d: float = None  # default is set in __post_init__
    table_step_size_volume_2d: float = None  # default is set in __post_init__

    # TODO --> https://github.com/nens/threedigrid-builder/issues/86
    manhole_storage_area: float = 0.0

    # obstacle detection could be a tool in QGis?
    # dem_obstacle_detection: bool
    # dem_obstacle_height: float

    ## from Groundwater
    groundwater_impervious_layer_level: float = None
    groundwater_impervious_layer_level_type: InitializationType = (
        InitializationType.NONE
    )
    phreatic_storage_capacity: float = None
    phreatic_storage_capacity_type: InitializationType = InitializationType.NONE
    equilibrium_infiltration_rate: float = None
    equilibrium_infiltration_rate_type: InitializationType = InitializationType.NONE
    initial_infiltration_rate: float = None
    initial_infiltration_rate_type: InitializationType = InitializationType.NONE
    infiltration_decay_period: float = None
    infiltration_decay_period_type: InitializationType = InitializationType.NONE
    groundwater_hydro_connectivity: float = None
    groundwater_hydro_connectivity_type: InitializationType = InitializationType.NONE

    ## from Interflow
    interflow_type: int = None  # InterflowType
    porosity: float = None
    porosity_type: InitializationType = InitializationType.NONE
    porosity_layer_thickness: float = None
    impervious_layer_elevation: float = None
    hydraulic_conductivity: float = None
    hydraulic_conductivity_type: InitializationType = InitializationType.NONE

    ## from SimpleInfiltration
    infiltration_rate: float = None
    infiltration_rate_type: InitializationType = InitializationType.NONE
    infiltration_surface_option: int = None  # InfiltrationSurfaceOption
    max_infiltration_capacity_type: InitializationType = InitializationType.NONE

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
