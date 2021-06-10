from dataclasses import dataclass, fields

__all__ = ["MakeGridSettings", "MakeTablesSettings"]


@dataclass
class MakeGridSettings:
    """Settings necessary for threedigrid-builder.
    """
    ## from GlobalSettings
    grid_space: float
    dist_calc_points: float
    kmax: int
    embedded_cutoff_threshold: float = 0.05
    max_angle_1d_advection: float = 90.0

    @classmethod
    def from_dict(cls, dct):
        """Construct skipping unknown fields and None values"""
        class_fields = {f.name for f in fields(cls)}
        return cls(**{k: v for k, v in dct.items() if k in class_fields and v is not None})


@dataclass
class MakeTablesSettings:
    """Settings necessary for threedi-tables.
    """
    ## from GlobalSettings
    table_step_size: float
    frict_type: int
    frict_coef: float
    frict_avg: float
    interception_global: float = None
    table_step_size_1d: float = None  # default is set in __post_init__
    table_step_size_volume_2d: float = None  # default is set in __post_init__

    # obstacle detection could be a tool in QGis?
    # dem_obstacle_detection: bool
    # dem_obstacle_height: float

    ## from Groundwater
    groundwater_impervious_layer_level: float = None
    groundwater_impervious_layer_level_type: int = None  # InitializationType
    phreatic_storage_capacity: float = None
    phreatic_storage_capacity_type: int = None  # InitializationType
    equilibrium_infiltration_rate: float = None
    equilibrium_infiltration_rate_type: int = None  # InitializationType
    initial_infiltration_rate: float = None
    initial_infiltration_rate_type: int = None  # InitializationType
    infiltration_decay_period: float = None
    infiltration_decay_period_type: int = None  # InitializationType
    groundwater_hydro_connectivity: float = None
    groundwater_hydro_connectivity_type: int = None  # InitializationType

    ## from Interflow
    interflow_type: int = None  # InterflowType
    porosity: float = None
    porosity_layer_thickness: float = None
    impervious_layer_elevation: float = None
    hydraulic_conductivity: float = None

    ## from SimpleInfiltration
    infiltration_rate: float = None
    infiltration_surface_option: int = None  # InfiltrationSurfaceOption

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
        return cls(**{k: v for k, v in dct.items() if k in class_fields and v is not None})
