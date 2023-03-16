from dataclasses import dataclass
from dataclasses import fields
from threedigrid_builder.base import is_int_enum
from threedigrid_builder.base import unpack_optional_type
from threedigrid_builder.constants import FrictionType
from threedigrid_builder.constants import InfiltrationSurfaceOption
from threedigrid_builder.constants import InitializationType
from threedigrid_builder.constants import InterflowType
from threedigrid_builder.exceptions import SchematisationError
from typing import Optional

import numpy as np


__all__ = ["GridSettings", "TablesSettings"]


def greater_zero_check(obj, attr):
    value = getattr(obj, attr, None)
    if value is not None and value <= 0:
        raise SchematisationError(f"'{attr}' must be greater than 0.")


@dataclass
class GridSettings:
    """Settings necessary for threedigrid-builder."""

    ## from GlobalSettings
    use_2d: bool
    use_1d_flow: bool
    use_2d_flow: bool
    use_0d_inflow: int
    grid_space: float
    dist_calc_points: float
    kmax: int
    embedded_cutoff_threshold: float = 0.05
    max_angle_1d_advection: float = 0.4 * np.pi

    @classmethod
    def from_dict(cls, dct):
        """Construct skipping unknown fields and None values"""
        class_fields = {f.name for f in fields(cls)}
        return cls(
            **{k: v for k, v in dct.items() if k in class_fields and v is not None}
        )

    def __post_init__(self):
        # validations
        if self.use_2d:
            for field in ["grid_space", "kmax"]:
                greater_zero_check(self, field)


@dataclass
class TablesSettings:
    """Settings necessary for threedi-tables."""

    ## from GlobalSettings
    table_step_size: float
    frict_coef: float
    frict_coef_type: InitializationType
    frict_type: FrictionType = FrictionType.MANNING
    interception_global: Optional[float] = None
    interception_type: Optional[InitializationType] = None
    table_step_size_1d: float = None  # actual default is set in __post_init__
    maximum_table_step_size: float = None  # actual default  is set in __post_init__

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
    interflow_type: InterflowType = InterflowType.NO_INTERLFOW
    porosity: Optional[float] = None
    porosity_type: Optional[InitializationType] = None
    porosity_layer_thickness: Optional[float] = None
    impervious_layer_elevation: Optional[float] = None
    hydraulic_conductivity: Optional[float] = None
    hydraulic_conductivity_type: Optional[InitializationType] = None

    ## from SimpleInfiltration
    infiltration_rate: Optional[float] = None
    infiltration_rate_type: Optional[InitializationType] = None
    infiltration_surface_option: InfiltrationSurfaceOption = (
        InfiltrationSurfaceOption.RAIN
    )
    max_infiltration_capacity: Optional[float] = None
    max_infiltration_capacity_type: Optional[InitializationType] = None

    ## from VegetationDrag
    vegetation_height: Optional[float] = None
    vegetation_height_type: Optional[InitializationType] = None
    vegetation_stem_count: Optional[float] = None
    vegetation_stem_count_type: Optional[InitializationType] = None
    vegetation_stem_diameter: Optional[float] = None
    vegetation_stem_diameter_type: Optional[InitializationType] = None
    vegetation_drag_coefficient: Optional[float] = None
    vegetation_drag_coefficient_type: Optional[InitializationType] = None

    def __post_init__(self):
        # defaults
        if self.table_step_size_1d is None:
            self.table_step_size_1d = self.table_step_size
        if self.maximum_table_step_size is None:
            self.maximum_table_step_size = 100 * self.table_step_size

        # validations
        for field in (
            "table_step_size",
            "table_step_size_1d",
            "maximum_table_step_size",
        ):
            greater_zero_check(self, field)

        if self.maximum_table_step_size < self.table_step_size:
            raise SchematisationError(
                f"'maximum_table_step_size' must not be less than 'table_step_size'."
            )

        # check enums
        for (name, elem_type) in self.__annotations__.items():
            elem_type = unpack_optional_type(elem_type)
            if is_int_enum(elem_type):
                value = getattr(self, name)
                if value is None:
                    continue
                try:
                    elem_type(value)
                except ValueError as e:
                    raise SchematisationError(str(e))

    @classmethod
    def from_dict(cls, dct):
        """Construct skipping unknown fields and None values"""
        class_fields = {f.name for f in fields(cls)}
        return cls(
            **{k: v for k, v in dct.items() if k in class_fields and v is not None}
        )
