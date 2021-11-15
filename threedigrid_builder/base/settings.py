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


def greater_zero_check(obj, attr, allow_zero=False):
    value = getattr(obj, attr, None)
    if value is None:
        return
    if allow_zero:
        if value < 0:
            raise SchematisationError(f"'{attr}' must be greater than or equal to 0.")
    else:
        if value <= 0:
            raise SchematisationError(f"'{attr}' must be greater than 0.")


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
        for field in ("grid_space", "dist_calc_points", "kmax"):
            greater_zero_check(self, field, allow_zero=False)
        for field in ("embedded_cutoff_threshold", "max_angle_1d_advection"):
            greater_zero_check(self, field, allow_zero=True)
        if self.max_angle_1d_advection > 0.5 * np.pi:
            raise SchematisationError(
                "'max_angle_1d_advection' must be less than 0.5 * pi."
            )


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
    infiltration_surface_option: Optional[InfiltrationSurfaceOption] = None
    max_infiltration_capacity_type: Optional[InitializationType] = None

    def __post_init__(self):
        # defaults
        if self.table_step_size_1d is None:
            self.table_step_size_1d = self.table_step_size
        if self.table_step_size_volume_2d is None:
            self.table_step_size_volume_2d = self.table_step_size

        # validations
        for field in (
            "table_step_size",
            "table_step_size_1d",
            "table_step_size_volume_2d",
            "equilibrium_infiltration_rate",
            "initial_infiltration_rate",
            "infiltration_decay_period",
            "groundwater_hydro_connectivity",
            "porosity_layer_thickness",
            "hydraulic_conductivity",
        ):
            greater_zero_check(self, field, allow_zero=False)
        for field in (
            "frict_coef",
            "interception_global",
            "manhole_storage_area",
            "phreatic_storage_capacity",
            "infiltration_rate",
            "porosity",
        ):
            greater_zero_check(self, field, allow_zero=True)
        if (
            self.phreatic_storage_capacity is not None
            and self.phreatic_storage_capacity > 1
        ):
            raise SchematisationError(
                "'phreatic_storage_capacity' must be less than or equal to 1."
            )
        if self.porosity is not None and self.porosity > 1:
            raise SchematisationError("'porosity' must be less than or equal to 1.")

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
