"""The interface layer (data conversion layer) of the project

A set of adapters that convert data from the external (HTTP, Database) to the internal
structure.

No code outside of this layer should contain details of a specific external data
structure.

The interface layer depends on the domain layer only.
"""

from .db import *  # NOQA
from .dict_out import *  # NOQA
from .geopackage import *  # NOQA
from .gridadmin import *  # NOQA
from .raster_gdal import *  # NOQA
