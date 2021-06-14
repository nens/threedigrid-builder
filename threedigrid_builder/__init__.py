from .application import *  # NOQA
from .exceptions import *  # NOQA

from pathlib import Path


_version_path = Path(__file__).parent.parent / "VERSION.rst"
try:
    __version__ = open(_version_path, "r").read().strip()
except FileNotFoundError:
    __version__ = "UNKNOWN"
del _version_path
