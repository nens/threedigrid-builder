
from setuptools import find_packages

import pathlib
import sys
import shutil
import os
import numpy

try:
    from skbuild import setup
except ImportError:
    if not any(x in sys.argv for x in {"sdist", "--version", "egg_info"}):
        sys.exit("ERROR: Numpy is required to build threedi-tables from source.")
    # stuff like "python setup.py --version" should be allowed without Numpy:
    from setuptools import Extension
    from setuptools import setup

ext_modules = []
build_type = "Release"

if "clean" in sys.argv:
    # delete any previously compiled files
    sys.argv.remove("clean")
    p = pathlib.Path(".")
    for pattern in [
        "libthreedigrid/*.c",
        "libthreedigrid/*.f",
        "libthreedigrid/*.mod",
        "threedigrid_builder/grid/fgrid/*.so",
    ]:
        for filename in p.glob(pattern):
            print("removing '{}'".format(filename))
            filename.unlink()
    print("removing build folder")
    if os.path.isdir(p / "_skbuild"):
        shutil.rmtree(p / "_skbuild")

if "debug" in sys.argv:
    # Set debug compiler flags of fortran lib when required.
    sys.argv.remove("debug")
    build_type = "Debug"
    macro = [("F2PY_REPORT_ON_ARRAY_COPY", "2")]
else:
    # Otherwise set normal compiler flags.
    macro = []

long_description = "\n\n".join([open("README.rst").read(), open("CHANGES.rst").read()])


def get_version():
    # Edited from https://packaging.python.org/guides/single-sourcing-package-version/
    init_path = pathlib.Path(__file__).parent / "threedigrid_builder/__init__.py"
    for line in init_path.open("r").readlines():
        if line.startswith("__version__"):
            delim = '"' if '"' in line else "'"
            return line.split(delim)[1]
    else:
        raise RuntimeError("Unable to find version string.")


install_requires = [
    "numpy>=1.15",
    "threedi-modelchecker>=0.12",
    "pygeos>=0.11.1",
    "pyproj>=3",
    "condenser[geo]>=0.1.1",
    "sqlalchemy",
]

test_requires = ["pytest"]

setup(
    name="threedigrid-builder",
    version=get_version(),
    description="Generate a 3Di simulation grid from a model schematisation.",
    long_description=long_description,
    url="https://docs.3di.lizard.net/",
    author="Martijn Siemerink",
    author_email="martijn.siemerink@nelen-schuurmans.nl",
    license="Proprietary",
    packages=find_packages(
        include=(
            "threedigrid_builder",
            "threedigrid_builder.*",
        ),
    ),
    install_requires=install_requires,
    extras_require={
        "test": test_requires,
        "gridadmin": ["h5py>=2.7"],
        "gpkg": ["geopandas"],
    },
    python_requires=">=3.7",
    include_package_data=True,
    classifiers=[
        "Programming Language :: Python :: 3",
        "Intended Audience :: Science/Research",
        "Intended Audience :: Developers",
        "Development Status :: 5 - Production/Stable",
        "Topic :: Scientific/Engineering",
        "Topic :: Software Development",
        "Operating System :: Unix",
        "Operating System :: MacOS",
        "Operating System :: Microsoft :: Windows",
        "License :: Other/Proprietary License",
    ],
    zip_safe=False,
)
