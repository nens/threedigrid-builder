from setuptools import find_packages

import os
import pathlib
import setuptools  # noqa
import shutil
import sys


try:
    from numpy.distutils.core import Extension
    from numpy.distutils.core import setup
except ImportError:
    if not any(x in sys.argv for x in {"sdist", "--version", "egg_info"}):
        sys.exit("ERROR: Numpy is required to build threedigrid-builder from source.")
    # stuff like "python setup.py --version" should be allowed without Numpy:
    from setuptools import Extension
    from setuptools import setup

def get_version():
    # Edited from https://packaging.python.org/guides/single-sourcing-package-version/
    init_path = pathlib.Path(__file__).parent / "threedigrid_builder/__init__.py"
    for line in init_path.open("r").readlines():
        if line.startswith("__version__"):
            delim = '"' if '"' in line else "'"
            return line.split(delim)[1]
    else:
        raise RuntimeError("Unable to find version string.")


long_description = "\n\n".join([open("README.rst").read(), open("CHANGES.rst").read()])


install_requires = [
    "numpy>=1.13",
    "threedi-modelchecker>=0.12",
    "pygeos>=0.10",
    "pyproj>=3",
    "condenser[geo]>=0.1.1",
    "sqlalchemy",
    "dataclasses ; python_version<'3.7'",
]

test_requires = ["pytest"]


if "clean" in sys.argv:
    # delete any previously compiled files
    sys.argv.remove("clean")
    p = pathlib.Path(".")
    for pattern in [
        "libthreedigrid/*.c",
        "libthreedigrid/*.f",
        "libthreedigrid/*.mod",
        "threedigrid_builder/*.so",
    ]:
        for filename in p.glob(pattern):
            print("removing '{}'".format(filename))
            filename.unlink()
    print("removing build folder")
    if os.path.isdir(p / "build"):
        shutil.rmtree(p / "build")

if "debug" in sys.argv:
    # Set debug compiler flags of fortran lib when required.
    sys.argv.remove("debug")
    comp_flags = ["-O0", "-g"]
    macro = [("F2PY_REPORT_ON_ARRAY_COPY", "2")]
else:
    # Otherwise set normal compiler flags.
    comp_flags = ["-O3"]
    macro = []

ext_modules = [
    Extension(
        name="threedigrid_builder.grid._fwrapper",
        sources=[
            "./libthreedigrid/wrapper.pyf",
            "./libthreedigrid/parameters.f90",
            "./libthreedigrid/array_utils.f90",
            "./libthreedigrid/geo_utils.f90",
            "./libthreedigrid/cells.f90",
            "./libthreedigrid/quadtree.f90",
        ],
        f2py_options=["f2cmap", "./libthreedigrid/.f2py_f2cmap"],
        extra_f90_compile_args=comp_flags,
        define_macros=macro + [("NPY_NO_DEPRECATED_API", "NPY_1_7_API_VERSION")],
    )
]


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
        "rasters": ["rasterio"],
        "gridadmin": ["h5py>=2.7"],
        "gpkg": ["geopandas"],
    },
    python_requires=">=3.7",
    include_package_data=True,
    ext_modules=ext_modules,
    classifiers=[
        "Programming Language :: Python :: 3",
        "Intended Audience :: Science/Research",
        "Intended Audience :: Developers",
        "Development Status :: 4 - Beta",
        "Topic :: Scientific/Engineering",
        "Topic :: Software Development",
        "Operating System :: Unix",
        "Operating System :: MacOS",
        "Operating System :: Microsoft :: Windows",
        "License :: Other/Proprietary License",
    ],
    zip_safe=False,
)
