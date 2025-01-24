import os
import pathlib
import shutil
import sys

from setuptools import find_packages

try:
    from skbuild import setup
except ImportError:
    from setuptools import setup
    if not any(x in sys.argv for x in {"sdist", "--version", "egg_info"}):
        sys.exit("ERROR: skbuild is required to build threedigrid-builder from source.")

ext_modules = []

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
    sys.exit()

if sys.platform == "win32":
    cmake_args = ["-G", "MSYS Makefiles", "-DCMAKE_GNUtoMS=ON"]
else:
    cmake_args = []

long_description = open("README.rst").read()


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
    "numpy>=1.15,<3.0",
    "threedi-schema==0.300.*",
    "shapely>=2",
    "pyproj>=3",
    "condenser[geo]>=0.1.1",
    "sqlalchemy>=1.4.1",
]

test_requires = ["pytest", "pytest-cov"]

setup(
    name="threedigrid-builder",
    version=get_version(),
    description="Generate a 3Di simulation grid from a model schematisation.",
    long_description=long_description,
    url="https://docs.3di.lizard.net/",
    author="Martijn Siemerink",
    author_email="martijn.siemerink@nelen-schuurmans.nl",
    license="GNU General Public License v3.0",
    packages=find_packages(
        include=(
            "threedigrid_builder",
            "threedigrid_builder.*",
        ),
    ),
    cmake_args=cmake_args,
    install_requires=install_requires,
    extras_require={
        "test": test_requires,
        "gridadmin": ["h5py>=2.7"],
        "gpkg": ["geopandas"],
        "cli": ["typer"],
    },
    python_requires=">=3.9",
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
    entry_points={
        "console_scripts": ["threedigrid-builder=threedigrid_builder.cli:run"],
    },
)
