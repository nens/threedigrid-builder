from setuptools import Extension
from setuptools import find_packages
from setuptools import setup
from setuptools.command.build_ext import build_ext as _build_ext

import builtins
import pathlib
import sys


# Skip Cython build if not available (for source distributions)
try:
    from Cython.Build import cythonize
except ImportError:
    cythonize = None


# Add numpy include dirs without importing numpy on module level.
# derived from scikit-hep:
# https://github.com/scikit-hep/root_numpy/pull/292
class build_ext(_build_ext):
    def finalize_options(self):
        _build_ext.finalize_options(self)

        # Prevent numpy from thinking it is still in its setup process:
        try:
            del builtins.__NUMPY_SETUP__
        except AttributeError:
            pass

        import numpy

        self.include_dirs.append(numpy.get_include())


ext_modules = []

if "clean" in sys.argv:
    # delete any previously Cythonized or compiled files
    p = pathlib.Path("threedigrid_builder")
    for pattern in ["*.c", "*.so", "*.pyd", "*.dll"]:
        for filename in p.glob("**/" + pattern):
            print("removing '{}'".format(filename))
            filename.unlink()
elif all(x not in sys.argv for x in {"sdist", "--version", "egg_info"}):
    # Cython is required (except for sdist or the commands used by zest.releaser)
    if not cythonize:
        sys.exit("ERROR: Cython is required to build threedigrid-builder from source.")
    if sys.platform == 'win32':
        libs = ["libthreedigrid"]
        runtime_lib_dirs = []
        include_dirs = ["libthreedigrid/include"]
    else:
        libs = ["threedigrid"]
        runtime_lib_dirs=["./libthreedigrid/lib"]
        include_dirs = []
    cython_opts = dict(
        libraries=libs,
        # We can enable this once Cython 0.3 is released:
        # define_macros=[("NPY_NO_DEPRECATED_API", "NPY_1_7_API_VERSION")],
        library_dirs=["./libthreedigrid/lib"],
        runtime_library_dirs=runtime_lib_dirs,
        include_dirs=include_dirs
    )

    cython_modules = [
        Extension(
            "*", sources=["threedigrid_builder/grid/fwrapper/*.pyx"], **cython_opts
        )
    ]

    ext_modules += cythonize(cython_modules, language_level=3)

long_description = "\n\n".join([open("README.rst").read(), open("CHANGES.rst").read()])

version = open("VERSION.rst").read()

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

setup(
    name="threedigrid-builder",
    version=version,
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
    cmdclass={"build_ext": build_ext},
    zip_safe=False,
    package_data={"threedigrid_builder.lib": ["../libthreedigrid/lib/*"]},
)
