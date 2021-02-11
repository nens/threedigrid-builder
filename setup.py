import builtins
from setuptools import setup, Extension
from setuptools.command.build_ext import build_ext as _build_ext
import numpy
import sys
import pathlib
import os

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


if "clean" in sys.argv:
    # delete any previously Cythonized or compiled files
    p = pathlib.Path("pygrid/pylibgrid")
    for pattern in ["*.c", "*.so", "*.pyd"]:
        for filename in p.glob(pattern):
            print("removing '{}'".format(filename))
            filename.unlink()
elif "sdist" not in sys.argv:
    ext_options = {
        "define_macros": [("NPY_NO_DEPRECATED_API", 0)],
        "libraries": ["gridgen3di"],
    }

    ext_modules = []  # TODO compile fortran library "gridgen3di" here

    # Cython is required
    if not cythonize:
        sys.exit("ERROR: Cython is required to build grid3di from source.")

    cython_modules = [
        Extension("grid3di.lib.lines", ["pygrid/pylibgrid/lines.pyx"], **ext_options),
        Extension("grid3di.lib.nodes", ["pygrid/pylibgrid/nodes.pyx"], **ext_options),
        Extension(
            "grid3di.lib.quadtree", ["pygrid/pylibgrid/quadtree.pyx"], **ext_options
        ),
    ]

    ext_modules += cythonize(
        cython_modules,
        compiler_directives={"language_level": "3"},
        # enable once Cython >= 0.3 is released
        # define_macros=[("NPY_NO_DEPRECATED_API", "NPY_1_7_API_VERSION")],
    )

long_description = "\n\n".join([open("README.rst").read(), open("CHANGES.rst").read()])

version = open("VERSION.rst").read()

install_requires = ["numpy>=1.13"]

test_requires = ["pytest"]

setup(
    name="grid3di",
    version=version,
    description="3Di Grid Generator",
    long_description=long_description,
    url="https://github.com/nens/grid3di",
    author="Martijn Siemerink",
    author_email="martijn.siemerink@nelen-schuurmans.nl",
    license="BSD 3-Clause",
    packages=["pygrid"],
    install_requires=install_requires,
    extras_require={"test": test_requires},
    python_requires=">=3",
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
    ],
    cmdclass={"build_ext": build_ext},
    zip_safe=False,
)
