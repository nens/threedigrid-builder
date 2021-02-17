from setuptools import Extension
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
elif "sdist" not in sys.argv:
    # Cython is required
    if not cythonize:
        sys.exit("ERROR: Cython is required to build threedigrid-builder from source.")

    cython_opts = dict(
        libraries=["threedigrid"],
        define_macros=[("NPY_NO_DEPRECATED_API", 0)],
        include_dirs=["./libthreedigrid/include"],
        library_dirs=["./libthreedigrid/lib"],
    )

    cython_modules = [
        Extension("*", sources=["threedigrid_builder/lib/*.pyx"], **cython_opts)
    ]

    ext_modules += cythonize(cython_modules, language_level=3)

long_description = "\n\n".join([open("README.rst").read(), open("CHANGES.rst").read()])

version = open("VERSION.rst").read()

install_requires = ["numpy>=1.13", "threedimodel-checker", "pygeos", "pyproj"]

test_requires = ["pytest"]

setup(
    name="threedigrid-builder",
    version=version,
    description="3Di Grid Generator",
    long_description=long_description,
    url="https://github.com/nens/grid3di",
    author="Martijn Siemerink",
    author_email="martijn.siemerink@nelen-schuurmans.nl",
    packages=["threedigrid_builder"],
    install_requires=install_requires,
    extras_require={"test": test_requires},
    python_requires=">=3.6",
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
