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


class CMakeExtension(Extension):
    def __init__(self, name):
        # don't invoke the original build_ext for this special extension
        super().__init__(name, sources=[])


class build_ext(_build_ext):
    def run(self):
        for ext in self.extensions:
            if isinstance(ext, CMakeExtension):
                self.build_cmake(ext)
        super().run()

    def build_cmake(self, ext):
        cwd = pathlib.Path().absolute()

        # these dirs will be created in build_py, so if you don't have
        # any python sources to bundle, the dirs will be missing
        build_temp = pathlib.Path(self.build_temp)
        build_temp.mkdir(parents=True, exist_ok=True)
        extdir = pathlib.Path(self.get_ext_fullpath(ext.name))
        extdir.mkdir(parents=True, exist_ok=True)

        # example of cmake args
        config = "Debug" if self.debug else "Release"
        cmake_args = [
            "-DCMAKE_LIBRARY_OUTPUT_DIRECTORY=" + str(extdir.parent.absolute()),
            "-DCMAKE_BUILD_TYPE=" + config,
        ]

        # example of build args
        build_args = ["--config", config, "--", "-j4"]

        os.chdir(str(build_temp))
        self.spawn(["cmake", str(cwd)] + cmake_args)
        if not self.dry_run:
            self.spawn(["cmake", "--build", "."] + build_args)
        # Troubleshooting: if fail on line above then delete all possible
        # temporary CMake files including "CMakeCache.txt" in top level dir.
        os.chdir(str(cwd))

    # Add numpy include dirs without importing numpy on module level.
    # derived from scikit-hep:
    # https://github.com/scikit-hep/root_numpy/pull/292
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
