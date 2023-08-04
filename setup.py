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

if sys.platform == "win32":
    cmake_args = ["-G", "MSYS Makefiles", "-DCMAKE_GNUtoMS=ON"]
else:
    cmake_args = []

setup(
    packages=find_packages(
        include=(
            "threedigrid_builder",
            "threedigrid_builder.*",
        ),
    ),
    cmake_args=cmake_args,
    include_package_data=True,
    zip_safe=False,
)
