from setuptools import setup, Extension
from Cython.Build import cythonize
import numpy as np

setup(
    name='grid3diwrap',
    ext_modules = cythonize(
        [
            Extension(
                "quadtree",
                sources=["quadtree.pyx", ],
                libraries=["gridgen3di"]
            ),
            Extension(
                "nodes",
                ["nodes.pyx"],
                libraries=["gridgen3di"]
            ),
            Extension(
                "lines",
                ["lines.pyx"],
                libraries=["gridgen3di"]
            ),
        ],
        language_level=3
    ),
    include_dirs=[np.get_include()],
    zip_safe=False,
)