from setuptools import setup, Extension
from Cython.Build import cythonize
import numpy as np

setup(
    name='gridwrap',
    ext_modules = cythonize(
        [Extension("gridwrap", ["quadtree.pyx"], libraries=["gridgen3di"])],
        compiler_directives={'language_level' : "3"}
    ),
    include_dirs=[np.get_include()],
    zip_safe=False,
)