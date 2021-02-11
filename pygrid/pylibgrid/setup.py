from setuptools import setup, Extension
from Cython.Build import cythonize
import numpy as np

numpy_includes = np.get_include()
print(numpy_includes)

setup(
    name='grid3diwrap',
    ext_modules = cythonize(
        Extension(
            "*",
            sources=["*.pyx", ],
            libraries=["gridgen3di"],
            define_macros=[("NPY_NO_DEPRECATED_API", 0)],
            include_dirs=[np.get_include()]
        ),
        include_path=[np.get_include()],
        language_level=3
    ),
    zip_safe=False,
)