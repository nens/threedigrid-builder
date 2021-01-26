from setuptools import setup
from Cython.Build import cythonize

setup(
    name='Test',
    ext_modules = cythonize([
    Extension("py3digrid", ["quadtree.pyx"],
              libraries=["3digridgen"])
    ])
    zip_safe=False,
)