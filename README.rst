threedigrid-builder
===================

.. image:: https://github.com/nens/threedigrid-builder/actions/workflows/test.yml/badge.svg
	:alt: Github Actions status
	:target: https://github.com/nens/threedigrid-builder/actions/workflows/test.yml

.. image:: https://img.shields.io/pypi/v/threedigrid-builder.svg
	:alt: PyPI
	:target: https://pypi.org/project/threedigrid-builder/


Generate a 3Di simulation grid from a model schematisation.


Usage
-----

This library converts a model schematisation to a 3Di simulation grid. This can be done
using a single function that reads data from an SQLite and TIFF and then outputs the
generated grid into a Geopackage or HDF5 file:

.. code:: python

  >>> from threedigrid_builder import make_gridadmin
  >>> sqlite_path = "/path/to/model.sqlite"
  >>> dem_path = "/path/to/dem.tiff"
  >>> out_path = "grid.gpkg"  # or "something.h5" for HDF5 output
  >>> make_gridadmin(sqlite_path, dem_path, out_path)


Alternatively, the generated grid can be output in-memory:

.. code:: python

  >>> make_gridadmin(sqlite_path, dem_path)
  {'nodes': {'id': array([   1,    2,    3, ..., 7903, 7904, 7905], dtype=int32), ...}


Installation
------------

This package is distributed as source and binary wheels on PyPI. The currently supported platforms are Windows, Linux, and OSX, all
64 bit versions only.

First install sqlite and spatialite libraries, e.g. on Ubuntu::

  $ sudo apt-get install sqlite3 libsqlite3-mod-spatialite

For raster input, GDAL is required to be present. We omitted these from the dependencies
because installation of GDAL depends on your platform an on your personal perference.
One option is to install gdal using apt-get, and then pygdal with a matching version:

  $ sudo apt-get libgdal-dev
  $ pip install pygdal=={your gdal version}.*

The current platform version of GDAL can be retrieved by:

  $ gdalinfo --version

Install the threedigrid-builder::

  $ pip install threedigrid-builder

For output into a file for the 3Di calculationcore, enable gridadmin output::

  $ pip install threedigrid-builder[gridadmin]

For output into Geopackage for display in e.g. QGis, enable gpkg output::

  $ pip install threedigrid-builder[gpkg]

The command line interface requires Typer::

  $ pip install threedigrid-builder[cli]

Development
-----------

Clone the repo and fetch the LFS objects::

  $ git lfs fetch origin refs/remotes/origin/master

Install the build tools::

  $ sudo apt install cmake gfortran

Install platform dependencies::
  
  $ sudo apt-get update && sudo apt-get install --yes --no-install-recommends libgdal-dev sqlite3 libsqlite3-mod-spatialite

Create and activate a virtual environment::

  $ python -m venv ./venv
  $ source ./venv/bin/activate

Install the dependencies. For your distribution, check the dependency matrix in .github/workflows/test.yml. For example, for Python 3.10 with numpy 1::

  $ pip install --upgrade pip wheel scikit-build
  $ pip install setuptools==63.*
  $ pip install numpy==1.23.*
  $ pip install -e .[test,gridadmin] --no-build-isolation h5py==3.7.* sqlalchemy==1.4.40 shapely==2.0.* pyproj==3.4.* "pygdal==$(gdal-config --version).*"

In case the Fortan code needs to be recompiled::

  $ python setup.py develop

Now you should be able to run the tests::

  $ pytest
  $ pytest integration_tests

For VSCode, optionally select the python interpreter corresponding to the virtual environment.