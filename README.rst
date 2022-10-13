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

Install the threedigrid-builder::

  $ pip install threedigrid-builder

For output into a file for the 3Di calculationcore, enable gridadmin output::

  $ pip install threedigrid-builder[gridadmin]

For output into Geopackage for display in e.g. QGis, enable gpkg output::

  $ pip install threedigrid-builder[gpkg]

The command line interface requires Typer::

  $ pip install threedigrid-builder[cli]
