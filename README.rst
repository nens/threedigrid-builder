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

  >>> from threedigrid_builder import make_grid
  >>> sqlite_path = "/path/to/model.sqlite"
  >>> dem_path = "/path/to/dem.tiff"
  >>> out_path = "grid.gpkg"  # or "something.h5" for HDF5 output
  >>> make_grid(sqlite_path, dem_path, out_path)


Installation
------------

This package is distributed as binary only, because its (Fortran) source code
is proprietary. The only currently supported platform is Linux.

First install sqlite and spatialite libraries, e.g. on Ubuntu::

  $ sudo apt-get install sqlite3 libsqlite3-mod-spatialite

Then install the threedigrid-builder::

  $ pip install threedigrid-builder

For output into a file for the 3Di calculationcore, enable gridadmin output::

  $ pip install threedigrid-builder[gridadmin]

For output into Geopackage for display in e.g. QGis, enable gpkg output::

  $ pip install threedigrid-builder[gpkg]
