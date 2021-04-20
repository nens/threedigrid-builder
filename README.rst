threedigrid-builder
===================

Generate threedi simulation grids from a model schematisation.


.. image:: hhttps://github.com/nens/threedigrid-builder/actions/workflows/test.yml/badge.svg
	:alt: Github Actions status
	:target: https://github.com/nens/threedigrid-builder/actions/workflows/test.yml

.. image:: https://img.shields.io/pypi/v/threedigrid-builder.svg
	:alt: PyPI
	:target: https://pypi.org/project/threedigrid-builder/


Installation
------------

Because part of its (Fortran) source code is proprietary, threedigrid-builder is
distributed as binary only. The current supported platform is Linux (manylinux1.).

First install sqlite and spatialite libraries, on Ubuntu::

  $ sudo apt-get install sqlite3 libsqlite3-mod-spatialite

Then install the threedigrid-builder:

    $ pip install threedigrid-builder

For output into a file for the 3Di calculationcore, enable gridadmin output:

    $ pip install threedigrid-builder[gridadmin]

For output into Geopackage for display in e.g. QGis, enable gpkg output:

    $ pip install threedigrid-builder[gpkg]


Usage
-----

The public API of threedigrid_builder can be found in ``threedigrid_builder.application``.
