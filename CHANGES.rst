Changelog of threedigrid-builder
================================

0.1.3 (unreleased)
------------------

- Added threedigrid_builder.grid.geo_utils with segmentize and line_substring functions.
  These are used to compute the Lines.line_geometries for channel lines.

- Fixed a bug in the refinement areas code (Fortran) on Ubuntu 20.04.

- Added the Pipes model that is able to compute Nodes & Lines from Pipes.
  Pipes are also included in the calculation_type and bottom_level computations.

- Added 1D-2D lines for connection nodes, manholes, and channels.


0.1.2 (2021-04-28)
------------------

- Added public API with 1 function: `threedigrid_builder.make_grid`.


0.1.1 (2021-04-20)
------------------

- Fixed automatic PyPI upload.


0.1.0 (2021-04-20)
------------------

- Partially ported functionality from inpy (generate 3di files, makegrid): 1D channel
  grid (including calculation_type and bottom_level), and 2D open water grid.

- Added gridadmin and geopackage output.

- Breaking change: the interpolation between cross section locations (channels)
  now also extrapolates for lines and nodes  that are not in between two
  connection nodes. This happens only if the channel has at least 2 cross section
  locations. When extrapolatic, the line.cross_weight is less than 0 or greater than 1.

- Breaking change: missing or empty values in float datasets in the output gridadmin are
  now denoted by NaN (not-a-number) instead of -9999.0.

- Breaking change: integers in the output gridadmin are now always 32-bit (instead of
  sometimes 32-bit and sometimes 64-bit).
