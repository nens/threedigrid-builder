Changelog of threedigrid-builder
================================

0.1 (unreleased)
----------------

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
