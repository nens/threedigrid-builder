Changelog of threedigrid-builder
================================

0.5.3 (unreleased)
------------------

- Fixed exchange_level (dpumax) for 1D2D lines attached to non-manhole connection nodes.
  The exchange_level is now derived from the bank_levels of attached channels.

- Fixed discharge_coefficients for culverts.


0.5.2 (2021-11-02)
------------------

- Consistently write NaN (and not -9999.0) in gridadmin float datasets.

- Fix tests with GEOS 3.10.0

- Make 'meta' group complete.


0.5.1 (2021-11-01)
------------------

- Add storage_area to calculation nodes. 

- Added ds1d_half to nodes.

- Added has_embedded to attrs.


0.5.0 (2021-10-21)
------------------

- Fixed nodes.is_manhole in the gridadmin output.

- Handle user-supplied 1D-2D lines (connected point / calculation point).

- Write initial_waterlevel for 1D nodes and add 'has_initial_waterlevels' to meta.


0.4.0 (2021-09-23)
------------------

- Added 1D boundary conditions.

- Added 2D boundary conditions.

- Enable compression in HDF5 output.

- Fixed 2D lines that connect a larger to a smaller cell in south east direction.


0.3.1 (2021-08-16)
------------------

- Handle embedded connection nodes. These are removed from the grid and written to a
  new dataset "nodes_embedded".

- Fixed bug with cross sections tables being None in Grid instance

- Handle embedded channels, pipes and culverts. Embedded objects result in
  embedded nodes and and lines with kcu LINE_1D_EMBEDDED between between 2D cells.

- Fixed a bug with lines that connect nodes to themselves in quadtree generation.

- Fixed a bug with wrong usage of lines.ds1d in bottom level and cross section weights
  computation. The added attribute lines.s1d is now used, and for clarity nodes.ds1d
  was renamed to nodes.s1d.

- Added invert_level_start_point and invert_level_end_point attributes to lines.

- Fixed coordinate order in lines.line_geometries field in gridadmin.h5.


0.3.0 (2021-07-28)
------------------

- Read and convert cross section definitions.

- Solve gridadmin off-by-one errors for pumps.

- Add 'dmax' to nodes output.

- Changed external API function name to "make_gridadmin".


0.2.1 (2021-07-20)
------------------

- Fixed issue when reprojecting 0 grid refinements with pyproj 2.*

- Fixed issue when writing 0 pumps with h5py 2.*

- Fixed missing transpose when writing pumps.coordinates to HDF5.

- Added obstacles.


0.2.0 (2021-07-15)
------------------

- Added threedigrid_builder.grid.geo_utils with segmentize and line_substring functions.
  These are used to compute the Lines.line_geometries for channel lines.

- Fixed a bug in the refinement areas code (Fortran) on Ubuntu 20.04.

- Added the Pipes model that is able to compute Nodes & Lines from Pipes.
  Pipes are also included in the calculation_type and bottom_level computations.

- Added 1D-2D lines for connection nodes, manholes, and channels.

- Added culverts, orifices, and weirs.

- Added pumps (pumpstations).

- Settings and metadata are read from the SQLite. Some metadata (like model_slug) can
  also be provided via the main (make_grid) function. The metadata is written to the
  root 'attrs' of the output gridadmin.h5. The settings are written into datasets inside
  newly addres groups "grid_settings" and "tables_settings".

- Fixes for models with no channels.

- Add an optional progress callback.


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
