Changelog of threedigrid-builder
================================


1.24.8 (2025-12-18)
-------------------

- Don't pin gfortran version for MacOS build.


1.24.7 (2025-12-15)
-------------------

- Match on all dimensions in cross section definitions when determining uniqueness.


1.24.6 (2025-08-20)
-------------------

- Test and build on Python 3.13 with Numpy 2.


1.24.5 (2025-06-05)
-------------------

- Stdout unix wrappers can now receive arbitrary length messages.


1.24.4 (2025-06-03)
-------------------

- Include the check_unassigned for 2D nodes for 1D2D groundwater lines.
- Exclude embedded nodes from creating 1D2D groundwater lines.
- Add preprocessor directives to only write to standard out on Linux. Fixes QGIS segfault for windows-2022 runner build.


1.24.3 (2025-04-28)
-------------------

- Use windows-2019 runner to fix QGIS segfault.


1.24.2 (2025-04-22)
-------------------

- Fix is_channel check by converting connection node ids to indexes.


1.24.1 (2025-04-09)
-------------------

- Add fixes for GDAL 3.10


1.24.0 (2025-01-24)
-------------------

- Bump schema version to 0.300


1.23.1 (2025-01-23)
-------------------

- Adapt to match some name changes in threedi-schema


1.23.0 (2025-01-16)
-------------------

- Remove transformations of geometries from schematisation


1.22.1 (2025-01-16)
-------------------

- Use pipe.geom, weir.geom and orifice.weir to build gridadmin


1.22.0 (2025-01-08)
-------------------

- Remove internal changes made for previous schema upgrades and limit changes to db interface


1.21.3 (2024-12-17)
-------------------

- Fix creating cross sections with schema 228


1.21.2 (2024-12-17)
-------------------

- Set MacOs release version to 13 to fix build


1.21.1 (2024-12-12)
-------------------

- Fix total discharge boundary type for groundwater not being properly defined.


1.21.0 (2024-11-25)
-------------------

- Set nodata value for neighbour nodes to 0 in quarters admin.

- Adapt for schema changes for 1D (schema 228)


1.20.1 (2024-10-23)
-------------------

- Fix interception file settings not being read properly in new schema version.


1.20.0 (2024-10-14)
-------------------

- Bump version for schema 227


1.19.1 (2024-09-30)
-------------------

- Define new boundary condition types.
- Handle open YZ profiles correctly when the middle is higher than the sides.


1.19.0 (2024-09-10)
-------------------

- Adapt for changes in schema upgrade for 2d and 1d2d


1.18.0 (2024-09-09)
-------------------

- Support 225 schema migration (boundary conditions and laterals).


1.17.1 (2024-09-02)
-------------------

- Fix for adding grid objects to include quarters attribute.
- Rename groundwater.equilibrium_infiltration_rate_type to equilibrium_infiltration_rate_aggregation


1.17.0 (2024-08-16)
-------------------

- Bump threedi-schema to 0.224 (updated structure control)


1.16.0 (2024-08-01)
-------------------

- Adapt grid-builder to work with schema upgrades for inflow (0.223)
- Create quarter administration.
- Add support for NumPy 2.


1.15.0 (2024-05-22)
-------------------

- Adapt gridbuilder to work with schema upgrades for model settings (0.222)



1.14.3 (2024-06-21)
-------------------

- Fix error for schematisation without any breaches


1.14.2 (2024-05-17)
-------------------

- Fixed bug where yz data wasn't tabulated correctly in case of an open
  profile with different heights for the left and right side.


1.14.1 (2024-04-26)
-------------------

- Fix 1D vegetation empty string bug


1.14.0 (2024-03-12)
-------------------

- Make package compatible with Python 3.12

- Make package compatible with GDAL 3.4 and 3.6. Note that in GDAL 3.6 refinements where the area exactly matches
  the grid may slightly differ from those in GDAL 3.4 and older.

- Add geopackage schematisation to test data

- Use geopackage, instead of spatialite, for tests

- Add geopackage to schematisations tested in integration test


1.13.1 (2024-02-19)
-------------------

- Remove python 3.12 from wheel build and set minimum python version to 3.7. 


1.13.0 (2024-02-16)
-------------------

- Add single vegetation and variable friction/vegetation for 1D elements. 

- Add a table for open YZ profile to support 1D variable friction and vegetation.

- Update unit tests and the integration test (bergemeer test).

- Remove the spatialie4 version of the integration test (bergemeer test).


1.12.2 (2023-12-07)
-------------------

- Add manhole_indicator field to gridadmin for future export.

- Reduce memory use by reducing DEM mask datatype to int8.


1.12.1 (2023-08-14)
-------------------

- Fixed failing Windows build.


1.12.0 (2023-08-14)
-------------------

- Add new friction types to support friction with conveyance.
- Pin numpy to 1.24.* on python 3.11 to prevent build errors 


1.11.4 (2023-05-17)
-------------------

- Build the release with the build package instead of setuptools.
- Rewrite release workflow to use a supported github action for github release.


1.11.3 (2023-05-08)
-------------------

- Change DB interface to match changed column names in threedi-schema.


1.11.2 (2023-05-03)
-------------------

- Added has_vegetation attribute to Meta() class.

- Add 1D boundary nodes to get_extent_1d.


1.11.1 (2023-05-01)
-------------------

- Fix MacOS build by bumping gfortran version.


1.11.0 (2023-05-01)
-------------------

- Added groundwater boundaries to line types.

- Refactored hydraulic conductivity information to be set on all
  Channels, Pipes, and Manholes in the form of hydraulic resistance.
  Previously a weighted average was taken and set on the 1D2D groundwater
  line.

- Filter boundry nodes from nodes which can have a 1D2D groundwater connection.


1.10.0 (2023-03-20)
-------------------

- Require schema version 216.

- Add vegetation drag settings to gridadmin.h5


1.9.0 (2023-03-06)
------------------

- Added 2D groundwater boundaries (boundary types 6: GROUNDWATERLEVEL
  and 7: GROUNDWATERDISCHARGE). These generate in boundary nodes connected
  to groundwater cells and lines having new kcu values
  600, 700, 800, 900.

- Added 1D-2D groundwater exchange lines (kcu 57). The generated lines have
  attributes 'cross_weight', 'frict_value1', 'frict_value2' set based on
  input 'exchange_thickness', 'hydraulic_conductivity_out' and
  'hydraulic_conductivity_in' on input Channels, Pipes and Manholes.

- Save memory by lazily creating empty columns.

- Add Python 3.11 and SQLAlchemy 2.0 support, drop SQLAlchemy 1.3.

- Raise comprehensive error for objects that connect to outside the 2D model domain.

- Set dpumax for 1D2D groundwater lines based on dmax of 1D node.

- Set dpumax for 1D2D open water lines based on the intersection of its line
  geometry with obstacles. The line geometry is the line from the 1D node to the
  2D cell center, except for potential_breaches, where it is geometry that
  was provided by the user.


1.8.0 (2023-01-19)
------------------

- Replaced threedi-modelchecker dependency with threedi-schema==214.*.

- Replaced pygeos with shapely 2.*.


1.7.1 (2023-01-18)
------------------

- Fixed breaches.line_id when there are boundary conditions in the model.

- Revert "1D-2D lines derived from exchange lines are also converted to breaches"
  from version 1.7.0.


1.7.0 (2023-01-11)
------------------

- Fix build with numpy >= 1.24

- Require schema version 214 (threedi-modelchecker >= 0.35).

- Set 1D-2D line dpumax based on v2_potential_breach.exchange_level >
  v2_exchange_line.exchange_level > highest intersected obstacle/levee
  > (existing logic) manhole/channel/pipe/culvert details.

- Associate breaches with 1D-2D lines: adapt the 2D side to the 2D side of
  the breach line. This overrides a possible exchange line. The content type
  is changed to v2_breach.

- Output breaches based on new breach lines.

- Adapt ds1d_half of 1D-2D lines to the spot where they cross a levee.

- Refactored connection node dpumax and calculation type assignment.

- Assign breaches to connection nodes according to the following priority:
  First, take the breach points of the first channel that has 2 breach.
  If there are no double breach points: take the breach points of the 
  first channel.

- Adapt 1D-2D lines generation for connection nodes to the exchange lines.
  A connection node derives its exchange lines from a particular channel.
  If the connection node has breaches assigned, take that channel. Else,
  take the first double connected channel. Else, take the first single
  connected channel.

- 1D-2D lines derived from exchange lines are also converted to breaches.
  These breaches have no properties.

- Draw breach points where the user-input linestring intersects the obstacle.


1.6.1 (2022-12-08)
------------------

- Fix setup.py (for sdist creation).


1.6.0 (2022-12-08)
------------------

- Adapt channel interpolated nodes based on the new v2_potential_breach table
  (only if the table exists).

- Adapt 1D-2D lines generation to v2_exchange_line table. Breaches are not implemented
  in that case. If there are no excange lines (or the table is missing),
  v2_connected_points are still used and breaches still work.

- Add TABULATED_YZ (7) and INVERTED_EGG (8) cross section shapes. Both are converted
  to tabulated trapezium.


1.5.1 (2022-11-30)
------------------

- Use the global 'max_infiltration_capacity', if present.

- Work around incompatibility of the system GDAL with the Fiona binary wheel
  distribution.


1.5.0 (2022-10-26)
------------------

- Add support for SQLITE schema migration.

- Added command-line interface (optionally installable via [cli]).

- Made quadtree building more efficient if refinement levels are not used.

- Fix error for models without CrossSectionDefinitions.


1.4.0 (2022-09-21)
------------------

- Add windshielding to Lines.

- Removed testbank action on local runner.

- Base table settings on modelchecker 0.28 (schema version 208). This
  includes the new 'maximum_table_step_size' and removes 
  'table_step_size_volume_2d'.


1.3.6 (2022-04-04)
------------------

- Write the literal WKT of the DEM projection into the gridadmin.


1.3.5 (2022-03-28)
------------------

- Refactored refinement level processing to reduce memory usage.


1.3.4 (2022-03-10)
------------------

- Bugfix: `half verhard` instead of `half_verhard` (without underscore).

- Include cross section information for weirs and orifices.

- Better error message when v2_connected_points are outside 2D grid domain. 

- 2D grid cannot contain uneven number of pixels in one grid cell.

- Bugfix in grid.sort() for models with no flowlines.


1.3.3 (2022-02-10)
------------------

- Use the bank_level for 1D2D lines to connection nodes without manhole but with
  storage area.


1.3.2 (2022-02-08)
------------------

- Fixed extracting EPSG code from outdated CRS WKT definitions.


1.3.1 (2022-02-07)
------------------

- Disable tests for CI Build Wheel, because GDAL is not included on build machine.


1.3.0 (2022-02-07)
------------------

- Write the WKT of the CRS (in addition to the epsg code) into the gridadmin.

- Drop rasterio as optional raster interface.

- Append / prepend coordinates to channel/culvert linestrings if they do not intersect
  the connection nodes they are attached to.

- Only give node_type 4 and kcu 52/54 to manholes with a not-NULL storage area.

- Use pygoes to calculate grid refinements.


1.2.1 (2022-01-27)
------------------

- Store epsg_code as string in GridMeta.


1.2.0 (2022-01-26)
------------------

- Interpret non-finite raster values (NaN, Inf, -Inf) as nodata.

- Use GDAL (instead of rasterio) for reading rasters, if present.


1.1.0 (2022-01-24)
------------------

- Write "grid_coordinate_attributes" also for pure 1D models.

- Make requesting spatial reference of GDAL dataset compatible with GDAL 2.x.

- Fix: do not ignore (Impervious)Surface records without geometries. These surfaces
  will get their location from their connection node.

- Do not ignore invalid geometries (surfaces, grid refinement areas, dem average areas)


1.0.2 (2022-01-17)
------------------

- Change in calculation_type type order of connection nodes. Embedded comes first.


1.0.1 (2022-01-13)
------------------

- Fixed the ordering of nodes and lines within node/line types.


1.0.0 (2022-01-12)
------------------

- Snap 2D boundary conditions to the closest edge if they are completely outside of the
  model domain.

- Raise SchematisationError instead of an internal error if the spatialite version is
  below 173.

- Raise FileNotFound instead of creating an empty file if spatialite does not exist.

- Added manhole fields (manhole_indicator, shape, width, surface_level) to nodes.

- Removed data from nodes.bottom_level for non-manhole nodes.

- Added dist_calc_points and material to lines.

- Added cross section width, height, shape to lines.

- Added sewerage_type (pipes) and sewerage (weirs/orifices) to lines.

- Added friction_type and friction_value (pipes/culverts/weirs/orifices) to lines.

- Fix: accept unknown sewerage types.


0.16.0 (2022-01-06)
-------------------

- Added crest level and crest type to to lines.

- Added connection node start and end id to lines.

- Handle non-ASCII characters in gridadmin HDF5 output.

- Fixed node ids in groundwater lines (they now connect groundwater cells instead of 
  open water cells).


0.15.0 (2022-01-05)
-------------------

- Small fix for use_2d_flow setting.

- Added zoom_category to nodes, lines and pumps.


0.14.0 (2022-01-04)
-------------------

- Add nodm and nodn for 2D boundary nodes.

- Handle use_2d_flow setting.

- Added display_name to nodes, lines and pumps.


0.13.0 (2021-12-28)
-------------------

- Enable groundwater and write dimp to nodes.


0.12.0 (2021-12-27)
-------------------

- Add drain_level of manholes to gridadmin.

- Bugfix: Set culvert calculation_type to isolated when not provided.

- Added display name to culverts, weirs, pipes, pumps, channels

- Added zoom category to pumps, pipes, culverts, orifices, weirs, manholes.


0.11.0 (2021-12-22)
-------------------

- Accept dist_calc_points <= 0; the effect is that there are no interpolated nodes.

- Ignore grid refinements with NULLs in their type or geometry fields.

- Ignore (impervious) surfaces, grid refinements, and dem averages areas with invalid
  geometries (mostly, polygons with self-intersections).

- Set ds1d of 1d2d lines to 2d cell_width.


0.10.0 (2021-12-21)
-------------------

- Bugfix: Added support for refinement geometries within smallest Grid cell.

- Reverse the order of coordinates in channel and culvert geometries if necessary.


0.9.2 (2021-12-17)
------------------

- Temporarily disable groundwater.

- Bugfix: Edge case with connected points.


0.9.1 (2021-12-16)
------------------

- Bugfix: use DEM epsg_code for 2D models.

- Bugfix: Small fix for lgrtot.

- Bugfix: Small fix adding groundwater cells.

- Bugfix: Fix pump.line remapping in case of embedded nodes.

- Bugfix: Remap surface_map.cci on grid.sort().

- Bugfix: also need to evaluate embedded nodes for connection node mapping for zero-d surface maps.

- Added pixel_width to groundwater nodes.


0.9.0 (2021-12-15)
------------------

- Add padding to area_mask for creating quadtree.

- Added groundwater 2D nodes, 2D vertical lines, and 2D groundwater lines.

- Bugfix: cross section tabulate used wrong width/height.


0.8.3 (2021-12-09)
------------------

- Only process cross section definitions that are actually used.

- Removed cross1 & cross2 in the in-memory or geopackage output and added cross_id1 and
  cross_id2.

- Fixed bug where writing a single line geometry or geometries of equal size would
  result in an incorrect cast to a numpy object dtype.

- Fixed bug with zero-d administration `cci` (index needs to be 1-based), removed `cid` field.


0.8.2 (2021-12-05)
------------------

- Fixed the Linux wheel distribution. These are now built with manylinux2014 instead of
  manylinux2010.


0.8.1 (2021-12-05)
------------------

- Added support for zero-d administration including surfaces and impervious surfaces.

- Dropped support for Python 3.6.

- Fixed __version__ attribute and  "threedigrid_builder_version" HDF5 attribute.

- Set the dpumax of a 1D line (channel, pipe, culvert, weir, orifice) always to the
  largest of its two invert levels. Previously, it was set to the largest of the two
  bottom_levels of the two adjacent nodes, which gave wrong results for lines attached
  to manholes.

- Disable extrapolation for channel node/line attributes that are derived from
  crosssection locations.

- Disable the SchematisationError when a Manhole has a bottom_level above a level
  of a connected object. Instead, emit a warning through the logger.


0.8.0 (2021-11-30)
------------------

- Added has_max_infiltration_capacity flag.

- Added breaches and levees.

- Implement GDAL as an alternative to RasterIO.

- Check the raster EPSG code and use it if the model does not have one.

- Removed 'model_area_path' feature from application.

- Added an in-memory output interface. Supply out_path=None to instead of writing the
  grid to a file, receive the grid as dictionaries of 1D ndarrays.

- Removed the "sqlalchemy<1.4" constraint, this library is compatible with SQLAlchemy 1.4


0.7.0 (2021-11-25)
------------------

- Raise SchematisationError on invalid settings.

- Removed SchematisationError on tabulated rectangle cross section definition with zero
  first "width" value.
  
- Add calculation_type for nodes to be Dem averaged.


0.6.1 (2021-11-10)
------------------

- Fixed l1dtot (exclude 1D boundaries).


0.6.0 (2021-11-09)
------------------

- Raise SchematisationError on embedding linear objects that begin/end outside of 2D
  domain. Added tests for edge cases.

- Fixed exchange_level (dpumax) for 1D2D lines attached to non-manhole connection nodes.
  The exchange_level is now derived from the bank_levels of attached channels.

- Add discharge_coefficients for structures.

- Swap the order in lines.line for 1D2D lines. The order is now (2D, 1D).

- Fixed kcu for lines attached to 1D boundary conditions.

- Copy crest_level from v2_levee if a v2_connected_point refers to one.


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
