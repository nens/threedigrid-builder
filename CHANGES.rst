Changelog of threedigrid-builder
================================

0.1 (unreleased)
----------------

- Partially ported functionality from inpy (generate 3di files, makegrid).

- Breaking change: the interpolation between cross section locations (channels)
  now also extrapolates for velocity points that are not in between 2 
  connection nodes if the channel has at least 2 connection nodes. This is
  reflected by the line.cross_weight less than 0 or greater than 1.
