from threedigrid_builder.base import array_of

import pygeos


class GridRefinement:
    id: int
    display_name: str
    code: str
    refinement_level: int
    the_geom: pygeos.Geometry


@array_of(GridRefinement)
class GridRefinements:
    pass
