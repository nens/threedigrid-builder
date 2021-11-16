from threedigrid_builder.base import array_of
from threedigrid_builder.constants import CalculationType

import pygeos


class DemAverageArea:
    id: int
    the_geom: pygeos.Geometry


@array_of(DemAverageArea)
class DemAverageAreas:
    pass


def apply_dem_average_areas(nodes, cell_tree, dem_average_areas):
    """Set calculation type to node being DEM_AVERAGED in tables. When the 
    dem_average_areas polygons intersect with a 2D cell the node will be dem_averaged
    during table preprocessing.
    
    Args:
        nodes (Nodes)
        cell_tree (STRTree)
        dem_average_areas (DemAverageAreas)
    """
    idx = cell_tree.query_bulk(dem_average_areas.the_geom, "intersects")
    nodes.calculation_type[idx[1, :]] = CalculationType.DEM_AVERAGED
