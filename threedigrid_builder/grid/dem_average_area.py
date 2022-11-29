import shapely

from threedigrid_builder.base import Array


class DemAverageArea:
    id: int
    the_geom: shapely.Geometry


class DemAverageAreas(Array[DemAverageArea]):
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
    idx = cell_tree.query(dem_average_areas.the_geom, "intersects")
    nodes.has_dem_averaged[idx[1, :]] = 1
