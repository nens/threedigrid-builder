import numpy as np
import shapely
from osgeo import gdal, ogr

from threedigrid_builder.base import Array
from threedigrid_builder.utils import Dataset

gdal.SetConfigOption("GDAL_MEM_ENABLE_OPEN", "YES")


class GridRefinement:
    id: int
    display_name: str
    code: str
    refinement_level: int
    the_geom: shapely.Geometry


class GridRefinements(Array[GridRefinement]):
    def rasterize(self, **kwargs):
        """Return a raster with refinement_level at cells where a refinement intersect.

        When multiple refinements intersect the same cell, the lowest refinement level
        is taken.
        """
        # sort by refinement_level (lowest values last, so that they are on top)
        sorter = np.argsort(self.refinement_level)[::-1]
        return rasterize(self.the_geom[sorter], self.refinement_level[sorter], **kwargs)


def rasterize(geoms, values, origin, width, height, cell_size, no_data_value):
    """Rasterize a set of geometries to a raster with values at cells where the
    geometries intersect.

    Args:
        arr (array of shapely.Geometry)
        values (array of int)
        origin (tuple of 2 floats): x, y order
        width (int)
        height (int)
        cell_size (float)
        no_data_value (int)

    Returns:
        array with shape (y, x) and dtype int32
    """
    assert np.issubdtype(values.dtype, np.integer)
    assert height > 0
    assert width > 0
    dtype = np.int32  # OGR knows int32 and int64, but GDAL only int32
    ogr_dtype = ogr.OFTInteger

    # initialize the array
    array = np.full((1, height, width), no_data_value, dtype=dtype)

    # drop empty geometries
    mask = shapely.is_geometry(geoms)
    geoms = geoms[mask]
    values = values[mask]

    # if there are no features, return directly
    if len(geoms) == 0:
        return array[0]

    # create an output datasource in memory
    driver = ogr.GetDriverByName("Memory")
    burn_attr = "BURN_IT"

    # prepare in-memory ogr layer
    ds_ogr = driver.CreateDataSource("")
    layer = ds_ogr.CreateLayer("")
    layer_definition = layer.GetLayerDefn()
    field_definition = ogr.FieldDefn(burn_attr, ogr_dtype)
    layer.CreateField(field_definition)

    for geometry, value in zip(geoms, values):
        feature = ogr.Feature(layer_definition)
        feature.SetGeometry(ogr.CreateGeometryFromWkb(shapely.to_wkb(geometry)))
        if ogr_dtype is not None:
            feature[burn_attr] = value
        layer.CreateFeature(feature)

    dataset_kwargs = {
        "no_data_value": no_data_value,
        "geo_transform": [
            float(x) for x in (origin[0], cell_size, 0, origin[1], 0, cell_size)
        ],
    }

    # ATTRIBUTE=BURN_ATTR burns the BURN_ATTR value of each feature
    # ALL_TOUCHED turns off the special GDAL algorithm for lines
    options = ["ATTRIBUTE=" + burn_attr, "ALL_TOUCHED=1"]

    with Dataset(array, **dataset_kwargs) as dataset:
        gdal.RasterizeLayer(dataset, (1,), layer, options=options)

    return array[0]
