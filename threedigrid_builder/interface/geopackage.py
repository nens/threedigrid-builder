try:
    import fiona  # NOQA, see Github issue 265
    import geopandas
except ImportError:
    geopandas = None
from threedigrid_builder.base import OutputInterface
from threedigrid_builder.grid import Grid

from .dict_out import DictOut

__all__ = ["GeopackageOut"]


class GeopackageOut(OutputInterface):
    def __init__(self, path):
        if not self.available():
            raise ImportError("Cannot write to GPKG if geopandas is not available.")
        if not path.suffix.lower() == ".gpkg":
            raise ValueError("Extension should be .gpkg")
        super().__init__(path)

    @staticmethod
    def available():
        return geopandas is not None

    def __enter__(self):
        return self

    def __exit__(self, *args, **kwargs):
        pass

    def write(self, grid: Grid):
        # use the DictOut writer to transform to dictionaries
        with DictOut(path=None) as out:
            layers = out.write(grid, geometry_format="native")

        # write the obtained layers using geopandas
        for name, layer in layers.items():
            if layer is None or len(layer.get("geometry", [])) == 0:
                continue
            df = geopandas.GeoDataFrame(layer, crs=grid.meta.crs_wkt)
            df.to_file(self.path, layer=name, driver="GPKG")
