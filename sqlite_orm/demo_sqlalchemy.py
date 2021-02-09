from threedi_modelchecker.threedi_database import ThreediDatabase
from threedi_modelchecker.threedi_model import models
from shapely.ops import transform
from shapely import wkb
import pyproj

transformer = pyproj.Transformer.from_crs(
    "epsg:4326", "epsg:28992", always_xy=True
)

def channel_points(path, out_path=None):
    sqlite_settings = {"db_path": path, "db_file": path}
    db = ThreediDatabase(
        connection_settings=sqlite_settings, db_type="spatialite", echo=False
    )
    session = db.get_session()

    nodes = []
    for channel in session.query(models.Channel).all():
        line = wkb.loads(channel.the_geom.data, hex=True)
        line_reproj = transform(transformer.transform, line)
        length = line_reproj.length
        dist = channel.dist_calc_points
        n_segments = max(round(line_reproj.length / dist), 1)
        dist /= n_segments
        nodes.extend([{
                "geometry": line_reproj.interpolate((i + 1) * dist),
                "ch_code": channel.code,
                "ch_name": channel.display_name,
                "ch_type": channel.calculation_type.value,
            } for i in range(1, n_segments)
        ])

    # transform to geopandas for IO
    if out_path:
        from geopandas import GeoDataFrame
        df = GeoDataFrame(nodes)
        df.set_geometry("geometry")
        df.to_file(out_path)

    return nodes


if __name__ == "__main__":
    channel_points(
        "/var/models/v2_bergermeer/3a58dfb2e95f10ca121e26e4259ed68f60e3268f/v2_bergermeer.sqlite",
        out_path="/var/models/noded_channels_shapely.shp",
    )
