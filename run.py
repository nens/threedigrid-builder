import argparse
import os
from threedigrid_builder.models import GlobalSettings
from threedigrid_builder.db import ModelFileDB
from threedigrid_builder.interface import SubgridArea
from geojson import Feature, LineString, Polygon, FeatureCollection
import json
import pyproj
import numpy as np
import time
from threedigrid_builder.lib.quadtree import QuadTree
from threedigrid_builder.lib.nodes import Nodes
from threedigrid_builder.lib.lines import Lines


def get_parser():
    """Return argument parser."""
    parser = argparse.ArgumentParser(
        description=__doc__,
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
    )
    parser.add_argument(
        "sqlite_path",
        help="Path to sqlite database.",
    )
    parser.add_argument(
        "-a",
        "--model_area_path",
        help="Path to model_area json database.",
    )
    parser.add_argument(
        "-d",
        "--debug",
        help="Option to start a debug session with run.py",
        action="store_true",
    )
    return parser


def command(**kwargs):

    """Run grid generation"""

    sqlite_path = kwargs.get("sqlite_path")
    model_root = os.path.dirname(sqlite_path)
    if not kwargs.get("model_area_path"):
        model_area = None
    else:
        model_area = kwargs.get("model_area_path")

    model_db = ModelFileDB(sqlite_path)
    dem_filepath = model_root + '/' + \
        model_db.session.query(GlobalSettings.dem_file).first()[0]

    subgrid = SubgridArea(dem_filepath, model_area)
    quadtree = QuadTree(
        subgrid=subgrid,
        max_ref_lvl=model_db.kmax,
        min_grid_size=model_db.grid_space
    )

    if kwargs.get("debug", False):
        import ipdb;ipdb.set_trace()
    refinements = model_db.get_grid_refinements()

    
    # start = time.time()
    # for refinement in refinements:
        # status = quadtree.set_refinement(refinement)
        # if status == 1:
            # quadtree.create_active_grid(subgrid.area_pix)
    # end = time.time()
    # print(
        # "==========Time needed for applying refinements is: ", end - start, "=========="
    # )
    # 
    start = time.time()
    #[quadtree.set_refinement(refinement) for refinement in refinements]
    quadtree.set_refinements(refinements)
    end = time.time()
    print(
        "==========Time needed for applying refinements is: ", end - start, "=========="
    )

    start = time.time()
    quadtree.create_active_grid(subgrid.area_pix)
    end = time.time()
    print(
        "==========Time needed for Quadtree init is: ", end - start, "=========="
    ) 

    start = time.time()
    nodes = Nodes(quadtree)
    end = time.time()
    print(
        "==========Time needed for Node init is: ", end - start, "=========="
    ) 
    if kwargs.get("debug", False):
        import ipdb;ipdb.set_trace()

    start = time.time()
    lines = Lines(nodes, subgrid.area_pix)
    end = time.time()
    print(
        "==========Time needed for Line init is: ", end - start, "=========="
    ) 

    if kwargs.get("debug", False):
        import ipdb;ipdb.set_trace()
    dump_node_grid_to_json(
        model_root,
        nodes.get_id(1, nodes.nodtot),
        nodes.get_bnds(1, nodes.nodtot),
        model_db.epsg_code
    )

    dump_line_grid_to_json(
        model_root,
        lines.get_id(1, lines.lintot),
        nodes.get_coordinates(1, nodes.nodtot),
        lines.get_line(1,lines.lintot),
        model_db.epsg_code
    )

    model_db.close_connection()


def dump_node_grid_to_json(folder, node_ids, bnds, srid):
    #xy = transform_xys(xcorner, ycorner, srid, srid)
    bnd_list = bnds.tolist()
    coords = []
    for i, bnd in enumerate(bnd_list):
        entry = [[
            (bnd[0], bnd[1]),
            (bnd[0], bnd[3]),
            (bnd[2], bnd[3]),
            (bnd[2], bnd[1]),
            (bnd[0], bnd[1])
        ]]
        polygon = Polygon(entry)
        coords.append(polygon)
    return dump_to_json(folder + '/nodes.geojson', coords, node_ids, srid)


def dump_line_grid_to_json(folder, line_ids, node_coords, line_node, srid):
    lines_list = line_node.tolist()
    #xy = transform_xys(node_x, node_y, srid, srid)
    xy = node_coords.tolist()
    coords_line = []
    for l_nodes in lines_list:
        entry = [
            (xy[l_nodes[0] - 1][0], xy[l_nodes[0] - 1][1]),
            (xy[l_nodes[1] - 1][0], xy[l_nodes[1] - 1][1]),
        ]
        linestring = LineString(entry)
        coords_line.append(linestring)

    return dump_to_json(folder + '/lines.geojson', coords_line, line_ids, srid, lines_list)


def dump_to_json(filepath, coords, ids, srid, node_ids=None):

    crs = {
        "type": "name",
        "properties": {
            "name":  "urn:ogc:def:crs:OGC::CRS84" + str(srid)
        }
    }

    features = []
    for i, coord in enumerate(coords):
        id = str(ids[i])
        if node_ids:
            node_id = node_ids[i]
            features.append(Feature(geometry=coord, properties={'id': id, 'node_a': node_id[0], 'node_b': node_id[1]}))
        else:
            features.append(Feature(geometry=coord, properties={'id': id}))
    feature_json = FeatureCollection(features) #

    with open(filepath, 'w') as f:
        json.dump(feature_json, f)

    return feature_json


def transform_xys(x_array, y_array, source_epsg, target_epsg):
    """
    Transform x_array, y_array from source_epsg_code to
    target_epsg code
    """
    if pyproj is None:
        raise_import_exception('pyproj')

    pyproj_version = pyproj.__version__.split('.')
    check_version = '2.2.0'.split('.')

    assert isinstance(x_array, np.ndarray)
    assert isinstance(y_array, np.ndarray)

    if x_array.size == 0 and y_array.size == 0:
        return np.array([[], []])
    # if threedigrid.orm.transform is None:
    #     raise ImportError('')
    projections = []
    for epsg_code in (source_epsg, target_epsg):
        if isinstance(epsg_code, bytes):
            epsg_code = epsg_code.decode('utf-8')
        epsg_str = u'epsg:{}'.format(epsg_code)

        if pyproj_version >= check_version:
            projection = pyproj.Proj(epsg_str)
        else:
            projection = pyproj.Proj(init=epsg_str)
        projections.append(projection)

    if pyproj_version >= check_version:
        reprojected = pyproj.transform(
            projections[0], projections[1], x_array, y_array, always_xy=True)
    else:
        reprojected = pyproj.transform(
            projections[0], projections[1], x_array, y_array)

    return np.array(reprojected)

def main():
    
    """Execute main program with multiprocessing."""
    try:
        return command(**vars(get_parser().parse_args()))
    except SystemExit:
        raise


if __name__ == "__main__":
    main()
