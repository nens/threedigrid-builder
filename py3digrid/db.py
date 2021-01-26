from sqlalchemy import create_engine
from sqlalchemy.orm import sessionmaker
from sqlalchemy.event import listen
from geoalchemy2.shape import to_shape
from pygrid.models import GridRefinement, GridRefinementArea, GlobalSettings
import numpy as np


def load_spatialite(dbapi_conn, connection_record):
    dbapi_conn.enable_load_extension(True)
    dbapi_conn.load_extension('/usr/lib/x86_64-linux-gnu/mod_spatialite.so')


class ModelFileDB:

    def __init__(self, filename):
        self.connect(filename)

    def connect(self, filename):
        engine = create_engine('sqlite:///' + filename) #, echo=True
        listen(engine, 'connect', load_spatialite)
        Session = sessionmaker(bind=engine)
        self.session = Session()

    def get_grid_refinements(self):
        refinements = []

        refinement_query = self.session.query(
            GridRefinement.id,
            GridRefinement.refinement_level,
            GridRefinement.the_geom.ST_TRANSFORM(self.epsg_code)
        ).all()

        for row in refinement_query:
            refinement = {
                'id': row.id,
                'refinement_level': row.refinement_level,
                'geometry_type': 'LINESTRING',
                'geometry': np.array(to_shape(row[2]).xy)
            }
            refinements.append(refinement)

        refinementarea_query = self.session.query(
            GridRefinementArea.id,
            GridRefinementArea.refinement_level,
            GridRefinementArea.the_geom.ST_TRANSFORM(self.epsg_code)
        ).all()
        for row in refinementarea_query:
            refinement = {
                'id': row.id,
                'refinement_level': row.refinement_level,
                'geometry_type': 'POLYGON',
                'geometry': np.array(to_shape(row[2]).exterior.xy)
            }
            refinements.append(refinement)

        return refinements

    def close_connection(self):
        self.session.close()

    @property
    def epsg_code(self):
        return self.session.query(GlobalSettings.epsg_code).first()[0]

    @property
    def kmax(self):
        return self.session.query(GlobalSettings.kmax).first()[0]

    @property
    def grid_space(self):
        return self.session.query(GlobalSettings.grid_space).first()[0]
