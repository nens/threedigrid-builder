from sqlalchemy import Column, String, Integer, Float
from sqlalchemy.ext.declarative import declarative_base
from sqlalchemy.orm import sessionmaker
from geoalchemy2 import Geometry


Base = declarative_base()


class GridRefinement(Base):
    """"""
    __tablename__ = "v2_grid_refinement"

    id = Column(Integer, primary_key=True)
    display_name = Column(String)
    code = Column(String)
    refinement_level = Column(Integer)
    the_geom = Column(
        Geometry(geometry_type='LINESTRING', from_text='ST_GeomFromEWKT')
    )


class GridRefinementArea(Base):
    """"""
    __tablename__ = "v2_grid_refinement_area"

    id = Column(Integer, primary_key=True)
    display_name = Column(String)
    code = Column(String)
    refinement_level = Column(Integer)
    the_geom = Column(
        Geometry(geometry_type='POLYGON', from_text='ST_GeomFromEWKT')
    )


class GlobalSettings(Base):
    """"""
    __tablename__ = "v2_global_settings"

    id = Column(Integer, primary_key=True)
    dem_file = Column(String)
    kmax = Column(Integer)
    grid_space = Column(Float)
    epsg_code = Column(Integer)
