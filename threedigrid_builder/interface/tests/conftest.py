from condenser import NumpyQuery
from condenser.utils import load_spatialite
from sqlalchemy import create_engine
from sqlalchemy.orm import scoped_session
from sqlalchemy.orm import sessionmaker
from sqlalchemy.event import listen
from sqlalchemy.sql import select, func
import pytest
from threedi_modelchecker.threedi_model import models
from threedi_modelchecker.threedi_model import constants


@pytest.fixture(scope="session")
def db_data():
    return {
        "global_settings": [
            models.GlobalSetting(
                use_2d_flow=True,
                use_1d_flow=True,
                sim_time_step=1,
                nr_timesteps=100,
                start_date="2010",
                grid_space=1,
                dist_calc_points=100,
                kmax=10,
                table_step_size=1,
                flooding_threshold=1,
                advection_1d=1,
                advection_2d=1,
                frict_coef=1,
                initial_waterlevel=1,
                dem_obstacle_detection=False,
                timestep_plus=1,
                use_0d_inflow=False,
                use_2d_rain=False,
                numerical_settings_id=1,
                epsg_code=4326,
            )
        ],
        "channels": [
            models.Channel(
                id=12,
                code="12",
                display_name="Channel 12",
                calculation_type=constants.CalculationType.CONNECTED,
                dist_calc_points=100.0,
                zoom_category=constants.ZoomCategories.LOW_VISIBILITY,
                the_geom="EPSG=4326; LINESTRING (5.2 52.1 5.21 52.11)",
                connection_node_start_id=1,
                connection_node_end_id=3,
            )
        ],
        "connection_nodes": [
            models.ConnectionNode(
                id=1,
                code="1",
                the_geom="EPSG=4326; POINT (5.2 52.1)",
                storage_area=15.0,
            ),
            models.ConnectionNode(
                id=3,
                code="3",
                the_geom="EPSG=4326; POINT (5.21 52.11)",
                storage_area=13.0,
            ),
        ],
    }


@pytest.fixture(scope="session")
def db_engine(db_data):
    """yields a SQLAlchemy engine which is suppressed after the test session"""
    engine = create_engine("sqlite://")

    # https://geoalchemy-2.readthedocs.io/en/latest/spatialite_tutorial.html
    listen(engine, "connect", load_spatialite)
    conn = engine.connect()
    conn.execute(select([func.InitSpatialMetaData()]))
    conn.close()

    models.Base.metadata.create_all(engine)

    # write the sample data
    session = scoped_session(sessionmaker(bind=engine))()
    for records in db_data.values():
        session.add_all(records)
    session.commit()
    session.close()

    yield engine
    engine.dispose()


@pytest.fixture(scope="session")
def db_session_factory(db_engine):
    """returns a SQLAlchemy scoped session factory"""
    return scoped_session(sessionmaker(bind=db_engine))


@pytest.fixture(scope="function")
def db_session(db_session_factory):
    """yields a SQLAlchemy connection which is rollbacked after the test"""
    session = db_session_factory(query_cls=NumpyQuery)

    yield session

    session.rollback()
    session.close()
