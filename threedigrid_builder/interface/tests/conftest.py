from condenser import NumpyQuery
from condenser.utils import load_spatialite
from sqlalchemy import create_engine
from sqlalchemy.event import listen
from sqlalchemy.orm import scoped_session
from sqlalchemy.orm import sessionmaker
from sqlalchemy.sql import func
from sqlalchemy.sql import select
from threedi_modelchecker.threedi_model import constants
from threedi_modelchecker.threedi_model import models
from threedigrid_builder.interface.db import SQLite

import os
import pathlib
import pytest


data_path = pathlib.Path(__file__).resolve().parents[0] / "data"


@pytest.fixture(scope="session")
def db():
    """Yields a threedigrid_builder.interfac.db.SQLite object with access
    to the test v2_bergermeer.sqlite."""
    return SQLite(data_path / "v2_bergermeer.sqlite")
