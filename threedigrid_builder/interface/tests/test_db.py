from unittest import mock
import pytest
import pygeos

from threedigrid_builder.grid1d import Channels
from threedigrid_builder.interface import SQLite

from numpy.testing import assert_equal


@pytest.fixture
def mocked_db():
    with mock.patch("threedigrid_builder.interface.db.ThreediDatabase") as db:
        yield db


@pytest.fixture
def mocked_sqlite(mocked_db, db_session):
    sqlite = SQLite("/some/path")
    with mock.patch.object(sqlite, "get_session") as get_session:
        get_session.return_value.__enter__.return_value = db_session
        yield sqlite


def test_init(mocked_db):
    path = "/some/path"
    sqlite = SQLite(path)

    mocked_db.assert_called_with(
        connection_settings={"db_path": path, "db_file": path}, db_type="spatialite"
    )

    assert sqlite.db is mocked_db.return_value


def test_get_channels(mocked_sqlite):
    channels = mocked_sqlite.get_channels()
    assert isinstance(channels, Channels)
