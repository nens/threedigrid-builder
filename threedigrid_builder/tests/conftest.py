from threedigrid_builder.interface.db import SQLite

import pathlib
import pytest


data_path = pathlib.Path(__file__).resolve().parents[0] / "data"


@pytest.fixture(scope="session")
def db():
    """Yields a threedigrid_builder.interfac.db.SQLite object with access
    to the test v2_bergermeer.sqlite."""
    return SQLite(data_path / "v2_bergermeer.sqlite")
