from threedigrid_builder.interface.db import SQLite

import os
import pathlib
import pytest


data_path = pathlib.Path(__file__).resolve().parents[0] / "data"


@pytest.fixture(scope="session")
def db():
    """Yields a threedigrid_builder.interface.db.SQLite object with access
    to the test v2_bergermeer.sqlite."""
    sqlite_path = data_path / "v2_bergermeer.sqlite"
    if not os.path.isfile(sqlite_path):
        pytest.skip("sample sqlite is not available", allow_module_level=True)
    return SQLite(sqlite_path)
