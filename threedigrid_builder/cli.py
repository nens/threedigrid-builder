from pathlib import Path

import typer

from .application import make_gridadmin


def _progress_callback(progress: float, message: str):
    typer.echo("Progress: %d, Message: %s" % (progress * 100, message))


def main(
    sqlite_path: Path,
    dem_path: Path,
    out_path: Path,
    upgrade: bool = False,
):
    make_gridadmin(
        sqlite_path=sqlite_path,
        dem_path=dem_path,
        out_path=out_path,
        progress_callback=_progress_callback,
        upgrade=upgrade,
    )


def run():
    typer.run(main)


if __name__ == "__main__":
    run()
