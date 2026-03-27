import pathlib
from typing import Optional

import typer

DATA_DIRECTORY: Optional[pathlib.Path] = typer.Option(
    None,
    "-d",
    "--data-dir",
    help=("Path to where files will be read/written"),
)
QUERY_EVENTID: str = typer.Argument(
    ...,
    help="An event ID used to identify an earthquake in the USGS Feeds or ComCat",
)

QUERY_SOURCE: str = typer.Option(
    "us", "-s", "--source", help="The source of the product to use"
)
