import pathlib

import typer

from ffm.acquisition.cmt import Cmt

app = typer.Typer(help="Get data to use in modelling")


@app.command(help="Generate a CMT file by an event ID")
def cmt(
    eventid: str = typer.Argument(
        ...,
        help="An event ID used to identify an earthquake in the USGS Feeds or ComCat",
    ),
    source: str = typer.Argument("us", help="The source of the origin product to use"),
    data_directory: pathlib.Path = typer.Option(
        pathlib.Path(),
        "-d",
        "--data-dir",
        help=(
            "Path to the data directory where the CMT file will be written, "
            "otherwise <directory>/data"
        ),
    ),
):
    """Get an event detail from the USGS feeds (automatically routes to ComCat if event does not exist)
    and use an origin product from the specified source (default is 'us') to generate a CMT file
    """
    cmt = Cmt.from_id(eventid=eventid, source=source)
    cmt.write(data_directory / f"{eventid}_cmt_CMT")
