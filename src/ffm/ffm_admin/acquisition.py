import json
import pathlib
from datetime import datetime
from typing import List, Optional

import typer

from ffm.acquisition.cmt import Cmt
from ffm.acquisition.irisquery import IrisQuery
from ffm.ffm_admin.common_args import (
    DATA_DIRECTORY,
    QUERY_EVENTID,
    QUERY_SOURCE,
)
from ffm.ffm_admin.fileutils import validate_data_directory

app = typer.Typer(help="Get data to use in modelling")


@app.command(help="Generate a CMT file by an event ID")
def cmt(
    eventid: str = QUERY_EVENTID,
    source: str = QUERY_SOURCE,
    data_directory: Optional[pathlib.Path] = DATA_DIRECTORY,
):
    """Get an event detail from the USGS feeds (automatically routes to ComCat if event does not exist)
    and use an origin product from the specified source (default is 'us') to generate a CMT file
    """
    cmt = Cmt.from_id(eventid=eventid, source=source)
    data_directory = validate_data_directory(data_directory)
    cmt.write(data_directory / f"{eventid}_cmt_CMT")


@app.command(help="Get teleseismic data from IRIS or USGS NEIC CWB")
def teleseismic(
    eventid: str = QUERY_EVENTID,
    source: str = QUERY_SOURCE,
    depth: float = typer.Option(
        None, "-dep", "--depth", help="Override origin's depth (km)"
    ),
    event_time: datetime = typer.Option(
        None,
        "-time",
        "--event-time",
        help=(
            "Override origin's event time "
            "(expects UTC isoformat like YYYY-MM-DDThh:mm:ss.sssZ)"
        ),
    ),
    latitude: float = typer.Option(
        None, "-lat", "--latitude", help="Override origin's latitude"
    ),
    longitude: float = typer.Option(
        None, "-lon", "--longitude", help="Override origin's longitude"
    ),
    min_distance: float = typer.Option(
        30, "-min", "--min-distance", help="The minimum distance in decimal degrees"
    ),
    max_distance: float = typer.Option(
        90, "-max", "--max-distance", help="The max distance in decimal degrees"
    ),
    seconds_before: int = typer.Option(
        default=0,
        help="How many seconds prior to the event time to search for waveforms",
    ),
    seconds_after: int = typer.Option(
        default=3000,
        help="How many seconds after the event time to search for waveforms",
    ),
    networks: List[str] = typer.Option(
        # ["IU", "II", "GE", "G", "US"],
        ["G"],
        "-n",
        "--networks",
        help=(
            "Networks to search. If not specified, all networks are accepted. "
            "Default stations file contents will always override networks specification"
        ),
    ),
    stations: pathlib.Path = typer.Option(
        None,
        "-st",
        "--stations",
        help=(
            "File containing stations to use (all other stations are ignored). "
            "File should be formatted like: {<Network>: {<station>: [<components>]}}. "
            "See example in ffm/acquisition/default_stations.json"
        ),
    ),
    data_directory: Optional[pathlib.Path] = DATA_DIRECTORY,
    debug: bool = typer.Option(
        False, "--debug", help="Run obspy clients with debug when running commands"
    ),
):
    """This first queries for all stations within the minimum/maximum distance
    in the list of networks, then queries for specific stations if they are
    provided in a station file"""
    data_directory = validate_data_directory(data_directory)

    waveform_query = IrisQuery.from_id(eventid=eventid, source=source)
    waveform_query.depth = depth or waveform_query.depth
    waveform_query.event_time = event_time or waveform_query.event_time
    waveform_query.latitude = latitude or waveform_query.latitude
    waveform_query.longitude = longitude or waveform_query.longitude
    waveform_query.min_distance = min_distance
    waveform_query.max_distance = max_distance
    waveform_query.seconds_before = seconds_before or waveform_query.seconds_before
    waveform_query.seconds_after = seconds_after or waveform_query.seconds_after
    if stations is not None:
        with open(stations) as s:
            stations_data = json.load(s)
    else:
        stations_data = None
    stream1 = waveform_query.get_data(
        networks=networks, stations=stations_data, debug=debug
    )
    print("HERE", stream1)
    stream2 = waveform_query.get_responses(stream=stream1, debug=debug)
    waveform_query.write_data(directory=data_directory, stream=stream2)
