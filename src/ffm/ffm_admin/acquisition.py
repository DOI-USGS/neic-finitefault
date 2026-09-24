import json
import pathlib
from datetime import datetime
from typing import List, Optional, Tuple, Union

import numpy as np
import typer

from ffm.acquisition.cmt import Cmt
from ffm.acquisition.cwbquery import CwbQuery
from ffm.acquisition.dartquery import DartQuery
from ffm.acquisition.irisquery import IrisQuery
from ffm.data_processing import process_dart
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
    allow_unreviewed: bool = typer.Option(
        False,
        "-u",
        "--allow-unreviewed",
        help="Allow use of unreviewed origin and moment-tensor products",
    ),
):
    """Get an event detail from the USGS feeds (automatically routes to ComCat if event does not exist)
    and use an origin product from the specified source (default is 'us') to generate a CMT file
    """
    cmt = Cmt.from_id(eventid=eventid, source=source, allow_unreviewed=allow_unreviewed)
    data_directory = validate_data_directory(data_directory)
    cmt.write(data_directory / f"{eventid}_cmt_CMT")


@app.command(
    help="Get DART data from NOAA. One of eventid, gcmt_tensor_file, or qcmt_tensor_file must be defined."
)
def dart(
    eventid: str = typer.Option(
        None, "-e", "--eventid",
        help="An event ID used to identify an earthquake in the USGS Feeds or ComCat",
    ),
    source: str = QUERY_SOURCE,
    gcmt_tensor_file: pathlib.Path = typer.Option(
        None, "-g", "--gcmt", help="Path to the GCMT moment tensor file"
    ),
    qcmt_tensor_file: pathlib.Path = typer.Option(
        None, "-q", "--qcmt", help="Path to the QuakeML moment tensor file"
    ),
    latitude: float = typer.Option(
        None, "-lat", "--latitude", help="Override origin or tensor's latitude"
    ),
    longitude: float = typer.Option(
        None, "-lon", "--longitude", help="Override origin or tensor's  longitude"
    ),
    dt: float = typer.Option(
        60.0,
        "-dt",
        "--sampling_rate",
        help="Desired data sampling rate, in seconds, default is 60s",
    ),
    max_distance: float = typer.Option(
        20,
        "-max",
        "--maximum_distance",
        help="Maximum station hypocentral distance, in kilometers, default is 20",
    ),
    hours_after: int = typer.Option(
        8,
        help="How many hours after the event to get data for",
    ),
    hours_before: int = typer.Option(
        1,
        help="How many hours before the event to get data for",
    ),
    stations_file: pathlib.Path = typer.Option(
        pathlib.Path(__file__).parent.parent
        / "acquisition"
        / "DART_station_locations.txt",
        "-sf",
        "--stations-file",
        help="A file of dart station locations.",
    ),
    plot: bool = typer.Option(
        False,
        "-p",
        "--plot_observations",
        help="Whether to save plots of collected observations, default is False",
    ),
    fcorner: Tuple[float, float] = typer.Option(
        [5, 120],
        "-f",
        "--frequency_bandpass",
        help="Corner frequencies for bandpass filter, in seconds, default is 5 to 120 s",
    ),
    # Although breaking convention with other workflows, the dart data should always
    # be processed to interpolated data right away, due to inconsistencies with sampling
    # In general processing should fall into a separate step from data acquisition
    raw: bool = typer.Option(
        False,
        "-r",
        "--raw",
        help="Get raw data. Default is to process with a bandpass filter.",
    ),
    data_directory: Optional[pathlib.Path] = DATA_DIRECTORY,
):
    data_directory = validate_data_directory(data_directory, make_directory=True)
    # initialize query
    if (
        len(
            [
                input_type
                for input_type in [eventid, gcmt_tensor_file, qcmt_tensor_file]
                if input_type is not None
            ]
        )
        > 1
    ):
        raise ValueError(
            "eventid, gcmt_tensor_file, and qcmt_tensor_file are mutually exclusive."
        )
    if gcmt_tensor_file is not None:
        dart_query = DartQuery.from_gcmt_tensor_file(gcmt_tensor_file)
    elif qcmt_tensor_file:
        dart_query = DartQuery.from_qcmt_tensor_file(qcmt_tensor_file)
    elif eventid is not None:
        dart_query = DartQuery.from_id(eventid=eventid, source=source)
    dart_query.latitude = latitude or dart_query.latitude
    dart_query.longitude = longitude or dart_query.longitude
    dart_query.dt = dt or dart_query.dt
    dart_query.max_distance = max_distance or dart_query.max_distance
    dart_query.hours_after = hours_after or dart_query.hours_after
    dart_query.hours_before = hours_before or dart_query.hours_before

    # get stations
    stations = dart_query.get_stations(stations_file)
    if len(stations) == 0:
        typer.echo(
            f"No stations within maximum distance boundary ({dart_query.max_distance} decimal degrees) for event at lat={dart_query.latitude} lon={dart_query.longitude}"
        )
        raise typer.Exit()

    # get data
    stations = dart_query.get_data(stations)

    # write data
    dart_query.write_data(
        stations=stations,
        directory=data_directory,
        frequency_corners=fcorner,
        dt=dt,
        raw=raw,
        plot=plot,
    )


@app.command(help="Get teleseismic data from EARTHSCOPE or USGS NEIC CWB")
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
        ["IU", "II", "GE", "G", "US"],
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
    edge_cwb_jar_path: Optional[pathlib.Path] = typer.Option(
        None,
        "-cwb",
        "--edge-cwb",
        help=(
            "Path to USGS Edge CWB Java client jar file. If not specified, "
            "teleseismic query is pointed to EARTHSCOPE (Requires java)"
        ),
    ),
    host: str = typer.Option(
        "edgecwb.usgs.gov",
        "-h",
        "--host",
        help="CWB host name. Defaults to public server.",
    ),
    java: pathlib.Path = typer.Option(
        pathlib.Path("/usr/bin/java"),
        "-j",
        "--java",
        help=("Path to Java"),
    ),
):
    """This first queries for all stations within the minimum/maximum distance
    in the list of networks, then queries for specific stations if they are
    provided in a station file"""
    data_directory = validate_data_directory(data_directory, make_directory=True)

    waveform_query: Union[CwbQuery, IrisQuery]
    if edge_cwb_jar_path is not None:
        edge_cwb_jar_path = edge_cwb_jar_path.resolve()
        java = java.resolve()
        waveform_query = CwbQuery(eventid=eventid, source=source)
    else:
        waveform_query = IrisQuery("teleseismic", eventid=eventid, source=source)
    waveform_query.query.depth = depth or waveform_query.query.depth
    waveform_query.query.event_time = event_time or waveform_query.query.event_time
    waveform_query.query.latitude = latitude or waveform_query.query.latitude
    waveform_query.query.longitude = longitude or waveform_query.query.longitude
    waveform_query.query.min_distance = min_distance
    waveform_query.query.max_distance = max_distance
    waveform_query.query.seconds_before = (
        seconds_before or waveform_query.query.seconds_before
    )
    waveform_query.query.seconds_after = (
        seconds_after or waveform_query.query.seconds_after
    )
    if stations is not None:
        with open(stations) as s:
            stations_data = json.load(s)
    else:
        stations_data = None

    if edge_cwb_jar_path is not None:
        waveform_query.get_data(  # type: ignore [union-attr]
            directory=data_directory,
            edge_cwb_jar_path=edge_cwb_jar_path,
            host=host,
            java=java,
            networks=networks,
            stations=stations_data,
        )
    else:
        stream1 = waveform_query.get_teleseismic_data(  # type: ignore [union-attr]
            networks=networks, stations=stations_data, debug=debug
        )
        stream2 = waveform_query.get_responses(stream=stream1, debug=debug)  # type: ignore [union-attr]
        waveform_query.write_data(directory=data_directory, stream=stream2)  # type: ignore [union-attr]


@app.command(help="Get strong motion data from EARTHSCOPE")
def strong_motion(
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
        0, "-min", "--min-distance", help="The minimum distance in decimal degrees"
    ),
    max_distance: float = typer.Option(
        5, "-max", "--max-distance", help="The max distance in decimal degrees"
    ),
    seconds_before: int = typer.Option(
        default=60,
        help="How many seconds prior to the event time to search for waveforms",
    ),
    seconds_after: int = typer.Option(
        default=5 * 60,
        help="How many seconds after the event time to search for waveforms",
    ),
    networks: List[str] = typer.Option(
        ["C", "C1", "II", "IU"],
        "-n",
        "--networks",
        help=("Networks to search. If not specified, default networks used."),
    ),
    include_gfz: bool = typer.Option(
        False,
        "-gfz",
        "--include-gfz",
        help="Include network 'CX' with 'HL*' station data.",
    ),
    data_directory: Optional[pathlib.Path] = DATA_DIRECTORY,
    debug: bool = typer.Option(
        False, "--debug", help="Run obspy clients with debug when running commands"
    ),
):
    """This first queries for all stations within the minimum/maximum distance
    in the list of networks, then queries for specific stations if they are
    provided in a station file"""
    data_directory = validate_data_directory(data_directory, make_directory=True)

    waveform_query = IrisQuery("strongmotion", eventid=eventid, source=source)
    waveform_query.query.depth = depth or waveform_query.query.depth
    waveform_query.query.event_time = event_time or waveform_query.query.event_time
    waveform_query.query.latitude = latitude or waveform_query.query.latitude
    waveform_query.query.longitude = longitude or waveform_query.query.longitude
    waveform_query.query.min_distance = min_distance
    waveform_query.query.max_distance = max_distance
    waveform_query.query.seconds_before = (
        seconds_before or waveform_query.query.seconds_before
    )
    waveform_query.query.seconds_after = (
        seconds_after or waveform_query.query.seconds_after
    )

    stream1 = waveform_query.get_strongmotion_data(
        networks=networks, debug=debug, include_gfz=include_gfz
    )
    print(stream1)
    stream2 = waveform_query.get_responses(stream=stream1, debug=debug)
    waveform_query.write_data(directory=data_directory, stream=stream2)
