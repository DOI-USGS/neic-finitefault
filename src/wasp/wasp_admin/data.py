import json
import pathlib
import time
from typing import List

import typer
from obspy.core.utcdatetime import UTCDateTime

from wasp.data_acquisition import acquisition  # type: ignore
from wasp.data_management import filling_data_dicts
from wasp.seismic_tensor import get_tensor
from wasp.wasp_admin.datautils import validate_data_types
from wasp.wasp_admin.fileutils import validate_files

app = typer.Typer(help="WASP data management")


DEFAULT_ACQUIRE_DATA = ["strong", "tele"]
DEFAULT_MANAGEMENT_DATA = [
    "cgps",
    "gps",
    "insar",
    "strong_motion",
    "surf_tele",
    "tele_body",
]


@app.command(help="Acquire strong motion and teleseismic data")
def acquire(
    directory: str = typer.Argument(..., help="Path to the waveform directory"),
    gcmt_tensor_file: str = typer.Argument(
        ..., help="Path to the GCMT moment tensor file"
    ),
    data_types: List[str] = typer.Option(
        [],
        "-d",
        "--data-type",
        help="Type to add to the data_types list, default is [strong, tele]",
    ),
):
    directory = pathlib.Path(directory)
    # get the tensor information
    tensor_info = get_tensor(cmt_file=gcmt_tensor_file)
    event_time = tensor_info["datetime"]
    event_time = UTCDateTime(event_time)
    lat_ep = tensor_info["lat"]
    lon_ep = tensor_info["lon"]
    depth = tensor_info["depth"]

    # set default data type
    if data_types == []:
        data_types = DEFAULT_ACQUIRE_DATA
    validate_data_types(data_types, DEFAULT_ACQUIRE_DATA)

    # aquire the data
    time0 = time.time()
    acquisition(
        event_time,
        lat_ep,
        lon_ep,
        depth,
        data_types,
        waveform_directory=directory,
    )
    typer.echo("time spent downloading metadata: ", time.time() - time0)


@app.command(help="Populate data dictionaries for managing data")
def fill_dicts(
    directory: pathlib.Path = typer.Argument(
        ..., help="Path to the waveform directory"
    ),
    gcmt_tensor_file: pathlib.Path = typer.Argument(
        ..., help="Path to the GCMT moment tensor file"
    ),
    data_types: List[str] = typer.Option(
        [],
        "-d",
        "--data-type",
        help=f"Type to add to the data_types list, default is []",
    ),
    insar_ascending: pathlib.Path = typer.Option(
        None, "-ina", "--insar-ascending", help="Path to an ascending insar file"
    ),
    insar_ascending_ramp: float = typer.Option(
        None, "-inda", "--ascending-ramp", help="Ascending insar ramp value"
    ),
    insar_descending: pathlib.Path = typer.Option(
        None, "-ind", "--insar-descending", help="Path to an descending insar file"
    ),
    insar_descending_ramp: float = typer.Option(
        None, "-indr", "--descending-ramp", help="Descending insar ramp value"
    ),
):
    # get the tensor information
    tensor_info = get_tensor(cmt_file=gcmt_tensor_file)
    event_time = tensor_info["datetime"]
    event_time = UTCDateTime(event_time)

    # set default data type
    if data_types == []:
        data_types = DEFAULT_MANAGEMENT_DATA
    validate_data_types(data_types, DEFAULT_MANAGEMENT_DATA)
    if (
        insar_ascending is not None or insar_descending is not None
    ) and "insar" not in data_types:
        data_types += ["insar"]

    # validate files
    files_to_validate = []
    for arg in [
        insar_ascending,
        insar_descending,
    ]:
        if arg is not None:
            files_to_validate += [arg]
    sampling_filtering_file = directory / "sampling_filter.json"
    files_to_validate += [sampling_filtering_file.resolve()]
    validate_files(files_to_validate)

    # get the sampling filtering properties
    with open(sampling_filtering_file) as sf:
        data_prop = json.load(sf)

    # fill data dictionaries
    filling_data_dicts(
        tensor_info,
        data_types,
        data_prop,
        directory,
        insar_asc=[insar_ascending],
        insar_desc=[insar_descending],
        ramp_asc=[insar_ascending_ramp],
        ramp_desc=[insar_descending_ramp],
        working_directory=directory,
    )
