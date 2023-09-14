import pathlib
import time
from typing import List

import typer
from obspy.core.utcdatetime import UTCDateTime

from wasp.data_acquisition import acquisition
from wasp.seismic_tensor import get_tensor

app = typer.Typer(help="WASP data management")

DEFAULT_ACQUIRE_DATA = ["strong", "tele"]


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
    for d in data_types:
        if d not in DEFAULT_ACQUIRE_DATA:
            typer.echo(
                f"'{d}' is not in the allowed data type list: {DEFAULT_ACQUIRE_DATA}."
            )
            typer.Exit(1)

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
