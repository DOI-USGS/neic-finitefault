import json
import pathlib
import time
from typing import Dict, List, Tuple

import typer
from obspy.core.utcdatetime import UTCDateTime  # type: ignore

from wasp.data_acquisition import acquisition
from wasp.data_management import filling_data_dicts
from wasp.management import default_dirs
from wasp.modify_jsons import modify_channels
from wasp.modify_sacs import correct_waveforms, plot_channels
from wasp.read_config import CONFIG_PATH
from wasp.seismic_tensor import get_tensor
from wasp.velocity_models import model2dict, select_velmodel, velmodel2json

from .datautils import (
    DEFAULT_MANAGEMENT_FILES,
    AcquireDataTypes,
    ManagedDataTypes,
    ModifiableDataTypes,
)
from .fileutils import validate_files

app = typer.Typer(help="Manage WASP data, faults, and property files")


def _get_correction(correction_string: str) -> Tuple[str, List[str], float]:
    "Convert the correction string to station, channel, correction value"
    split_station_string = correction_string.split(":")
    station = split_station_string[0]
    correct_split = split_station_string[-1].split("=")
    correction = float(correct_split[-1])
    channels = correct_split[0].split(",")
    return station, channels, correction


@app.command(help="Acquire strong motion and teleseismic bodywave data")
def acquire(
    directory: pathlib.Path = typer.Argument(..., help="Path to the data directory"),
    gcmt_tensor_file: str = typer.Argument(
        ..., help="Path to the GCMT moment tensor file"
    ),
    data_types: List[AcquireDataTypes] = typer.Option(
        [],
        "-t",
        "--data-type",
        help="Type to add to the data_types list, default is [strong, body]",
    ),
):
    # get the tensor information
    tensor_info = get_tensor(cmt_file=gcmt_tensor_file)
    event_time = tensor_info["datetime"]
    event_time = UTCDateTime(event_time)
    lat_ep = tensor_info["lat"]
    lon_ep = tensor_info["lon"]
    depth = tensor_info["depth"]

    # set default data type
    chosen_data_types: List[str]
    if data_types == []:
        chosen_data_types = [d.value for d in AcquireDataTypes]
    else:
        chosen_data_types = [d.value for d in data_types]

    # acquire the data
    time0 = time.time()
    acquisition(
        event_time,
        lat_ep,
        lon_ep,
        depth,
        chosen_data_types,
        waveform_directory=directory,
    )
    typer.echo(f"Time spent downloading metadata: {time.time() - time0}")


@app.command(help="Populate data dictionaries used for managing data")
def fill_dicts(
    directory: pathlib.Path = typer.Argument(..., help="Path to the data directory"),
    gcmt_tensor_file: pathlib.Path = typer.Argument(
        ..., help="Path to the GCMT moment tensor file"
    ),
    data_types: List[ManagedDataTypes] = typer.Option(
        [],
        "-t",
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
    chosen_data_types: List[str]
    if data_types == []:
        chosen_data_types = [d.value for d in ManagedDataTypes]
    else:
        chosen_data_types = [d.value for d in data_types]
    if (
        insar_ascending is not None or insar_descending is not None
    ) and "insar" not in data_types:
        chosen_data_types += ["insar"]

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
        chosen_data_types,
        data_prop,
        directory,
        insar_asc=[insar_ascending],
        insar_desc=[insar_descending],
        ramp_asc=[insar_ascending_ramp],
        ramp_desc=[insar_descending_ramp],
        working_directory=directory,
    )


@app.command(help="Modify data dictionaries used for managing data")
def modify_dicts(
    directory: pathlib.Path = typer.Argument(..., help="Path to the data directory"),
    modification: str = typer.Argument(
        ...,
        help="The type of modification to apply to the selected channels (delete or downweight)",
    ),
    data_type: ModifiableDataTypes = typer.Argument(
        ...,
        help="Type of data being modified",
    ),
    input: bool = typer.Option(
        False,
        "-i",
        "--input-channels",
        help=f"Whether to manually input channels",
    ),
    station_channels: List[str] = typer.Option(
        [],
        "-sc",
        "--station-channel",
        help=f'Station and channels to modify in format "STATION:CHANNEL1,CHANNEL2"',
    ),
):
    # validate files
    management_file = directory / DEFAULT_MANAGEMENT_FILES[data_type]
    validate_files([management_file])

    if input:
        # Modify dict
        if modification.lower() == "delete":
            modify_channels(management_file, method="delete", input=True)
        if modification.lower() == "downweight":
            modify_channels(management_file, method="downweight", input=True)
    else:
        # construct dict
        modify_dict: dict = {}
        for station in station_channels:
            split_station_string = station.split(":")
            station = split_station_string[0]
            channels = split_station_string[-1].split(",")
            if station not in modify_dict:
                modify_dict[station] = []
            for channel in channels:
                if channel not in modify_dict[station]:
                    modify_dict[station] += [channel]

        # Modify dict
        if modification.lower() == "delete":
            modify_channels(
                management_file, method="delete", input=False, station_dict=modify_dict
            )
        if modification.lower() == "downweight":
            modify_channels(
                management_file,
                method="downweight",
                input=False,
                station_dict=modify_dict,
            )


@app.command(help="Modify data in sac files")
def modify_sacs(
    directory: pathlib.Path = typer.Argument(..., help="Path to the data directory"),
    data_type: ModifiableDataTypes = typer.Argument(
        ...,
        help="Type of data being modified",
    ),
    input: bool = typer.Option(
        False,
        "-i",
        "--input-channels",
        help=f"Whether to manually input channels",
    ),
    baseline_corrections: List[str] = typer.Option(
        [],
        "-b",
        "--baseline-correct",
        help=f'Station/channels to baseline correct in format "STATION:CHANNEL1,CHANNEL2=-10"',
    ),
    time_corrections: List[str] = typer.Option(
        [],
        "-t",
        "--time-correct",
        help=f'Station/channels to time correct in format "STATION:CHANNEL1,CHANNEL2=-10"',
    ),
    plot: bool = typer.Option(
        False,
        "-p",
        "--plot",
        help=f"Create more detailed plots of waveforms",
    ),
    plot_directory: str = typer.Option(
        None,
        "-pd",
        "--plot-directory",
        help=f"Path to where plots should be written, default is the data directory",
    ),
):
    # validate files
    management_file = directory / DEFAULT_MANAGEMENT_FILES[data_type]
    validate_files([management_file])

    if input:
        # modify waveforms
        correct_waveforms(management_file, input=True)
    else:
        # construct dict
        modify_dict: dict = {}
        # for baseline correction
        for correction_string in baseline_corrections:
            station, channels, correction = _get_correction(correction_string)
            if station not in modify_dict:
                modify_dict[station] = {}
            for channel in channels:
                if channel not in modify_dict[station]:
                    modify_dict[station][channel] = {}
                modify_dict[station][channel]["baseline_correction"] = correction
        for correction_string in time_corrections:
            station, channels, correction = _get_correction(correction_string)
            if station not in modify_dict:
                modify_dict[station] = {}
            for channel in channels:
                if channel not in modify_dict[station]:
                    modify_dict[station][channel] = {}
                modify_dict[station][channel]["time_correction"] = correction
        # modify dict
        correct_waveforms(management_file, input=False, station_dict=modify_dict)

    # optional plot
    if plot:
        plot_channels(
            management_file,
            plot_directory=plot_directory or directory,
        )


@app.command(help="Write a velocity model")
def velmodel_from_tensor(
    gcmt_tensor_file: str = typer.Argument(
        ..., help="Path to the GCMT moment tensor file"
    ),
    vel_model_file: pathlib.Path = typer.Argument(
        ..., help="Path to the output velocity model file"
    ),
    config_file: pathlib.Path = typer.Option(
        CONFIG_PATH, "-c", "--config-file", help="Path to config file"
    ),
):
    # get tensor information
    tensor_info = get_tensor(cmt_file=gcmt_tensor_file)

    # get default directories
    default_directories = default_dirs(config_path=config_file)

    # get velocity model
    velmodel = select_velmodel(
        tensor_info=tensor_info,
        default_dirs=default_directories,
        directory=vel_model_file.parent,
    )
    with open(vel_model_file, "w") as outf:
        nlen = len(velmodel["p_vel"])
        outf.write("{}\n".format(nlen))
        zipped = zip(
            velmodel["p_vel"],
            velmodel["s_vel"],
            velmodel["thick"],
            velmodel["dens"],
            velmodel["qa"],
            velmodel["qb"],
        )
        for p, s, thick, dens, qa, qb in zipped:
            outf.write("{} {} {} {} {} {}\n".format(p, s, dens, thick, qa, qb))


@app.command(help="Convert the velocity model to json")
def velmodel_to_json(
    directory: pathlib.Path = typer.Argument(..., help="Path to the data directory"),
    vel_model_file: pathlib.Path = typer.Argument(
        ..., help="Path to the velocity model file"
    ),
):
    velmodel = model2dict(vel_model_file)
    velmodel2json(velmodel, directory=directory)
