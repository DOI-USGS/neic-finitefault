import datetime
import json
import pathlib
import time
from typing import List, Tuple

import numpy as np
import typer
from obspy.core.utcdatetime import UTCDateTime  # type: ignore

from wasp.data_acquisition import acquisition
from wasp.data_management import filling_data_dicts
from wasp.fault_plane import create_finite_fault, event_mult_in_to_json
from wasp.get_outputs import read_solution_static_format
from wasp.management import default_dirs
from wasp.modify_jsons import modify_channels
from wasp.modify_sacs import correct_waveforms, plot_channels
from wasp.read_config import CONFIG_PATH
from wasp.seismic_tensor import get_tensor, modify_tensor, write_tensor
from wasp.static2fsp import static_to_fsp as convert_static_to_fsp
from wasp.static2srf import static_to_srf as convert_static_to_srf
from wasp.traces_properties import properties_json
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


@app.command(help="Create the finite fault from the tensor and plane information")
def create_ff(
    directory: pathlib.Path = typer.Argument(..., help="Path to read/write from"),
    gcmt_tensor_file: pathlib.Path = typer.Argument(
        ..., help="Path to the GCMT moment tensor file"
    ),
    nodal_plane: Tuple[float, float, float] = typer.Argument(
        ..., help="The nodal plane (strike, dip, rake)"
    ),
    data_types: List[ManagedDataTypes] = typer.Option(
        [],
        "-t",
        "--data-type",
        help="Type to add to the data_types list, default is []",
    ),
    rupture_velocity: float = typer.Option(
        None,
        "-v",
        "--rupture-velocity",
        help="The rupture velocity to use in the finite fault creation",
    ),
    water_level: float = typer.Option(
        0,
        "-w",
        "--water-level",
        help="The water level to use in the finite fault creation",
    ),
):
    # set default data type
    chosen_data_types: List[str]
    chosen_data_types = [d.value for d in data_types]

    # get tensor information
    tensor_info = get_tensor(cmt_file=gcmt_tensor_file)

    event_mult_in_to_json(directory=directory)
    create_finite_fault(
        tensor_info=tensor_info,
        np_plane_info={
            "strike": nodal_plane[0],
            "dip": nodal_plane[1],
            "rake": nodal_plane[2],
        },
        data_type=chosen_data_types,
        water_level=water_level,
        rupture_vel=rupture_velocity,
        directory=directory,
    )


@app.command(help="Convert Event_mult.in to json")
def eventmult_to_json(
    directory: pathlib.Path = typer.Argument(..., help="Path to read/write from"),
):
    event_mult_in_to_json(directory=directory)


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


@app.command(help="Write trace properties (sampling_filter.json)")
def sampling_filtering(
    directory: pathlib.Path = typer.Argument(..., help="Path to the data directory"),
    gcmt_tensor_file: str = typer.Argument(
        ..., help="Path to the GCMT moment tensor file"
    ),
    cgps_dt: float = typer.Option(
        None,
        "-cdt",
        "--cgps-dt",
        help="CGPS dt to be used in in USGS routines",
    ),
):
    # get tensor information
    tensor_info = get_tensor(cmt_file=gcmt_tensor_file)
    properties_json(tensor_info, dt_cgps=cgps_dt, data_directory=directory)


@app.command(help="Convert the static solution to SRF format")
def static_to_srf(
    directory: pathlib.Path = typer.Argument(
        ..., help="Path to the directory to read/write from"
    ),
    gcmt_tensor_file: pathlib.Path = typer.Argument(
        ..., help="Path to the GCMT moment tensor file"
    ),
    data_types: List[ManagedDataTypes] = typer.Option(
        [],
        "-t",
        "--data-type",
        help="Type to add to the data_types list, default is []",
    ),
    segments_file: pathlib.Path = typer.Option(
        None,
        "-s",
        "--segments",
        help="Path to the segments file (otherwise assumes in specified directory)",
    ),
    velocity_model_file: pathlib.Path = typer.Option(
        None,
        "-v",
        "--velocity-model",
        help="Path to the velocity model file (otherwise assumes in specified directory)",
    ),
):
    # set default data type
    chosen_data_types: List[str]
    chosen_data_types = [d.value for d in data_types]
    for d in data_types:
        if d not in ["insar", "gps"]:
            validate_files([directory / DEFAULT_MANAGEMENT_FILES[d]])

    # get tensor information
    tensor_info = get_tensor(cmt_file=gcmt_tensor_file)

    # get velocity model
    if velocity_model_file:
        with open(velocity_model_file) as v:
            vel_model = json.load(v)
    else:
        with open(directory / "velmodel_data.json") as v:
            vel_model = json.load(v)

    # get velocity model
    if segments_file:
        with open(segments_file) as s:
            segments_data = json.load(s)
    else:
        with open(directory / "segments_data.json") as s:
            segments_data = json.load(s)

    # get the solution
    solution = read_solution_static_format(
        segments=segments_data["segments"], data_dir=directory
    )

    # convert static to srf
    convert_static_to_srf(
        tensor_info=tensor_info,
        segments_data=segments_data,
        used_data=chosen_data_types,
        vel_model=vel_model,
        solution=solution,
        directory=directory,
    )


@app.command(help="Convert the static solution to FSP format")
def static_to_fsp(
    directory: pathlib.Path = typer.Argument(
        ..., help="Path to the directory to read/write from"
    ),
    gcmt_tensor_file: pathlib.Path = typer.Argument(
        ..., help="Path to the GCMT moment tensor file"
    ),
    data_types: List[ManagedDataTypes] = typer.Option(
        [],
        "-t",
        "--data-type",
        help="Type to add to the data_types list, default is []",
    ),
    segments_file: pathlib.Path = typer.Option(
        None,
        "-s",
        "--segments",
        help="Path to the segments file (otherwise assumes in specified directory)",
    ),
    velocity_model_file: pathlib.Path = typer.Option(
        None,
        "-v",
        "--velocity-model",
        help="Path to the velocity model file (otherwise assumes in specified directory)",
    ),
):
    # set default data type
    chosen_data_types: List[str]
    chosen_data_types = [d.value for d in data_types]
    for d in data_types:
        if d not in ["insar", "gps"]:
            validate_files([directory / DEFAULT_MANAGEMENT_FILES[d]])

    # get tensor information
    tensor_info = get_tensor(cmt_file=gcmt_tensor_file)

    # get velocity model
    if velocity_model_file:
        with open(velocity_model_file) as v:
            vel_model = json.load(v)
    else:
        with open(directory / "velmodel_data.json") as v:
            vel_model = json.load(v)

    # get velocity model
    if segments_file:
        with open(segments_file) as s:
            segments_data = json.load(s)
    else:
        with open(directory / "segments_data.json") as s:
            segments_data = json.load(s)

    # get the solution
    solution = read_solution_static_format(
        segments=segments_data["segments"], data_dir=directory
    )

    # convert static to srf
    convert_static_to_fsp(
        tensor_info=tensor_info,
        segments_data=segments_data,
        used_data=chosen_data_types,
        vel_model=vel_model,
        solution=solution,
        directory=directory,
    )


@app.command(help="Write the tensor file from a GCMT moment tensor")
def tensor_from_gcmt(
    gcmt_tensor_file: str = typer.Argument(
        ..., help="Path to the GCMT moment tensor file"
    ),
    directory: pathlib.Path = typer.Option(
        pathlib.Path(),
        "-d",
        "--directory",
        help="Directory where to write the tensor file. Default is working directory",
    ),
):
    # get tensor information
    tensor_info = get_tensor(cmt_file=gcmt_tensor_file)

    tensor_info = get_tensor(cmt_file=gcmt_tensor_file)
    tensor_info = modify_tensor(tensor_info)
    date_origin: UTCDateTime = tensor_info["date_origin"]
    delta = datetime.datetime.utcnow() - date_origin.datetime
    tensor_info["timedelta"] = delta.seconds
    moment_mag = tensor_info["moment_mag"]
    moment_mag = 2 * np.log10(moment_mag) / 3 - 10.7
    write_tensor(tensor_info, directory=directory)


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
