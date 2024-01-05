import datetime
import json
import pathlib
import time
from enum import Enum
from typing import List, Optional, Tuple, Union

import numpy as np
import typer
from obspy.core.utcdatetime import UTCDateTime  # type: ignore

from wasp.data_acquisition import acquisition
from wasp.data_management import filling_data_dicts
from wasp.fault_plane import create_finite_fault, event_mult_in_to_json
from wasp.get_outputs import read_solution_static_format
from wasp.input_files import (
    input_chen_insar,
    input_chen_near_field,
    input_chen_static,
    input_chen_tele_body,
    input_chen_tele_surf,
    inputs_simmulated_annealing,
)
from wasp.input_files import model_space as model_space_update
from wasp.input_files import plane_for_chen
from wasp.management import default_dirs
from wasp.many_events import (
    get_model_space_events,
    get_moment_events,
    get_segments_events,
    get_waveforms_events,
)
from wasp.modelling_parameters import modelling_prop
from wasp.modify_jsons import modify_channels
from wasp.modify_sacs import correct_waveforms, plot_channels
from wasp.read_config import CONFIG_PATH, PROJECT_DIRECTORY
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


VALID_INSAR_RAMPS = ["bilinear", "linear", "quadratic"]


def _get_correction(correction_string: str) -> Tuple[str, List[str], float]:
    "Convert the correction string to station, channel, correction value"
    split_station_string = correction_string.split(":")
    station = split_station_string[0]
    correct_split = split_station_string[-1].split("=")
    correction = float(correct_split[-1])
    channels = correct_split[0].split(",")
    return station, channels, correction


def _parse_insar(value: str) -> Tuple[pathlib.Path, Optional[str]]:
    """Parse an insar file string in format <path>:<ramp>"""
    if ":" not in value:
        return pathlib.Path(value), None
    else:
        parts = value.split(":")
        filepath = pathlib.Path(parts[0])
        validate_files([filepath])
        ramp = parts[-1]
        if ramp not in VALID_INSAR_RAMPS:
            raise ValueError(
                f"The insar ramp provided ({ramp}) is not valid. "
                f"Must be one of {VALID_INSAR_RAMPS}."
            )
        return filepath, ramp


@app.command(help="Acquire strong motion and teleseismic bodywave data")
def acquire(
    directory: pathlib.Path = typer.Argument(..., help="Path to the data directory"),
    gcmt_tensor_file: pathlib.Path = typer.Argument(
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
    data_types: List[ManagedDataTypes] = typer.Option(
        [],
        "-t",
        "--data-type",
        help=f"Type to add to the data_types list, default is []",
    ),
    insar_ascending: List[str] = typer.Option(
        [],
        "-ina",
        "--insar-ascending",
        help=("Path and ramp ascending insar file. " "Example: -ina <path>:<ramp>"),
    ),
    insar_descending: List[str] = typer.Option(
        [],
        "-ind",
        "--insar-descending",
        help=("Path and ramp descending insar file. " "Example: -ind <path>:<ramp>"),
    ),
):
    # set default data type
    chosen_data_types: List[str]
    chosen_data_types = [d.value for d in data_types]
    if (
        insar_ascending is not None or insar_descending is not None
    ) and "insar" not in data_types:
        chosen_data_types += ["insar"]

    # Parse files and ramps
    insar_ascending_files: Optional[List[Union[str, pathlib.Path]]]
    insar_descending_files: Optional[List[Union[str, pathlib.Path]]]
    if not len(insar_ascending):
        insar_ascending_files = None
        insar_ascending_ramps = None
    else:
        insar_ascending_files = []
        insar_ascending_ramps = []
        for ia in insar_ascending:
            filepath, ramp = _parse_insar(ia)
            insar_ascending_files += [filepath]
            insar_ascending_ramps += [ramp]
    if not len(insar_descending):
        insar_descending_files = None
        insar_descending_ramps = None
    else:
        insar_descending_files = []
        insar_descending_ramps = []
        for id in insar_descending:
            filepath, ramp = _parse_insar(id)
            insar_descending_files += [filepath]
            insar_descending_ramps += [ramp]

    # validate files
    files_to_validate = []
    sampling_filtering_file = directory / "sampling_filter.json"
    files_to_validate += [sampling_filtering_file.resolve()]
    tensor_file = directory / "tensor_info.json"
    files_to_validate += [tensor_file]
    validate_files(files_to_validate)

    # get the tensor information
    with open(tensor_file) as tf:
        tensor_info = json.load(tf)
    event_time = tensor_info["datetime"]
    event_time = UTCDateTime(event_time)

    # get the sampling filtering properties
    with open(sampling_filtering_file) as sf:
        data_prop = json.load(sf)

    # fill data dictionaries
    filling_data_dicts(
        tensor_info,
        chosen_data_types,
        data_prop,
        directory,
        insar_asc=insar_ascending_files,
        insar_desc=insar_descending_files,
        ramp_asc=insar_ascending_ramps,
        ramp_desc=insar_descending_ramps,
        working_directory=directory,
    )


@app.command(help="Manage the data for many events")
def many_events(
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
    event_folders: List[pathlib.Path] = typer.Option(
        [],
        "-e",
        "--event-folder",
        help=f"The path to the folders where different events are stored",
    ),
):
    # set data types
    chosen_data_types: List[str]
    if data_types == []:
        chosen_data_types = [d.value for d in ManagedDataTypes]
    else:
        chosen_data_types = [d.value for d in data_types]

    for event_folder in event_folders:
        tensors: list = []
        segments_events: list = []
        model_events: list = []
        annealing_events: list = []
        strong_motion_events: list = []
        tele_events: list = []
        surf_events: list = []
        cgps_events: list = []
        static_events: list = []

        # validate files
        files_to_validate = []
        for d in chosen_data_types:
            if d not in ["insar"]:
                files_to_validate += [event_folder / DEFAULT_MANAGEMENT_FILES[d]]
        sampling_filtering_file = event_folder / "sampling_filter.json"
        files_to_validate += [sampling_filtering_file.resolve()]
        annealing_file = event_folder / "annealing_prop.json"
        files_to_validate += [annealing_file.resolve()]
        model_file = event_folder / "model_space.json"
        files_to_validate += [model_file.resolve()]
        segments_file = event_folder / "segments_data.json"
        files_to_validate += [segments_file.resolve()]
        files_to_validate += [gcmt_tensor_file]
        validate_files(files_to_validate)

        # get the tensor information
        tensor_info = get_tensor(cmt_file=gcmt_tensor_file)

        # get the sampling filtering properties
        with open(sampling_filtering_file) as sf:
            data_prop = json.load(sf)

        # get the annealing properties
        with open(annealing_file) as af:
            annealing_props = json.load(af)

        # get the segment properties
        with open(segments_file) as sf:
            segments_data = json.load(sf)

        # get the model properties
        with open(model_file) as mf:
            model_data = json.load(mf)

        # fill data dictionaries
        if "strong" in chosen_data_types:
            with open(event_folder / "strong_motion_waves.json") as smw:
                traces_strong = json.load(smw)
            strong_motion_events = strong_motion_events + [traces_strong]
        if "body" in chosen_data_types:
            with open(event_folder / "tele_waves.json") as tw:
                traces_tele = json.load(tw)
            tele_events = tele_events + [traces_tele]
        if "surf" in chosen_data_types:
            with open(event_folder / "surf_waves.json") as sw:
                traces_surf = json.load(sw)
            surf_events = surf_events + [traces_surf]
        if "cgps" in chosen_data_types:
            with open(event_folder / "cgps_waves.json") as cw:
                traces_cgps = json.load(cw)
            cgps_events = cgps_events + [traces_cgps]
        if "gps" in chosen_data_types:
            with open(event_folder / "static_data.json") as sd:
                static_data = json.load(sd)
            static_events = static_events + [static_data]
        segments_events += [segments_data]
        annealing_events += [annealing_props]
        model_events += [model_data]
        tensors += [tensor_info]
    get_moment_events(annealing_events, directory=directory)
    get_model_space_events(model_events, directory=directory)
    get_segments_events(segments_events, tensors, directory=directory)
    if "strong" in chosen_data_types:
        get_waveforms_events(strong_motion_events, "strong", directory=directory)
    if "cgps" in chosen_data_types:
        get_waveforms_events(cgps_events, "cgps", directory=directory)
    if "body" in chosen_data_types:
        get_waveforms_events(tele_events, "body", directory=directory)
    if "surf" in chosen_data_types:
        get_waveforms_events(surf_events, "surf", directory=directory)
    if "gps" in chosen_data_types:
        get_waveforms_events(static_events, "gps", directory=directory)


@app.command(help="Write modelling properties file")
def model_props(
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
):
    # set default data type
    chosen_data_types: List[str]
    chosen_data_types = [d.value for d in data_types]

    # validate required file(s)
    validate_files([gcmt_tensor_file])

    # get tensor information
    tensor_info = get_tensor(cmt_file=gcmt_tensor_file)

    # get velocity model
    if segments_file:
        with open(segments_file) as s:
            segments_data = json.load(s)
    else:
        with open(directory / "segments_data.json") as s:
            segments_data = json.load(s)

    # convert write modelling props
    modelling_prop(
        tensor_info=tensor_info,
        segments_data=segments_data,
        data_type=chosen_data_types,
        directory=directory,
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
    gcmt_tensor_file: pathlib.Path = typer.Argument(
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


@app.command(help="Display the configuration file")
def show_config(
    config_file: pathlib.Path = typer.Option(
        CONFIG_PATH, "-c", "--config-file", help="Path to config file"
    ),
):
    validate_files([config_file])
    with open(config_file) as f:
        typer.echo(f.read())


@app.command(help="Convert static solution to FSP format")
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
    gcmt_tensor_file: pathlib.Path = typer.Argument(
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


@app.command(help="Update the input data dict after changes")
def update_inputs(
    annealing: bool = typer.Option(
        False, "-a", "--annealing", help="Compute files for annealing"
    ),
    config_file: pathlib.Path = typer.Option(
        CONFIG_PATH, "-c", "--config-file", help="Path to config file"
    ),
    data_types: List[ManagedDataTypes] = typer.Option(
        [],
        "-t",
        "--data-type",
        help="Type to add to the data_types list, default is []",
    ),
    directory: pathlib.Path = typer.Option(
        pathlib.Path(),
        "-d",
        "--directory",
        help="Directory where to write the tensor file. Default is working directory",
    ),
    model_space: bool = typer.Option(
        False, "-m", "--model_space", help="compute files for model space"
    ),
    plane: bool = typer.Option(
        False, "-p", "--plane", help="Compute Fault.pos, Fault.time, Niu_model"
    ),
):
    # set data types
    chosen_data_types = [d.value for d in data_types]

    # get sampling filtering data
    sampling_file = directory / "sampling_filter.json"
    tensor_file = directory / "tensor_info.json"
    validate_files([sampling_file, tensor_file])
    with open(sampling_file) as sf:
        data_prop = json.load(sf)

    # get tensor information
    with open(tensor_file) as tf:
        tensor_info = json.load(tf)

    if plane:
        # validate files
        velocity__file = directory / "velmodel_data.json"
        segments_file = directory / "segments_data.json"
        validate_files([velocity__file, segments_file, sampling_file])

        with open(velocity__file) as vm:
            velmodel = json.load(vm)
        with open(segments_file) as sd:
            segments_data = json.load(sd)
        segments = segments_data["segments"]
        min_vel = segments[0]["min_vel"]
        max_vel = segments[0]["max_vel"]
        plane_for_chen(
            tensor_info=tensor_info,
            segments_data=segments_data,
            min_vel=min_vel,
            max_vel=max_vel,
            velmodel=velmodel,
            directory=directory,
        )
    if "body" in chosen_data_types:
        validate_files([directory / "tele_waves.json"])
        input_chen_tele_body(
            tensor_info=tensor_info,
            data_prop=data_prop,
            directory=directory,
        )
    if "surf" in chosen_data_types:
        validate_files([directory / "surf_waves.json"])
        input_chen_tele_surf(
            tensor_info=tensor_info,
            data_prop=data_prop,
            directory=directory,
            config_path=config_file,
        )
    if "strong" in chosen_data_types:
        validate_files([directory / "strong_motion_waves.json"])
        input_chen_near_field(
            tensor_info=tensor_info,
            data_prop=data_prop,
            data_type="strong",
            directory=directory,
        )
    if "cgps" in chosen_data_types:
        validate_files([directory / "cgps_waves.json"])
        input_chen_near_field(
            tensor_info=tensor_info,
            data_prop=data_prop,
            data_type="cgps",
            directory=directory,
        )
    if "gps" in chosen_data_types:
        validate_files([directory / "static_data.json"])
        input_chen_static(directory=directory)
    if "insar" in chosen_data_types:
        validate_files([directory / "insar_data.json"])
        input_chen_insar(directory=directory)
    if model_space:
        validate_files([directory / "model_space.json"])
        with open(directory / "model_space.json") as mf:
            model_dict = json.load(mf)
        model_space_update(model_dict, directory=directory)
    if annealing:
        validate_files([directory / "annealing_prop.json"])
        with open(directory / "annealing_prop.json") as af:
            annealing_dict = json.load(af)
        inputs_simmulated_annealing(annealing_dict, directory=directory)


@app.command(help="Write a velocity model")
def velmodel_from_tensor(
    gcmt_tensor_file: pathlib.Path = typer.Argument(
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


@app.command(help="Display the configuration")
def write_config(
    config_file: pathlib.Path = typer.Option(
        CONFIG_PATH, "-c", "--config-file", help="Path to config file"
    ),
    code_path: pathlib.Path = typer.Option(
        PROJECT_DIRECTORY, "-code", "--code-path", help="The path to the code directory"
    ),
):
    if config_file.exists():
        raise FileExistsError(
            f"The config file '{config_file}' already exists. Exiting."
        )
    with open(config_file, "w") as cf:
        cf.write("[PATHS]\n")
        cf.write(f"code_path = {str(code_path.resolve())}\n")
        cf.write("surf_gf_bank = %(code_path)s/fortran_code/gfs_nm/long/low.in\n")
        cf.write("modelling = %(code_path)s/fortran_code/bin_inversion_gfortran_f95\n")
        cf.write("get_near_gf = %(code_path)s/fortran_code/bin_str_f95\n")
        cf.write("compute_near_gf = %(code_path)s/fortran_code/src_dc_f95\n")
        cf.write("info = %(code_path)s/fortran_code/info\n")
        cf.write("cartopy_files = %(code_path)s/fortran_code/tectonicplates")
