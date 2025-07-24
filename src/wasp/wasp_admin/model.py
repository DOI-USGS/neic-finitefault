import json
import os
import pathlib
from enum import Enum
from shutil import copy2
from typing import List

import typer

from wasp.inversion_chen_new import (
    automatic_usgs,
    checkerboard,
    forward_modelling,
    manual_modelling,
    modelling_new_data,
    set_directory_structure,
)
from wasp.management import default_dirs
from wasp.multiple_solutions import multiple_solutions
from wasp.read_config import CONFIG_PATH
from wasp.seismic_tensor import get_tensor
from wasp.velocity_models import model2dict
from wasp.wasp_admin.fileutils import validate_files

from .datautils import ManagedDataTypes

app = typer.Typer(help="WASP modelling scripts")


class ModellingOption(str, Enum):
    checkerboard = "Checkerboard"
    faultrisetime = "fault&rise_time.txt"
    patches = "Patches"
    pointsource = "point_source"
    solucion = "Solution.txt"


class ModellingOption2(str, Enum):
    ffmmodel = "FFM modelling"
    forward = "forward"


class ModellingRoutine(str, Enum):
    checkerboard_model = "checkerboard_model"
    forward_model = "forward_model"
    manual_model = "manual_model"
    manual_model_add_data = "manual_model_add_data"
    auto_model = "auto_model"


@app.command(help="Run a modelling routine")
def run(
    directory: pathlib.Path = typer.Argument(
        ..., help="Path to the modelling directory"
    ),
    modelling_routine: ModellingRoutine = typer.Argument(
        ..., help="Chose the modelling routine to use"
    ),
    add_error: bool = typer.Option(
        False,
        "-e",
        "--add-error",
        help="Whether to add error in checkerboard and forward routines",
    ),
    cgps_dt: float = typer.Option(
        None,
        "-cdt",
        "--cgps-dt",
        help="CGPS dt to be used in in USGS routines",
    ),
    config_file: pathlib.Path = typer.Option(
        CONFIG_PATH, "-c", "--config-file", help="Path to config file"
    ),
    data_directory: pathlib.Path = typer.Option(
        None,
        "-d",
        "--data-dir",
        help="Path to the data directory, otherwise <directory>/data",
    ),
    data_types: List[ManagedDataTypes] = typer.Option(
        [],
        "-t",
        "--data-type",
        help="Type to add to the data_types list, default is []",
    ),
    gcmt_tensor_file: pathlib.Path = typer.Option(
        None, "-g", "--gcmt", help="Path to the GCMT moment tensor file"
    ),
    qcmt_tensor_file: pathlib.Path = typer.Option(
        None, "-q", "--qcmt", help="Path to the QuakeML moment tensor file"
    ),
    insar_ascending: pathlib.Path = typer.Option(
        None, "-ina", "--insar-ascending", help="Path to an ascending insar file"
    ),
    insar_descending: pathlib.Path = typer.Option(
        None, "-ind", "--insar-descending", help="Path to an descending insar file"
    ),
    max_slip: float = typer.Option(
        None,
        "-m",
        "--max-slip",
        help="The maximum slip to be used in checkerboard and forward routines",
    ),
    option: ModellingOption = typer.Option(
        None,
        "-o",
        "--option",
        help="The model output option to be used in checkerboard and forward routines",
    ),
    option2: ModellingOption2 = typer.Option(
        None,
        "-o2",
        "--option2",
        help="The second model output option to be used in checkerboard routines",
    ),
    remove_response: bool = typer.Option(
        True,
        "-r",
        "--remove-response",
        help=(
            "Whether to remove response for processing strong motion data "
            "used by new data and usgs routines"
        ),
    ),
    velocity_model_file: pathlib.Path = typer.Option(
        None, "-v", "--velocity-model", help="Path to the velocity model file"
    ),
):
    # set default data type
    chosen_data_types: List[str]
    chosen_data_types = [d.value for d in data_types]
    if (
        insar_ascending is not None or insar_descending is not None
    ) and "insar" not in data_types:
        chosen_data_types += ["insar"]

    # Set and validate files
    paths_to_validate = []
    directory = directory.resolve()
    data_directory = (data_directory or directory.parent / "data").resolve()
    paths_to_validate += [directory, data_directory]
    if modelling_routine != ModellingRoutine.auto_model:
        tensor_file = directory / "tensor_info.json"
    else:
        if (gcmt_tensor_file is None and qcmt_tensor_file is None) or (
            gcmt_tensor_file is not None and qcmt_tensor_file is not None
        ):
            raise ValueError(
                "Either gcmt_tensor_file or qcmt_tensor_file must be defined."
            )
        tensor_file = gcmt_tensor_file or qcmt_tensor_file
    paths_to_validate += [tensor_file]
    if insar_ascending is not None:
        paths_to_validate += [insar_ascending]
    if insar_descending is not None:
        paths_to_validate += [insar_descending]
    if velocity_model_file is not None:
        paths_to_validate += [velocity_model_file]
    if modelling_routine != ModellingRoutine.auto_model:
        segments_file = directory / "segments_data.json"
        paths_to_validate += [segments_file]
    validate_files(paths_to_validate)

    # get default directories
    default_directories = default_dirs(config_path=config_file)

    # get the tensor information
    if modelling_routine != ModellingRoutine.auto_model:
        with open(directory / tensor_file) as tf:
            tensor_info = json.load(tf)
    else:
        tensor_info = get_tensor(cmt_file=gcmt_tensor_file, quake_file=qcmt_tensor_file)

    # get the segments data
    if modelling_routine != ModellingRoutine.auto_model:
        with open(segments_file) as sf:
            segments_data = json.load(sf)

    if modelling_routine == ModellingRoutine.checkerboard_model:
        checkerboard(
            tensor_info=tensor_info,
            data_type=chosen_data_types,
            default_dirs=default_directories,
            segments_data=segments_data,
            max_slip=max_slip,
            add_error=add_error,
            option=option,
            option2=option2,
            directory=directory,
        )
    if modelling_routine == ModellingRoutine.forward_model:
        forward_modelling(
            tensor_info=tensor_info,
            data_type=chosen_data_types,
            default_dirs=default_directories,
            segments_data=segments_data,
            option=option,
            max_slip=max_slip,
            directory=directory,
        )
    if modelling_routine == ModellingRoutine.manual_model:
        manual_modelling(
            tensor_info=tensor_info,
            data_type=chosen_data_types,
            default_dirs=default_directories,
            segments_data=segments_data,
            directory=directory,
        )

    if modelling_routine == ModellingRoutine.manual_model_add_data:
        modelling_new_data(
            tensor_info=tensor_info,
            data_type=chosen_data_types,
            default_dirs=default_directories,
            data_folder=data_directory,
            segments_data=segments_data,
            st_response=remove_response,
            directory=directory,
        )
    if modelling_routine == ModellingRoutine.auto_model:
        solution_folder = set_directory_structure(tensor_info, directory=directory)
        if data_directory:
            for file in os.listdir(data_directory):
                if os.path.isfile(os.path.join(data_directory, file)):
                    copy2(
                        os.path.join(data_directory, file),
                        solution_folder / "data",
                    )
        if velocity_model_file:
            velmodel = model2dict(velocity_model_file)
        else:
            velmodel = None
        automatic_usgs(
            tensor_info=tensor_info,
            data_type=data_types,
            default_dirs=default_directories,
            velmodel=velmodel,
            dt_cgps=cgps_dt,
            st_response=remove_response,
            config_path=config_file,
            directory=solution_folder,
        )


@app.command(help="Run multiple models varying strike, dip, and/or rupture velocity")
def run_multiple(
    directory: pathlib.Path = typer.Argument(
        ..., help="Path to the modelling directory"
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
    dips: List[float] = typer.Option(
        [],
        "-d",
        "--dips",
        help="Dip values to model (use flag repeatedly for multiple dips)",
    ),
    gcmt_tensor_file: pathlib.Path = typer.Option(
        None, "-g", "--gcmt", help="Path to the GCMT moment tensor file"
    ),
    qcmt_tensor_file: pathlib.Path = typer.Option(
        None, "-q", "--qcmt", help="Path to the QuakeML moment tensor file"
    ),
    insar_ascending: pathlib.Path = typer.Option(
        None, "-ina", "--insar-ascending", help="Path to an ascending insar file"
    ),
    insar_descending: pathlib.Path = typer.Option(
        None, "-ind", "--insar-descending", help="Path to an descending insar file"
    ),
    rupture_velocities: List[float] = typer.Option(
        [],
        "-vr",
        "--rupture_velocities",
        help="Rupture velocity values to model (use flag repeatedly for multiple rupture velocities)",
    ),
    solution_folders: str = typer.Option(
        "NP3",
        "-sf",
        "--solution-folder",
        help="The general directory name for output (e.g., NP3). Each solution will be appended by solution number (e.g., NP3.0, NP3.1,... NP3.n)",
    ),
    strikes: List[float] = typer.Option(
        [],
        "-s",
        "--strikes",
        help="Strike values to model (use flag repeatedly for multiple strikes)",
    ),
    velocity_model_file: pathlib.Path = typer.Option(
        None, "-v", "--velocity-model", help="Path to the velocity model file"
    ),
):

    # set default data type
    chosen_data_types: List[str]
    chosen_data_types = [d.value for d in data_types]
    if (
        insar_ascending is not None or insar_descending is not None
    ) and "insar" not in data_types:
        chosen_data_types += ["insar"]

    # Set and validate files
    paths_to_validate: list[pathlib.Path] = []
    if (gcmt_tensor_file is None and qcmt_tensor_file is None) or (
        gcmt_tensor_file is not None and qcmt_tensor_file is not None
    ):
        raise ValueError("Either gcmt_tensor_file or qcmt_tensor_file must be defined.")
    tensor_file = gcmt_tensor_file or qcmt_tensor_file
    paths_to_validate += [tensor_file]
    if insar_ascending is not None:
        paths_to_validate += [insar_ascending]
    if insar_descending is not None:
        paths_to_validate += [insar_descending]
    if velocity_model_file is not None:
        paths_to_validate += velocity_model_file
    tensor_info = get_tensor(cmt_file=gcmt_tensor_file or qcmt_tensor_file)
    segments_file = directory / "segments_data.json"
    paths_to_validate += [segments_file]
    validate_files(paths_to_validate)

    # get default directories
    default_directories = default_dirs(config_path=config_file)

    # get the tensor information
    tensor_info = get_tensor(cmt_file=gcmt_tensor_file, quake_file=qcmt_tensor_file)
    multiple_solutions(
        tensor_info=tensor_info,
        data_type=chosen_data_types,
        default_dirs=default_directories,
        folders=solution_folders,
        strike=strikes,
        dip=dips,
        rupt_vel=rupture_velocities,
        directory=directory,
    )
