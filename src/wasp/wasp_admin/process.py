import json
import os
import pathlib
import time
from glob import glob
from typing import List

import typer

from wasp.data_processing import (
    select_process_cgps,
    select_process_strong,
    select_process_surf_tele,
    select_process_tele_body,
)
from wasp.green_functions import fk_green_fun1, gf_retrieve
from wasp.management import default_dirs
from wasp.read_config import CONFIG_PATH
from wasp.seismic_tensor import get_tensor
from wasp.shift_match import (
    print_arrival,
    save_waveforms,
    shift_match2,
    shift_match_regional,
)
from wasp.wang_baseline_removal_v1 import wang_process

from .datautils import (
    DEFAULT_MANAGEMENT_FILES,
    ManagedDataTypes,
    ProcessDataTypes,
    ShiftMatchDataTypes,
)
from .fileutils import validate_files

app = typer.Typer(help="WASP data processing")


@app.command(help="Recalculate Green's functions")
def greens(
    directory: pathlib.Path = typer.Argument(..., help="Path to the data directory"),
    config_file: pathlib.Path = typer.Option(
        CONFIG_PATH, "-c", "--config-file", help="Path to config file"
    ),
    data_types: List[ManagedDataTypes] = typer.Option(
        [],
        "-t",
        "--data-type",
        help="Data types to plot misfit for, default is []",
    ),
    max_depth: float = typer.Option(
        None, "-m", "--max-depth", help="The maximum depth"
    ),
):
    # set data type
    chosen_data_types = [d.value for d in data_types]

    # validate files
    files_to_validate = []
    sampling_filtering_file = directory / "sampling_filter.json"
    tensor_file = directory / "tensor_info.json"
    files_to_validate += [sampling_filtering_file.resolve(), tensor_file]
    if "cgps" in chosen_data_types:
        gf_bank_cgps = directory / "GF_cgps"
        files_to_validate += [gf_bank_cgps]
    if "strong" in chosen_data_types:
        gf_bank_strong = directory / "GF_strong"
        files_to_validate += [gf_bank_strong]
    validate_files(files_to_validate)

    # get the tensor information
    with open(tensor_file) as tf:
        tensor_info = json.load(tf)

    # get the sampling filtering properties
    with open(sampling_filtering_file) as sf:
        data_prop = json.load(sf)

    # set default data type
    chosen_data_types = [d.value for d in data_types]

    # get default directories
    default_directories = default_dirs(config_path=config_file)

    if "strong" in chosen_data_types:
        green_dict = fk_green_fun1(
            data_prop=data_prop,
            tensor_info=tensor_info,
            location=gf_bank_strong,
            max_depth=max_depth,
            directory=directory,
        )
    if "cgps" in chosen_data_types:
        green_dict = fk_green_fun1(
            data_prop=data_prop,
            tensor_info=tensor_info,
            location=gf_bank_cgps,
            cgps=True,
            max_depth=max_depth,
            directory=directory,
        )
    gf_retrieve(chosen_data_types, default_directories, directory=directory)


@app.command(help="Base data processing routine")
def process_all(
    directory: pathlib.Path = typer.Argument(..., help="Path to the data directory"),
    data_type: ProcessDataTypes = typer.Argument(
        ...,
        help="Type to add to the data_types list, default is []",
    ),
    qcmt_tensor_file: pathlib.Path = typer.Option(
        None, "-q", "--qcmt", help="Path to the QuakeML moment tensor file"
    ),
    remove_response: bool = typer.Option(
        True,
        "-r",
        "--remove-response",
        help="Whether to remove response for processing strong motion data",
    ),
):
    # validate files
    management_file = directory / DEFAULT_MANAGEMENT_FILES[data_type]
    sampling_filtering_file = directory / "sampling_filter.json"
    tensor_file = directory / "tensor_info.json"
    validate_files([management_file, sampling_filtering_file, tensor_file])

    # get the sampling filtering properties
    with open(sampling_filtering_file) as sf:
        data_prop = json.load(sf)

    # get the tensor information
    with open(tensor_file) as tf:
        tensor_info = json.load(tf)

    dir_str = str(directory.resolve())
    if data_type == "body":
        tele_files = glob(dir_str + "/*BH*SAC") + glob(dir_str + "/*BH*sac")
        select_process_tele_body(
            tele_files, tensor_info, data_prop, directory=directory
        )
    if data_type == "surf":
        tele_files = glob(dir_str + "/*BH*SAC") + glob(dir_str + "/*BH*sac")
        select_process_surf_tele(
            tele_files, tensor_info, data_prop, directory=directory
        )
    if data_type == "strong":
        strong_files = (
            glob(dir_str + "/**.HN*SAC")
            + glob(dir_str + "/**.HL*SAC")
            + glob(dir_str + "/**.HN*sac")
            + glob(dir_str + "/**.HL*sac")
            + glob(dir_str + "/**.AH?.*")
            + glob(dir_str + "/**.AH?.*")
            + glob(dir_str + "/**_HN*sac")
            + glob(dir_str + "/**_HL*sac")
        )
        select_process_strong(
            strong_files,
            tensor_info,
            data_prop,
            remove_response=remove_response,
            directory=directory,
        )
    if data_type == "cgps":
        cgps_files = glob(dir_str + "/**.L[HXY]*SAC") + glob(dir_str + "/**.L[HXY]*sac")
        select_process_cgps(cgps_files, tensor_info, data_prop, directory=directory)


@app.command(help="Perform wang baseline removal")
def remove_baseline(
    directory: pathlib.Path = typer.Argument(..., help="Path to the data directory"),
    plot: bool = typer.Option(
        False,
        "-p",
        "--plot",
        help=f"Plot results for baseline removal procedure",
    ),
):
    plot_dir = directory / "plots"
    if not plot_dir.exists():
        os.mkdir(plot_dir)
    int_dir = directory.parent.parent / "int_STR"
    if not int_dir.exists():
        os.mkdir(int_dir)
    time0 = time.time()
    files = glob(str(directory) + "/acc*")
    for file in files:
        wang_process(file, plot=plot, directory=directory)
    print("Time spent: ", time.time() - time0)


@app.command(help="Shift data to match pick")
def shift_match(
    directory: pathlib.Path = typer.Argument(..., help="Path to the data directory"),
    data_type: ShiftMatchDataTypes = typer.Argument(
        ...,
        help="Type of data being processed",
    ),
    option: str = typer.Option(
        "auto",
        "-o",
        "--option",
        help=(
            "whether to shift by cross-correlation or plot to pick manually (auto or manual)"
        ),
    ),
    plot: bool = typer.Option(
        False,
        "-p",
        "--plot",
        help=f"Plot results of the shift",
    ),
    multiple_events: int = typer.Option(
        None,
        "-me",
        "--multiple-events",
        help=f"Create plots for multiple events, defined by the specified number",
    ),
):
    # validate files
    validate_files([directory / DEFAULT_MANAGEMENT_FILES[data_type]])

    if option == "match":
        if multiple_events is None:
            if data_type in ["body", "surf"]:
                files = shift_match2(data_type, plot=plot, directory=directory)
            else:
                files = shift_match_regional(data_type, plot=plot, directory=directory)
        else:
            if data_type in ["body", "surf"]:
                files = []
                for i in range(multiple_events + 1):
                    shift_files = shift_match2(
                        data_type, plot=plot, event=i + 1, directory=directory
                    )
                    files += shift_files
            else:
                files = []
                for i in range(multiple_events + 1):
                    shift_files = shift_match_regional(
                        data_type, plot=plot, event=i + 1, directory=directory
                    )
                    files += shift_files
        save_waveforms(data_type, files)
    else:
        tensor_file = directory / "tensor_info.json"
        validate_files([tensor_file])
        with open(tensor_file) as tf:
            tensor_info = json.load(tf)
        print_arrival(tensor_info, directory=directory)
