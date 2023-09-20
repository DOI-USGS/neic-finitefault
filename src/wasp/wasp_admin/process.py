import json
import os
import pathlib
import time
from glob import glob

import typer

from wasp.data_processing import (
    select_process_cgps,
    select_process_strong,
    select_process_surf_tele,
    select_process_tele_body,
)
from wasp.seismic_tensor import get_tensor
from wasp.shift_match import (
    print_arrival,
    save_waveforms,
    shift_match2,
    shift_match_regional,
)
from wasp.wang_baseline_removal_v1 import wang_process

from .datautils import DEFAULT_MANAGEMENT_FILES, ProcessDataTypes, ShiftMatchDataTypes
from .fileutils import validate_files

app = typer.Typer(help="WASP data processing")

DATA_TYPES = [
    "cgps",
    "gps",
    "insar",
    "strong",
    "surf",
    "body",
]


@app.command(help="Base data processing routine")
def process_all(
    directory: pathlib.Path = typer.Argument(..., help="Path to the data directory"),
    data_type: ProcessDataTypes = typer.Argument(
        ...,
        help="Type to add to the data_types list, default is []",
    ),
    gcmt_tensor_file: pathlib.Path = typer.Option(
        None, "-g", "--gcmt", help="Path to the GCMT moment tensor file"
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
    # get tensor information
    if (gcmt_tensor_file is None and qcmt_tensor_file is None) or (
        gcmt_tensor_file is not None and qcmt_tensor_file is not None
    ):
        raise ValueError("Either gcmt_tensor_file or qcmt_tensor_file must be defined.")
    tensor_info = get_tensor(cmt_file=gcmt_tensor_file or qcmt_tensor_file)

    # validate files
    management_file = directory / DEFAULT_MANAGEMENT_FILES[data_type]
    sampling_filtering_file = directory / "sampling_filter.json"
    validate_files([management_file, sampling_filtering_file])

    # get the sampling filtering properties
    with open(sampling_filtering_file) as sf:
        data_prop = json.load(sf)

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
    gcmt_tensor_file: str = typer.Argument(
        ..., help="Path to the GCMT moment tensor file"
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
    # get the tensor information
    tensor_info = get_tensor(cmt_file=gcmt_tensor_file)

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
        print_arrival(tensor_info, directory=directory)
