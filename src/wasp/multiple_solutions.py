import glob
import json
import os
import pathlib
import shutil
from itertools import product
from multiprocessing import Pool, cpu_count
from typing import List, Optional, Union

import numpy as np
import pandas as pd  # type:ignore

import wasp.fault_plane as pf
import wasp.inversion_chen_new as inv
import wasp.management as mng
import wasp.seismic_tensor as tensor


def multiple_solutions(
    tensor_info: dict,
    data_type: List[str],
    default_dirs: dict,
    folders: list,
    strike: Optional[float] = None,
    dip: Optional[float] = None,
    rupt_vel: Optional[float] = None,
    directory: Union[pathlib.Path, str] = pathlib.Path(),
    test: bool = False,
):
    """Run multiple solutions, each with specified parameters

    :param tensor_info: The tensor information
    :type tensor_info: dict
    :param data_type: The list of available data types
    :type data_type: List[str]
    :param default_dirs: The locations of default directories
    :type default_dirs: dict
    :param folders: The list of solution folders
    :type folders: list
    :param strike: The list of strike values to use, defaults to None
    :type strike: Optional[float], optional
    :param dip: The list of dip values to use, defaults to None
    :type dip: Optional[float], optional
    :param rupt_vel: The list of rupture velocity values to use, defaults to None
    :type rupt_vel: Optional[float], optional
    :param test: Skip multiprocessing for tests, defaults to False
    :type test: bool, optional
    :param directory: The base directory to read/write to, defaults to pathlib.Path()
    :type directory: Union[pathlib.Path, str], optional
    """
    if not len(folders):
        print("No folders defined. Exiting.")
        return
    directory = pathlib.Path(directory).resolve()
    event_folder = pathlib.Path(os.path.dirname(directory)).resolve()
    with open(directory / "segments_data.json") as s:
        segments_data = json.load(s)
    segments = segments_data["segments"]
    strike = strike if len(strike) > 0 else [segments[0]["strike"]]  # type:ignore
    dip = dip if len(dip) > 0 else [segments[0]["dip"]]  # type:ignore
    rupt_vel = (
        rupt_vel if len(rupt_vel) > 0 else [segments[0]["rupture_vel"]]  # type:ignore
    )
    new_iter = product(strike, dip, rupt_vel)  # type:ignore
    subfolders: list = []
    for i, (strike1, dip1, rupt_vel1) in enumerate(new_iter):
        subfolder = event_folder / f"{folders}.{i}"
        shutil.copytree(directory, subfolder)
        subfolders += [subfolder]
        new_segments_data = segments_data.copy()
        new_segments_data["segments"][0]["strike"] = strike1
        new_segments_data["segments"][0]["dip"] = dip1
        new_segments_data["segments"][0]["rupture_vel"] = rupt_vel1
        segments2 = new_segments_data["segments"]
        rise_time = new_segments_data["rise_time"]
        force_plane_above_ground(tensor_info, segments2)
        with open(subfolder / "segments_data.json", "w") as f:
            json.dump(
                new_segments_data,
                f,
                sort_keys=True,
                indent=4,
                separators=(",", ": "),
                ensure_ascii=False,
            )
    cpus = cpu_count()
    processes = int(cpus / 3)
    if test:
        return
    with Pool(processes=processes) as pool:
        results = pool.starmap(
            __worker,
            [
                [tensor_info, data_type, default_dirs, subfolder]
                for subfolder in subfolders
            ],
        )
    for subfolder in subfolders:
        inv.execute_plot(
            tensor_info, data_type, segments_data, default_dirs, directory=subfolder
        )
    df = get_summary_all_models(subfolders)


def __worker(
    tensor_info: dict, data_type: List[str], default_dirs: dict, subfolder: str
):
    """Run the manual inversion

    :param tensor_info: The tensor information
    :type tensor_info: dict
    :param data_type: The data types available
    :type data_type: List[str]
    :param default_dirs: The locations of default directories
    :type default_dirs: dict
    :param subfolder: The subfolder for the specific solution
    :type subfolder: str
    """
    with open(pathlib.Path(subfolder) / "segments_data.json") as s:
        segments_data = json.load(s)
    inv.manual_modelling(
        tensor_info, data_type, default_dirs, segments_data, subfolder, plot_sol=False
    )
    return


def force_plane_above_ground(tensor_info: dict, segments: List[dict]) -> List[dict]:
    """Update the plane if it intersects with the surface

    :param tensor_info: The tensor information
    :type tensor_info: dict
    :param segments: The segments' properties
    :type segments: List[dict]
    :return: The updated segments' properties
    :rtype: List[dict]
    """
    segment = segments[0]
    correct = pf.is_fault_correct(tensor_info, segment)
    if correct:
        return segments
    else:
        depth = tensor_info["depth"]
        dip = segment["dip"]
        hyp_dip = segment["hyp_dip"]
        delta_dip = 0.99 * depth / np.sin(dip * np.pi / 180) / hyp_dip
        segments[0]["delta_dip"] = delta_dip
        return segments


def get_summary_all_models(subfolders: list) -> pd.DataFrame:
    """Get a table with misfits and parameters in all subfolders

    :param subfolders: List of folders where solutions exist
    :type subfolders: list
    :return: The model summaries
    :rtype: pd.DataFrame
    """
    strike: list[float] = []
    dip: list[float] = []
    rupt_vel: list[float] = []
    objective_error: list[float] = []
    misfit_error: list[float] = []
    short_subfolder: list[str] = []
    for subfolder in subfolders:
        subfolder = pathlib.Path(subfolder)
        short_subfolder = short_subfolder + [str(subfolder).split("/")[-1]]
        with open(subfolder / "segments_data.json") as f:
            segments_data = json.load(f)
        segments = segments_data["segments"]
        strike = strike + [segments[0]["strike"]]
        dip = dip + [segments[0]["dip"]]
        rupt_vel = rupt_vel + [segments[0]["rupture_vel"]]
        errors = get_summary(directory=subfolder)
        objective_error = objective_error + [errors["objective_error"]]
        misfit_error = misfit_error + [errors["misfit_error"]]
    df = pd.DataFrame(
        {
            "subfolders": short_subfolder,
            "strike": strike,
            "dip": dip,
            "rupt_vel": rupt_vel,
            "objective_error": objective_error,
            "misfit_error": misfit_error,
        }
    )
    print(df)
    return df


def get_summary(directory: Union[pathlib.Path, str] = pathlib.Path()) -> dict:
    """Get summary from FFM model

    :param directory: The directory to read files from, defaults to pathlib.Path()
    :type directory: Union[pathlib.Path, str], optional
    :return: The summary dictionary
    :rtype: dict
    """
    directory = pathlib.Path(directory)
    with open(directory / "modelling_summary.txt", "r") as infile:
        lines = [line.split() for line in infile]

    error1: Union[float, int] = 0
    error2: Union[float, int] = 0
    for line in lines:
        if "averaged" in line:
            error1 = float(line[-1])
        if "objective" in line:
            error2 = float(line[-1])
    errors = {"misfit_error": error1, "objective_error": error2}
    return errors
