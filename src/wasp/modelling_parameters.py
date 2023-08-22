#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""Module to find and store modelling parameters, such as number of iterations
of simmulated annealing and its cooling rate, weight of regularization
constraints, and boundaries of model space.
"""


import argparse
import json
import os
import pathlib
from typing import Optional, Tuple, Union

import numpy as np

import wasp.management as mng
import wasp.seismic_tensor as tensor


def modelling_prop(
    tensor_info: dict,
    segments_data: dict,
    data_type: list,
    moment_mag: Optional[float] = None,
    directory: Union[pathlib.Path, str] = pathlib.Path(),
) -> Tuple[dict, list]:
    """Write a file with input for performing FFM modelling using simulated annealing

    :param tensor_info: The tensor information
    :type tensor_info: dict
    :param segments_data: The segments information
    :type segments_data: dict
    :param data_type: The types of data available, defaults to None
    :type data_type: list
    :param moment_mag: The moment magnitude, defaults to None
    :type moment_mag: Optional[float], optional
    :param directory: The location of files to read/write, defaults to pathlib.Path()
    :type directory: Union[pathlib.Path, str], optional
    :return: The property dictionary, the updated segments information
    :rtype: Tuple[dict, list]
    """
    directory = pathlib.Path(directory)
    moment_mag = tensor_info["moment_mag"] if not moment_mag else moment_mag
    time_shift = tensor_info["time_shift"]
    data_type2 = set(data_type) - {"gps"}

    syn_len = int(2.5 * time_shift)
    #    factor = 1 / len(data_type2)
    moment_weight = 0.1
    slip_weight = 0.15
    time_weight = 0.15
    dictionary = {
        "initial_temperature": 0.01,
        "seismic_moment": moment_mag,
        "moment_weight": moment_weight,
        "slip_weight": slip_weight,
        "time_weight": time_weight,
        "max_source_dur": syn_len,
        "iterations": 250,
        "cooling_rate": 0.970,
    }

    with open(directory / "annealing_prop.json", "w") as f:
        json.dump(
            dictionary,
            f,
            sort_keys=True,
            indent=4,
            separators=(",", ": "),
            ensure_ascii=False,
        )

    moment_mag = tensor_info["moment_mag"]
    segments = segments_data["segments"]
    rise_time = segments_data["rise_time"]
    #
    # After Blaser et al (2010)
    #
    seismic_moment = moment_mag
    moment_mag = 2.0 * (np.log10(seismic_moment)) / 3.0 - 10.7  # type:ignore
    length = 10 ** (-2.31 + 0.57 * moment_mag - 0.2)  # type:ignore
    width = 10 ** (-1.56 + 0.41 * moment_mag - 0.17)  # type:ignore
    avg_slip = seismic_moment / length / width / (3 * 10**21)  # type:ignore
    peak_slip = 2 * avg_slip
    peak_slip = int(peak_slip / 100) * 100
    step = 20
    nstep = min(int(peak_slip / step) + 1, 51)

    dictionary2 = {
        "min_slip": 0,
        "max_upper_slip": peak_slip,
        "max_lower_slip": peak_slip,
        "max_left_slip": peak_slip,
        "max_right_slip": peak_slip,
        "max_center_slip": peak_slip,
        "max_slip_delta": peak_slip,
        "slip_step": nstep,
    }
    regularization = {
        "regularization": {
            "neighbour_up": None,
            "neighbour_down": None,
            "neighbour_left": None,
            "neighbour_right": None,
        }
    }

    segments2: list = []
    for segment in segments:
        rake = segment["rake"]
        rake_min = rake - 20
        rake_max = rake + 20
        rstep = 21
        if rake_min < 0:
            rake_min = rake_min + 360
            rake_max = rake_max + 360
        dictionary3 = {"rake_min": rake_min, "rake_max": rake_max, "rake_step": rstep}
        dictionary3.update(dictionary2)
        dictionary3.update(regularization)
        segments2 = segments2 + [dictionary3]

    with open(directory / "model_space.json", "w") as f:
        json.dump(
            segments2,
            f,
            sort_keys=True,
            indent=4,
            separators=(",", ": "),
            ensure_ascii=False,
        )
    return dictionary, segments2


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "-f", "--folder", default=os.getcwd(), help="folder where there are input files"
    )
    parser.add_argument(
        "-gcmt", "--gcmt_tensor", help="location of GCMT moment tensor file"
    )
    parser.add_argument(
        "-t", "--tele", action="store_true", help="compute files with teleseismic data"
    )
    parser.add_argument(
        "-su",
        "--surface",
        action="store_true",
        help="compute files with surface waves data",
    )
    parser.add_argument(
        "-st",
        "--strong",
        action="store_true",
        help="compute files with strong motion data",
    )
    parser.add_argument(
        "--cgps", action="store_true", help="compute files with cGPS data"
    )
    args = parser.parse_args()

    os.chdir(args.folder)
    if args.gcmt_tensor:
        cmt_file = args.gcmt_tensor
        tensor_info = tensor.get_tensor(cmt_file=cmt_file)
    else:
        tensor_info = tensor.get_tensor()
    data_type: list = []
    if args.tele:
        data_type = data_type + ["tele_body"]
    if args.surface:
        data_type = data_type + ["surf_tele"]
    if args.strong:
        data_type = data_type + ["strong_motion"]
    if args.cgps:
        data_type = data_type + ["cgps"]
    modelling_prop(tensor_info, data_type=data_type)  # type:ignore
