# -*- coding: utf-8 -*-
"""
"""


import json
import logging
import os
import pathlib
import subprocess
from typing import List, Optional, Union

import wasp.modulo_logs as ml


def gf_retrieve(
    used_data_type: List[str],
    default_dirs: dict,
    directory: Union[pathlib.Path, str] = pathlib.Path(),
):
    """Compute and store (in binary files) Green functions for each station, both
    for teleseismic body waves, as for strong motion, cGPS and static data

    :param used_data_type: List of data types to process
    :type used_data_type: List[str]
    :param default_dirs: Dictionary of default directories
    :type default_dirs: dict
    :param directory: The directory to read/write from, defaults to pathlib.Path()
    :type directory: Union[pathlib.Path, str], optional
    :raises RuntimeError: If the fortran code produces any error
    """
    directory = pathlib.Path(directory)
    green_fun_tele = default_dirs["tele_gf"]
    green_fun_str = default_dirs["strong_motion_gf"]
    green_fun_gps = default_dirs["gps_gf"]

    processes: List[subprocess.Popen] = []
    loggers: List[logging.Logger] = []
    data_types: List[str] = []
    ch = logging.StreamHandler()
    formatter = logging.Formatter("%(levelname)s - %(message)s")
    ch.setFormatter(formatter)
    ch.setLevel(logging.ERROR)

    if "body" in used_data_type:
        print("Computing teleseismic GFs")
        logger1 = ml.create_log("body_wave_GF", directory / "logs" / "green_tele_log")
        logger1.addHandler(ch)
        p1 = subprocess.Popen(
            [green_fun_tele, f"{str(directory)}/"],
            stdout=subprocess.PIPE,
            stderr=subprocess.PIPE,
        )
        processes = processes + [p1]
        loggers = loggers + [logger1]
        data_types = data_types + ["body waves"]
    if "strong" in used_data_type:
        print("Computing strong motion GFs")
        logger2 = ml.create_log(
            "get_strong_motion_GF", directory / "logs" / "green_str_log"
        )
        logger2.addHandler(ch)
        p2 = subprocess.Popen(
            [green_fun_str, "str", f"{str(directory)}/"],
            stdout=subprocess.PIPE,
            stderr=subprocess.PIPE,
        )
        processes = processes + [p2]
        loggers = loggers + [logger2]
        data_types = data_types + ["strong motion"]
    if "cgps" in used_data_type:
        print("Computing cGPS GFs")
        logger3 = ml.create_log("get_cgps_GF", directory / "logs" / "green_cgps_log")
        logger3.addHandler(ch)
        p3 = subprocess.Popen(
            [green_fun_str, "cgps", f"{str(directory)}/"],
            stdout=subprocess.PIPE,
            stderr=subprocess.PIPE,
        )
        processes = processes + [p3]
        loggers = loggers + [logger3]
        data_types = data_types + ["cgps"]
    if "gps" in used_data_type:
        print("Computing static GPS GFs")
        logger4 = ml.create_log("GPS_GF", directory / "logs" / "green_gps_log")
        logger4.addHandler(ch)
        p4 = subprocess.Popen(
            [green_fun_gps, "gps", f"{str(directory)}/"],
            stdout=subprocess.PIPE,
            stderr=subprocess.PIPE,
        )
        processes = processes + [p4]
        loggers = loggers + [logger4]
        data_types = data_types + ["gps"]
    if "insar" in used_data_type:
        print("Computing InSAR GFs")
        logger5 = ml.create_log("INSAR_GF", directory / "logs" / "green_insar_log")
        logger5.addHandler(ch)
        p5 = subprocess.Popen(
            [green_fun_gps, "insar", f"{str(directory)}/"],
            stdout=subprocess.PIPE,
            stderr=subprocess.PIPE,
        )
        processes = processes + [p5]
        loggers = loggers + [logger5]
        data_types = data_types + ["insar"]

    for p, log, data_type in zip(processes, loggers, data_types):
        out, err = p.communicate(timeout=100000)
        log.info(out.decode("utf-8"))
        if err:
            log.error(err.decode("utf-8", "ignore"))
        ml.close_log(log)
        if err:
            raise RuntimeError(
                "Got following error while retrieving GF "
                "for {}:\n{}".format(data_type, err)
            )


def fk_green_fun1(
    data_prop: dict,
    tensor_info: dict,
    location: Union[pathlib.Path, str],
    cgps: bool = False,
    max_depth: Optional[Union[float, int]] = None,
    directory: Union[pathlib.Path, str] = pathlib.Path(),
) -> dict:
    """Write data for computing/retrieving strong motion Green's functions

    :param data_prop: The sampling filtering properties
    :type data_prop: dict
    :param tensor_info: The moment tensor information
    :type tensor_info: dict
    :param location: The file location
    :type location: Union[pathlib.Path, str]
    :param cgps: Whether cgps data is included, defaults to False
    :type cgps: bool, optional
    :param max_depth: The max depth, defaults to None
    :type max_depth: Optional[Union[float,int]], optional
    :param directory: The directory to read/write at, defaults to pathlib.Path()
    :type directory: Union[pathlib.Path, str], optional
    :return: the green dictionary
    :rtype: dict
    """
    directory = pathlib.Path(directory)
    sampling = data_prop["sampling"]
    dt = sampling["dt_strong"] if not cgps else sampling["dt_cgps"]
    depth = tensor_info["depth"]
    time_shift = tensor_info["time_shift"]
    min_depth = max(1, depth - 100)
    if max_depth is None:
        max_depth = max(30, 2 * depth)
        max_depth = min(max_depth, depth + 60)
    max_depth = max_depth + 5  # type:ignore
    min_dist = 0
    max_dist = 600 if time_shift < 40 else 1000
    time_corr = 10 if not cgps else 25

    green_dict = {
        "location": str(location),
        "min_depth": min_depth,
        "max_depth": max_depth,
        "min_dist": min_dist,
        "max_dist": max_dist,
        "dt": dt,
        "time_corr": time_corr,
    }

    name = "strong_motion_gf.json" if not cgps else "cgps_gf.json"
    with open(directory / name, "w") as f:
        json.dump(
            green_dict,
            f,
            sort_keys=True,
            indent=4,
            separators=(",", ": "),
            ensure_ascii=False,
        )
    return green_dict


if __name__ == "__main__":
    import argparse

    import wasp.manage_parser as mp
    import wasp.management as mng
    import wasp.seismic_tensor as tensor

    parser = argparse.ArgumentParser()
    parser.add_argument(
        "-f", "--folder", default=os.getcwd(), help="folder where there are input files"
    )
    parser = mp.parser_add_tensor(parser)
    parser = mp.parser_add_gf(parser)
    parser.add_argument(
        "-dt", type=float, default=0.2, help="sampling step of strong motion data"
    )
    args = parser.parse_args()
    os.chdir(args.folder)
    if args.gcmt_tensor:
        cmt_file = args.gcmt_tensor
        tensor_info = tensor.get_tensor(cmt_file=cmt_file)
    else:
        tensor_info = tensor.get_tensor()
    tensor_info["timedelta"] = 81 * 90
    used_data = mp.get_used_data(args)
    default_dirs = mng.default_dirs()
    # fk_green_fun0 removed since it isn't used anywhere
    # if "strong" in used_data and not os.path.isfile("strong_motion_gf.json"):
    #     green_dict = fk_green_fun0(args.dt, tensor_info, default_dirs)
    # #        write_green_file(green_dict)
    # if "cgps" in used_data and not os.path.isfile("cgps_gf.json"):
    #     green_dict = fk_green_fun0(1.0, tensor_info)
    # #        write_green_file(green_dict)
    gf_retrieve(used_data, default_dirs)
