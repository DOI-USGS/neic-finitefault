# -*- coding: utf-8 -*-
"""Script for performing FFM modelling and forward, using the method of
Chen-Ji.
"""


import errno
import glob
import json
import logging
import os
import pathlib
import subprocess
import time
import warnings
from multiprocessing import Process
from shutil import copy2, move
from typing import List, Optional, Tuple, Union

import numpy as np

import wasp.data_management as dm
import wasp.data_processing as proc
import wasp.fault_plane as pf
import wasp.green_functions as gf
import wasp.management as mng
import wasp.modelling_parameters as mp
import wasp.modulo_logs as ml
import wasp.plot_graphic as plot
import wasp.seismic_tensor as tensor
import wasp.traces_properties as tp
import wasp.velocity_models as mv
from wasp import get_outputs, input_files
from wasp.load_ffm_model import load_ffm_model
from wasp.static2fsp import static_to_fsp


def automatic_usgs(
    tensor_info: dict,
    data_type: List[str],
    default_dirs: dict,
    velmodel: Optional[dict] = None,
    dt_cgps: Optional[float] = 1.0,
    st_response: bool = True,
    config_path: Optional[Union[str, pathlib.Path]] = None,
    directory: Union[pathlib.Path, str] = pathlib.Path(),
):
    """Routine for automatically running the FFM modelling

    :param tensor_info: The moment tensor information
    :type tensor_info: dict
    :param data_type: The data types available
    :type data_type: List[str]
    :param default_dirs: The location of default directories
    :type default_dirs: dict
    :param velmodel: The velocity model, defaults to None
    :type velmodel: Optional[dict], optional
    :param dt_cgps: The time delta for cgps data, defaults to 1.0
    :type dt_cgps: Optional[float], optional
    :param st_response: Whether to remove paz response of strong motion, defaults to True
    :type st_response: bool, optional
    :param config_path: The path to the config file, defaults to None
    :type config_path: Optional[Union[str, pathlib.Path]], optional
    :param directory: Where the file(s) should be read/written, defaults to pathlib.Path()
    :type directory: Union[pathlib.Path, str], optional
    """
    directory = pathlib.Path(directory)
    logger = ml.create_log(
        "automatic_ffm", os.path.join(directory, "logs", "automatic_ffm.log")
    )
    logger = ml.add_console_handler(logger)
    logger.info("Starting fff program")
    time0 = time.time()
    if "gps" in data_type:
        if os.path.isfile(os.path.join(directory, "data", "gps_data")):
            copy2(os.path.join(directory, "data", "gps_data"), directory)
    if "insar" in data_type:
        insar_files = glob.glob(os.path.join(directory, "data", "insar_a*txt"))
        insar_files = insar_files + glob.glob(
            os.path.join(directory, "data", "insar_d*txt")
        )
        for file in insar_files:
            if os.path.isfile(file):
                copy2(file, directory)
    data_dir = directory / "data"
    data_prop = tp.properties_json(
        tensor_info, dt_cgps=dt_cgps, data_directory=directory
    )
    time2 = time.time()
    logger.info("Process data")
    processing(
        tensor_info, data_type, data_prop, st_response=st_response, directory=data_dir
    )
    time2 = time.time() - time2
    logger.info("Time spent processing traces: {}".format(time2))
    data_folder = os.path.join(directory, "data")
    insar_asc = glob.glob(str(directory) + "/insar_asc*txt")
    insar_asc = None if len(insar_asc) == 0 else insar_asc  # type:ignore
    insar_desc = glob.glob(str(directory) + "/insar_desc*txt")
    insar_desc = None if len(insar_desc) == 0 else insar_desc  # type:ignore
    dm.filling_data_dicts(
        tensor_info,
        data_type,
        data_prop,
        data_folder,
        insar_asc=insar_asc,  # type:ignore
        insar_desc=insar_desc,  # type:ignore
        working_directory=directory,
    )
    writing_inputs0(
        tensor_info, data_type, config_path=config_path, directory=directory
    )
    logger.info("Compute GF bank")
    if not velmodel:
        velmodel = mv.select_velmodel(tensor_info, default_dirs, directory=directory)
    input_files.write_velmodel(velmodel, directory=directory)
    gf_bank_str = os.path.join(directory, "GF_strong")
    gf_bank_cgps = os.path.join(directory, "GF_cgps")
    get_gf_bank = default_dirs["strong_motion_gf_bank2"]
    if "cgps" in data_type:
        logger.info("Compute cGPS GF bank")
        green_dict = gf.fk_green_fun1(
            data_prop, tensor_info, gf_bank_cgps, cgps=True, directory=directory
        )
        input_files.write_green_file(green_dict, cgps=True, directory=directory)
        with open(os.path.join(directory, "logs", "GF_cgps_log"), "w") as out_gf_cgps:
            p1 = subprocess.Popen(
                [get_gf_bank, "cgps", f"{(directory)}/"], stdout=out_gf_cgps
            )
        p1.wait()
    if "strong_motion" in data_type:
        logger.info("Compute strong motion GF bank")
        green_dict = gf.fk_green_fun1(
            data_prop, tensor_info, gf_bank_str, directory=directory
        )
        input_files.write_green_file(green_dict, directory=directory)
        with open(
            os.path.join(directory, "logs", "GF_strong_log"), "w"
        ) as out_gf_strong:
            p2 = subprocess.Popen(
                [get_gf_bank, "strong", f"{(directory)}/"],
                stdout=out_gf_strong,
            )
        p2.wait()
    files = [
        directory / "Green_strong.txt",
        directory / "Green_cgps.txt",
        directory / "modelling_stats.json",
        directory / "gps_data",
        directory / "strong_motion_gf.json",
        directory / "cgps_gf.json",
        directory / "sampling_filter.json",
    ]
    files2 = glob.glob(str(directory) + "/channels_*txt")
    files3 = glob.glob(str(directory) + "/wavelets_*txt")
    files4 = glob.glob(str(directory) + "/waveforms_*txt")
    files5 = glob.glob(str(directory) + "/*waves.json")
    files6 = glob.glob(str(directory) + "/static*")
    files7 = glob.glob(str(directory) + "/filtro*") + glob.glob(
        str(directory) + "/surf_filter*"
    )
    files8 = [
        directory / "instrumental_response.txt",
        directory / "body_wave_weight.txt",
    ]
    files9 = glob.glob(str(directory) + "/insar*")
    files = (
        files
        + files2
        + files3
        + files4
        + files5
        + files6
        + files7
        + files8
        + files9  # type:ignore
    )
    folders = [directory / "NP1", directory / "NP2"]
    for folder in folders:
        for file in files:  # type:ignore
            if os.path.isfile(file):
                copy2(file, folder)
    info_np1, info_np2 = tensor.planes_from_tensor(tensor_info)
    plane1_folder = directory / "NP1"
    keywords = {"velmodel": velmodel, "directory": plane1_folder}
    p1 = Process(  # type:ignore
        target=_automatic2,
        args=(
            tensor_info,
            info_np1,
            data_type,
            data_prop,
            default_dirs,
            logger,
        ),
        kwargs=keywords,
    )
    p1.start()  # type:ignore
    plane2_folder = directory / "NP2"
    keywords["directory"] = plane2_folder
    p2 = Process(  # type:ignore
        target=_automatic2,
        args=(
            tensor_info,
            info_np2,
            data_type,
            data_prop,
            default_dirs,
            logger,
        ),
        kwargs=keywords,
    )
    p2.start()  # type:ignore
    [p.join() for p in [p1, p2]]  # type:ignore
    logger.info("Time spent: {}".format(time.time() - time0))
    ml.close_log(logger)
    return


def _automatic2(
    tensor_info: dict,
    plane_data: dict,
    data_type: List[str],
    data_prop: dict,
    default_dirs: dict,
    logger: logging.Logger,
    velmodel: Optional[dict] = None,
    directory: pathlib.Path = pathlib.Path(),
):
    """Routine for automatic FFM modelling for each nodal plane

    :param tensor_info: The moment tensor information
    :type tensor_info: dict
    :param plane_data: The properties of the plane
    :type plane_data: dict
    :param data_type: The data types available
    :type data_type: List[str]
    :param data_prop: The data properties
    :type data_prop: dict
    :param default_dirs: The location of default directories
    :type default_dirs: dict
    :param logger: The logger used to log information
    :type logger: logging.Logger
    :param velmodel: The velocity model, defaults to None
    :type velmodel: Optional[dict], optional
    :param directory: Where the file(s) should be read/written, defaults to pathlib.Path()
    :type directory: Union[pathlib.Path, str], optional
    """
    #
    # Create JSON files
    #
    logger.info("Create input files for Fortran scripts")
    logger.info("Create automatic JSON")
    tensor.write_tensor(tensor_info, directory=directory)
    if velmodel:
        mv.velmodel2json(velmodel, directory=directory)
    if not velmodel:
        velmodel = mv.select_velmodel(tensor_info, default_dirs, directory=directory)
    np_plane_info = plane_data["plane_info"]
    data_folder = os.path.join(directory.parent.parent, "data")
    insar_asc = glob.glob(str(directory) + "/insar_asc*txt")
    insar_asc = None if len(insar_asc) == 0 else insar_asc  # type:ignore
    insar_desc = glob.glob(str(directory) + "/insar_desc*txt")
    insar_desc = None if len(insar_desc) == 0 else insar_desc  # type:ignore
    dm.filling_data_dicts(
        tensor_info,
        data_type,
        data_prop,
        data_folder,
        insar_asc=insar_asc,  # type:ignore
        insar_desc=insar_desc,  # type:ignore
        working_directory=directory,
    )
    segments_data = pf.create_finite_fault(
        tensor_info, np_plane_info, data_type, directory=directory
    )
    segments = segments_data["segments"]
    rise_time = segments_data["rise_time"]
    connections = None
    if "connections" in segments_data:
        connections = segments_data["connections"]
    point_sources = pf.point_sources_param(
        segments, tensor_info, rise_time, connections=connections
    )
    data_type = _check_surf_GF(point_sources, data_type, logger=logger)
    mp.modelling_prop(
        tensor_info, segments_data, data_type=data_type, directory=directory
    )
    #
    # write text files from JSONs
    #
    rupt_vel = segments[0]["rupture_vel"]
    lambda_min = 0.4
    lambda_max = 1.25
    min_vel, max_vel = [lambda_min * rupt_vel, lambda_max * rupt_vel]
    logger.info("Write input files")
    writing_inputs(
        tensor_info, data_type, segments_data, min_vel, max_vel, directory=directory
    )
    #
    # Modelling and plotting results
    #
    inversion(data_type, default_dirs, logger, directory=directory)
    logger.info("Plot data in folder {}".format(directory))
    execute_plot(
        tensor_info,
        data_type,
        segments_data,
        default_dirs,
        velmodel=velmodel,
        directory=directory,
    )
    #
    # write solution in FSP format
    #
    logging.info
    solution = get_outputs.read_solution_static_format(segments, data_dir=directory)
    static_to_fsp(
        tensor_info, segments_data, data_type, velmodel, solution, directory=directory
    )
    for file in glob.glob(str(directory) + "/*png"):
        if os.path.isfile(os.path.join(directory, file)):
            copy2(os.path.join(directory, file), os.path.join(directory, "plots"))


def _check_surf_GF(
    point_sources: list, used_data: List[str], logger: Optional[logging.Logger] = None
) -> List[str]:
    """Check maximum fault depth to determine whether the surface wave data can be used (No surface wave GFs below 125km)
    :param point_sources: The location of point sources
    :type point_sources: list
    :param data_type: The data types available
    :type used_data: List[str]
    :param logger: The logger used to log information, defaults to None
    :type logger: logging.Logger, optional
    :return: The updated data types list
    :rtype: List[str]
    """
    new_used_data = used_data.copy()
    depths = [ps[:, :, :, :, 2] for ps in point_sources]
    depths = [np.max(depths1) for depths1 in depths]
    max_depth = np.max(depths)
    is_surf = "surf_tele" in used_data
    if max_depth > 125 and is_surf:
        warnings.warn("Fault plane extends below 125 km depth limit for surface wave Green's functions. Surface waves won't be used.")
        new_used_data.remove("surf_tele")
        if logger:
            logger.info("Maximum depth larger than 125 km.")
    return new_used_data


def modelling_new_data(
    tensor_info: dict,
    data_type: List[str],
    default_dirs: dict,
    data_folder: Union[pathlib.Path, str],
    segments_data: dict,
    st_response: bool = True,
    directory: Union[pathlib.Path, str] = pathlib.Path(),
):
    """Routine for manual finite fault modelling with new data types

    :param tensor_info: The moment tensor information
    :type tensor_info: dict
    :param data_type: The data types available
    :type data_type: List[str]
    :param default_dirs: The location of default directories
    :type default_dirs: dict
    :param data_folder: The location of the data used in modelling
    :type data_folder: Union[pathlib.Path, str]
    :param segments_data: The segments properties
    :type segments_data: dict
    :param st_response: Whether to remove paz response of strong motion, defaults to True
    :type st_response: bool, optional
    :param directory: Where the file(s) should be read/written, defaults to pathlib.Path()
    :type directory: Union[pathlib.Path, str], optional
    """
    directory = pathlib.Path(directory)
    if os.path.isfile(os.path.join(data_folder, "gps_data")):
        copy2(os.path.join(data_folder, "gps_data"), directory)
    insar_asc = glob.glob(os.path.join(data_folder, "insar_a*txt"))
    insar_desc = glob.glob(os.path.join(data_folder, "insar_d*txt"))
    insar_files = insar_asc + insar_desc
    for file in insar_files:
        if os.path.isfile(file):
            copy2(file, directory)
    with open(directory / "sampling_filter.json") as sf:
        data_prop = json.load(sf)
    time2 = time.time()
    processing(
        tensor_info, data_type, data_prop, st_response=st_response, directory=directory
    )
    dm.filling_data_dicts(
        tensor_info,
        data_type,
        data_prop,
        data_folder,
        insar_asc=insar_asc,  # type:ignore
        insar_desc=insar_desc,  # type:ignore
        working_directory=directory,
    )
    gf_bank_str = os.path.join(directory, "GF_strong")
    gf_bank_cgps = os.path.join(directory, "GF_cgps")
    get_gf_bank = default_dirs["strong_motion_gf_bank2"]
    segments = segments_data["segments"]
    rise_time = segments_data["rise_time"]
    connections = None
    if "connections" in segments_data:
        connections = segments_data["connections"]
    point_sources = pf.point_sources_param(
        segments, tensor_info, rise_time, connections=connections
    )
    depths = [ps[:, :, :, :, 2] for ps in point_sources]
    max_depths = [np.max(depth1.flatten()) for depth1 in depths]
    max_depth = np.max(max_depths)
    if "cgps" in data_type:
        green_dict = gf.fk_green_fun1(
            data_prop,
            tensor_info,
            gf_bank_cgps,
            max_depth=max_depth,
            cgps=True,
            directory=directory,
        )
        input_files.write_green_file(green_dict, cgps=True, directory=directory)
        with open(os.path.join(directory / "logs", "GF_cgps_log"), "w") as out_gf_cgps:
            p1 = subprocess.Popen(
                [get_gf_bank, "cgps", f"{(directory)}/"], stdout=out_gf_cgps
            )
        p1.wait()
    if "strong_motion" in data_type:
        green_dict = gf.fk_green_fun1(
            data_prop,
            tensor_info,
            gf_bank_str,
            max_depth=max_depth,
            directory=directory,
        )
        input_files.write_green_file(green_dict, directory=directory)
        with open(
            os.path.join(directory / "logs", "GF_strong_log"), "w"
        ) as out_gf_strong:
            p2 = subprocess.Popen(
                [get_gf_bank, "strong", f"{(directory)}/"],
                stdout=out_gf_strong,
            )
        p2.wait()
    data_type2: list = []
    if os.path.isfile(directory / "tele_waves.json"):
        data_type2 = data_type2 + ["tele_body"]
    if os.path.isfile(directory / "surf_waves.json"):
        data_type2 = data_type2 + ["surf_tele"]
    if os.path.isfile(directory / "strong_motion_waves.json"):
        data_type2 = data_type2 + ["strong_motion"]
    if os.path.isfile(directory / "cgps_waves.json"):
        data_type2 = data_type2 + ["cgps"]
    if os.path.isfile(directory / "static_data.json"):
        data_type2 = data_type2 + ["gps"]
    if os.path.isfile(directory / "insar_data.json"):
        data_type2 = data_type2 + ["insar"]
    if os.path.isfile(directory / "dart_waves.json"):
        data_type2 = data_type2 + ["dart"]
    manual_modelling(
        tensor_info, data_type2, default_dirs, segments_data, directory=directory
    )
    return


def manual_modelling(
    tensor_info: dict,
    data_type: List[str],
    default_dirs: dict,
    segments_data: dict,
    directory: Union[pathlib.Path, str] = pathlib.Path(),
):
    """Routine for manual finite fault modelling

    :param tensor_info: The moment tensor information
    :type tensor_info: dict
    :param data_type: The data types available
    :type data_type: List[str]
    :param default_dirs: The location of default directories
    :type default_dirs: dict
    :param segments_data: The segments properties
    :type segments_data: dict
    :param directory: Where the file(s) should be read/written, defaults to pathlib.Path()
    :type directory: Union[pathlib.Path, str], optional
    """
    directory = pathlib.Path(directory)
    if not os.path.isdir(directory / "logs"):
        os.mkdir(directory / "logs")
    if not os.path.isdir(directory / "plots"):
        os.mkdir(directory / "plots")
    min_vel, max_vel = __ask_velrange(directory=directory)
    logger = ml.create_log(
        "manual_ffm", os.path.join(directory, "logs", "manual_ffm.log")
    )
    logger.info("Write input files")
    tensor.write_tensor(tensor_info, directory=directory)
    writing_inputs(
        tensor_info, data_type, segments_data, min_vel, max_vel, directory=directory
    )
    writing_inputs0(tensor_info, data_type, directory=directory)
    inversion(data_type, default_dirs, logger, directory=directory)
    logger.info("Plot data in folder {}".format(directory))
    execute_plot(
        tensor_info, data_type, segments_data, default_dirs, directory=directory
    )
    ml.close_log(logger)


def forward_modelling(
    tensor_info: dict,
    data_type: List[str],
    default_dirs: dict,
    segments_data: dict,
    option: str = "Solucion.txt",
    max_slip: Union[float, int] = 200,
    directory: Union[pathlib.Path, str] = pathlib.Path(),
):
    """Routine for forward modelling

    :param tensor_info: The moment tensor information
    :type tensor_info: dict
    :param data_type: The data types available
    :type data_type: List[str]
    :param default_dirs: The location of default directories
    :type default_dirs: dict
    :param segments_data: The segments properties
    :type segments_data: dict
    :param option: Name of the file with the kinematic model, defaults to "Solucion.txt"
    :type option: str, optional
    :param max_slip: Maximum slip in case of checkerboard test, defaults to 200
    :type max_slip: Union[float, int], optional
    :param directory: Where the file(s) should be read/written, defaults to pathlib.Path()
    :type directory: Union[pathlib.Path, str], optional
    :raises FileNotFoundError: If the velocity model file is not found
    """
    directory = pathlib.Path(directory)
    tensor.write_tensor(tensor_info, directory=directory)
    if not os.path.isdir(directory / "logs"):
        os.mkdir(directory / "logs")
    if not os.path.isdir(directory / "plots"):
        os.mkdir(directory / "plots")
    len_stk = 5 if not option == "point_source" else 8
    len_dip = 5 if not option == "point_source" else 1
    segments = segments_data["segments"]
    rise_time = segments_data["rise_time"]
    connections = None
    if "connections" in segments_data:
        connections = segments_data["connections"]
    point_sources = pf.point_sources_param(
        segments, tensor_info, rise_time, connections=connections
    )
    #
    # Get input model
    #
    model = load_ffm_model(
        segments_data,
        point_sources,
        option=option,
        max_slip=max_slip,
        len_stk=len_stk,
        len_dip=len_dip,
        directory=directory,
    )
    if not os.path.isfile(directory / "velmodel_data.json"):
        raise FileNotFoundError(
            errno.ENOENT, os.strerror(errno.ENOENT), "velmodel_data.json"
        )
    with open(directory / "velmodel_data.json") as vm:
        velmodel = json.load(vm)
    min_vel, max_vel = __ask_velrange(directory=directory)

    logger = ml.create_log(
        "forward_model", os.path.join(directory, "logs", "forward_model.log")
    )
    logger.info("Write input files")
    shear = pf.shear_modulous(point_sources, velmodel=velmodel)
    dx = segments[0]["delta_strike"]
    dy = segments[0]["delta_dip"]
    slip = model["slip"]
    zipped = zip(slip, shear)
    moment_sub = [dx * dy * slip_seg * shear_seg for slip_seg, shear_seg in zipped]
    moment = np.sum([np.sum(moment_seg.flatten()) for moment_seg in moment_sub])
    moment = 10**10 * moment
    writing_inputs(
        tensor_info,
        data_type,
        segments_data,
        min_vel,
        max_vel,
        moment_mag=moment,
        forward_model=model,
        directory=directory,
    )
    inversion(data_type, default_dirs, logger, forward=True, directory=directory)
    logger.info("Plot data in folder {}".format(directory))
    execute_plot(
        tensor_info,
        data_type,
        segments_data,
        default_dirs,
        velmodel=velmodel,
        directory=directory,
    )
    ml.close_log(logger)


def checkerboard(
    tensor_info: dict,
    data_type: List[str],
    default_dirs: dict,
    segments_data: dict,
    max_slip: Union[float, int] = 200,
    add_error: bool = False,
    option: str = "Checkerboard",
    option2: str = "FFM modelling",
    directory: Union[pathlib.Path, str] = pathlib.Path(),
):
    """Routine for running checkerboard tests

    :param tensor_info: The moment tensor information
    :type tensor_info: dict
    :param data_type: The data types available
    :type data_type: List[str]
    :param default_dirs: The location of default directories
    :type default_dirs: dict
    :param segments_data: The segments properties
    :type segments_data: dict
    :param max_slip: Maximum slip in case of checkerboard test, defaults to 200
    :type max_slip: Union[float, int], optional
    :param add_error: Whether to add noise to synthetic waveforms, defaults to False
    :type add_error: bool, optional
    :param option: Name of the file with the kinematic model, defaults to "Checkerboard"
    :type option: str, optional
    :param option2: Whether to invert the checkerboard model or not, defaults to "FFM modelling"
    :type option2: str, optional
    :param directory: Where the file(s) should be read/written, defaults to pathlib.Path()
    :type directory: Union[pathlib.Path, str], optional
    """
    directory = pathlib.Path(directory)
    if max_slip > 0:
        folder_name = directory / "checkerboard_resolution"
    else:
        folder_name = directory / "checkerboard_noise"
    if not option == "Checkerboard":
        folder_name = directory / option
    if not os.path.isdir(folder_name):
        os.mkdir(folder_name)
    if not os.path.isdir(directory / "logs"):
        os.mkdir(directory / "logs")
    if not os.path.isdir(directory / "plots"):
        os.mkdir(directory / "plots")
    for file in os.listdir():
        if os.path.isfile(file):
            copy2(file, folder_name)
    segments = segments_data["segments"]
    rise_time = segments_data["rise_time"]
    connections = None
    if "connections" in segments_data:
        connections = segments_data["connections"]
    point_sources = pf.point_sources_param(
        segments, tensor_info, rise_time, connections=connections
    )
    forward_modelling(
        tensor_info,
        data_type,
        default_dirs,
        segments_data,
        option=option,
        max_slip=max_slip,
        directory=folder_name,
    )
    with open(folder_name / "sampling_filter.json") as sf:
        data_prop = json.load(sf)
    for data_type0 in data_type:
        if data_type0 == "tele_body":
            json_dict = "tele_waves.json"
        if data_type0 == "surf_tele":
            json_dict = "surf_waves.json"
        if data_type0 == "strong_motion":
            json_dict = "strong_motion_waves.json"
        if data_type0 == "cgps":
            json_dict = "cgps_waves.json"
        if data_type0 == "gps":
            json_dict = "static_data.json"
        if data_type0 == "dart":
            json_dict = "dart_waves.json"
        with open(folder_name / json_dict) as jd:
            files = json.load(jd)
        input_files.from_synthetic_to_obs(
            files,
            data_type0,
            data_prop,
            add_error=add_error,
        )
    logger = ml.create_log(
        "checkerboard_ffm", os.path.join(folder_name, "logs", "checkerboard_ffm.log")
    )
    if option2 == "FFM modelling":
        inversion(data_type, default_dirs, logger, directory=folder_name)
        execute_plot(
            tensor_info,
            data_type,
            segments_data,
            default_dirs,
            plot_input=True,
            directory=folder_name,
        )
    ml.close_log(logger)


def set_directory_structure(
    tensor_info: dict, directory: Union[pathlib.Path, str] = pathlib.Path()
):
    """Create directory structure

    :param tensor_info: The moment tensor information
    :type tensor_info: dict
    :param directory: Where the file(s) should be read/written, defaults to pathlib.Path()
    :type directory: Union[pathlib.Path, str], optional
    """
    directory = pathlib.Path(directory)
    sol_folder = directory / mng.start_time_id(tensor_info)
    if not os.path.isdir(sol_folder):
        os.mkdir(sol_folder)
    version = len(glob.glob(os.path.join(sol_folder, "ffm*")))
    sol_folder2 = os.path.join(sol_folder, "ffm.{}".format(version))
    os.mkdir(sol_folder2)
    folders = [
        "data",
        os.path.join("data", "cGPS"),
        os.path.join("data", "STR"),
        os.path.join("data", "P"),
        os.path.join("data", "SH"),
        os.path.join("data", "LONG"),
        os.path.join("data", "final"),
        os.path.join("data", "final_r"),
        "NP1",
        os.path.join("NP1", "logs"),
        os.path.join("NP1", "plots"),
        "NP2",
        os.path.join("NP2", "logs"),
        os.path.join("NP2", "plots"),
        "logs",
        "plots",
        os.path.join("plots", "NP1"),
        os.path.join("plots", "NP2"),
    ]
    for folder in folders:
        new_folder = os.path.join(sol_folder2, folder)
        if os.path.exists(new_folder):
            continue
        os.mkdir(new_folder)
    return


def processing(
    tensor_info: dict,
    data_type: List[str],
    data_prop: dict,
    st_response: bool = True,
    directory: Union[pathlib.Path, str] = pathlib.Path(),
):
    """Run all waveform data processing

    :param tensor_info: The moment tensor information
    :type tensor_info: dict
    :param data_type: The data types available
    :type data_type: List[str]
    :param data_prop: The data properties
    :type data_prop: dict
    :param st_response: Whether to remove paz response of strong motion, defaults to True
    :type st_response: bool, optional
    :param directory: Where the file(s) should be read/written, defaults to pathlib.Path()
    :type directory: Union[pathlib.Path, str], optional
    """
    directory = pathlib.Path(directory)
    tele_files = (
        glob.glob(str(directory) + "/*.BH*SAC")
        + glob.glob(str(directory) + "/*.BH*sac")
        + glob.glob(str(directory) + "/*_BH*sac")
        + glob.glob(str(directory) + "/*_BH*sac")
    )
    strong_files = (
        glob.glob(str(directory) + "/*.HN*SAC")
        + glob.glob(str(directory) + "/*.HL*SAC")
        + glob.glob(str(directory) + "/*.HN*sac")
        + glob.glob(str(directory) + "/*.HL*sac")
        + glob.glob(str(directory) + "/*.AH?.*")
        + glob.glob(str(directory) + "/*_HN*sac")
        + glob.glob(str(directory) + "/*_HL*sac")
        + glob.glob(str(directory) + "/*HG*sac")
    )
    cgps_files = glob.glob(str(directory) + "/*L[HXY]*sac") + glob.glob(
        str(directory) + "/*L[HXY]*SAC"
    )
    if "tele_body" in data_type:
        proc.select_process_tele_body(
            tele_files, tensor_info, data_prop, directory=directory
        )
    if "surf_tele" in data_type:
        proc.select_process_surf_tele(
            tele_files, tensor_info, data_prop, directory=directory
        )
    if "strong_motion" in data_type:
        proc.select_process_strong(
            strong_files,
            tensor_info,
            data_prop,
            remove_response=st_response,
            directory=directory,
        )
    if "cgps" in data_type:
        proc.select_process_cgps(
            cgps_files, tensor_info, data_prop, directory=directory
        )


def writing_inputs0(
    tensor_info: dict,
    data_type: List[str],
    config_path: Optional[Union[str, pathlib.Path]] = None,
    directory: Union[pathlib.Path, str] = pathlib.Path(),
):
    """Write the input files

    :param tensor_info: The moment tensor information
    :type tensor_info: dict
    :param data_type: The data types available
    :type data_type: List[str]
    :param config_path: The path to the config file, defaults to None
    :type config_path: Optional[Union[str, pathlib.Path]], optional
    :param directory: Where the file(s) should be read/written, defaults to pathlib.Path()
    :type directory: Union[pathlib.Path, str], optional
    :raises FileNotFoundError: If sampling_filter.json cannot be found
    """
    directory = pathlib.Path(directory)
    if not os.path.isfile(directory / "sampling_filter.json"):
        raise FileNotFoundError(
            errno.ENOENT, os.strerror(errno.ENOENT), "sampling_filter.json"
        )
    with open(directory / "sampling_filter.json") as dp:
        data_prop = json.load(dp)
    if "tele_body" in data_type:
        input_files.input_chen_tele_body(tensor_info, data_prop, directory=directory)
    if "surf_tele" in data_type:
        input_files.input_chen_tele_surf(
            tensor_info, data_prop, config_path=config_path, directory=directory
        )
    if "strong_motion" in data_type:
        input_files.input_chen_near_field(
            tensor_info, data_prop, "strong_motion", directory=directory
        )
    if "cgps" in data_type:
        input_files.input_chen_near_field(
            tensor_info, data_prop, "cgps", directory=directory
        )
    if "gps" in data_type:
        input_files.input_chen_static(directory=directory)
    if "insar" in data_type:
        input_files.input_chen_insar(directory=directory)
    if "dart" in data_type:
        input_files.input_chen_dart(tensor_info, data_prop, directory=directory)


def writing_inputs(
    tensor_info: dict,
    data_type: List[str],
    segments_data: dict,
    min_vel: float,
    max_vel: float,
    moment_mag: Optional[float] = None,
    forward_model: Optional[dict] = None,
    directory: Union[pathlib.Path, str] = pathlib.Path(),
):
    """Write all required text files from the information found in the JSONs

    :param tensor_info: The moment tensor information
    :type tensor_info: dict
    :param data_type: The data types available
    :type data_type: List[str]
    :param segments_data: The segments properties
    :type segments_data: dict
    :param min_vel: The minimum rupture velocity
    :type min_vel: float
    :param max_vel: The maximum rupture velocity
    :type max_vel: float
    :param moment_mag: The input seismic moment, defaults to None
    :type moment_mag: Optional[float], optional
    :param forward_model: The input kinematic model, defaults to None
    :type forward_model: Optional[np.ndarray], optional
    :param directory: Where the file(s) should be read/written, defaults to pathlib.Path()
    :type directory: Union[pathlib.Path,str], optional
    :raises FileNotFoundError: Whether a property json file cannot be found
    """
    directory = pathlib.Path(directory)
    if not os.path.isfile(directory / "velmodel_data.json"):
        raise FileNotFoundError(
            errno.ENOENT, os.strerror(errno.ENOENT), "velmodel_data.json"
        )
    with open(directory / "velmodel_data.json") as vm:
        velmodel = json.load(vm)
    input_files.write_velmodel(velmodel, directory=directory)
    input_files.plane_for_chen(
        tensor_info, segments_data, min_vel, max_vel, velmodel, directory=directory
    )
    if forward_model:
        input_files.forward_model(
            tensor_info,
            segments_data,
            forward_model,
            min_vel,
            max_vel,
            directory=directory,
        )
    if not os.path.isfile(directory / "annealing_prop.json"):
        raise FileNotFoundError(
            errno.ENOENT, os.strerror(errno.ENOENT), "annealing_prop.json"
        )
    with open(directory / "annealing_prop.json") as ap:
        dictionary = json.load(ap)
    if moment_mag:
        dictionary["seismic_moment"] = moment_mag
    input_files.inputs_simmulated_annealing(dictionary, directory=directory)
    if "cgps" in data_type:
        if not os.path.isfile(directory / "cgps_gf.json"):
            raise FileNotFoundError(
                errno.ENOENT, os.strerror(errno.ENOENT), "cgps_gf.json"
            )
        with open(directory / "cgps_gf.json") as cg:
            green_dict = json.load(cg)
        input_files.write_green_file(green_dict, cgps=True, directory=directory)
    if "strong_motion" in data_type:
        if not os.path.isfile(directory / "strong_motion_gf.json"):
            raise FileNotFoundError(
                errno.ENOENT, os.strerror(errno.ENOENT), "strong_motion_gf.json"
            )
        with open(directory / "strong_motion_gf.json") as sm:
            green_dict = json.load(sm)
        input_files.write_green_file(green_dict, directory=directory)
    if not os.path.isfile(directory / "model_space.json"):
        raise FileNotFoundError(
            errno.ENOENT, os.strerror(errno.ENOENT), "model_space.json"
        )
    with open(directory / "model_space.json") as ms:
        segments2 = json.load(ms)
    input_files.model_space(segments2, directory=directory)
    return


def inversion(
    data_type: List[str],
    default_dirs: dict,
    logger: logging.Logger,
    forward: bool = False,
    directory: Union[pathlib.Path, str] = pathlib.Path(),
) -> dict:
    """Get the binaries with gf for each station, run the ffm code, and
    proceed to plot the results

    :param data_type: The data types available
    :type data_type: List[str]
    :param default_dirs: The location of default directories
    :type default_dirs: dict
    :param logger: The logger to log to
    :type logger: logging.Logger
    :param forward: Whether to forward model, defaults to False
    :type forward: bool, optional
    :param directory:  Where the file(s) should be read/written, defaults to pathlib.Path()
    :type directory: Union[pathlib.Path,str], optional
    :return: The run's statistics
    :rtype: dict
    """
    directory = pathlib.Path(directory)
    logger.info("Green_functions")
    time1 = time.time()
    gf.gf_retrieve(data_type, default_dirs, directory=directory)
    time1 = time.time() - time1
    run_stats = dict()
    run_stats["gf_time"] = time1
    logger.info("Elapsed time of green_fun: {}".format(time1))
    time3 = time.time()
    args = ["auto"]
    args = args + ["strong"] if "strong_motion" in data_type else args
    args = args + ["cgps"] if "cgps" in data_type else args
    args = args + ["body"] if "tele_body" in data_type else args
    args = args + ["surf"] if "surf_tele" in data_type else args
    args = args + ["gps"] if "gps" in data_type else args
    args = args + ["dart"] if "dart" in data_type else args
    args = args + ["insar"] if "insar" in data_type else args
    if not forward:
        logger.info("Inversion at folder {}".format(directory))
        finite_fault = default_dirs["finite_fault"]
    else:
        logger.info("Forward at folder {}".format(directory))
        finite_fault = default_dirs["forward"]
    working_dir = str(pathlib.Path().resolve())
    os.chdir(directory)
    p1 = subprocess.Popen([finite_fault, *args], stderr=subprocess.PIPE)
    os.chdir(working_dir)
    #
    # need to wait until FFM modelling is finished.
    #
    outs, errs = p1.communicate(timeout=40 * 60)
    if errs:
        logger.error(errs.decode("utf-8"))
        raise Exception(errs)
    time3 = time.time() - time3
    logger.info("Elapsed time of finite fault modelling: {}".format(time3))
    run_stats["ffm_time"] = time3
    delete_binaries(directory=directory)
    return run_stats


def execute_plot(
    tensor_info: dict,
    data_type: List[str],
    segments_data: dict,
    default_dirs: dict,
    velmodel: Optional[dict] = None,
    plot_input: bool = False,
    directory: Union[pathlib.Path, str] = pathlib.Path(),
):
    """Plot modelling results

    :param tensor_info: The moment tensor information
    :type tensor_info: dict
    :param data_type: The data types available
    :type data_type: List[str]
    :param segments_data: The segments properties
    :type segments_data: dict
    :param default_dirs: The location of default directories
    :type default_dirs: dict
    :param velmodel: The velocity model, defaults to None
    :type velmodel: Optional[dict], optional
    :param plot_input: Whether to plot initial kinematic model as well, defaults to False
    :type plot_input: bool, optional
    :param directory: Where the file(s) should be read/written, defaults to pathlib.Path()
    :type directory: Union[pathlib.Path,str], optional
    """
    directory = pathlib.Path(directory)
    print("Plot results")
    segments = segments_data["segments"]
    rise_time = segments_data["rise_time"]
    connections = None
    if "connections" in segments_data:
        connections = segments_data["connections"]
    solution = get_outputs.read_solution_static_format(segments, data_dir=directory)
    if not velmodel:
        velmodel = mv.select_velmodel(tensor_info, default_dirs, directory=directory)
    point_sources = pf.point_sources_param(
        segments, tensor_info, rise_time, connections=connections
    )
    shear = pf.shear_modulous(point_sources, velmodel=velmodel)
    use_waveforms = mng.use_waveforms(data_type)
    plot.plot_ffm_sol(
        tensor_info,
        segments_data,
        point_sources,
        shear,
        solution,
        velmodel,
        default_dirs,
        use_waveforms=use_waveforms,
        directory=directory,
    )
    plot.plot_misfit(data_type, directory=directory)
    traces_info, stations_gps = [None, None]
    if "strong_motion" in data_type:
        with open(directory / "strong_motion_waves.json") as smw:
            traces_info = json.load(smw)
    if "gps" in data_type:
        names, lats, lons, observed, synthetic, error = get_outputs.retrieve_gps(
            directory=directory
        )
        stations_gps = zip(names, lats, lons, observed, synthetic, error)
    if "strong_motion" in data_type or "gps" in data_type:
        plot._PlotMap(
            tensor_info,
            segments,
            point_sources,
            solution,
            default_dirs,
            files_str=traces_info,
            stations_gps=stations_gps,
            directory=directory,
        )
    if "insar" in data_type:
        insar_data = get_outputs.get_insar(data_dir=directory)
        if "ascending" in insar_data:
            asc_properties = insar_data["ascending"]
            for i, asc_property in enumerate(asc_properties):
                insar_points = asc_property["points"]
                plot._PlotInsar(
                    tensor_info,
                    segments,
                    point_sources,
                    default_dirs,
                    insar_points,
                    los="ascending{}".format(i),
                    directory=directory,
                )
        if "descending" in insar_data:
            desc_properties = insar_data["descending"]
            for i, desc_property in enumerate(desc_properties):
                insar_points = desc_property["points"]
                plot._PlotInsar(
                    tensor_info,
                    segments,
                    point_sources,
                    default_dirs,
                    insar_points,
                    los="descending{}".format(i),
                    directory=directory,
                )
    if plot_input:
        input_model = load_ffm_model(
            segments_data,
            point_sources,
            option="fault&rise_time.txt",
            directory=directory,
        )
        plot._PlotSlipDist_Compare(
            segments, point_sources, input_model, solution, directory=directory
        )
        plot._PlotComparisonMap(
            tensor_info,
            segments,
            point_sources,
            input_model,
            solution,
            directory=directory,
        )
    plot_files = glob.glob(os.path.join(directory, "plots", "*png"))
    for plot_file in plot_files:
        os.remove(plot_file)
    plot_files = glob.glob(str(directory) + "/*png")
    for plot_file in plot_files:
        move(plot_file, directory / "plots")


def delete_binaries(directory: Union[pathlib.Path, str] = pathlib.Path()):
    """Remove the files with Green function data

    :param directory: Where the file(s) should be read/written, defaults to pathlib.Path()
    :type directory: Union[pathlib.Path,str], optional
    """
    directory = pathlib.Path(directory)
    deletables = glob.glob("*.GRE") + glob.glob("*.TDE") + glob.glob("*[1-2-3]")
    deletables = deletables + glob.glob("*.H[LN][E-N-Z]")
    deletables = deletables + glob.glob("*.L[HXY][E-N-Z]")
    for file in deletables:
        if os.path.isfile(file):
            os.remove(file)


def __ask_velrange(
    directory: Union[pathlib.Path, str] = pathlib.Path()
) -> Tuple[float, float]:
    """Get the velocity range from fault&rise_time.txt

    :param directory: Where the file(s) should be read/written, defaults to pathlib.Path()
    :type directory: Union[pathlib.Path,str], optional
    :return: The minimum and maximum velocity
    :rtype: Tuple[float,float]
    """
    directory = pathlib.Path(directory)
    with open(directory / "fault&rise_time.txt", "r") as infile:
        lines = [line.split() for line in infile]

    min_vel = float(lines[1][5])
    max_vel = float(lines[1][6])
    return min_vel, max_vel


if __name__ == "__main__":
    import argparse

    import wasp.manage_parser as mpar

    parser = argparse.ArgumentParser()
    parser.add_argument(
        "-f", "--folder", default=os.getcwd(), help="folder where there are input files"
    )
    parser = mpar.parser_add_tensor(parser)
    parser = mpar.parser_ffm_data(parser)
    parser.add_argument(
        "-o",
        "--option",
        choices=[
            "auto",
            "manual",
            "forward",
            "point_source",
            "checker_mod",
            "checker_noise",
            "point_source_err",
            "forward_patch",
            "add_data",
        ],
        required=True,
        help="which method to run",
    )
    parser.add_argument("-d", "--data", help="direction of folder with seismic data")
    parser.add_argument("-v", "--velmodel", help="direction of velocity model file")
    parser.add_argument(
        "-rst",
        "--st_response",
        action="store_false",
        help="whether to remove response of strong_motion or not",
    )
    args = parser.parse_args()
    if args.gcmt_tensor:
        args.gcmt_tensor = os.path.abspath(args.gcmt_tensor)
    if args.qcmt_tensor:
        args.qcmt_tensor = os.path.abspath(args.qcmt_tensor)
    if args.data:
        args.data = os.path.abspath(args.data)
    velmodel = args.velmodel if args.velmodel else None
    if velmodel:
        velmodel = mv.model2dict(velmodel)
    os.chdir(args.folder)
    data_type = mpar.get_used_data(args)
    data_type = data_type + ["insar"] if args.insar else data_type

    default_dirs = mng.default_dirs()
    if args.option not in ["auto"]:
        segments_data = json.load(open("segments_data.json"))
    if args.option == "auto":
        if not args.gcmt_tensor and not args.qcmt_tensor:
            raise RuntimeError("You must select direction of input GCMT file")
        if args.gcmt_tensor:
            tensor_info = tensor.get_tensor(
                cmt_file=args.gcmt_tensor, directory=pathlib.Path()
            )
        if args.qcmt_tensor:
            tensor_info = tensor.get_tensor(
                quake_file=args.qcmt_tensor, directory=pathlib.Path()
            )
        set_directory_structure(tensor_info, directory=pathlib.Path())
        if args.data:
            for file in os.listdir(args.data):
                if os.path.isfile(os.path.join(args.data, file)):
                    copy2(os.path.join(args.data, file), "data")
        data_type = data_type if len(data_type) >= 1 else ["tele_body"]
        automatic_usgs(
            tensor_info,
            data_type,
            default_dirs,
            velmodel=velmodel,
            dt_cgps=None,
            st_response=args.st_response,
            directory=pathlib.Path(),
        )
    if args.option == "add_data":
        if args.gcmt_tensor:
            cmt_file = args.gcmt_tensor
            tensor_info = tensor.get_tensor(cmt_file=cmt_file, directory=pathlib.Path())
        else:
            tensor_info = tensor.get_tensor(directory=pathlib.Path())
        if len(data_type) == 0:
            raise RuntimeError("You must input at least one data type")
        data_folder = args.data if args.data else None
        modelling_new_data(
            tensor_info,
            data_type,
            default_dirs,
            data_folder,
            segments_data,
            st_response=args.st_response,
            directory=pathlib.Path(),
        )
    if args.option == "manual":
        if args.gcmt_tensor:
            cmt_file = args.gcmt_tensor
            tensor_info = tensor.get_tensor(cmt_file=cmt_file, directory=pathlib.Path())
        else:
            tensor_info = tensor.get_tensor(directory=pathlib.Path())
        if len(data_type) == 0:
            raise RuntimeError("You must input at least one data type")
        data_folder = args.data if args.data else None
        manual_modelling(
            tensor_info,
            data_type,
            default_dirs,
            segments_data,
            directory=pathlib.Path(),
        )
    if args.option == "forward":
        if args.gcmt_tensor:
            cmt_file = args.gcmt_tensor
            tensor_info = tensor.get_tensor(cmt_file=cmt_file, directory=pathlib.Path())
        else:
            tensor_info = tensor.get_tensor(directory=pathlib.Path())
        if len(data_type) == 0:
            raise RuntimeError("You must input at least one data type")
        data_folder = args.data if args.data else None
        forward_modelling(
            tensor_info,
            data_type,
            default_dirs,
            segments_data,
            option="Solucion.txt",
            directory=pathlib.Path(),
        )
    if args.option == "forward_patch":
        if args.gcmt_tensor:
            cmt_file = args.gcmt_tensor
            tensor_info = tensor.get_tensor(cmt_file=cmt_file, directory=pathlib.Path())
        else:
            tensor_info = tensor.get_tensor(directory=pathlib.Path())
        if len(data_type) == 0:
            raise RuntimeError("You must input at least one data type")
        data_folder = args.data if args.data else None
        checkerboard(
            tensor_info,
            data_type,
            default_dirs,
            segments_data,
            max_slip=400,
            option="Patches",
            option2="forward",
            directory=pathlib.Path(),
        )
    if args.option == "checker_mod":
        if args.gcmt_tensor:
            cmt_file = args.gcmt_tensor
            tensor_info = tensor.get_tensor(cmt_file=cmt_file, directory=pathlib.Path())
        else:
            tensor_info = tensor.get_tensor(directory=pathlib.Path())
        if len(data_type) == 0:
            raise RuntimeError("You must input at least one data type")
        data_folder = args.data if args.data else None
        checkerboard(
            tensor_info,
            data_type,
            default_dirs,
            segments_data,
            add_error=True,
            directory=pathlib.Path(),
        )
    if args.option == "checker_noise":
        if args.gcmt_tensor:
            cmt_file = args.gcmt_tensor
            tensor_info = tensor.get_tensor(cmt_file=cmt_file, directory=pathlib.Path())
        else:
            tensor_info = tensor.get_tensor(directory=pathlib.Path())
        if len(data_type) == 0:
            raise RuntimeError("You must input at least one data type")
        data_folder = args.data if args.data else None
        checkerboard(
            tensor_info,
            data_type,
            default_dirs,
            segments_data,
            max_slip=0,
            add_error=True,
            directory=pathlib.Path(),
        )
    if args.option == "point_source":
        if args.gcmt_tensor:
            cmt_file = args.gcmt_tensor
            tensor_info = tensor.get_tensor(cmt_file=cmt_file, directory=pathlib.Path())
        else:
            tensor_info = tensor.get_tensor(directory=pathlib.Path())
        if len(data_type) == 0:
            raise RuntimeError("You must input at least one data type")
        data_folder = args.data if args.data else None
        checkerboard(
            tensor_info,
            data_type,
            default_dirs,
            segments_data,
            max_slip=500,
            option="Patches",
            directory=pathlib.Path(),
        )
    if args.option == "point_source_err":
        if args.gcmt_tensor:
            cmt_file = args.gcmt_tensor
            tensor_info = tensor.get_tensor(cmt_file=cmt_file, directory=pathlib.Path())
        else:
            tensor_info = tensor.get_tensor(directory=pathlib.Path())
        if len(data_type) == 0:
            raise RuntimeError("You must input at least one data type")
        data_folder = args.data if args.data else None
        checkerboard(
            tensor_info,
            data_type,
            default_dirs,
            segments_data,
            max_slip=300,
            option="Patches",
            add_error=True,
            directory=pathlib.Path(),
        )
