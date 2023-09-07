# -*- coding: utf-8 -*-
"""Script for performing FFM modelling and forward, using the method of
Chen-Ji.
"""


import errno
import glob
import json
import os
import subprocess
import time
import warnings
from multiprocessing import Process
from shutil import copy2, move

import numpy as np

import wasp.data_management as dm
import wasp.data_processing as proc
import wasp.fault_plane as pf
import wasp.green_functions as gf
import wasp.management as mng
import wasp.modelling_parameters as mp
import wasp.modulo_logs as ml
import wasp.plane_management as pl_mng
import wasp.plot_graphic as plot
import wasp.seismic_tensor as tensor
import wasp.traces_properties as tp
import wasp.velocity_models as mv
from wasp import get_outputs, input_files
from wasp.data_acquisition import acquisition
from wasp.load_ffm_model import load_ffm_model
from wasp.static2fsp import static_to_fsp


def automatic_usgs(
    tensor_info, data_type, default_dirs, velmodel=None, dt_cgps=1.0, st_response=True
):
    """Routine for automatic FFM modelling

    :param tensor_info: dictionary with moment tensor properties
    :param data_type: list with data types to be used in modelling
    :param default_dirs: dictionary with default directories to be used
    :param velmodel: dictionary with velocity model
    :param dt_cgps: sampling interval for cgps data
    :param st_response: whether to remove paz response of strong motion
    :type tensor_info: dict
    :type data_type: list
    :type default_dirs: dict
    :type velmodel: dict, optional
    :type dt_cgps: float, optional
    :type st_response: bool, optional
    """
    logger = ml.create_log("automatic_ffm", os.path.join("logs", "automatic_ffm.log"))
    logger = ml.add_console_handler(logger)
    logger.info("Starting fff program")
    sol_folder = os.getcwd()
    sol_folder = os.path.abspath(sol_folder)
    time0 = time.time()
    if "gps" in data_type:
        if os.path.isfile(os.path.join("data", "gps_data")):
            copy2(os.path.join("data", "gps_data"), sol_folder)
    if "insar" in data_type:
        insar_files = glob.glob(os.path.join("data", "insar_a*txt"))
        insar_files = insar_files + glob.glob(os.path.join("data", "insar_d*txt"))
        for file in insar_files:
            if os.path.isfile(file):
                copy2(file, sol_folder)
    data_prop = tp.properties_json(tensor_info, dt_cgps=dt_cgps)
    os.chdir(os.path.join(sol_folder, "data"))
    time2 = time.time()
    logger.info("Process data")
    processing(tensor_info, data_type, data_prop, st_response=st_response)
    time2 = time.time() - time2
    logger.info("Time spent processing traces: {}".format(time2))
    os.chdir(sol_folder)
    data_folder = os.path.join(sol_folder, "data")
    insar_asc = glob.glob("insar_asc*txt")
    insar_asc = None if len(insar_asc) == 0 else insar_asc
    insar_desc = glob.glob("insar_desc*txt")
    insar_desc = None if len(insar_desc) == 0 else insar_desc
    dm.filling_data_dicts(
        tensor_info,
        data_type,
        data_prop,
        data_folder,
        insar_asc=insar_asc,
        insar_desc=insar_desc,
    )
    writing_inputs0(tensor_info, data_type)
    logger.info("Compute GF bank")
    if not velmodel:
        velmodel = mv.select_velmodel(tensor_info, default_dirs)
    input_files.write_velmodel(velmodel)
    gf_bank_str = os.path.join(sol_folder, "GF_strong")
    gf_bank_cgps = os.path.join(sol_folder, "GF_cgps")
    get_gf_bank = default_dirs["strong_motion_gf_bank2"]
    if "cgps" in data_type:
        logger.info("Compute cGPS GF bank")
        green_dict = gf.fk_green_fun1(data_prop, tensor_info, gf_bank_cgps, cgps=True)
        input_files.write_green_file(green_dict, cgps=True)
        with open(os.path.join("logs", "GF_cgps_log"), "w") as out_gf_cgps:
            p1 = subprocess.Popen([get_gf_bank, "cgps"], stdout=out_gf_cgps)
        p1.wait()
    if "strong_motion" in data_type:
        logger.info("Compute strong motion GF bank")
        green_dict = gf.fk_green_fun1(data_prop, tensor_info, gf_bank_str)
        input_files.write_green_file(green_dict)
        with open(os.path.join("logs", "GF_strong_log"), "w") as out_gf_strong:
            p2 = subprocess.Popen(
                [
                    get_gf_bank,
                ],
                stdout=out_gf_strong,
            )
        p2.wait()

    files = [
        "Green_strong.txt",
        "Green_cgps.txt",
        "modelling_stats.json",
        "gps_data",
        "strong_motion_gf.json",
        "cgps_gf.json",
        "sampling_filter.json",
    ]
    files2 = glob.glob("channels_*txt")
    files3 = glob.glob("wavelets_*txt")
    files4 = glob.glob("waveforms_*txt")
    files5 = glob.glob("*waves.json")
    files6 = glob.glob("static*")
    files7 = glob.glob("filtro*") + glob.glob("surf_filter*")
    files8 = ["instrumental_response.txt", "body_wave_weight.txt"]
    files9 = glob.glob("insar*")
    files = (
        files + files2 + files3 + files4 + files5 + files6 + files7 + files8 + files9
    )
    folders = ["NP1", "NP2"]
    for folder in folders:
        for file in files:
            if os.path.isfile(file):
                copy2(file, folder)
    info_np1, info_np2 = tensor.planes_from_tensor(tensor_info)
    keywords = {"velmodel": velmodel}
    os.chdir(os.path.join(sol_folder, "NP1"))
    p1 = Process(
        target=_automatic2,
        args=(tensor_info, info_np1, data_type, data_prop, default_dirs, logger),
        kwargs=keywords,
    )
    p1.start()
    os.chdir(os.path.join(sol_folder, "NP2"))
    p2 = Process(
        target=_automatic2,
        args=(tensor_info, info_np2, data_type, data_prop, default_dirs, logger),
        kwargs=keywords,
    )
    p2.start()
    [p.join() for p in [p1, p2]]
    logger.info("Time spent: {}".format(time.time() - time0))
    ml.close_log(logger)
    return


def _automatic2(
    tensor_info,
    plane_data,
    data_type,
    data_prop,
    default_dirs,
    logger,
    velmodel=None,
    check_surf=True,
):
    """Routine for automatic FFM modelling for each nodal plane

    :param tensor_info: dictionary with moment tensor properties
    :param plane_data: dictionary with fault plane mechanism
    :param data_type: list with data types to be used in modelling
    :param default_dirs: dictionary with default directories to be used
    :param data_prop: dictionary with properties for different waveform types
    :param logger: logging object
    :param velmodel: dictionary with velocity model
    :param check_surf: check whether surface waves can be used
    :type default_dirs: dict
    :type tensor_info: dict
    :type plane_data: dict
    :type data_type: list
    :type data_prop: dict
    :type logger: Logger
    :type velmodel: dict, optional
    :type check_surf: bool, optional
    """
    #
    # Create JSON files
    #
    logger.info("Create input files for Fortran scripts")
    logger.info("Create automatic JSON")
    tensor.write_tensor(tensor_info)
    if velmodel:
        mv.velmodel2json(velmodel)
    if not velmodel:
        velmodel = mv.select_velmodel(tensor_info, default_dirs)
    np_plane_info = plane_data["plane_info"]
    data_folder = os.path.join("..", "data")
    insar_asc = glob.glob("insar_asc*txt")
    insar_asc = None if len(insar_asc) == 0 else insar_asc
    insar_desc = glob.glob("insar_desc*txt")
    insar_desc = None if len(insar_desc) == 0 else insar_desc
    dm.filling_data_dicts(
        tensor_info,
        data_type,
        data_prop,
        data_folder,
        insar_asc=insar_asc,
        insar_desc=insar_desc,
    )
    segments_data = pf.create_finite_fault(tensor_info, np_plane_info, data_type)
    segments = segments_data["segments"]
    rise_time = segments_data["rise_time"]
    connections = None
    if "connections" in segments_data:
        connections = segments_data["connections"]
    point_sources = pf.point_sources_param(
        segments, tensor_info, rise_time, connections=connections
    )
    data_type = _check_surf_GF(point_sources, data_type, logger=logger)
    mp.modelling_prop(tensor_info, segments_data, data_type=data_type)
    #
    # write text files from JSONs
    #
    rupt_vel = segments[0]["rupture_vel"]
    lambda_min = 0.4
    lambda_max = 1.25
    min_vel, max_vel = [lambda_min * rupt_vel, lambda_max * rupt_vel]
    logger.info("Write input files")
    writing_inputs(tensor_info, data_type, segments_data, min_vel, max_vel)
    #
    # Modelling and plotting results
    #
    inversion(tensor_info, data_type, default_dirs, logger)
    logger.info("Plot data in folder {}".format(os.getcwd()))
    execute_plot(tensor_info, data_type, segments_data, default_dirs, velmodel=velmodel)
    base = os.path.basename(os.getcwd())
    dirname = os.path.abspath(os.getcwd())
    #
    # write solution in FSP format
    #
    solution = get_outputs.read_solution_static_format(segments)
    static_to_fsp(tensor_info, segments_data, data_type, velmodel, solution)
    for file in glob.glob("*png"):
        if os.path.isfile(os.path.join(dirname, base, file)):
            copy2(os.path.join(dirname, base, file), os.path.join(dirname, "plots"))


def _check_surf_GF(point_sources, used_data, logger=None):
    """ """
    new_used_data = used_data.copy()
    depths = [ps[:, :, :, :, 2] for ps in point_sources]
    depths = [np.max(depths1) for depths1 in depths]
    max_depth = np.max(depths)
    is_surf = "surf_tele" in used_data
    if max_depth > 125 and is_surf:
        warnings.warn(
            "Maximum depth larger than 125 km. " "Surface waves won't be used"
        )
        new_used_data.remove("surf_tele")
        if logger:
            logger.info("Maximum depth larger than 125 km.")
    return new_used_data


def modelling_new_data(
    tensor_info, data_type, default_dirs, data_folder, segments_data, st_response=True
):
    """Routine for manual finite fault modelling with new data types.

    :param tensor_info: dictionary with moment tensor properties
    :param data_type: list with data types to be used in modelling
    :param default_dirs: dictionary with default directories to be used
    :param data_folder: location of data used for modelling
    :param segments_data: properties of fault segments and rise time
    :param st_response: whether to remove paz response of strong motion
    :type default_dirs: dict
    :type tensor_info: dict
    :type data_type: list
    :type data_folder: string
    :type segments_data: dict
    :type st_response: bool, optional
    """
    sol_folder = os.getcwd()
    sol_folder = os.path.abspath(sol_folder)
    if os.path.isfile(os.path.join(data_folder, "gps_data")):
        copy2(os.path.join(data_folder, "gps_data"), sol_folder)
    insar_asc = glob.glob(os.path.join(data_folder, "insar_a*txt"))
    insar_desc = glob.glob(os.path.join(data_folder, "insar_d*txt"))
    insar_files = insar_asc + insar_desc
    for file in insar_files:
        if os.path.isfile(file):
            copy2(file, sol_folder)
    data_prop = json.load(open("sampling_filter.json"))
    os.chdir(os.path.join(data_folder))
    time2 = time.time()
    processing(tensor_info, data_type, data_prop, st_response=st_response)
    os.chdir(sol_folder)
    dm.filling_data_dicts(
        tensor_info,
        data_type,
        data_prop,
        data_folder,
        insar_asc=insar_asc,
        insar_desc=insar_desc,
    )
    gf_bank_str = os.path.join(sol_folder, "GF_strong")
    gf_bank_cgps = os.path.join(sol_folder, "GF_cgps")
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
            data_prop, tensor_info, gf_bank_cgps, max_depth=max_depth, cgps=True
        )
        input_files.write_green_file(green_dict, cgps=True)
        with open(os.path.join("logs", "GF_cgps_log"), "w") as out_gf_cgps:
            p1 = subprocess.Popen([get_gf_bank, "cgps"], stdout=out_gf_cgps)
        p1.wait()
    if "strong_motion" in data_type:
        green_dict = gf.fk_green_fun1(
            data_prop, tensor_info, gf_bank_str, max_depth=max_depth
        )
        input_files.write_green_file(green_dict)
        with open(os.path.join("logs", "GF_strong_log"), "w") as out_gf_strong:
            p2 = subprocess.Popen(
                [
                    get_gf_bank,
                ],
                stdout=out_gf_strong,
            )
        p2.wait()
    data_type2 = []
    if os.path.isfile("tele_waves.json"):
        data_type2 = data_type2 + ["tele_body"]
    if os.path.isfile("surf_waves.json"):
        data_type2 = data_type2 + ["surf_tele"]
    if os.path.isfile("strong_motion_waves.json"):
        data_type2 = data_type2 + ["strong_motion"]
    if os.path.isfile("cgps_waves.json"):
        data_type2 = data_type2 + ["cgps"]
    if os.path.isfile("static_data.json"):
        data_type2 = data_type2 + ["gps"]
    if os.path.isfile("insar_data.json"):
        data_type2 = data_type2 + ["insar"]
    if os.path.isfile("dart_waves.json"):
        data_type2 = data_type2 + ["dart"]
    manual_modelling(tensor_info, data_type2, default_dirs, segments_data)
    return


def manual_modelling(tensor_info, data_type, default_dirs, segments_data):
    """Routine for manual finite fault modelling.

    :param tensor_info: dictionary with moment tensor properties
    :param data_type: list with data types to be used in modelling
    :param default_dirs: dictionary with default directories to be used
    :param segments_data: properties of fault segments and rise time
    :type default_dirs: dict
    :type tensor_info: dict
    :type data_type: list
    :type segments_data: dict
    """
    if not os.path.isdir("logs"):
        os.mkdir("logs")
    if not os.path.isdir("plots"):
        os.mkdir("plots")
    # segments_data = json.load(open('segments_data.json'))
    min_vel, max_vel = __ask_velrange()
    logger = ml.create_log("manual_ffm", os.path.join("logs", "manual_ffm.log"))
    logger.info("Write input files")
    tensor.write_tensor(tensor_info)
    writing_inputs(tensor_info, data_type, segments_data, min_vel, max_vel)
    writing_inputs0(tensor_info, data_type)
    inversion(tensor_info, data_type, default_dirs, logger)
    logger.info("Plot data in folder {}".format(os.getcwd()))
    execute_plot(tensor_info, data_type, segments_data, default_dirs)
    ml.close_log(logger)


def forward_modelling(
    tensor_info,
    data_type,
    default_dirs,
    segments_data,
    option="Solucion.txt",
    max_slip=200,
):
    """Routine for forward modelling.

    :param tensor_info: dictionary with moment tensor properties
    :param data_type: set with data types to be used in modelling
    :param option: string with location of input file with kinematic model to
     use
    :param max_slip: maximum slip in case of checkerboard test
    :param default_dirs: dictionary with default directories to be used
    :param segments_data: properties of fault segments and rise time
    :type default_dirs: dict
    :type tensor_info: dict
    :type data_type: set, optional
    :type option: string, optional
    :type max_slip: float, optional
    :type segments_data: dict
    """
    tensor.write_tensor(tensor_info)
    if not os.path.isdir("logs"):
        os.mkdir("logs")
    if not os.path.isdir("plots"):
        os.mkdir("plots")
    len_stk = 5 if not option == "point_source" else 8
    len_dip = 5 if not option == "point_source" else 1
    # segments_data = json.load(open('segments_data.json'))
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
    )
    if not os.path.isfile("velmodel_data.json"):
        raise FileNotFoundError(
            errno.ENOENT, os.strerror(errno.ENOENT), "velmodel_data.json"
        )
    velmodel = json.load(open("velmodel_data.json"))
    min_vel, max_vel = __ask_velrange()

    logger = ml.create_log("forward_model", os.path.join("logs", "forward_model.log"))
    logger.info("Write input files")
    # segments, rise_time, point_sources = pl_mng.__read_planes_info()
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
    )
    inversion(tensor_info, data_type, default_dirs, logger, forward=True)
    logger.info("Plot data in folder {}".format(os.getcwd()))
    execute_plot(tensor_info, data_type, segments_data, default_dirs, velmodel=velmodel)
    ml.close_log(logger)


def checkerboard(
    tensor_info,
    data_type,
    default_dirs,
    segments_data,
    max_slip=200,
    add_error=False,
    option="Checkerboard",
    option2="FFM modelling",
):
    """Routine for running checkerboard tests.

    :param tensor_info: dictionary with moment tensor properties
    :param data_type: set with data types to be used in modelling
    :param option: string with location of input file with kinematic model to
     use
    :param max_slip: maximum slip in case of checkerboard test
    :param add_error: whether we add noise to synthetic waveforms
    :param option2: whether we invert the checkerboard model or not
    :param default_dirs: dictionary with default directories to be used
    :param segments_data: properties of fault segments and rise time
    :type default_dirs: dict
    :type tensor_info: dict
    :type data_type: set
    :type option: string, optional
    :type max_slip: float, optional
    :type add_error: bool, optional
    :type option2: string, optional
    :type segments_data: dict
    """
    if max_slip > 0:
        folder_name = "checkerboard_resolution"
    else:
        folder_name = "checkerboard_noise"
    if not option == "Checkerboard":
        folder_name = option
    if not os.path.isdir(folder_name):
        os.mkdir(folder_name)
    if not os.path.isdir("logs"):
        os.mkdir("logs")
    if not os.path.isdir("plots"):
        os.mkdir("plots")
    for file in os.listdir():
        if os.path.isfile(file):
            copy2(file, folder_name)
    os.chdir(folder_name)
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
    )
    data_prop = json.load(open("sampling_filter.json"))
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
        files = json.load(open(json_dict))
        input_files.from_synthetic_to_obs(
            files, data_type0, tensor_info, data_prop, add_error=add_error
        )
    logger = ml.create_log(
        "checkerboard_ffm", os.path.join("logs", "checkerboard_ffm.log")
    )
    if option2 == "FFM modelling":
        inversion(tensor_info, data_type, default_dirs, logger)
        execute_plot(
            tensor_info, data_type, segments_data, default_dirs, plot_input=True
        )
    ml.close_log(logger)


def set_directory_structure(tensor_info):
    """Create directory structure

    :param tensor_info: dictionary with moment tensor properties
    :type tensor_info: dict
    """
    sol_folder = mng.start_time_id(tensor_info)
    if not os.path.isdir(sol_folder):
        os.mkdir(sol_folder)
    sol_folder = os.path.abspath(sol_folder)
    version = len(glob.glob(os.path.join(sol_folder, "ffm*")))
    sol_folder2 = os.path.join(sol_folder, "ffm.{}".format(version))
    os.mkdir(sol_folder2)
    os.mkdir(os.path.join(sol_folder2, "data"))
    os.mkdir(os.path.join(sol_folder2, "data", "cGPS"))
    os.mkdir(os.path.join(sol_folder2, "data", "STR"))
    os.mkdir(os.path.join(sol_folder2, "data", "P"))
    os.mkdir(os.path.join(sol_folder2, "data", "SH"))
    os.mkdir(os.path.join(sol_folder2, "data", "LONG"))
    os.mkdir(os.path.join(sol_folder2, "data", "final"))
    os.mkdir(os.path.join(sol_folder2, "data", "final_r"))
    os.mkdir(os.path.join(sol_folder2, "NP1"))
    os.mkdir(os.path.join(sol_folder2, "NP1", "logs"))
    os.mkdir(os.path.join(sol_folder2, "NP1", "plots"))
    os.mkdir(os.path.join(sol_folder2, "NP2"))
    os.mkdir(os.path.join(sol_folder2, "NP2", "logs"))
    os.mkdir(os.path.join(sol_folder2, "NP2", "plots"))
    os.mkdir(os.path.join(sol_folder2, "logs"))
    os.mkdir(os.path.join(sol_folder2, "plots"))
    os.mkdir(os.path.join(sol_folder2, "plots", "NP1"))
    os.mkdir(os.path.join(sol_folder2, "plots", "NP2"))
    os.chdir(sol_folder2)
    return


def processing(tensor_info, data_type, data_prop, st_response=True):
    """Run all waveform data processing.

    :param tensor_info: dictionary with moment tensor properties
    :param data_type: set with data types to be used in modelling
    :param data_prop: dictionary with properties for different waveform types
    :param st_response: whether to remove paz response of strong motion
    :type tensor_info: dict
    :type data_type: set
    :type data_prop: dict
    :type st_response: bool, optional
    """
    tele_files = (
        glob.glob("*.BH*SAC")
        + glob.glob("*.BH*sac")
        + glob.glob("*_BH*sac")
        + glob.glob("*_BH*sac")
    )
    strong_files = (
        glob.glob("*.HN*SAC")
        + glob.glob("*.HL*SAC")
        + glob.glob("*.HN*sac")
        + glob.glob("*.HL*sac")
        + glob.glob("*.AH?.*")
        + glob.glob("*_HN*sac")
        + glob.glob("*_HL*sac")
        + glob.glob("*HG*sac")
    )
    cgps_files = glob.glob("*L[HXY]*sac") + glob.glob("*L[HXY]*SAC")
    if "tele_body" in data_type:
        proc.select_process_tele_body(tele_files, tensor_info, data_prop)
    if "surf_tele" in data_type:
        proc.select_process_surf_tele(tele_files, tensor_info, data_prop)
    if "strong_motion" in data_type:
        proc.select_process_strong(
            strong_files, tensor_info, data_prop, remove_response=st_response
        )
    if "cgps" in data_type:
        proc.select_process_cgps(cgps_files, tensor_info, data_prop)


def writing_inputs0(tensor_info, data_type):
    """ """
    if not os.path.isfile("sampling_filter.json"):
        raise FileNotFoundError(
            errno.ENOENT, os.strerror(errno.ENOENT), "sampling_filter.json"
        )
    data_prop = json.load(open("sampling_filter.json"))
    if "tele_body" in data_type:
        input_files.input_chen_tele_body(tensor_info, data_prop)
    if "surf_tele" in data_type:
        input_files.input_chen_tele_surf(tensor_info, data_prop)
    if "strong_motion" in data_type:
        input_files.input_chen_near_field(tensor_info, data_prop, "strong_motion")
    if "cgps" in data_type:
        input_files.input_chen_near_field(tensor_info, data_prop, "cgps")
    if "gps" in data_type:
        input_files.input_chen_static(tensor_info)
    if "insar" in data_type:
        input_files.input_chen_insar()
    if "dart" in data_type:
        input_files.input_chen_dart(tensor_info, data_prop)


def writing_inputs(
    tensor_info,
    data_type,
    segments_data,
    min_vel,
    max_vel,
    moment_mag=None,
    forward_model=None,
):
    """Write all required text files from the information found in the JSONs.

    :param tensor_info: dictionary with moment tensor properties
    :param data_type: set with data types to be used in modelling
    :param segments_data: properties of fault segments and rise time
    :param min_vel: minimum rupture velocity
    :param max_vel: maximum rupture velocity
    :param moment_mag: input seismic moment
    :param forward_model: input kinematic model
    :type tensor_info: dict
    :type data_type: set
    :type min_vel: float
    :type max_vel: float
    :type segments_data: dict
    :type moment_mag: float, optional
    :type forward_model: dict, optional
    """
    if not os.path.isfile("velmodel_data.json"):
        raise FileNotFoundError(
            errno.ENOENT, os.strerror(errno.ENOENT), "velmodel_data.json"
        )
    velmodel = json.load(open("velmodel_data.json"))
    input_files.write_velmodel(velmodel)
    input_files.plane_for_chen(tensor_info, segments_data, min_vel, max_vel, velmodel)
    if forward_model:
        input_files.forward_model(
            tensor_info, segments_data, forward_model, min_vel, max_vel
        )
    if not os.path.isfile("annealing_prop.json"):
        raise FileNotFoundError(
            errno.ENOENT, os.strerror(errno.ENOENT), "annealing_prop.json"
        )
    dictionary = json.load(open("annealing_prop.json"))
    if moment_mag:
        dictionary["seismic_moment"] = moment_mag
    input_files.inputs_simmulated_annealing(dictionary, data_type)
    if "cgps" in data_type:
        if not os.path.isfile("cgps_gf.json"):
            raise FileNotFoundError(
                errno.ENOENT, os.strerror(errno.ENOENT), "cgps_gf.json"
            )
        green_dict = json.load(open("cgps_gf.json"))
        input_files.write_green_file(green_dict, cgps=True)
    if "strong_motion" in data_type:
        if not os.path.isfile("strong_motion_gf.json"):
            raise FileNotFoundError(
                errno.ENOENT, os.strerror(errno.ENOENT), "strong_motion_gf.json"
            )
        green_dict = json.load(open("strong_motion_gf.json"))
        input_files.write_green_file(green_dict)
    if not os.path.isfile("model_space.json"):
        raise FileNotFoundError(
            errno.ENOENT, os.strerror(errno.ENOENT), "model_space.json"
        )
    segments2 = json.load(open("model_space.json"))
    input_files.model_space(segments2)
    return


def inversion(tensor_info, data_type, default_dirs, logger, forward=False):
    """We get the binaries with gf for each station, run the ffm code, and
    proceed to plot the results.

    :param tensor_info: dictionary with moment tensor properties
    :param data_type: set with data types to be used in modelling
    :param logger: name of logfile to write
    :param forward: whether we solve the inverse problem or the forward problem
    :param default_dirs: dictionary with default directories to be used
    :type default_dirs: dict
    :type tensor_info: dict
    :type data_type: set
    :type logger: Logger
    :type forward: bool, optional
    """
    # print('Retrieve GF for modelling')
    logger.info("Green_functions")
    time1 = time.time()
    gf.gf_retrieve(data_type, default_dirs)
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
        # print('Perform kinematic modelling')
        logger.info("Inversion at folder {}".format(os.getcwd()))
        finite_fault = default_dirs["finite_fault"]
    else:
        logger.info("Forward at folder {}".format(os.getcwd()))
        finite_fault = default_dirs["forward"]
    p1 = subprocess.Popen([finite_fault, *args], stderr=subprocess.PIPE)
    #
    # need to wait until FFM modelling is finished.
    #
    outs, errs = p1.communicate(timeout=40 * 60)
    if errs:
        logger.error(errs.decode("utf-8"))
    # p1.wait()
    time3 = time.time() - time3
    logger.info("Elapsed time of finite fault modelling: {}".format(time3))
    run_stats["ffm_time"] = time3
    delete_binaries()
    return run_stats


def execute_plot(
    tensor_info, data_type, segments_data, default_dirs, velmodel=None, plot_input=False
):
    """We plot modelling results

    :param tensor_info: dictionary with moment tensor properties
    :param data_type: set with data types to be used in modelling
    :param plot_input: choose whether to plot initial kinematic model as well
    :param default_dirs: dictionary with default directories to be used
    :param velmodel: dictionary with velocity model
    :param segments_data: properties of fault segments and rise time
    :type velmodel: dict, optional
    :type default_dirs: dict
    :type tensor_info: dict
    :type data_type: set
    :type segments:data: dict
    :type plot_input: bool, optional
    """

    print("Plot results")
    segments = segments_data["segments"]
    rise_time = segments_data["rise_time"]
    connections = None
    if "connections" in segments_data:
        connections = segments_data["connections"]
    solution = get_outputs.read_solution_static_format(segments)
    if not velmodel:
        velmodel = mv.select_velmodel(tensor_info, default_dirs)
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
    )
    plot.plot_misfit(data_type)
    # plot.plot_beachballs(tensor_info, data_type)
    traces_info, stations_gps = [None, None]
    if "strong_motion" in data_type:
        traces_info = json.load(open("strong_motion_waves.json"))
    if "gps" in data_type:
        names, lats, lons, observed, synthetic, error = get_outputs.retrieve_gps()
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
        )
    if "insar" in data_type:
        insar_data = get_outputs.get_insar()
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
                )
    if plot_input:
        input_model = load_ffm_model(
            segments_data, point_sources, option="fault&rise_time.txt"
        )
        plot._PlotSlipDist_Compare(segments, point_sources, input_model, solution)
        plot._PlotComparisonMap(
            tensor_info, segments, point_sources, input_model, solution
        )
    plot_files = glob.glob(os.path.join("plots", "*png"))
    for plot_file in plot_files:
        os.remove(plot_file)
    plot_files = glob.glob("*png")
    for plot_file in plot_files:
        move(plot_file, "plots")


def delete_binaries():
    """to remove the files with Green function data."""
    deletables = glob.glob("*.GRE") + glob.glob("*.TDE") + glob.glob("*[1-2-3]")
    deletables = deletables + glob.glob("*.H[LN][E-N-Z]")
    deletables = deletables + glob.glob("*.L[HXY][E-N-Z]")
    for file in deletables:
        if os.path.isfile(file):
            os.remove(file)


def __ask_velrange():
    with open("fault&rise_time.txt", "r") as infile:
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
            tensor_info = tensor.get_tensor(cmt_file=args.gcmt_tensor)
        if args.qcmt_tensor:
            tensor_info = tensor.get_tensor(quake_file=args.qcmt_tensor)
        # if not args.gcmt_tensor:
        #     raise RuntimeError('You must select direction of input GCMT file')
        # tensor_info = tensor.get_tensor(cmt_file=args.gcmt_tensor)
        set_directory_structure(tensor_info)
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
        )
    if args.option == "add_data":
        if args.gcmt_tensor:
            cmt_file = args.gcmt_tensor
            tensor_info = tensor.get_tensor(cmt_file=cmt_file)
        else:
            tensor_info = tensor.get_tensor()
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
        )
    if args.option == "manual":
        if args.gcmt_tensor:
            cmt_file = args.gcmt_tensor
            tensor_info = tensor.get_tensor(cmt_file=cmt_file)
        else:
            tensor_info = tensor.get_tensor()
        if len(data_type) == 0:
            raise RuntimeError("You must input at least one data type")
        data_folder = args.data if args.data else None
        manual_modelling(tensor_info, data_type, default_dirs, segments_data)
    if args.option == "forward":
        if args.gcmt_tensor:
            cmt_file = args.gcmt_tensor
            tensor_info = tensor.get_tensor(cmt_file=cmt_file)
        else:
            tensor_info = tensor.get_tensor()
        if len(data_type) == 0:
            raise RuntimeError("You must input at least one data type")
        data_folder = args.data if args.data else None
        forward_modelling(
            tensor_info, data_type, default_dirs, segments_data, option="Solucion.txt"
        )
    if args.option == "forward_patch":
        if args.gcmt_tensor:
            cmt_file = args.gcmt_tensor
            tensor_info = tensor.get_tensor(cmt_file=cmt_file)
        else:
            tensor_info = tensor.get_tensor()
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
        )
    if args.option == "checker_mod":
        if args.gcmt_tensor:
            cmt_file = args.gcmt_tensor
            tensor_info = tensor.get_tensor(cmt_file=cmt_file)
        else:
            tensor_info = tensor.get_tensor()
        if len(data_type) == 0:
            raise RuntimeError("You must input at least one data type")
        data_folder = args.data if args.data else None
        checkerboard(
            tensor_info, data_type, default_dirs, segments_data, add_error=True
        )
    if args.option == "checker_noise":
        if args.gcmt_tensor:
            cmt_file = args.gcmt_tensor
            tensor_info = tensor.get_tensor(cmt_file=cmt_file)
        else:
            tensor_info = tensor.get_tensor()
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
        )
    if args.option == "point_source":
        if args.gcmt_tensor:
            cmt_file = args.gcmt_tensor
            tensor_info = tensor.get_tensor(cmt_file=cmt_file)
        else:
            tensor_info = tensor.get_tensor()
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
        )
    if args.option == "point_source_err":
        if args.gcmt_tensor:
            cmt_file = args.gcmt_tensor
            tensor_info = tensor.get_tensor(cmt_file=cmt_file)
        else:
            tensor_info = tensor.get_tensor()
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
        )
