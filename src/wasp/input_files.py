# -*- coding: utf-8 -*-
"""Module with routines for writing text files which are necessary inputs for
Chen's fortran scripts.
"""


import errno
import io
import json
import os
import pathlib
from typing import Any, List, Literal, Optional, Union

import numpy as np
from obspy import UTCDateTime, read  # type: ignore
from obspy.core.stream import Stream  # type: ignore
from obspy.geodetics import kilometers2degrees  # type: ignore
from obspy.io.sac.util import SacIOError  # type: ignore
from obspy.taup import TauPyModel  # type: ignore
from scipy.signal import butter, filtfilt  # type: ignore

import wasp.fault_plane as pf
import wasp.get_outputs as get_outputs
import wasp.management as mng
import wasp.seismic_tensor as tensor

################################
# velocity models
################################


def write_velmodel(
    velmodel: dict, directory: Union[pathlib.Path, str] = pathlib.Path()
):
    """Write velocity model file for fortran scripts

    :param velmodel: The velocity model properties (from velmodel_data.json)
    :type velmodel: dict
    :param directory: The directory to write to, defaults to pathlib.Path()
    :type directory: Union[pathlib.Path, str], optional
    """
    directory = pathlib.Path(directory)
    p_vel = velmodel["p_vel"]
    s_vel = velmodel["s_vel"]
    dens = velmodel["dens"]
    thick = velmodel["thick"]
    qa = velmodel["qa"]
    qb = velmodel["qb"]
    zipped = zip(p_vel, s_vel, dens, thick, qa, qb)
    with open(directory / "vel_model.txt", "w") as outfile:
        outfile.write("{}\n".format(len(thick)))
        for pv, sv, rho, th, qaa, qbb in zipped:
            outfile.write("{} {} {} {} {} {}\n".format(pv, sv, rho, th, qaa, qbb))
    return


################################
# fault plane
################################


def forward_model(
    tensor_info: dict,
    segments_data: dict,
    model: dict,
    vel0: float,
    vel1: float,
    directory: Union[pathlib.Path, str] = pathlib.Path(),
):
    """Rewrite the Fault.time file with an input model

    :param tensor_info: The tensor information (from tensor_info.json)
    :type tensor_info: dict
    :param segments_data: The segments information (from segments_data.json)
    :type segments_data: dict
    :param model: The model
    :type model: np.ndarray
    :param vel0: The minimum velocity
    :type vel0: float
    :param vel1: The maximum velocity
    :type vel1: float
    :param directory: The directory to read/write to, defaults to pathlib.Path()
    :type directory: Union[pathlib.Path, str], optional
    """
    directory = pathlib.Path(directory)
    slip_segs = model["slip"]
    rake_segs = model["rake"]
    trup_segs = model["trup"]
    tris_segs = model["trise"]
    tfall_segs = model["tfall"]

    segments = segments_data["segments"]
    rise_time = segments_data["rise_time"]
    connections = None
    if "connections" in segments_data:
        connections = segments_data["connections"]

    hyp_stk = segments[0]["hyp_stk"]
    hyp_dip = segments[0]["hyp_dip"]
    delta_strike = segments[0]["delta_strike"]
    delta_dip = segments[0]["delta_dip"]
    rupture_vel = segments[0]["rupture_vel"]
    subfaults = {"delta_strike": delta_strike, "delta_dip": delta_dip}

    subfaults2 = pf._point_sources_def(rise_time, rupture_vel, subfaults)
    strike_ps = subfaults2["strike_ps"]
    dip_ps = subfaults2["dip_ps"]
    t1 = rise_time["min_rise"]
    t2 = rise_time["delta_rise"]
    windows = rise_time["windows"]

    depth = tensor_info["depth"]

    disp_or_vel = 0
    string = "{} {} {} {} {}\n"

    point_sources0 = pf.point_sources_param(
        segments, tensor_info, rise_time, connections=connections
    )
    ny = int(dip_ps / 2)
    nx = int(strike_ps / 2)
    times = [point_sources[:, :, ny, nx, 4] for point_sources in point_sources0]
    trup_segs2 = [rupt_seg - time for time, rupt_seg in zip(times, trup_segs)]

    zipped = zip(segments, slip_segs, rake_segs, trup_segs2, tris_segs, tfall_segs)
    with open(directory / "fault&rise_time.txt", "w") as outfile:
        outfile.write("{} {} {} 10\n".format(hyp_stk, hyp_dip, depth))
        outfile.write(
            "{} {} {} {} {} {} {} {} {}\n".format(
                len(segments),
                delta_strike,
                delta_dip,
                strike_ps,
                dip_ps,
                vel0,
                vel1,
                -100,
                100,
            )
        )
        outfile.write(
            "{} {} {} {} {}\n".format(t1, t2, windows, rupture_vel, disp_or_vel)
        )
        for i_segment, (
            segment,
            slip_seg,
            rake_seg,
            trup_seg,
            tris_seg,
            tfall_seg,
        ) in enumerate(zipped):
            dip = segment["dip"]
            strike = segment["strike"]
            n_stk = segment["stk_subfaults"]
            n_dip = segment["dip_subfaults"]
            outfile.write("{} {} {}\n".format(i_segment + 1, dip, strike))
            outfile.write("{} {} 0\n".format(n_stk, n_dip))
            for i in range(n_dip):
                for j in range(n_stk):
                    outfile.write(
                        string.format(
                            slip_seg[i, j],
                            rake_seg[i, j],
                            trup_seg[i, j],
                            tris_seg[i, j],
                            tfall_seg[i, j],
                        )
                    )
    return


def plane_for_chen(
    tensor_info: dict,
    segments_data: dict,
    min_vel: float,
    max_vel: float,
    velmodel: dict,
    directory: Union[pathlib.Path, str] = pathlib.Path(),
):
    """Write Fault.time, Fault.pos, and Niu_model files

    :param tensor_info: The tensor information (from tensor_info.json)
    :type tensor_info: dict
    :param segments_data: The segments information (from segments_data.json)
    :type segments_data: dict
    :param min_vel: Minimum velocity
    :type min_vel: float
    :param max_vel: Maximum velocity
    :type max_vel: float
    :param velmodel: The velocity model (from velmodel_data.json)
    :type velmodel: dict
    :param directory: The directory to read/write to, defaults to pathlib.Path()
    :type directory: Union[pathlib.Path, str], optional
    """
    directory = pathlib.Path(directory)
    segments = segments_data["segments"]
    rise_time = segments_data["rise_time"]
    connections = None
    if "connections" in segments_data:
        connections = segments_data["connections"]
    delta_strike = segments[0]["delta_strike"]
    delta_dip = segments[0]["delta_dip"]
    rupture_vel = segments[0]["rupture_vel"]
    subfaults = {"delta_strike": delta_strike, "delta_dip": delta_dip}
    subfaults2 = pf._point_sources_def(rise_time, rupture_vel, subfaults)
    strike_ps = subfaults2["strike_ps"]
    dip_ps = subfaults2["dip_ps"]
    t1 = rise_time["min_rise"]
    t2 = rise_time["delta_rise"]
    windows = rise_time["windows"]

    hyp_stk = segments[0]["hyp_stk"]
    hyp_dip = segments[0]["hyp_dip"]
    delta_strike = segments[0]["delta_strike"]
    delta_dip = segments[0]["delta_dip"]

    depth = tensor_info["depth"]
    point_sources = pf.point_sources_param(
        segments, tensor_info, rise_time, connections=connections
    )
    shear = pf.shear_modulous(point_sources, velmodel=velmodel)

    disp_or_vel = 0
    string = "{} {} {} {} {}\n"

    with open(directory / "fault&rise_time.txt", "w") as outfile:
        outfile.write("{} {} {} 10\n".format(hyp_stk, hyp_dip, depth))
        outfile.write(
            "{} {} {} {} {} {} {} {} {}\n".format(
                len(segments),
                delta_strike,
                delta_dip,
                strike_ps,
                dip_ps,
                min_vel,
                max_vel,
                -100,
                100,
            )
        )
        outfile.write(
            "{} {} {} {} {}\n".format(t1, t2, windows, rupture_vel, disp_or_vel)
        )
        for i_segment, segment in enumerate(segments):
            dip = segment["dip"]
            strike = segment["strike"]
            rake = segment["rake"]
            n_stk = segment["stk_subfaults"]
            n_dip = segment["dip_subfaults"]
            delay = 0
            if "delay_segment" in segment:
                delay = segment["delay_segment"]
            hyp_stk = segment["hyp_stk"]
            hyp_dip = segment["hyp_dip"]
            outfile.write("{} {} {}\n".format(i_segment + 1, dip, strike))
            outfile.write("{} {} {}\n".format(n_stk, n_dip, delay))
            for i in range(n_dip):
                for j in range(n_stk):
                    slip = 300 if j == hyp_stk - 1 and i == hyp_dip - 1 else 0
                    outfile.write(string.format(slip, rake, 0, t1, t1))

    with open(directory / "point_sources.txt", "w") as outfile:
        for i_segment, (ps_seg, segment) in enumerate(zip(point_sources, segments)):
            dip = segment["dip"]
            strike = segment["strike"]
            n_stk = segment["stk_subfaults"]
            n_dip = segment["dip_subfaults"]
            outfile.write("{} {} {}\n".format(i_segment + 1, dip, strike))
            for j1 in range(n_dip):
                for i1 in range(n_stk):
                    for j2 in range(dip_ps):
                        for i2 in range(strike_ps):
                            outfile.write(
                                "{} {} {} {} {} {} {}\n".format(*ps_seg[j1, i1, j2, i2])
                            )

    with open(directory / "shear_model.txt", "w") as outfile:
        outfile.write("{}\n".format(len(shear)))
        for i_segment, (shear_seg, segment) in enumerate(zip(shear, segments)):
            n_stk = segment["stk_subfaults"]
            n_dip = segment["dip_subfaults"]
            outfile.write("{} {} {}\n".format(i_segment + 1, n_stk, n_dip))
            ratio = len(shear_seg[0, :]) // 5
            format_str = ("{} " * 5 + "\n") * ratio
            remain = len(shear_seg[0, :]) % 5
            format_str = (
                format_str if remain == 0 else format_str + ("{} " * remain + "\n")
            )
            for i in range(n_dip):
                outfile.write(format_str.format(*shear_seg[i, :]))
    return


################################
# station data
################################


def input_chen_tele_body(
    tensor_info: dict,
    data_prop: dict,
    directory: Union[pathlib.Path, str] = pathlib.Path(),
) -> Optional[str]:
    """Write text files, which are inputs for Chen's scripts, with information
       about the teleseismic body wave data

    :param tensor_info: The tensor information (from tensor_info.json)
    :type tensor_info: dict
    :param data_prop: The sampling/filtering properties (from sampling_filter.json)
    :type data_prop: dict
    :param directory: The directory to read/write to, defaults to pathlib.Path()
    :type directory: Union[pathlib.Path, str], optional
    :return: The data type
    :rtype: Optional[str]
    """
    directory = pathlib.Path(directory)
    if not os.path.isfile(directory / "tele_waves.json"):
        return None
    with open(directory / "tele_waves.json", "r") as t:
        traces_info = json.load(t)

    date_origin = UTCDateTime(tensor_info["datetime"])
    dt = traces_info[0]["dt"]
    dt = round(dt, 1)
    filtro = data_prop["tele_filter"]
    low_freq = filtro["low_freq"]
    high_freq = filtro["high_freq"]

    with open(directory / "filtro_tele.txt", "w") as outfile:
        outfile.write("Corners: {} {}\n".format(low_freq, high_freq))
        outfile.write("dt: {}".format(dt))

    nsta = len(traces_info)
    model = TauPyModel(model="ak135f_no_mud")
    depth = tensor_info["depth"]
    event_lat = tensor_info["lat"]
    event_lon = tensor_info["lon"]

    string = (
        "{0:2d}   FAR GDSN {1:>6} {1:>6}BHZ.DAT {2:5.2f} {3:6.2f} "
        "{4:5.2f} {5:6.2f} {6:6.2f} {7} 0 {8} 0 {9}  1 0\n"
    )
    sin_fun = lambda p: p * 3.6 / 111.12
    angle_fun = (
        lambda p: np.arctan2(sin_fun(p), np.sqrt(1 - sin_fun(p) ** 2)) * 180.0 / np.pi
    )
    string_fun1 = (
        lambda i, name, dist, az, lat, lon, p_slowness, disp_or_vel: string.format(
            i, name, dist, az, lat, lon, angle_fun(p_slowness), 1.0, disp_or_vel, 0
        )
    )
    string_fun2 = (
        lambda i, name, dist, az, lat, lon, s_slowness, disp_or_vel: string.format(
            i, name, dist, az, lat, lon, angle_fun(s_slowness), 4.0, disp_or_vel, 2
        )
    )

    with open(directory / "channels_body.txt", "w") as outfile:
        outfile.write("30 30 30 0 0 0 0 0 0 1.1e+20\n")
        outfile.write(
            "3 10 {}\n{}{}{}{}{}{}.{}\n{}\n".format(
                dt,
                date_origin.year,
                date_origin.month,
                date_origin.day,
                date_origin.hour,
                date_origin.minute,
                date_origin.second,
                date_origin.microsecond,
                nsta,
            )
        )
        i = 0
        for file in traces_info:  # header in headers:
            name = file["name"]
            channel = file["component"]
            lat, lon = file["location"]
            dist, az, back_azimuth = mng._distazbaz(lat, lon, event_lat, event_lon)
            dist = kilometers2degrees(dist)
            derivative = False if not "derivative" in file else file["derivative"]
            derivative = int(derivative)
            arrivals = mng.theoretic_arrivals(model, dist, depth)
            p_slowness = arrivals["p_slowness"]
            s_slowness = arrivals["s_slowness"]
            slowness = p_slowness if channel == "BHZ" else s_slowness
            string_fun3 = string_fun1 if channel == "BHZ" else string_fun2
            outfile.write(
                string_fun3(i + 1, name, dist, az, lat, lon, slowness, derivative)
            )
            i = i + 1

    with (
        open(directory / "wavelets_body.txt", "w") as file1,
        open(directory / "waveforms_body.txt", "w") as file2,
    ):
        write_files_wavelet_observed(file1, file2, dt, data_prop, traces_info)
    #
    # instrumental response common to all body waves
    #
    string2 = (
        "\n3\n" + "0. 0.\n" * 3 + "4\n-6.17E-03  6.17E-03\n"
        "-6.17E-03 -6.17E-03\n-39.18    49.12\n-39.18   "
        "-49.12\n3948\n"
    )
    with open(directory / "instrumental_response.txt", "w") as outfile:
        outfile.write("{}\n".format(nsta))
        outfile.write(string2 * len(traces_info))

    write_wavelet_freqs(dt, directory / "Wavelets_tele_body.txt")

    with open(directory / "body_wave_weight.txt", "w") as outfile:
        for info in traces_info:
            sta = info["name"]
            channel = info["component"]
            weight = info["trace_weight"]
            outfile.write("{} {} {}\n".format(weight, sta, channel))
    return "body"


def input_chen_tele_surf(
    tensor_info: dict,
    data_prop: dict,
    config_path: Optional[Union[str, pathlib.Path]] = None,
    directory: Union[pathlib.Path, str] = pathlib.Path(),
):
    """Write text files, which are inputs for Chen's scripts, with information
       about the teleseismic surface wave data

    :param tensor_info: The tensor information (from tensor_info.json)
    :type tensor_info: dict
    :param data_prop: The sampling/filtering properties (from sampling_filter.json)
    :type data_prop: dict
    :param config_path: The path to the config file, defaults to None
    :type config_path: Optional[Union[str, pathlib.Path]], optional
    :param directory: The directory to read/write to, defaults to pathlib.Path()
    :type directory: Union[pathlib.Path, str], optional
    :param config_directory: The directory where the config file is located, defaults to pathlib.Path()
    :type config_directory: U Optional[Union[pathlib.Path,str]], optional
    """
    directory = pathlib.Path(directory)
    if not os.path.isfile(directory / "surf_waves.json"):
        return
    if config_path is not None:
        dirs = mng.default_dirs(config_path=config_path)
    else:
        dirs = mng.default_dirs()
    gf_bank = dirs["long_gf_bank"]
    with open(directory / "surf_waves.json", "r") as t:
        traces_info = json.load(t)

    depth = tensor_info["depth"]
    event_lat = tensor_info["lat"]
    event_lon = tensor_info["lon"]

    nsta = len(traces_info)

    filtro = data_prop["surf_filter"]
    freq1 = filtro["freq1"]
    freq2 = filtro["freq2"]
    freq3 = filtro["freq3"]
    freq4 = filtro["freq4"]
    with open(directory / "surf_filter.txt", "w") as outfile:
        outfile.write("{} {} {} {}".format(freq1, freq2, freq3, freq4))
    date_origin = UTCDateTime(tensor_info["datetime"])
    string = (
        "{:3d} {:>6} {:>8.3f} {:>9.3f} 31"
        + 3 * "  {:>1}"
        + 2 * " {:>7.2f}"
        + "  1"
        + 3 * "  {:>7.2f}"
        + " 0\n"
    )
    string_fun = lambda i, name, lat, lon, a, b, c, d, e, weight: string.format(
        i, name, lat, lon, a, b, c, d, e, weight, weight, weight
    )

    with open(directory / "channels_surf.txt", "w") as outfile:
        outfile.write(
            "{}{}{}{}{}{}.{}\n".format(
                date_origin.year,
                date_origin.month,
                date_origin.day,
                date_origin.hour,
                date_origin.minute,
                date_origin.second,
                date_origin.microsecond,
            )
        )
        outfile.write(
            "{} {} {} {} {} {} {} {} {}\n".format(
                event_lat,
                event_lon,
                depth,
                date_origin.year,
                date_origin.julday,
                date_origin.hour,
                date_origin.minute,
                date_origin.second,
                date_origin.microsecond,
            )
        )
        outfile.write("0.0 90.0 0.0 10 4.0 1.0e+26\n")
        outfile.write("4.0 4.0 10 1.0 {}\n".format(0))
        outfile.write("{} {}\n".format(nsta, nsta))
        outfile.write("No STA Lat Lon M V H1 H2 Angle1 Angle2 Io_s Weight\n")
        i = 0
        for file in traces_info:
            weight = file["trace_weight"]
            name = file["name"]
            channel = file["component"]
            lat, lon = file["location"]
            if channel == "BHZ":  # Rayleigh
                outfile.write(string_fun(i + 1, name, lat, lon, 1, 0, 0, 0, 0, weight))
            else:  # Love
                outfile.write(string_fun(i + 1, name, lat, lon, 0, 1, 0, 90, 0, weight))
            i = i + 1

    with (
        open(directory / "wavelets_surf.txt", "w") as file1,
        open(directory / "waveforms_surf.txt", "w") as file2,
    ):
        write_files_wavelet_observed(
            file1, file2, 4.0, data_prop, traces_info, gf_bank=gf_bank
        )

    write_wavelet_freqs(4.0, directory / "Wavelets_surf_tele.txt")


def input_chen_near_field(
    tensor_info: dict,
    data_prop: dict,
    data_type: Literal["cgps", "strong"],
    directory: Union[pathlib.Path, str] = pathlib.Path(),
) -> Optional[str]:
    """Write text files, which are inputs for Chen's scripts, with information
       about the strong motion data

    :param tensor_info: The tensor information (from tensor_info.json)
    :type tensor_info: dict
    :param data_prop: The sampling/filtering properties (from sampling_filter.json)
    :type data_prop: dict
    :param data_type: The data type
    :type data_type: Literal[&quot;cgps&quot;, &quot;strong_motion&quot;]
    :param directory: The directory to read/write to, defaults to pathlib.Path()
    :type directory: Union[pathlib.Path, str], optional
    :return: The data type
    :rtype: Optional[str]

    .. warning::

        Make sure the filters of strong motion data agree with the values in
        sampling_filter.json!
    """
    directory = pathlib.Path(directory)
    if data_type == "strong":
        dict1 = directory / "strong_motion_waves.json"
        filename1 = directory / "filtro_strong.txt"
        filtro = data_prop["strong_filter"]
        filename2 = directory / "channels_strong.txt"
        filename3 = directory / "wavelets_strong.txt"
        filename4 = directory / "waveforms_strong.txt"
        filename5 = directory / "Wavelets_strong_motion.txt"
    elif data_type == "cgps":
        dict1 = directory / "cgps_waves.json"
        filename1 = directory / "filtro_cgps.txt"
        filtro = data_prop["cgps_filter"]
        filename2 = directory / "channels_cgps.txt"
        filename3 = directory / "wavelets_cgps.txt"
        filename4 = directory / "waveforms_cgps.txt"
        filename5 = directory / "Wavelets_cgps.txt"
    else:
        return

    if not os.path.isfile(dict1):
        return None

    with open(dict1, "r") as t:
        traces_info = json.load(t)
    date_origin = UTCDateTime(tensor_info["datetime"])
    moment_mag = tensor_info["moment_mag"]
    event_lat = tensor_info["lat"]
    event_lon = tensor_info["lon"]
    depth = tensor_info["depth"]
    dt_strong = traces_info[0]["dt"]
    dt_strong = round(dt_strong, 2)
    low_freq = filtro["low_freq"]
    high_freq = filtro["high_freq"]

    nsta = len(traces_info)

    with open(filename1, "w") as outfile:
        outfile.write("Corners: {} {}".format(low_freq, high_freq))

    disp_or_vel = 0
    string = "{0:3d} {1:>5}{2:>9.3f}{3:>10.3f} 31{4:>5} {5} 0\n"
    string_fun = lambda i, name, lat, lon, a, w: string.format(
        i + 1, name, lat, lon, a, w
    )

    with open(filename2, "w") as outfile:
        outfile.write(
            "{}{}{}{}{}{}{}\n".format(
                date_origin.year,
                date_origin.month,
                date_origin.day,
                date_origin.hour,
                date_origin.minute,
                date_origin.second,
                date_origin.microsecond,
            )
        )
        outfile.write("{} {} {}\n".format(event_lat, event_lon, depth))
        outfile.write("10 {} {}\n".format(dt_strong, moment_mag))
        outfile.write("{}\n".format(disp_or_vel))
        outfile.write("{} {}\n".format(nsta, nsta))
        outfile.write("No STA Lat Lon M Comp Weight\n")
        for i, file in enumerate(traces_info):
            weight = file["trace_weight"]
            name = file["name"]
            channel = file["component"]
            lat, lon = file["location"]
            outfile.write(string_fun(i, name, lat, lon, channel, weight))

    with open(filename3, "w") as file1, open(filename4, "w") as file2:
        write_files_wavelet_observed(file1, file2, dt_strong, data_prop, traces_info)

    write_wavelet_freqs(dt_strong, directory / "Wavelets_strong_motion")
    return "strong"


def input_chen_static(
    directory: Union[pathlib.Path, str] = pathlib.Path(),
):
    """Write text files, which are inputs for Chen's scripts, with information
       about the static gps data

    :param directory: The directory to read/write to, defaults to pathlib.Path()
    :type directory: Union[pathlib.Path, str], optional
    """
    directory = pathlib.Path(directory)
    if not os.path.isfile(directory / "static_data.json"):
        return

    with open(directory / "static_data.json", "r") as s:
        static_info = json.load(s)
    string = "{0:3d} {1:>5}{2:>10.3f}{3:>10.3f} {4} {5} {6} {7} {8} {9}\n"
    string_fun = lambda i, name, lat, lon, a, b, c, d, e, f: string.format(
        i, name, lat, lon, a, b, c, d, e, f
    )
    with open(directory / "static_data.txt", "w") as outfile:
        outfile.write("{}\n\n".format(len(static_info)))
        for i, info in enumerate(static_info):
            name = info["name"]
            lat, lon = info["location"]
            ud_disp, ns_disp, ew_disp = info["observed"]
            ud_weight, ns_weight, ew_weight = info["trace_weight"]
            outfile.write(
                string_fun(
                    i,
                    name,
                    lat,
                    lon,
                    ud_disp,
                    ns_disp,
                    ew_disp,
                    ud_weight,
                    ns_weight,
                    ew_weight,
                )
            )


def input_chen_dart(
    tensor_info: dict,
    data_prop: dict,
    directory: Union[pathlib.Path, str] = pathlib.Path(),
) -> Optional[str]:
    """Write text files, which are inputs for Chen's scripts, with information
       about the dart data

    :param tensor_info: The tensor information (from tensor_info.json)
    :type tensor_info: dict
    :param data_prop: The sampling/filtering properties (from sampling_filter.json)
    :type data_prop: dict
    :param directory: The directory to read/write to, defaults to pathlib.Path()
    :type directory: Union[pathlib.Path, str], optional
    :return: The data type
    :rtype: Optional[str]
    """
    directory = pathlib.Path(directory)
    if not os.path.isfile(directory / "dart_waves.json"):
        return None
    with open(directory / "dart_waves.json", "r") as t:
        traces_info = json.load(t)
    date_origin = UTCDateTime(tensor_info["datetime"])
    moment_mag = tensor_info["moment_mag"]
    event_lat = tensor_info["lat"]
    event_lon = tensor_info["lon"]
    depth = tensor_info["depth"]
    dt_dart = traces_info[0]["dt"]
    dt_dart = round(dt_dart, 2)

    nsta = len(traces_info)

    io_vd = 0
    string = "{0:3d} {1:>5}{2:>9.3f}{3:>10.3f} 31{4:>5} {5} 0\n"
    string_fun = lambda i, name, lat, lon, a, w: string.format(
        i + 1, name, lat, lon, a, w
    )

    with open(directory / "channels_dart.txt", "w") as outfile:
        outfile.write(
            "{}{}{}{}{}{}{}\n".format(
                date_origin.year,
                date_origin.month,
                date_origin.day,
                date_origin.hour,
                date_origin.minute,
                date_origin.second,
                date_origin.microsecond,
            )
        )
        outfile.write("{} {} {}\n".format(event_lat, event_lon, depth))
        outfile.write("10 {} {}\n".format(dt_dart, moment_mag))
        outfile.write("{}\n".format(io_vd))
        outfile.write("{} {}\n".format(nsta, nsta))
        outfile.write("No STA Lat Lon M V H1 H2 Weight\n")
        for i, file in enumerate(traces_info):
            name = file["name"]
            channel = file["component"]
            lat, lon = file["location"]
            weight = file["trace_weight"]
            outfile.write(string_fun(i, name, lat, lon, channel, weight))

    with (
        open(directory / "wavelets_dart.txt", "w") as file1,
        open(directory / "waveforms_dart.txt", "w") as file2,
    ):
        write_files_wavelet_observed(
            file1, file2, dt_dart, data_prop, traces_info, dart=True, zero_start=True
        )
    return "cgps"


def input_chen_insar(
    directory: Union[pathlib.Path, str] = pathlib.Path(),
):
    """Write text files, which are inputs for Chen's scripts, with information
       about the insar data

    :param directory: The directory to read/write to, defaults to pathlib.Path()
    :type directory: Union[pathlib.Path, str], optional
    """
    directory = pathlib.Path(directory)
    if not os.path.isfile(directory / "insar_data.json"):
        return

    with open(directory / "insar_data.json", "r") as i:
        insar_info = json.load(i)
    lines: List[str] = []
    weights: List[float] = []
    sizes: List[int] = []
    ramps: List[float] = []
    if "ascending" in insar_info:
        properties = insar_info["ascending"]
        for asc_property in properties:
            track = asc_property["name"]
            weights = weights + [asc_property["weight"]]
            ramps = ramps + [asc_property["ramp"]]
            new_lines: List[str] = []
            with open(track, "r") as infile:
                for line in infile:
                    if line.startswith("#"):
                        continue
                    else:
                        new_lines.append(line.split())  # type: ignore
            lines = lines + [new_lines]  # type: ignore
            sizes = sizes + [len(new_lines)]  # type: ignore
        lines_asc = len(lines)
    if "descending" in insar_info:
        properties = insar_info["descending"]
        for desc_property in properties:
            track = desc_property["name"]
            weights = weights + [desc_property["weight"]]
            ramps = ramps + [desc_property["ramp"]]
            new_lines = []
            with open(track, "r") as infile:
                for line in infile:
                    if line.startswith("#"):
                        continue
                    else:
                        new_lines.append(line.split())  # type: ignore
            lines = lines + [new_lines]  # type: ignore
            sizes = sizes + [len(new_lines)]  # type: ignore
    points = np.sum(sizes)

    string = (
        "{0:3d} {1:>5} {2:>12.6f} {3:>12.6f} {4:>12.6f} {5:>12.6f}"
        " {6:>12.6f} {7:>12.6f}\n"
    )
    string_fun = lambda i, name, lat, lon, a, b, c, d: string.format(
        i, name, lat, lon, a, b, c, d
    )
    points2 = 0
    with open(directory / "insar_data.txt", "w") as outfile:
        outfile.write("{}\n\n".format(points))
        for new_lines in lines:  # type: ignore
            for i, line in enumerate(new_lines):  # type: ignore
                points2 = points2 + 1
                lat = float(line[1])
                lon = float(line[0])
                observed = 100 * float(line[2])
                look_ew = float(line[3])  # look_ew
                look_ns = float(line[4])  # loow_ns
                look_ud = float(line[5])  # look_ud
                outfile.write(
                    string_fun(
                        points2, i, lat, lon, observed, look_ud, look_ns, look_ew
                    )
                )

    with open(directory / "insar_weights.txt", "w") as outfile:
        outfile.write("{}\n".format(len(weights)))
        zipped = zip(sizes, weights)
        for length, weight in zipped:
            outfile.write("{} {}\n".format(length, weight))

    if not any(ramps):
        if os.path.isfile(directory / "ramp_gf.txt"):
            os.remove(directory / "ramp_gf.txt")
        return

    latitudes: List[float] = []
    longitudes: List[float] = []
    for new_lines in lines:  # type: ignore
        latitudes = latitudes + [float(line[1]) for line in new_lines]
        longitudes = longitudes + [float(line[0]) for line in new_lines]
    ref_lon = longitudes[0]
    utm_coords = [
        mng.coords2utm(lat, lon, ref_lon) for lat, lon in zip(latitudes, longitudes)
    ]
    eastings = [easting for easting, northing in utm_coords]
    northings = [northing for easting, northing in utm_coords]
    min_northing = np.min(northings)
    min_easting = np.min(eastings)

    block_matrix = None
    start = 0
    for ramp, length in zip(ramps, sizes):
        if ramp is not None:
            size1 = 3
            if ramp == "bilinear":
                size1 = 6
            elif ramp == "quadratic":
                size1 = 5
            east1 = np.array(eastings[start : start + length]) - min_easting
            north1 = np.array(northings[start : start + length]) - min_northing
            east1 = east1 / np.max(np.abs(east1))
            north1 = north1 / np.max(np.abs(north1))
            east2 = east1**2
            north2 = north1**2
            east2 = east2 / np.max(east2)
            north2 = north2 / np.max(north2)
            east_north = east1 * north1
            east_north = east_north / np.max(east_north)
            if ramp == "linear":
                size1 = 3
                zipped = zip(east1, north1)
                gf_ramp2: Union[np.ndarray[Any, np.dtype[Any]], List[List[float]]] = [
                    [east, north, 1] for east, north in zipped
                ]
                gf_ramp2 = np.array(gf_ramp2)
            elif ramp == "bilinear":
                size1 = 3
                gf_ramp2 = [
                    [e1, n1, 1, en, e2, n2]
                    for e1, n1, en, e2, n2 in zip(
                        east1, north1, east_north, east2, north2
                    )
                ]
                gf_ramp2 = np.array(gf_ramp2)
            elif ramp == "quadratic":
                gf_ramp2 = [
                    [e1, n1, 1, en, en**2]
                    for e1, n1, en in zip(east1, north1, east_north)
                ]
                gf_ramp2 = np.array(gf_ramp2)
            if block_matrix is None:
                block_matrix = np.block(gf_ramp2)
            elif len(block_matrix) == 0:
                block_matrix = np.block(gf_ramp2)
            else:
                shape1 = block_matrix.shape
                rows1, cols1 = shape1
                block_matrix = np.block(
                    [
                        [block_matrix, np.zeros((rows1, size1))],
                        [np.zeros((length, cols1)), gf_ramp2],
                    ]
                )
        start = start + length

    with open(directory / "ramp_gf.txt", "w") as outf:
        string = " ".join(str(v) for v in ramps) + "\n"
        outf.write(string)
        shape = block_matrix.shape  # type: ignore
        rows, cols = shape
        for row in range(0, rows):
            new_row = [str(a) for a in block_matrix[row]]  # type: ignore
            string = " ".join(new_row) + "\n"
            outf.write(string)


def write_files_wavelet_observed(
    wavelet_file: io.TextIOWrapper,
    obse_file: io.TextIOWrapper,
    dt: float,
    data_prop: dict,
    traces_info: List[dict],
    gf_bank: Optional[Union[pathlib.Path, str]] = None,
    zero_start: bool = True,
    dart: bool = False,
):
    """Write files with observed waveforms and properties of wavelet for all
    selected stations and channels


    :param wavelet_file: The open/writeable wavelet file
    :type wavelet_file: io.TextIOWrapper
    :param obse_file: The open/writeable observation file
    :type obse_file: io.TextIOWrapper
    :param dt: The dt of the waveform
    :type dt: float
    :param data_prop: The sampling/filtering properties (from sampling_filter.json)
    :type data_prop: dict
    :param traces_info: The trace information to write out
    :type traces_info: List[dict]
    :param gf_bank: The path to the gf_bank file, defaults to None
    :type gf_bank: Optional[Union[pathlib.Path, str]], optional
    :param zero_start: Whether to consider a zero start, defaults to True
    :type zero_start: bool, optional
    :param dart: Whether the data is dart data, defaults to False
    :type dart: bool, optional
    """
    string = lambda name, channel, a, b: "{} {}\n{}{}".format(name, channel, a, b)

    n_begin, n_end = data_prop["wavelet_scales"]
    input_length = 256
    wavelet_file.write("{} {} {}\n".format(n_begin, n_end, input_length))
    if not gf_bank:
        wavelet_file.write("\n")
    else:
        wavelet_file.write("{}\n".format(gf_bank))
    wavelet_file.write("{}\n".format(len(traces_info)))
    error_norm = "3 " * (n_end - n_begin) + "3\n"
    for file in traces_info:
        name = file["name"]
        channel = file["component"]
        ffm_duration = file["duration"]
        error_norm = "3 " * (n_end - n_begin) + "3\n"
        derivative = False if not "derivative" in file else file["derivative"]
        if derivative:
            error_norm = "3 " * (5 - n_begin) + "4 " * (n_end - 5) + "4\n"
        wavelet_file.write(string(name, channel, error_norm, file["wavelet_weight"]))
        if file["file"]:
            start = file["start_signal"]
            stream = __get_stream(file)
            waveform = stream[0].data[start:]
            if zero_start:
                stream[0].data = stream[0].data - waveform[0]
                stream.write(file["file"], format="SAC", byteorder=0)
                waveform = waveform - waveform[0]
            waveform = np.gradient(waveform, dt) if derivative else waveform
            del stream
        else:
            waveform = [5000 * np.random.randn() for i in range(ffm_duration)]
        write_observed_file(file, dt, obse_file, waveform, dart=dart)
    return


def __get_stream(file: dict) -> Stream:
    """Read data into an obspy Stream

    :param file: The dictionary with the data file's properties
    :type file: dict
    :return: The stream
    :rtype: Stream
    """
    stream = None
    i = 8
    expected_error = True
    while expected_error and i >= 0:
        try:
            stream = read(file["file"], format="SAC")
            expected_error = False
        except IndexError:
            print(i)
            print(
                "Obspy bug when reading the file {}. "
                "Waveform set to random".format(file["file"])
            )
            i -= 1
        except SacIOError:
            print(i)
            print(
                "Obspy bug when reading the file {}. "
                "Waveform set to random".format(file["file"])
            )
            i -= 1
    return stream


def write_observed_file(
    file: dict,
    dt: float,
    data_file: io.TextIOWrapper,
    waveform: np.ndarray[Any, np.dtype[np.signedinteger]],
    dart: bool = False,
):
    """Routine for writing the computing Obser.x file

    :param file: The dictionary with the data file's properties
    :type file: dict
    :param dt: The dt of the timeseries data
    :type dt: float
    :param data_file: The open/writeable data file
    :type data_file: io.TextIOWrapper
    :param waveform: The waveform data
    :type waveform: np.ndarray[Any, np.signedinteger]
    :param dart: Whether the data is dart data, defaults to False
    :type dart: bool, optional
    """
    ffm_duration = file["duration"]
    channel = file["component"] if not dart else "dart"
    name = file["name"]
    length = len(waveform)
    trace = ["{}\n".format(val) for val in waveform]
    trace_str = "".join(trace)
    data_file.write(
        "name: {}\nchannel: {}\ndt: {}\nlength: {}\nffm_duration: {}\n"
        "data: \n".format(name, channel, dt, length, ffm_duration)
    )
    data_file.write(trace_str)
    return


def write_wavelet_freqs(dt: float, name: Union[pathlib.Path, str]):
    """Write a range of frequencies modelled by different wavelet coefficients

    :param dt: The dt of the timeseries data
    :type dt: float
    :param name: The path to the file
    :type name: Union[pathlib.Path,str]
    """
    with open(name, "w") as outfile:
        for j in range(1, 9):
            min_freq = float(2**j) / float(3 * 2**10 * dt)
            outfile.write(
                "j :{}\nFrequency range for these wavelet coefficients is"
                ": {:.4f} {:.4f} Hz\n".format(j, min_freq, 4 * min_freq)
            )


def from_synthetic_to_obs(
    files: List[dict],
    data_type: str,
    data_prop: dict,
    add_error: bool = False,
    directory: Union[pathlib.Path, str] = pathlib.Path(),
):
    """Write synthetic waveforms to files

    :param files: List of dictionaries with data file properties
    :type files: List[dict]
    :param data_type: The data type
    :type data_type: Literal[&quot;dart&quot;, &quot;cgps&quot;, &quot;gps&quot;, &quot;strong_motion&quot;, &quot;surf_tele&quot;, &quot;tele_body&quot;]
    :param data_prop: The sampling/filtering properties (from sampling_filter.json)
    :type data_prop: dict
    :param add_error: Whether to add error, defaults to False
    :type add_error: bool, optional
    :param directory: The directory to read/write to, defaults to pathlib.Path()
    :type directory: Union[pathlib.Path, str], optional
    """
    directory = pathlib.Path(directory)
    dt = files[0]["dt"]
    filtro_strong = data_prop["strong_filter"]
    filtro_cgps = data_prop["strong_filter"]
    if "cgps_filter" in data_prop:
        filtro_cgps = data_prop["cgps_filter"]
    filtro_tele = data_prop["tele_filter"]
    nyq = 0.5 / dt if not data_type == "gps" else 10000
    std_shift: Union[float, int] = 0
    if data_type == "strong":
        max_val = 0.1
        syn_file = directory / "synthetics_strong.txt"
        obser_file = directory / "waveforms_strong.txt"
        std_shift = 2
        low_freq = filtro_strong["low_freq"]
        high_freq = filtro_strong["high_freq"]
        corners = [low_freq / nyq, high_freq / nyq]
        filters = ["highpass", "lowpass"]
        orders = [4]
    if data_type == "cgps":
        max_val = 0.1
        syn_file = directory / "synthetics_cgps.txt"
        obser_file = directory / "waveforms_cgps.txt"
        std_shift = 0.5
        low_freq = 0
        high_freq = filtro_cgps["high_freq"]
        corners = [high_freq / nyq]
        filters = ["lowpass"]
        orders = [4]
    if data_type == "body":
        max_val = 10
        syn_file = directory / "synthetics_body.txt"
        obser_file = directory / "waveforms_body.txt"
        std_shift = 0.5
        high_freq = 1.0
        low_freq = filtro_tele["low_freq"]
        high_freq = filtro_tele["high_freq"]
        corners = [low_freq / nyq, high_freq / nyq]
        filters = ["highpass", "lowpass"]
        orders = [2, 2]
    if data_type == "surf":
        max_val = 0.01
        syn_file = directory / "synthetics_surf.txt"
        obser_file = directory / "waveforms_surf.txt"
        std_shift = 2
        low_freq = 0.004
        high_freq = 0.006
        corners = [[low_freq / nyq, high_freq / nyq]]
        filters = ["bandpass"]
        orders = [2]
    if data_type == "dart":
        max_val = 0.005
        syn_file = directory / "synthetics_dart.txt"
        obser_file = directory / "waveforms_dart.txt"
        std_shift = 30
        corners = []
        filters = []
        orders = []

    dart = "dart" in data_type
    string = "{0:3d} {1:>5}{2:>10.3f}{3:>10.3f} {4} {5} {6} {7} {8} {9}\n"
    string_fun = lambda i, name, lat, lon, a, b, c, d, e, f: string.format(
        i, name, lat, lon, a, b, c, d, e, f
    )
    if not data_type == "gps":
        files = get_outputs.get_data_dict(files, syn_file=syn_file)
        with open(obser_file, "w") as outfile:
            for file in files:
                dt = file["dt"]
                channel = file["component"]
                max_val0 = max_val
                if data_type == "cgps" and channel[-1] == "Z":
                    max_val0 = 5 * max_val
                waveform = 0 * np.array(file["synthetic"])
                shift = np.random.randn(1) * std_shift
                shift = int(shift / dt)  # type:ignore
                length = len(waveform)
                error = np.zeros(length)
                if shift > 0:
                    waveform[shift:] = file["synthetic"][:-shift]
                elif shift < 0:
                    waveform[:shift] = file["synthetic"][-shift:]
                else:
                    waveform = file["synthetic"]
                if add_error:
                    error = max_val0 * np.random.randn(length)
                for order, filter, corner in zip(orders, filters, corners):
                    b, a = butter(order, corner, btype=filter)
                    error = filtfilt(b, a, error)
                waveform = waveform + error  # type: ignore
                write_observed_file(
                    file, dt, outfile, waveform, dart=dart
                )  # type:ignore
    else:
        (
            names,
            lats,
            lons,
            observed,
            synthetic,
            error,
        ) = get_outputs.retrieve_gps()  # type:ignore
        with open(directory / "static_data.txt", "r") as infile:
            orig_lines = [line.split() for line in infile]
        with open(directory / "static_data.txt", "w") as outfile:
            outfile.write("{}\n\n".format(orig_lines[0][0]))
            for i, line in enumerate(orig_lines[2:]):
                name = line[1]
                lat = float(line[2])
                lon = float(line[3])
                new_obs = next(
                    syn for name2, syn in zip(names, synthetic) if name2 == name
                )
                weight1 = float(line[7])
                weight2 = float(line[8])
                weight3 = float(line[9])
                error = np.random.randn(3)
                new_obs[0] = float(new_obs[0]) + 0.1 * error[0]
                new_obs[1] = float(new_obs[1]) + 0.1 * error[1]
                new_obs[2] = float(new_obs[2]) + 0.5 * error[2]
                outfile.write(
                    string_fun(
                        i,
                        name,
                        lat,
                        lon,
                        new_obs[0],
                        new_obs[1],
                        new_obs[2],
                        weight1,
                        weight2,
                        weight3,
                    )
                )


###########################
# inputs annealing
###########################


def inputs_simmulated_annealing(
    dictionary: dict,
    directory: Union[pathlib.Path, str] = pathlib.Path(),
):
    """Write the HEAT.IN file with modelling information

    :param dictionary: The modelling dictionary
    :type dictionary: dict
    :param directory: The directory to read/write to, defaults to pathlib.Path()
    :type directory: Union[pathlib.Path, str], optional
    """
    directory = pathlib.Path(directory)
    moment_mag = dictionary["seismic_moment"]
    weight0 = dictionary["moment_weight"]
    weight1 = dictionary["slip_weight"]
    weight2 = dictionary["time_weight"]
    source_dur = dictionary["max_source_dur"]
    iters = dictionary["iterations"]
    cooling_rate = dictionary["cooling_rate"]
    initial_temp = dictionary["initial_temperature"]

    type_of_inversion = 1
    with open(directory / "annealing.txt", "w") as filewrite:
        filewrite.write("{} -7 {} {} 90\n".format(iters, type_of_inversion, moment_mag))
        filewrite.write(
            "{} {} {} 0.1 {} {} {}\n".format(
                initial_temp, cooling_rate, 4 * 10**-6, weight0, weight1, weight2
            )
        )
        filewrite.write("0 {} 0 {}\n".format(10**-4, source_dur))
        filewrite.write("1\n")
    return


def model_space(
    segments: dict,
    directory: Union[pathlib.Path, str] = pathlib.Path(),
):
    """Write a file describing the file of feasible FFM models

    :param segments: Segment data (from segments_data.json)
    :type segments: dict
    :param directory: The directory to read/write to, defaults to pathlib.Path()
    :type directory: Union[pathlib.Path, str], optional
    """
    directory = pathlib.Path(directory)
    with open(directory / "model_space.txt", "w") as filewrite:
        filewrite.write("0\n")
        for i, segment in enumerate(segments):
            peak_upper_slip = segment["max_upper_slip"]
            peak_lower_slip = segment["max_lower_slip"]
            peak_left_slip = segment["max_left_slip"]
            peak_right_slip = segment["max_right_slip"]
            peak_slip = segment["max_center_slip"]
            peak_slip_delta = segment["max_slip_delta"]
            nstep = segment["slip_step"]
            rake_max = segment["rake_max"]
            rake_min = segment["rake_min"]
            rstep = segment["rake_step"]
            filewrite.write("{}\n".format(i + 1))
            filewrite.write("{} {}\n".format(peak_slip_delta, peak_slip_delta))
            filewrite.write("1 1 1 1\n")
            filewrite.write("{} 0.0 {}\n".format(peak_slip, nstep))
            filewrite.write("{} 0.0 {}\n".format(peak_left_slip, nstep))
            filewrite.write("{} 0.0 {}\n".format(peak_right_slip, nstep))
            filewrite.write("{} 0.0 {}\n".format(peak_upper_slip, nstep))
            filewrite.write("{} 0.0 {}\n".format(peak_lower_slip, nstep))
            filewrite.write("{} {} {}\n".format(rake_max, rake_min, rstep))
            filewrite.write("2.6 2.4 3\n")
            filewrite.write("5 8\n")

    with open(directory / "regularization_borders.txt", "w") as filewrite:
        for i, segment in enumerate(segments):
            filewrite.write("{}\n".format(i + 1))

            if "regularization" not in segment:
                filewrite.write("0 0\n" * 4)
            else:
                reg_borders = segment["regularization"]
                neighbour_up = reg_borders["neighbour_up"]
                if neighbour_up:
                    up_segment = neighbour_up["segment"]
                    subfault_up_segment = neighbour_up["subfault"]
                    filewrite.write("{} {}\n".format(up_segment, subfault_up_segment))
                else:
                    filewrite.write("0 0\n")

                neighbour_down = reg_borders["neighbour_down"]
                if neighbour_down:
                    down_segment = neighbour_down["segment"]
                    subfault_down_segment = neighbour_down["subfault"]
                    filewrite.write(
                        "{} {}\n".format(down_segment, subfault_down_segment)
                    )
                else:
                    filewrite.write("0 0\n")

                neighbour_left = reg_borders["neighbour_left"]
                if neighbour_left:
                    left_segment = neighbour_left["segment"]
                    subfault_left_segment = neighbour_left["subfault"]
                    filewrite.write(
                        "{} {}\n".format(left_segment, subfault_left_segment)
                    )
                else:
                    filewrite.write("0 0\n")

                neighbour_right = reg_borders["neighbour_right"]
                if neighbour_right:
                    right_segment = neighbour_right["segment"]
                    subfault_right_segment = neighbour_right["subfault"]
                    filewrite.write(
                        "{} {}\n".format(right_segment, subfault_right_segment)
                    )
                else:
                    filewrite.write("0 0\n")
    with open(directory / "special_model_space.txt", "w") as filewrite:
        filewrite.write("0\n")
    with open(directory / "special_regularization_borders.txt", "w") as file:
        file.write("0\n")
    return


################################
# Green functions
################################


def write_green_file(
    green_dict: dict,
    cgps: bool = False,
    directory: Union[pathlib.Path, str] = pathlib.Path(),
):
    """Write the file needed to run fortran code and retrieve strong motion GF

    :param green_dict: Dictionary with Greens Function (GF) properties
    :type green_dict: dict
    :param cgps: Whether the data is cGPS data (otherwise strong motion),
                defaults to False
    :type cgps: bool, optional
    :param directory: The directory to read/write to, defaults to pathlib.Path()
    :type directory: Union[pathlib.Path, str], optional
    """
    directory = pathlib.Path(directory)
    dt = green_dict["dt"]
    min_depth = green_dict["min_depth"]
    max_depth = green_dict["max_depth"]
    min_dist = green_dict["min_dist"]
    max_dist = green_dict["max_dist"]
    location = green_dict["location"]
    time_corr = green_dict["time_corr"]
    name = directory / "Green_strong.txt" if not cgps else directory / "Green_cgps.txt"
    with open(name, "w") as green_file:
        green_file.write(
            "vel_model.txt\n{} {} 1\n{} {} 1\n".format(
                max_depth, min_depth, max_dist, min_dist
            )
        )
        green_file.write("10 {} 50000 {}\n".format(dt, time_corr))
        green_file.write(location)
