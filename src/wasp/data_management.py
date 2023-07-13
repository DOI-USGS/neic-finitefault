# -*- coding: utf-8 -*-
"""Module with routines to create json files for the different data sets to be
used in kinematic modelling, from the station and channel metadata stored in
the files with the waveforms-

.. warning::

    Currently, only waveform files in sac format are supported.
"""


import argparse
import errno
import glob
import json
import os
import pathlib
from typing import Dict, List, Literal, Optional, Tuple, Union

import numpy as np
from obspy import read
from obspy.geodetics import kilometers2degrees
from obspy.io.sac import SACTrace
from obspy.taup import TauPyModel
from scipy.stats import norm

import wasp.management as mng


def _dict_trace(
    file: Optional[Union[pathlib.Path, str]],
    name_station: str,
    channel: Optional[Union[List[str], str]],
    azimuth: Optional[float],
    distance: Optional[float],
    dt: Optional[float],
    duration: float,
    n_start_obs: int,
    trace_weight: float,
    wavelet_weight: Optional[float],
    synthetic_trace: Optional[List[float]] = None,
    observed_trace: Optional[List[float]] = None,
    location: Optional[Tuple[float, float]] = None,
    derivative: bool = False,
) -> Dict[str, Union[List[int], Tuple[float, float], float, int, str]]:
    """Organize trace information into a dictionary

    :param file: The file associated with the waveform
    :type file: Optional[Union[pathlib.Path, str]]
    :param name_station: The name of the station
    :type name_station: str
    :param channel: The channel
    :type channel: Optional[Union[List[str], str]]
    :param azimuth: The azimuth of the channel
    :type azimuth: Optional[float]
    :param distance: The distance from the station to the event
    :type distance: Optional[float]
    :param dt: Sampling distance in seconds
    :type dt: Optional[float]
    :param duration: The duration
    :type duration: float
    :param n_start_obs: Number of observations before the arrival
    :type n_start_obs: int
    :param trace_weight: The weight of the trace
    :type trace_weight: float
    :param wavelet_weight: The weight of the wavelet
    :type wavelet_weight: Optional[float]
    :param synthetic_trace: The synthetic trace, defaults to None
    :type synthetic_trace:  Optional[List[float]]
    :param observed_trace: The observed trace, defaults to None
    :type observed_trace: Optional[List[float]]
    :param location: The location (latitude, longitude), defaults to None
    :type location: Optional[Tuple[float, float]], optional
    :param derivative: Whether a derivative, defaults to False
    :type derivative: bool, optional
    :return: The information dictionary
    :rtype: Dict[str, Union[List[int], Tuple[float, float], float, int, str]]
    """
    info = {
        "file": str(file),
        "dt": dt,
        "duration": duration,
        "wavelet_weight": wavelet_weight,
        "start_signal": n_start_obs,
        "trace_weight": trace_weight,
        "synthetic": synthetic_trace,
        "observed": observed_trace,
        "name": name_station,
        "component": channel,
        "azimuth": azimuth,
        "distance": distance,
        "location": location,
        "derivative": derivative,
    }
    return info


def tele_body_traces(
    files: List[Union[pathlib.Path, str]],
    tensor_info: dict,
    data_prop: dict,
    directory: Union[pathlib.Path, str] = pathlib.Path(),
) -> dict:
    """Write json dictionary with properties for teleseismic data

    :param files: List of files holding teleseismic data
    :type files: List[Union[pathlib.Path, str]]
    :param tensor_info: The tensor information
    :type tensor_info: dict
    :param data_prop: The data properties (from sampling_filter.json)
    :type data_prop: dict
    :param directory: Where the file should be written, defaults to pathlib.Path()
    :type directory: Union[pathlib.Path, str], optional
    :return: The teleseismic properties written to tele_waves.json
    :rtype: dict

    .. warning::

        Make sure the filters of teleseismic data agree with the values in
        sampling_filter.json!
    """
    directory = pathlib.Path(directory)
    if len(files) == 0:
        return
    p_files = [file for file in files if "_BHZ" in file]
    sh_files = [file for file in files if "_SH" in file]
    p_files = select_tele_stations(p_files, "P", tensor_info)
    sh_files = select_tele_stations(sh_files, "SH", tensor_info)
    files = p_files + sh_files
    origin_time = tensor_info["date_origin"]
    headers = [SACTrace.read(file) for file in files]
    dt = headers[0].delta
    dt = round(dt, 1)
    n0, n1 = data_prop["wavelet_scales"]
    filter0 = data_prop["tele_filter"]
    duration = duration_tele_waves(tensor_info, dt)
    wavelet_weight0, wavelet_weight1 = wavelets_body_waves(
        duration, filter0, dt, n0, n1
    )
    info_traces = []
    event_lat = tensor_info["lat"]
    event_lon = tensor_info["lon"]
    depth = tensor_info["depth"]
    model = TauPyModel(model="ak135f_no_mud")
    header = (header for header in headers)

    for file, header in zip(files, headers):
        stream = read(file)
        __failsafe(filter0, header)
        distance, azimuth, back_azimuth = mng._distazbaz(
            header.stla, header.stlo, event_lat, event_lon
        )
        arrivals = mng.theoretic_arrivals(model, distance / 111.11, depth)
        if header.kcmpnm == "BHZ":
            arrival = arrivals["p_arrival"][0].time
            weight = 1.0
            wavelet_weight = wavelet_weight0
        elif header.kcmpnm == "SH":
            arrival = arrivals["s_arrival"][0].time
            weight = 0.5
            wavelet_weight = wavelet_weight1
        starttime = stream[0].stats.starttime
        begin = starttime - origin_time
        n_start_obs = int((arrival - begin) / dt + 0.5)
        info = _dict_trace(
            file=file,
            name_station=header.kstnm,
            channel=header.kcmpnm,
            azimuth=azimuth,
            distance=distance / 111.11,
            dt=dt,
            duration=duration,
            n_start_obs=n_start_obs,
            trace_weight=weight,
            wavelet_weight=wavelet_weight,
            synthetic_trace=[],
            location=[header.stla, header.stlo],
            derivative=False,
        )
        info_traces.append(info)
    with open(directory / "tele_waves.json", "w") as f:
        json.dump(
            info_traces,
            f,
            sort_keys=True,
            indent=4,
            separators=(",", ": "),
            ensure_ascii=False,
        )
    return info_traces


def tele_surf_traces(
    files: List[Union[pathlib.Path, str]],
    tensor_info: dict,
    data_prop: dict,
    directory: Union[pathlib.Path, str] = pathlib.Path(),
):
    """Write json dictionary with properties for surface wave data

    :param files: List of files holding surface wave  data
    :type files: List[Union[pathlib.Path, str]]
    :param tensor_info: The tensor information
    :type tensor_info: dict
    :param data_prop: The data properties (from sampling_filter.json)
    :type data_prop: dict
    :param directory: Where the file should be written, defaults to pathlib.Path()
    :type directory: Union[pathlib.Path, str], optional
    :return: The surface wave properties written to surf_waves.json
    :rtype: dict

    .. warning::

        Make sure the filters of strong motion data agree with the values in
        sampling_filter.json!
    """
    directory = pathlib.Path(directory)
    if len(files) == 0:
        return
    p_files = [file for file in files if "_BHZ." in file]
    sh_files = [file for file in files if "_SH." in file]
    p_files = select_tele_stations(p_files, "Rayleigh", tensor_info)
    sh_files = select_tele_stations(sh_files, "Love", tensor_info)
    files = p_files + sh_files
    origin_time = tensor_info["date_origin"]
    n0, n1 = data_prop["wavelet_scales"]
    surf_filter = data_prop["surf_filter"]
    wavelet_weight = wavelets_surf_tele(surf_filter, n0, n1)
    info_traces = []
    event_lat = tensor_info["lat"]
    event_lon = tensor_info["lon"]
    for file in files:
        header = SACTrace.read(file)
        stream = read(file)
        npts = header.npts
        starttime = stream[0].stats.starttime
        begin = origin_time - starttime
        n_start = int(begin / 4.0)
        if n_start < 0:
            continue
        if npts < n_start:
            continue
        if npts + n_start < 0:
            continue
        length = 900
        if 900 >= (npts - n_start):
            length = npts - n_start if n_start > 0 else npts
        distance, azimuth, back_azimuth = mng._distazbaz(
            header.stla, header.stlo, event_lat, event_lon
        )
        info = _dict_trace(
            file=file,
            name_station=header.kstnm,
            channel=header.kcmpnm,
            azimuth=azimuth,
            distance=distance / 111.11,
            dt=4.0,
            duration=length,
            n_start_obs=n_start,
            trace_weight=1.0,
            wavelet_weight=wavelet_weight,
            synthetic_trace=[],
            location=[header.stla, header.stlo],
        )
        info_traces.append(info)
    with open(directory / "surf_waves.json", "w") as f:
        json.dump(
            info_traces,
            f,
            sort_keys=True,
            indent=4,
            separators=(",", ": "),
            ensure_ascii=False,
        )
    return info_traces


def strong_motion_traces(
    files: List[Union[pathlib.Path, str]],
    tensor_info: dict,
    data_prop: dict,
    directory: Union[pathlib.Path, str] = pathlib.Path(),
) -> dict:
    """Write json dictionary with properties for strong motion data

    :param files: List of files holding strong motion data
    :type files: List[Union[pathlib.Path, str]]
    :param tensor_info: The tensor information
    :type tensor_info: dict
    :param data_prop: The data properties (from sampling_filter.json)
    :type data_prop: dict
    :param directory: Where the file should be written, defaults to pathlib.Path()
    :type directory: Union[pathlib.Path, str], optional
    :return: The strong motion properties written to strong_motion_waves.json
    :rtype: dict

    .. warning::

        Make sure the filters of strong_motion data agree with the values in
        sampling_filter.json!
    """
    directory = pathlib.Path(directory)
    if len(files) == 0:
        return
    event_lat = tensor_info["lat"]
    event_lon = tensor_info["lon"]
    depth = tensor_info["depth"]
    origin_time = tensor_info["date_origin"]
    headers = [SACTrace.read(file) for file in files]
    dt_strong = headers[0].delta
    dt_strong = round(dt_strong, 2)
    fun1 = lambda header: mng._distazbaz(header.stla, header.stlo, event_lat, event_lon)
    values = map(fun1, headers)
    distances = [value[0] for value in values]
    zipped = zip(distances, headers)
    arrivals = [np.sqrt(dist**2 + depth**2) / 5 for dist in distances]
    filter0 = data_prop["strong_filter"]
    duration = duration_strong_motion(distances, arrivals, tensor_info, dt_strong)
    n0, n1 = data_prop["wavelet_scales"]
    wavelet_weight = wavelets_strong_motion(duration, filter0, dt_strong, n0, n1)
    black_list = {"PB02": ["HNE", "HNN", "HLE", "HLN"], "PX02": ["HNE", "HLE"]}

    info_traces = []
    outlier_traces = []
    streams = [read(file) for file in files]
    weights = [1.0 for st, file in zip(streams, files)]
    fun2 = (
        lambda header: header.kstnm in black_list
        and header.kcmpnm in black_list[header.kstnm]
    )
    zipped = zip(weights, headers)
    weights = (0 if fun2(header) else weight for weight, header in zipped)
    streams = (st for st in streams)
    headers = (header for header in headers)

    zipped = zip(files, headers, weights, streams, arrivals)

    for file, header, weight, stream, arrival in zipped:
        start = origin_time - stream[0].stats.starttime
        distance, azimuth, back_azimuth = mng._distazbaz(
            header.stla, header.stlo, event_lat, event_lon
        )
        info = _dict_trace(
            file=file,
            name_station=header.kstnm,
            channel=header.kcmpnm,
            azimuth=azimuth,
            distance=distance / 111.11,
            dt=dt_strong,
            duration=duration,
            n_start_obs=int(start / dt_strong),
            trace_weight=float(weight),
            wavelet_weight=wavelet_weight,
            synthetic_trace=[],
            location=[header.stla, header.stlo],
        )
        info_traces.append(info)
    with open(directory / "strong_motion_waves.json", "w") as f:
        json.dump(
            info_traces,
            f,
            sort_keys=True,
            indent=4,
            separators=(",", ": "),
            ensure_ascii=False,
        )
    with open(directory / "outlier_strong_motion_waves.json", "w") as f:
        json.dump(
            outlier_traces,
            f,
            sort_keys=True,
            indent=4,
            separators=(",", ": "),
            ensure_ascii=False,
        )
    return info_traces


def cgps_traces(
    files: List[Union[pathlib.Path, str]],
    tensor_info: dict,
    data_prop: dict,
    directory: Union[pathlib.Path, str] = pathlib.Path(),
) -> dict:
    """Write json dictionary with properties for cGPS data

    :param files: List of files holding cGPS data
    :type files: List[Union[pathlib.Path, str]]
    :param tensor_info: The tensor information
    :type tensor_info: dict
    :param data_prop: The data properties (from sampling_filter.json)
    :type data_prop: dict
    :param directory: Where the file should be written, defaults to pathlib.Path()
    :type directory: Union[pathlib.Path, str], optional
    :return: The cGPS properties written to cgps_waves.json
    :rtype: dict

    .. warning::

        Make sure the filters of cGPS data agree with the values in
        sampling_filter.json!
    """
    directory = pathlib.Path(directory)
    if len(files) == 0:
        return
    event_lat = tensor_info["lat"]
    event_lon = tensor_info["lon"]
    depth = tensor_info["depth"]
    origin_time = tensor_info["date_origin"]
    headers = [SACTrace.read(file) for file in files]
    dt_cgps = headers[0].delta
    dt_cgps = round(dt_cgps, 2)
    fun1 = lambda header: mng._distazbaz(header.stla, header.stlo, event_lat, event_lon)
    values = map(fun1, headers)
    distances = [value[0] for value in values]
    zipped = zip(distances, headers)
    arrivals = [np.sqrt(dist**2 + depth**2) / 5 for dist in distances]
    duration = duration_strong_motion(distances, arrivals, tensor_info, dt_cgps)
    filter0 = data_prop["strong_filter"]
    if "cgps_filter" in data_prop:
        filter0 = data_prop["cgps_filter"]
    n0, n1 = data_prop["wavelet_scales"]
    wavelet_weight = wavelets_strong_motion(
        duration, filter0, dt_cgps, n0, n1, cgps=True
    )
    info_traces = []
    vertical = ["LXZ", "LHZ", "LYZ"]
    headers = (header for header in headers)
    streams = [read(file) for file in files]
    channels = (st[0].stats.channel for st in streams)
    is_horiz = lambda channel: channel not in vertical
    zipped = zip(files, channels)
    weights = (1.2 if is_horiz(channel) else 0.6 for file, channel in zipped)
    streams = (st for st in streams)

    for file, header, stream, weight in zip(files, headers, streams, weights):
        start = origin_time - stream[0].stats.starttime
        distance, azimuth, back_azimuth = mng._distazbaz(
            header.stla, header.stlo, event_lat, event_lon
        )
        info = _dict_trace(
            file=file,
            name_station=header.kstnm,
            channel=header.kcmpnm,
            azimuth=azimuth,
            distance=distance / 111.11,
            dt=dt_cgps,
            duration=duration,
            n_start_obs=int(start // dt_cgps),
            trace_weight=weight,
            wavelet_weight=wavelet_weight,
            synthetic_trace=[],
            location=[header.stla, header.stlo],
        )
        info_traces.append(info)
    with open(directory / "cgps_waves.json", "w") as f:
        json.dump(
            info_traces,
            f,
            sort_keys=True,
            indent=4,
            separators=(",", ": "),
            ensure_ascii=False,
        )
    return info_traces


def static_data(
    tensor_info: dict,
    unit: Literal["cm", "m", "mm"] = "m",
    directory: Union[pathlib.Path, str] = pathlib.Path(),
) -> dict:
    """Write json dictionary with properties for gps data

    :param tensor_info: The tensor information
    :type tensor_info: dict
    :param unit: The unit, defaults to "m"
    :type unit: Literal[&#39;cm&#39;,&#39;m&#39;,&#39;mm&#39;], optional
    :param directory: Where the file should be written, defaults to pathlib.Path()
    :type directory: Union[pathlib.Path, str], optional
    :return: The gps properties written to static_data.json
    :rtype: dict
    """
    directory = pathlib.Path(directory)
    event_lat = tensor_info["lat"]
    event_lon = tensor_info["lon"]
    info_traces = []
    if not os.path.isfile(directory / "cgps_waves.json") and not os.path.isfile(
        directory / "gps_data"
    ):
        return

    if os.path.isfile(directory / "cgps_waves.json") and not os.path.isfile(
        directory / "gps_data"
    ):
        with open(directory / "cgps_waves.json") as cgps_f:
            cgps_data = json.load(cgps_f)

        names = [file["name"] for file in cgps_data]
        names = list(set(names))

        for name in names:
            observed = [0, 0, 0]
            weights = [0, 0, 0]
            error = [0, 0, 0]
            for file in cgps_data:
                name2 = file["name"]
                if not name2 == name:
                    continue
                comp = file["component"]
                lat, lon = file["location"]
                duration = file["duration"]
                stream = read(file["file"])
                trace = stream[0].data
                # if len(trace) < 200:
                # continue
                baseline0 = np.mean(trace[:20])
                baseline1 = np.mean(trace[-70:-20])
                # baseline1 = np.mean(trace[-100:])
                baseline = baseline1 - baseline0
                variance = np.std(trace[-70:-20])
                # variance = np.std(trace[-100:])
                weight = 1.0 if not comp[-1] == "Z" else 0.5
                index = 0 if comp[-1] == "Z" else 1 if comp[-1] == "N" else 2
                observed[index] = str(baseline)
                weights[index] = str(weight)
                error[index] = str(variance)
            info = _dict_trace(
                file=None,
                name_station=name,
                channel=None,
                azimuth=None,
                distance=None,
                dt=None,
                duration=1,
                n_start_obs=10,
                trace_weight=weights,
                wavelet_weight=None,
                synthetic_trace=None,
                observed_trace=observed,
                location=[lat, lon],
            )
            info["data_error"] = error
            info_traces.append(info)

    if os.path.isfile(directory / "gps_data"):
        if unit == "cm":
            factor = 1
        elif unit == "m":
            factor = 100
        elif unit == "mm":
            factor = 1 / 10
        with open(directory / "gps_data", "r") as infile:
            lines = [line.split() for line in infile]

        for line in lines[2:]:
            if len(line) < 9:
                continue
            name = line[0]
            lon = float(line[1])
            lat = float(line[2])
            observed = [0, 0, 0]
            observed[2] = float(line[3]) * factor
            observed[1] = float(line[4]) * factor
            observed[0] = float(line[5]) * factor
            weights = [0, 0, 0]
            error = [0, 0, 0]
            variance = float(line[6]) * factor
            error[2] = variance
            weights[2] = 1.0
            variance = float(line[7]) * factor
            error[1] = variance
            weights[1] = 1.0
            variance = float(line[8]) * factor
            error[0] = variance
            weights[0] = 0.5
            observed = [str(obs) for obs in observed]
            weights = [str(w) for w in weights]
            distance, azimuth, back_azimuth = mng._distazbaz(
                lat, lon, event_lat, event_lon
            )
            info = _dict_trace(
                file=None,
                name_station=name,
                channel=["Z,N,E"],
                azimuth=azimuth,
                distance=distance / 111.11,
                dt=None,
                duration=1,
                n_start_obs=10,
                trace_weight=weights,
                wavelet_weight=None,
                synthetic_trace=None,
                observed_trace=observed,
                location=[lat, lon],
            )
            info["data_error"] = error
            info_traces.append(info)

    if len(info_traces) >= 250:
        observed_data = [trace["observed"] for trace in info_traces]
        observed_data = [[float(v) for v in obs] for obs in observed_data]
        pgd = [np.max(obs) for obs in observed_data]
        zipped = zip(info_traces, pgd)
        new_zipped = sorted(zipped, key=lambda val: val[1], reverse=True)
        info_traces = [a for a, b in new_zipped]
        info_traces = info_traces[:250]

    with open(directory / "static_data.json", "w") as f:
        json.dump(
            info_traces,
            f,
            sort_keys=True,
            indent=4,
            separators=(",", ": "),
            ensure_ascii=False,
        )
    return info_traces


def insar_data(
    insar_asc: Optional[List[Union[pathlib.Path, str]]] = None,
    insar_desc: Optional[List[Union[pathlib.Path, str]]] = None,
    ramp_asc: Optional[List[float]] = None,
    ramp_desc: Optional[List[float]] = None,
    directory: Union[pathlib.Path, str] = pathlib.Path(),
) -> dict:
    """Write json dictionary with properties for InSar data

    :param insar_asc: The paths to ascending InSar files, defaults to None
    :type insar_asc: Optional[List[Union[pathlib.Path, str]]], optional
    :param insar_desc: The paths to descending InSar files, defaults to None
    :type insar_desc: Optional[List[Union[pathlib.Path, str]]], optional
    :param ramp_asc: The ascending ramp values, defaults to None
    :type ramp_asc: Optional[List[float]], optional
    :param ramp_desc: The descending ramp values, defaults to None
    :type ramp_desc: Optional[List[float]], optional
    :param directory: _description_, defaults to pathlib.Path()
    :type directory: Union[pathlib.Path, str], optional
    :raises ValueError: If the length of ascending ramps and files are different
    :raises ValueError: If the length of descending ramps and files are different
    :return: The InSar properties written to cgps_waves.json
    :rtype: dict
    """
    directory = pathlib.Path(directory)
    print("InSAR data")
    print(f"insar_asc: {insar_asc}")
    print(f"insar_desc: {insar_desc}")
    if not insar_asc and not insar_desc:
        return

    insar_dict = {}
    if insar_asc:
        if not ramp_asc or len(ramp_asc) == 0 or ramp_asc == None:
            ramp_asc = [None] * len(insar_asc)
        if not len(ramp_asc) == len(insar_asc):
            raise ValueError(
                "You need to input the same amount of input to"
                + " both the amount of ascending tracks and their ramp-types"
            )
        zipped = zip(insar_asc, ramp_asc)
        properties = []
        for track, ramp in zipped:
            new_dict = {"name": track, "weight": 1.0, "ramp": ramp}
            properties = properties + [new_dict]
        insar_dict["ascending"] = properties

    if insar_desc:
        if not ramp_desc or len(ramp_desc) == 0 or ramp_desc == None:
            ramp_desc = [None] * len(insar_desc)
        if not len(ramp_desc) == len(insar_desc):
            raise ValueError(
                "You need to input the same amount of input to"
                + " both the amount of descending tracks and their ramp-types"
            )
        zipped = zip(insar_desc, ramp_desc)
        properties = []
        for track, ramp in zipped:
            new_dict = {"name": track, "weight": 1.0, "ramp": ramp}
            properties = properties + [new_dict]
        insar_dict["descending"] = properties

    with open(directory / "insar_data.json", "w") as f:
        json.dump(
            insar_dict,
            f,
            sort_keys=True,
            indent=4,
            separators=(",", ": "),
            ensure_ascii=False,
        )
    return insar_dict


def select_tele_stations(
    files: List[Union[pathlib.Path, str]],
    phase: Literal["Love", "P", "Rayleigh", "SH"],
    tensor_info: dict,
) -> List[Union[pathlib.Path, str]]:
    """Select body and/or surface waves to use in finite fault modelling.

    :param files: The list of files
    :type files: List[Union[pathlib.Path, str]]
    :param phase: The phase to choose
    :type phase: Literal[&quot;Love&quot;, &quot;P&quot;, &quot;Rayleigh&quot;, &quot;SH&quot;]
    :param tensor_info: The tensor information
    :type tensor_info: dict
    :return: The updated/filtered list
    :rtype: List[Union[pathlib.Path, str]]
    """
    event_lat = tensor_info["lat"]
    event_lon = tensor_info["lon"]
    weight = 1.0 if not phase == "SH" else 0.5
    if phase in ["P", "Rayleigh"]:
        # min_snr = 6.0 if phase == 'P' else 4.0
        min_snr = 5.0 if phase == "P" else 4.0
        window = 100 if phase == "P" else 1500
        jump = 2
        total = __used_stations(jump, files, tensor_info)
        if total > 50:
            jump = 4
            total = __used_stations(jump, files, tensor_info)
        if total > 50:
            jump = 8
    if phase in ["SH", "Love"]:
        min_snr = 2.5  # 3.0#4.0
        window = 200 if phase == "SH" else 1500
        jump = 8
        total = __used_stations(jump, files, tensor_info)
        if total > 30:
            jump = 10
            total = __used_stations(jump, files, tensor_info)
        if total > 30:
            jump = 12
    sacheaders = [SACTrace.read(sac) for sac in files]
    signal2noise = [__s2nr(sacfile, phase, window) for sacfile in files]
    min_signal2noise = len([ratio for ratio in signal2noise if ratio > min_snr])
    if min_signal2noise < 4:
        min_snr = min_snr / 2.0
    score = (
        lambda x, y, z, y0, y1, z0: x
        - (y - y0) * (y - y1) / ((y - y0) ** 2.0 + (y - y1) ** 2.0)
        + (z - z0) / (20 - z0)
    )

    # Limiting the number of data. Choose one station for every degree in
    # azimuth

    phase = "LONG" if phase in ["Rayleigh", "Love"] else phase
    new_files = []
    for az0 in range(0, 360, jump):
        az1 = az0 + jump
        best = ""
        best_sta = ""
        best_score = -1.0
        add_channel = False
        for sacheader, snr, sac in zip(sacheaders, signal2noise, files):
            dis, az, baz = mng._distazbaz(
                sacheader.stla, sacheader.stlo, event_lat, event_lon
            )
            if az0 > az or az >= az1 or snr <= min_snr:
                continue
            value = 1 if kilometers2degrees(dis) >= 45 else 0
            value = value if phase in ["P", "SH"] else 1.0
            value = score(value, az, snr, az0, az1, min_snr)
            if value > best_score:
                add_channel = True
                best_score = value
                best_sta = sac
        new_files = new_files + [best_sta] if add_channel else new_files

    return new_files


def __s2nr(sacfile: Union[pathlib.Path, str], phase: str, signal_length: int) -> float:
    """Signal to noise ratio for data

    :param sacfile: The path to the sac file
    :type sacfile: Union[pathlib.Path,str]
    :param phase: The phase
    :type phase: str
    :param signal_length: The window of the signal
    :type signal_length: int
    :return: The SNR
    :rtype: float
    """
    sacheader = SACTrace.read(sacfile)
    stream = read(sacfile)
    trace = stream[0].data
    begin = sacheader.b
    dt = sacheader.delta
    arrival = sacheader.t1 if phase == "P" else sacheader.t5
    error_window = 5.0 if phase == "P" else 10.0
    if phase in ["Rayleigh", "Love"]:
        arrival = sacheader.t1
        error_window = 100.0
    arrival = int((arrival - begin) / dt)
    length = int(signal_length / dt + 0.1)
    error_window = int(error_window / dt)
    length = min(length, arrival - error_window)
    begin_signal = arrival + error_window - 1
    signal = trace[begin_signal : begin_signal + length]
    begin_noise = arrival - error_window - length + 1
    noise = trace[begin_noise : begin_noise + length]
    signal_std = np.std(signal)
    noise_std = np.std(noise)
    return min(20, signal_std / noise_std)


def __used_stations(
    jump: int, sacfiles: Union[pathlib.Path, str], tensor_info: dict
) -> int:
    """TODO: Define what this does and what jump is"""
    sacheaders = (SACTrace.read(sac) for sac in sacfiles)
    event_lat = tensor_info["lat"]
    event_lon = tensor_info["lon"]
    fun1 = lambda header: mng._distazbaz(header.stla, header.stlo, event_lat, event_lon)
    values = map(fun1, sacheaders)
    azimuths = [az for dis, az, baz in values]
    total = 0
    for az0 in range(0, 360, jump):
        list_ = [az for az in azimuths if az0 <= az < az0 + jump]
        if list_:
            total = total + 1
    return total


def __failsafe(filtro: dict, header: SACTrace, cgps: bool = False):
    """Check that the selected filter matches the filter applied to the data

    :param filtro: The filter properties
    :type filtro: dict
    :param header: The sac class
    :type header: SACTrace
    :param cgps: Whether it is cgps data, defaults to False
    :type cgps: bool, optional
    :raises ValueError: If the filter selected does not match that applied to the data
    """
    low_freq = filtro["low_freq"]
    high_freq = filtro["high_freq"]
    low_freq2 = header.t8
    high_freq2 = header.t9
    if not cgps:
        if not abs(low_freq - low_freq2) + abs(high_freq - high_freq2) < 0.0001:
            print([low_freq, high_freq])
            print([low_freq2, high_freq2])
            raise ValueError("Selected filter doesn't match filter applied to data")
    else:
        if not abs(high_freq - high_freq2) < 0.0001:
            print([high_freq, high_freq2])
            raise ValueError("Selected filter doesn't match filter applied to data")


def filling_data_dicts(
    tensor_info,
    data_type,
    data_prop,
    data_folder,
    insar_asc=None,
    insar_desc=None,
    ramp_asc=None,
    ramp_desc=None,
):
    """Routine to fill JSON dictionaries containing data properties, for all
    data types selected.

    :param tensor_info: dictionary with moment tensor information
    :param data_type: list with data types to use in modelling.
    :param data_prop: dictionary with moment tensor information
    :param data_folder: string with folder where data is located.
    :param insar_asc: name of ascending insar track
    :param insar_desc: name of descending insar track
    :param ramp_asc: type of ramp to invert for ascending track
    :param ramp_desc: type of ramp to invert for descending track
    :type insar_asc: string, optional
    :type insar_desc: string, optional
    :type ramp_asc: string, optional
    :type ramp_desc: string, optional
    :type tensor_info: dict
    :type data_type: list
    :type data_prop: dict
    :type data_folder: string

    .. rubric:: Example:

    >>> data_prop = {
            "sampling": {
                "dt_strong": 0.4,
                "dt_tele": 0.2
            },
            "strong_filter": {
                "high_freq": 0.125,
                "low_freq": 0.01
            },
            "tele_filter": {
                "freq0": 0.002,
                "freq3": 1.2,
                "high_freq": 1.0,
                "low_freq": 0.004
            },
            "wavelet_scales": [
                1,
                8
            ]
        }
    >>> tensor_info = {
            'moment_mag': 7 * 10 ** 27,
            'date_origin': UTCDateTime(2014, 04, 01, 23, 46, 47)
            'lat': -19.5,
            'lon': -70.5,
            'depth': 25,
            'time_shift': 44,
            'half_duration': 40,
            'centroid_lat': -21,
            'centroid_lon': -70,
            'centroid_depth': 35
        }
    >>> data_folder = '/path/to/data_folder'
    >>> filling_data_dicts(tensor_info, data_type, data_prop, data_folder)
    """
    folder = os.path.abspath(os.getcwd())
    if "tele_body" in data_type:
        os.chdir(data_folder)
        tele_traces = get_traces_files("tele_body")
        os.chdir(folder)
        if not os.path.isfile("tele_waves.json"):
            tele_body_traces(tele_traces, tensor_info, data_prop)
    if "surf_tele" in data_type:
        os.chdir(data_folder)
        surf_traces = get_traces_files("surf_tele")
        os.chdir(folder)
        if not os.path.isfile("surf_waves.json"):
            tele_surf_traces(surf_traces, tensor_info, data_prop)
    if "strong_motion" in data_type:
        os.chdir(data_folder)
        strong_traces = get_traces_files("strong_motion")
        os.chdir(folder)
        if not os.path.isfile("strong_motion_waves.json"):
            strong_motion_traces(strong_traces, tensor_info, data_prop)
    if "cgps" in data_type:
        os.chdir(data_folder)
        cgps_data = get_traces_files("cgps")
        os.chdir(folder)
        if not os.path.isfile("cgps_waves.json"):
            cgps_traces(cgps_data, tensor_info, data_prop)
    if "gps" in data_type:
        static_data(tensor_info, unit="m")
    if "insar" in data_type:
        insar_data(
            insar_asc=insar_asc,
            insar_desc=insar_desc,
            ramp_asc=ramp_asc,
            ramp_desc=ramp_desc,
        )


def get_traces_files(
    data_type: Literal["cgps", "strong_motion", "surf_tele", "tele_body"],
    directory: Union[pathlib.Path, str] = pathlib.Path(),
) -> List[str]:
    """Get list with waveform files (in sac format) for stations and
    channels selected for modelling

    :param data_type: The type of data
    :type data_type: Literal[&quot;cgps&quot;, &quot;strong_motion&quot;, &quot;surf_tele&quot;, &quot;tele_body&quot;]
    :param directory: The directory where files are located, defaults to pathlib.Path()
    :type directory: Union[pathlib.Path, str], optional
    :return: The list of sac files
    :rtype: List[str]
    """
    directory = pathlib.Path(directory)
    if data_type == "tele_body":
        p_traces_files = glob.glob(os.path.join(directory / "P", "final*"))
        sh_traces_files = glob.glob(os.path.join(directory / "SH", "final*"))
        traces_files = p_traces_files + sh_traces_files
    if data_type == "surf_tele":
        traces_files = glob.glob(os.path.join(directory / "LONG", "final*"))
    if data_type == "strong_motion":
        traces_files = glob.glob(os.path.join(directory / "STR", "final*"))
    if data_type == "cgps":
        traces_files = glob.glob(os.path.join(directory / "cGPS", "final*"))
    traces_files = [os.path.abspath(file) for file in traces_files]
    return traces_files


def wavelets_body_waves(duration, filtro_tele, dt_tele, n_begin, n_end):
    """Automatic determination of weight of wavelet scales

    :param duration: length of signal to be considered for modelling
    :param filtro_tele: filtering properties for body waves
    :param dt_tele: sampling interval for body waves
    :param n_begin: minimum wavelet scale
    :param n_end: maximum wavelet scale
    :type duration: float
    :type filtro_tele: float
    :type dt_tele: float
    :type n_begin: int
    :type n_end: int
    """
    low_freq = filtro_tele["low_freq"]
    min_wavelet = int(np.log2(3 * 2**8 * dt_tele * low_freq)) + 1
    min_index = int(duration / dt_tele)
    min_wavelet = max(min_wavelet, 10 - int(np.log2(min_index)))
    p_wavelet_weight = ["0"] + ["1"] * 7
    p_wavelet_weight[:min_wavelet] = ["0"] * min_wavelet
    p_wavelet_weight = " ".join(p_wavelet_weight[n_begin - 1 : n_end])
    p_wavelet_weight = "".join([p_wavelet_weight, "\n"])
    #
    # the maximum frequency to be modelled for SH waves is 0.5 Hz
    #
    s_wavelet_weight = ["0"] + ["1"] * 6 + ["0"] if dt_tele <= 0.2 else ["1"] * 8
    s_wavelet_weight[:min_wavelet] = ["0"] * min_wavelet
    s_wavelet_weight = " ".join(s_wavelet_weight[n_begin - 1 : n_end])
    s_wavelet_weight = "".join([s_wavelet_weight, "\n"])
    return p_wavelet_weight, s_wavelet_weight


def wavelets_surf_tele(surf_filter, n_begin, n_end):
    """Automatic determination of weight of wavelet scales

    :param surf_filter: filtering properties of surface wave data
    :param n_begin: minimum wavelet scale
    :param n_end: maximum wavelet scale
    :type surf_filter: dict
    :type n_begin: int
    :type n_end: int
    """
    low_freq = surf_filter["freq2"]
    high_freq = surf_filter["freq3"]
    min_wavelet = int(np.log2(3 * 2**8 * 4.0 * low_freq)) + 1
    max_wavelet = max(int(np.log2(3 * 2**10 * 4.0 * high_freq)), 1)
    #
    # largest frequency we can model with wavelets
    #
    wavelet_weight = ["2"] * 10
    wavelet_weight[:min_wavelet] = ["0"] * min_wavelet
    wavelet_weight = wavelet_weight[:n_end]
    wavelet_weight = wavelet_weight[:max_wavelet] + ["0"] * (n_end - max_wavelet)
    wavelet_weight = " ".join(wavelet_weight[n_begin - 1 : n_end])
    wavelet_weight = "".join([wavelet_weight, "\n"])
    return wavelet_weight


def __wavelets_dart(n_begin, n_end):
    """Automatic determination of weight of wavelet scales

    :param n_begin: minimum wavelet scale
    :param n_end: maximum wavelet scale
    :type n_begin: int
    :type n_end: int
    """
    wavelet_weight = ["0"] * 2 + ["2"] * 6
    wavelet_weight = " ".join(wavelet_weight[n_begin - 1 : n_end])
    wavelet_weight = "".join([wavelet_weight, "\n"])
    return wavelet_weight


def wavelets_strong_motion(
    duration, filtro_strong, dt_strong, n_begin, n_end, cgps=False
):
    """Automatic determination of weight of wavelet scales

    :param duration: length of signal to be considered for modelling
    :param filtro_strong: filtering properties for strong motion data
    :param dt_strong: sampling interval for strong motion data
    :param n_begin: minimum wavelet scale
    :param n_end: maximum wavelet scale
    :param cgps: whether wavelet coefficients are for cGPS or strong motion data
    :type duration: float
    :type filtro_tele: float
    :type dt_tele: float
    :type n_begin: int
    :type n_end: int
    :type cgps: bool, optional
    """
    low_freq = filtro_strong["low_freq"]
    high_freq = filtro_strong["high_freq"]
    min_wavelet = int(np.log2(3 * 2**8 * dt_strong * low_freq)) + 1
    if cgps:
        min_wavelet = 3  # 4
    min_index = int(duration / dt_strong)
    min_wavelet = max(min_wavelet, 10 - int(np.log2(min_index)))
    max_wavelet = max(int(np.log2(3 * 2**10 * dt_strong * high_freq)), 1)
    #
    # largest frequency we can model with wavelets
    #
    wavelet_weight = ["2"] * 10
    wavelet_weight[:min_wavelet] = ["0"] * min_wavelet
    wavelet_weight = wavelet_weight[:n_end]
    wavelet_weight = wavelet_weight[:max_wavelet] + ["0"] * (n_end - max_wavelet)
    wavelet_weight = " ".join(wavelet_weight[n_begin - 1 : n_end])
    wavelet_weight = "".join([wavelet_weight, "\n"])
    return wavelet_weight


def duration_tele_waves(tensor_info, dt_tele):
    """Source duration plus depth phases

    :param tensor_info: dictionary with properties of strong motion data
    :param dt_tele: sampling interval for body waves
    :type tensor_info: dict
    :type dt_tele: float
    """
    depth = tensor_info["depth"]
    time_shift = tensor_info["time_shift"]
    if depth < 33.0:
        rdur = 2 * (time_shift + depth / 3.5)
    else:
        rdur = 2 * (time_shift + 33.0 / 3.5 + (depth - 33) / 5.0)
    max_length = 180 if depth < 300.0 else 300
    tele_syn_len = min(int(rdur * 1.5 / dt_tele), int(max_length / dt_tele))
    tele_syn_len = max(int(60 / dt_tele), tele_syn_len)
    return tele_syn_len


def duration_strong_motion(distances, arrivals, tensor_info, dt_strong):
    """Maximum duration estimation for strong motion records

    :param distances: distances from stations to hypocenter
    :param arrivals: estimation of first arrival from hypocenter to station
    :param dt_strong: sampling interval for strong motion data
    :param tensor_info: dictionary with moment tensor properties
    :type distances: array
    :type arrivals: array
    :type tensor_info: dict
    :type dt_strong: float
    """
    depth = float(tensor_info["depth"])
    time_shift = float(tensor_info["time_shift"])
    max_dist = max(distances)
    max_len = 350
    if dt_strong == 0.2:
        max_len = 200
    elif dt_strong == 0.4:
        max_len = 350
    elif dt_strong >= 0.8:
        max_len = 600
    max_arrival = max([arrival for arrival in arrivals])
    syn_len = 4 * time_shift + 7 * depth / 50
    syn_len = min(int((min(syn_len, max_len) + 2 * max_arrival) / dt_strong), 950)
    return syn_len


def duration_dart(distances, arrivals, tensor_info, dt_dart):
    """Maximum duration estimation for DART records

    :param distances: distances from stations to hypocenter
    :param arrivals: estimation of first arrival from hypocenter to station
    :param dt_dart: sampling interval for dart data
    :param tensor_info: dictionary with moment tensor properties
    :type distances: array
    :type arrivals: array
    :type tensor_info: dict
    :type dt_dart: float
    """
    depth = float(tensor_info["depth"])
    time_shift = float(tensor_info["time_shift"])
    max_dist = max(distances)
    max_arrival = 0
    max_len = 300
    syn_len = 4 * time_shift + 15 * max_dist / 111.11 + 7 * depth / 50
    syn_len = min(int((min(syn_len, max_len) + max_arrival) / dt_dart), 950)
    return syn_len


def __is_number(value: Union[str, float]) -> bool:
    """Determine if the provided value is a number

    :param value: The value
    :type value: str
    :return: Whether the provided value is a number
    :rtype: bool
    """
    try:
        float(value)
        return True
    except ValueError:
        return False


if __name__ == "__main__":
    import wasp.manage_parser as mp
    import wasp.seismic_tensor as tensor

    parser = argparse.ArgumentParser()
    parser.add_argument(
        "-f", "--folder", default=os.getcwd(), help="folder where there are input files"
    )
    parser.add_argument(
        "-d", "--data_folder", default=os.getcwd(), help="folder with waveform data"
    )
    parser = mp.parser_add_tensor(parser)
    parser = mp.parser_data_dict(parser)
    args = parser.parse_args()
    data_folder = os.path.abspath(args.data_folder)
    if not os.path.isfile("sampling_filter.json"):
        raise FileNotFoundError(
            errno.ENOENT, os.strerror(errno.ENOENT), "sampling_filter.json"
        )
    data_prop = json.load(open("sampling_filter.json"))
    os.chdir(args.folder)
    if args.gcmt_tensor:
        cmt_file = args.gcmt_tensor
        tensor_info = tensor.get_tensor(cmt_file=cmt_file)
    else:
        tensor_info = tensor.get_tensor()
    data_type = mp.get_used_data(args)
    filling_data_dicts(
        tensor_info,
        data_type,
        data_prop,
        data_folder,
        insar_asc=args.insar_asc,
        insar_desc=args.insar_desc,
        ramp_asc=args.insar_asc_ramp,
        ramp_desc=args.insar_desc_ramp,
    )
