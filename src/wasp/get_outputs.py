#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""Script with routines for retrieving both the kinematic model which solves
the inverse problem, and the synthetic waveforms produced by such model,
"""


import errno
import json
import os
import pathlib
from typing import Dict, List, Optional, Tuple, Union

import numpy as np
from obspy import read  # type: ignore

import wasp.plane_management as pl_mng

##########################
# Get FFM model
##########################


def read_solution_static_format(
    segments: List[Dict[str, str]], data_dir: Union[str, pathlib.Path] = pathlib.Path()
) -> dict:
    """Read the solution file in static format

    :param segments: List of dictionaries with properties of fault segments
    :type segments: List[Dict[str, str]]
    :param data_dir: The path to the data directory, defaults to pathlib.Path()
    :type data_dir: Union[str, pathlib.Path], optional
    :raises FileNotFoundError: If Solution.txt cannot be found
    :raises RuntimeError: If the length of segments does not match the number of
                            of fault segments
    :raises RuntimeError: If the number of subfaults does not match
                            (in segments and in the solution file)
    :return: The parsed information from Solucion.txt
    :rtype: dict
    """
    lat: List[float] = []
    lon: List[float] = []
    depth: List[float] = []
    slip: List[float] = []
    rake: List[float] = []
    trup: List[float] = []
    trise: List[float] = []
    tfall: List[float] = []
    moment: List[float] = []

    data_dir = pathlib.Path(data_dir)
    if not os.path.isfile(data_dir / "Solucion.txt"):
        raise FileNotFoundError(
            errno.ENOENT, os.strerror(errno.ENOENT), str(data_dir / "Solucion.txt")
        )
    with open(data_dir / "Solucion.txt", "r") as input_file:
        jk = [line.split() for line in input_file]

    faults_data = [
        index + 1
        for index, line in enumerate(jk)
        if set(["#Lat.", "Lon.", "depth", "slip"]) <= set(line)
    ]
    headers = [
        index
        for index, line in enumerate(jk)
        if set(["#Fault_segment", "nx(Along-strike)=", "ny(downdip)="]) <= set(line)
    ]
    headers = headers[1:] + [len(jk)]
    if not len(headers) == len(segments):
        raise RuntimeError(
            "Inconsistency between Fault.time and Solucion.txt."
            " Different amount of fault segments"
        )
    for segment, start, end in zip(segments, faults_data, headers):
        #
        # load FFM solution model and coordinates
        #
        lat_fault = np.array([float(line[0]) for line in jk[start:end]])
        lon_fault = np.array([float(line[1]) for line in jk[start:end]])
        depth_fault = np.array([float(line[2]) for line in jk[start:end]])
        slip_fault = np.array([float(line[3]) for line in jk[start:end]])
        rake_fault = np.array([float(line[4]) for line in jk[start:end]])
        trup_fault = np.array([float(line[7]) for line in jk[start:end]])
        trise_fault = np.array([float(line[8]) for line in jk[start:end]])
        tfall_fault = np.array([float(line[9]) for line in jk[start:end]])
        moment_fault = np.array([float(line[-1]) for line in jk[start:end]])
        (
            stk_subfaults,
            dip_subfaults,
            delta_strike,
            delta_dip,
            hyp_stk,
            hyp_dip,
        ) = pl_mng.__unpack_plane_data(segment)
        #
        # Reshape the rupture process
        #
        if not slip_fault.size == stk_subfaults * dip_subfaults:
            raise RuntimeError(
                "Inconsistency between Fault.time and Solucion.txt."
                " Different size of fault segment"
            )
        lat_fault.shape = dip_subfaults, stk_subfaults  # type: ignore
        lon_fault.shape = dip_subfaults, stk_subfaults  # type: ignore
        depth_fault.shape = dip_subfaults, stk_subfaults  # type: ignore
        slip_fault.shape = dip_subfaults, stk_subfaults  # type: ignore
        rake_fault.shape = dip_subfaults, stk_subfaults  # type: ignore
        trup_fault.shape = dip_subfaults, stk_subfaults  # type: ignore
        trise_fault.shape = dip_subfaults, stk_subfaults  # type: ignore
        tfall_fault.shape = dip_subfaults, stk_subfaults  # type: ignore
        moment_fault.shape = dip_subfaults, stk_subfaults  # type: ignore
        lat = lat + [lat_fault]  # type: ignore
        lon = lon + [lon_fault]  # type: ignore
        depth = depth + [depth_fault]  # type: ignore
        slip = slip + [slip_fault]  # type: ignore
        rake = rake + [rake_fault]  # type: ignore
        trup = trup + [trup_fault]  # type: ignore
        trise = trise + [trise_fault]  # type: ignore
        tfall = tfall + [tfall_fault]  # type: ignore
        moment = moment + [moment_fault]  # type: ignore

    solution = {
        "slip": slip,
        "rake": rake,
        "rupture_time": trup,
        "trise": trise,
        "tfall": tfall,
        "lat": lat,
        "lon": lon,
        "depth": depth,
        "moment": moment,
    }
    return solution


def read_solution_fsp_format(
    file_name: Union[str, pathlib.Path] = pathlib.Path(),
    custom: bool = False,
) -> Tuple[dict, dict, dict]:
    """Read the solution file in fsp format

    :param file_name: The path to the solution file, defaults to pathlib.Path()
    :type file_name: Union[str, pathlib.Path], optional
    :param custom: If the data has a custom start point, defaults to False
    :type custom: bool, optional
    :raises FileNotFoundError: If the the solution file cannot be found
    :return: The parsed tensor info, solution, and subfault data
    :rtype: Tuple[dict, dict, dict]
    """
    lat: List[float] = []
    lon: List[float] = []
    depth: List[float] = []
    slip: List[float] = []
    rake: List[float] = []
    trup: List[float] = []
    trise: List[float] = []
    width: List[float] = []

    if not os.path.isfile(file_name):
        raise FileNotFoundError(errno.ENOENT, os.strerror(errno.ENOENT), file_name)
    with open(file_name, "r") as input_file:
        jk = [line.split() for line in input_file]

    tensor_info = {
        "lat": float(jk[5][5]),
        "lon": float(jk[5][8]),
        "depth": float(jk[5][11]),
    }
    n_segments = int(jk[14][8])
    subfaults_data = {
        "stk_subfaults": int(jk[12][5]),
        "dip_subfaults": int(jk[12][8]),
        "delta_strike": float(jk[13][5]),
        "delta_dip": float(jk[13][9]),
        "strike": float(jk[7][5]),
    }

    line0: List[int] = [
        i for i, line in enumerate(jk) if {"SOURCE", "MODEL", "PARAMETERS"} < set(line)
    ]
    line0 = line0[0] + 9  # type:ignore
    if n_segments == 1:
        for line in jk[line0:]:  # type:ignore
            lat0 = float(line[0])
            lon0 = float(line[1])
            depth0 = float(line[4])
            slip0 = float(line[5])
            rake0 = float(line[6]) if len(line) >= 7 else 0
            trup0 = float(line[7]) if len(line) >= 8 else 0
            trise0 = float(line[8]) if len(line) >= 9 else 0
            lat = lat + [lat0]
            lon = lon + [lon0]
            depth = depth + [depth0]
            slip = slip + [slip0]
            rake = rake + [rake0]
            trup = trup + [trup0]
            trise = trise + [trise0]
            width = width + [0]
    else:
        for i_segment in range(n_segments):
            width0 = 0 if not custom else float(jk[line0 + 2][7])  # type:ignore
            subfaults_seg: List[int] = [
                int(jk[line0 + 6][3]) if not custom else int(jk[line0 + 7][3])  # type: ignore
            ]
            line0 = line0 + 10 if not custom else line0 + 11  # type:ignore
            for line in jk[line0 : line0 + subfaults_seg]:  # type:ignore
                lat0 = float(line[0])
                lon0 = float(line[1])
                depth0 = float(line[4])
                slip0 = float(line[5])
                rake0 = float(line[6]) if len(line) >= 7 else 0
                trup0 = float(line[7]) if len(line) >= 8 else 0
                trise0 = float(line[8]) if len(line) >= 9 else 0
                lat = lat + [lat0]
                lon = lon + [lon0]
                depth = depth + [depth0]
                slip = slip + [slip0]
                rake = rake + [rake0]
                trup = trup + [trup0]
                trise = trise + [trise0]
                width = width + [width0]
            line0 = line0 + subfaults_seg + 1  # type:ignore

    solution = {
        "slip": slip,
        "rake": rake,
        "trup": trup,
        "trise": trise,
        "lat": lat,
        "lon": lon,
        "depth": depth,
        "new_width": width,
    }
    return tensor_info, solution, subfaults_data


##########################
# Get data
##########################


def get_data_dict(
    traces_info: List[dict],
    syn_file: Optional[Union[pathlib.Path, str]] = None,
    observed: bool = True,
    obs_file: Optional[Union[pathlib.Path, str]] = None,
    margin: int = 10,
    directory: Union[pathlib.Path, str] = pathlib.Path(),
) -> List[dict]:
    """Fill dictionary with synthetic data at station and channel

    :param traces_info: The properties of each stations/channel
    :type traces_info: List[dict]
    :param syn_file: The path to the synthetic file, defaults to None
    :type syn_file: Optional[Union[pathlib.Path, str]], optional
    :param observed: Whether or not to include observed, defaults to True
    :type observed: bool, optional
    :param obs_file: The path to the observed file, defaults to None
    :type obs_file: Optional[Union[pathlib.Path, str]], optional
    :param margin: The file margin, defaults to 10
    :type margin: int, optional
    :param directory: Where the files should be read from, defaults to pathlib.Path()
    :type directory: Union[pathlib.Path, str], optional
    :return: The updated traces_info dictionaries
    :rtype: List[dict]
    """
    directory = pathlib.Path(directory)
    used_channels: List[List[str]] = []
    if syn_file:
        with open(directory / syn_file, "r") as infile:
            lines = [line.split() for line in infile]

        lines0 = [i for i, line in enumerate(lines) if len(line) > 2]
        for file in traces_info:
            name = file["name"]
            channel = file["component"]
            channel = __get_channel(channel)
            indexes = [
                i for i in lines0 if lines[i][2] == name and lines[i][3] in channel
            ]
            this_channel = [name, channel]
            if not this_channel in used_channels:
                used_channels = used_channels + [this_channel]
                index = indexes[0]
            else:
                index = indexes[1]
            npts = int(lines[index][0])
            synthetic = [float(real) for real, imag in lines[index + 1 : index + npts]]
            file["synthetic"] = np.array(synthetic)
    if observed:
        if obs_file:
            for file in traces_info:
                file["observed"] = _get_observed_from_chen(file, obs_file)
        else:
            for file in traces_info:
                file["observed"] = _get_observed_from_chen2(file, margin=margin)
                derivative = False if not "derivative" in file else file["derivative"]
                dt = file["dt"]
                file["observed"] = (
                    np.gradient(file["observed"], dt, edge_order=2)
                    if derivative
                    else file["observed"]
                )

    else:
        for file in traces_info:
            file["observed"] = [0 for i in range(1024)]
    return traces_info


def _get_observed_from_chen(
    file: dict, obse_file: Union[pathlib.Path, str]
) -> np.ndarray:
    """Get the observed waveform from a observed waveform file

    :param file: The waveform file properties
    :type file: dict
    :param syn_file: The path to the synthetic file
    :type syn_file: Union[pathlib.Path, str]
    :return: The observed data
    :rtype: np.ndarray
    """
    name = file["name"]
    channel = file["component"]
    channel = __get_channel(channel)
    with open(obse_file, "r") as infile:
        lines = [line.split() for line in infile]

    lines0 = [i for i, line in enumerate(lines) if not __is_number(line[0])]
    indexes = [i for i in lines0 if (len(lines[i]) > 1 and lines[i][1] == name)]
    index = next(i for i in indexes if lines[i + 1][1] in channel)
    npts = int(lines[index + 3][1])
    npts = min(file["duration"], npts)
    observed = [float(real[0]) for real in lines[index + 6 : index + 5 + npts]]
    return np.array(observed)


def _get_observed_from_chen2(file: dict, margin: int = 10) -> np.ndarray:
    """Get the observed waveform from a observed waveform file

    :param file: The waveform file properties
    :type file: dict
    :param margin: The margin of the file, defaults to 10
    :type margin: int, optional
    :return: The observed waveform data
    :rtype: np.ndarray
    """
    file2 = file["file"]
    st = read(file2)
    data = st[0].data
    dt = file["dt"]
    index0 = file["start_signal"] - int(margin // dt)
    # index0 = max(index0, 0)
    index1 = file["start_signal"] + file["duration"]
    index1 = max(index1, index0 + 1)
    if index0 >= 0:
        data = data[index0:index1]
    else:
        data = data[index0:]
    return data


def __get_channel(channel: str) -> List[str]:
    """Determine the channels associated with the provided channel

    :param channel: The channel
    :type channel: str
    :return: Channels
    :rtype: List[str]
    """
    if channel in ["P", "BHZ"]:
        channels = ["P", "BHZ"]
    if channel in ["SH", None, "BH1", "BH2", "BHE", "BHN"]:
        channels = ["SH"]
    if channel in ["HNZ", "HLZ", "BNZ"]:
        channels = ["HNZ", "HLZ", "BNZ"]
    if channel in ["HNE", "HLE", "BNE"]:
        channels = ["HNE", "HLE", "BNE"]
    if channel in ["HNN", "HLN", "BNN"]:
        channels = ["HNN", "HLN", "BNN"]
    if channel in ["LXZ", "LHZ", "LYZ"]:
        channels = ["LXZ", "LHZ", "LYZ"]
    if channel in ["LXE", "LHE", "LYE"]:
        channels = ["LXE", "LHE", "LYE"]
    if channel in ["LXN", "LHN", "LYN"]:
        channels = ["LXN", "LHN", "LYN"]
    if channel in ["dart"]:
        channels = ["dart"]
    return channels


def retrieve_gps(
    syn_name="static_synthetics.txt",
    directory: Union[pathlib.Path, str] = pathlib.Path(),
) -> Tuple[
    List[str],
    List[float],
    List[float],
    List[List[float]],
    List[List[float]],
    List[List[float]],
]:
    """Get the inverted and observed GPS data and station location

    :param syn_name: The name of the synthetic data file, defaults to "static_synthetics.txt"
    :type syn_name: str, optional
    :param directory: The directory where files are located, defaults to pathlib.Path()
    :type directory: Union[pathlib.Path, str], optional
    :return: The station names, the station latitudes, the station longitudes,
                the observed data, the synthetic data, the error
    :rtype: Tuple[ List[str], List[float], List[float], List[List[float]], List[List[float]], List[List[float]], ]
    """
    directory = pathlib.Path(directory)
    with open(directory / "static_data.json") as obs_file:
        obs_data = json.load(obs_file)

    lats = [data["location"][0] for data in obs_data]
    lons = [data["location"][1] for data in obs_data]
    names = [data["name"] for data in obs_data]
    observed = [data["observed"] for data in obs_data]
    error = [data["data_error"] for data in obs_data]

    with open(directory / syn_name, "r") as inf:
        lines = [line.split() for line in inf]

    name_synthetic = [[line[1], line[4:]] for line in lines[1:]]
    synthetic: list = []
    for index, name in enumerate(names):
        synthetic_gps = [
            gps_syn for gps_name, gps_syn in name_synthetic if gps_name == name
        ]
        synthetic = synthetic + synthetic_gps
    return names, lats, lons, observed, synthetic, error


def get_insar(data_dir: Union[str, pathlib.Path] = pathlib.Path()) -> dict:
    """_summary_

    :param data_dir: The path to the data directory, defaults to pathlib.Path()
    :type data_dir: Union[str, pathlib.Path], optional
    :return: The insar data
    :rtype: dict
    """
    data_dir = pathlib.Path(data_dir)
    insar_file = data_dir / "insar_data.json"
    with open(insar_file, "r") as f:
        insar_data = json.load(f)
    insar_synthetics_file = data_dir / "insar_synthetics.txt"
    with open(insar_synthetics_file, "r") as syn_file:
        lines_syn = [line.split() for line in syn_file]

    lines_ramp: list = []
    ramps: list = []
    lines0 = 0
    lines1 = 1
    if "ascending" in insar_data:
        asc_properties = insar_data["ascending"]
        ramps = [asc_property["ramp"] for asc_property in asc_properties]
    if "descending" in insar_data:
        desc_properties = insar_data["descending"]
        ramps = ramps + [desc_property["ramp"] for desc_property in desc_properties]
    lines_ramp = []
    if any(ramps):
        with open(data_dir / "insar_ramp.txt", "r") as ramp_file:
            lines_ramp = [line.split() for line in ramp_file]
    if "ascending" in insar_data:
        asc_properties = insar_data["ascending"]
        for asc_property in asc_properties:
            insar_points: list = []
            insar_asc = asc_property["name"]
            ramp_asc = asc_property["ramp"]
            with open(insar_asc, "r") as asc_file:
                lines_asc = [line.split() for line in asc_file]
            lines_asc = [line for line in lines_asc if not "#" in "".join(line)]
            lines1 = lines0 + len(lines_asc)
            ramp_track = [[0] * 5] * len(lines_asc)
            if ramp_asc:
                ramp_track = lines_ramp[lines0 + 1 : lines1]
            zipped = zip(lines_asc[:], lines_syn[lines0 + 1 : lines1], ramp_track[:])
            for line1, line2, line3 in zipped:
                lat = float(line1[1])
                lon = float(line1[0])
                observed = 100 * float(line1[2])
                synthetic = float(line2[4])
                ramp = float(line3[4])
                new_dict = {
                    "lat": lat,
                    "lon": lon,
                    "observed": observed,
                    "synthetic": synthetic,
                    "ramp": ramp,
                }
                insar_points = insar_points + [new_dict]
            asc_property["points"] = insar_points
            lines0 = lines0 + len(lines_asc)
    if "descending" in insar_data:
        desc_properties = insar_data["descending"]
        for desc_property in desc_properties:
            insar_points = []
            insar_desc = desc_property["name"]
            ramp_desc = desc_property["ramp"]
            with open(insar_desc, "r") as desc_file:
                lines_desc = [line.split() for line in desc_file]
            lines_desc = [line for line in lines_desc if not "#" in "".join(line)]
            lines1 = lines0 + len(lines_desc)
            ramp_track = [[0] * 5] * len(lines_desc)
            if ramp_desc:
                ramp_track = lines_ramp[lines0 + 1 : lines1]
            zipped = zip(lines_desc[:], lines_syn[lines0 + 1 : lines1], ramp_track[:])
            for line1, line2, line3 in zipped:
                lat = float(line1[1])
                lon = float(line1[0])
                observed = 100 * float(line1[2])
                synthetic = float(line2[4])
                ramp = float(line3[4])
                new_dict = {
                    "lat": lat,
                    "lon": lon,
                    "observed": observed,
                    "synthetic": synthetic,
                    "ramp": ramp,
                }
                insar_points = insar_points + [new_dict]
            desc_property["points"] = insar_points
            lines0 = lines0 + len(lines_desc)
    return insar_data


def __is_number(value: str) -> bool:
    """Test if the provided value is a number

    :param value: The provided value string
    :type value: str
    :return: Whether or not the value is a number
    :rtype: bool
    """
    try:
        float(value)
        return True
    except ValueError:
        return False
