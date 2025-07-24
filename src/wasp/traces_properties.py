# -*- coding: utf-8 -*-
"""Module with automatic settings for sampling, length and frequency filtering
for teleseismic and near field data
"""


import json
import pathlib
from typing import Dict, List, Optional, Union

from obspy import UTCDateTime  # type: ignore

###################################
# automatic settings
###################################


def nyquist_frequency(dt: float) -> float:
    """Calculate the nyquist frequency 1/(2*dt)

    :param dt: Delta t
    :type dt: float
    :return: The nyquist frequency
    :rtype: float
    """
    return 1 / 2 / dt


def sampling(
    tensor_info: Dict[str, Union[UTCDateTime, float, str]],
    dt_cgps: Optional[float] = None,
) -> Dict[str, float]:
    """Set data sampling properties

    :param tensor_info: The tensor information
    :type tensor_info: Dict[str, Union[UTCDateTime, float, str]]
    :param dt_cgps: The dt of cgps, defaults to None
    :type dt_cgps: Optional[float], optional
    :return: The sampling properties
    :rtype: Dict[str, float]
    """
    moment_mag = float(tensor_info["moment_mag"])
    time_shift = float(tensor_info["time_shift"])
    depth = float(tensor_info["depth"])
    dt_tele = 0.2
    if moment_mag > 10**29:
        dt_tele = 0.4
    if 200 < depth < 400:
        dt_tele = 0.4
    elif depth >= 400:
        dt_tele = 0.5
    #
    # I had some issues modelling deep earthquakes with smaller sampling.
    #
    dt_strong = 0.2
    if time_shift <= 40:
        dt_strong = 0.2
    elif time_shift <= 80:
        dt_strong = 0.4
    else:
        dt_strong = 0.8
    dt_strong = dt_strong if depth < 500 else 0.8
    dt_cgps = dt_strong if not dt_cgps else dt_cgps
    return {"dt_tele": dt_tele, "dt_strong": dt_strong, "dt_cgps": dt_cgps}


def filter_tele(
    tensor_info: Dict[str, Union[UTCDateTime, float, str]],
) -> Dict[str, float]:
    """Set filter properties for teleseismic data

    :param tensor_info: The tensor information
    :type tensor_info: Dict[str, Union[UTCDateTime, float, str]]
    :return: The filter properties
    :rtype: Dict[str,float]
    """
    time_shift = float(tensor_info["time_shift"])
    depth = float(tensor_info["depth"])
    freq0 = 0.003
    freq2 = 1.0
    freq3 = 1.2
    if time_shift < 10:
        freq1 = 0.01
    elif time_shift < 30:
        freq1 = 0.006
    elif time_shift < 80:
        freq0 = 0.002
        freq1 = 0.004
    if time_shift >= 80 or depth >= 200:
        freq0 = 0.001
        freq1 = 0.002
        freq2 = 0.8
        freq3 = 0.9
    filter_info = {
        "freq0": freq0,
        "low_freq": freq1,
        "high_freq": freq2,
        "freq3": freq3,
    }
    return filter_info


def filter_surf() -> Dict[str, float]:
    """Sef filter properties for surface waves

    :return: The filter properties
    :rtype: Dict[str,float]
    """
    filter_info = {"freq1": 0.003, "freq2": 0.004, "freq3": 0.006, "freq4": 0.007}
    return filter_info


def filter_strong(
    tensor_info: Dict[str, Union[UTCDateTime, float, str]], cgps: bool = False
) -> Dict[str, float]:
    """Set filter properties for strong motion data

    :param tensor_info: The tensor information
    :type tensor_info: Dict[str, Union[UTCDateTime, float, str]]
    :param cgps: Whether cgps, defaults to False
    :type cgps: bool, optional
    :return: The filter properties
    :rtype: Dict[str, float]
    """
    time_shift = float(tensor_info["time_shift"])

    if time_shift <= 10:
        min_freq = 0.02
    elif time_shift < 25:
        min_freq = 0.01
    elif time_shift < 50:
        min_freq = 0.005
    elif time_shift < 100:
        min_freq = 0.002
    else:
        min_freq = 0.001
    min_freq = max(min_freq, 1 / 100)

    max_freq = 0.125  # if tensor_info['time_shift'] > 10 else 0.25
    if cgps:
        max_freq = 0.3
    max_freq = max_freq if float(tensor_info["depth"]) < 300 else 0.05
    filter_info = {"low_freq": min_freq, "high_freq": max_freq}
    return filter_info


def properties_json(
    tensor_info: Dict[str, Union[UTCDateTime, float, str]],
    dt_cgps: Optional[float],
    data_directory: Union[pathlib.Path, str] = pathlib.Path(),
) -> Dict[str, Union[List[int], Dict[str, float]]]:
    """Set sampling and filtering properties in sampling_filtering.json

    :param tensor_info: The tensor information
    :type tensor_info: Dict[str, Union[UTCDateTime, float, str]]
    :param dt_cgps: The dt for cgps data
    :type dt_cgps: Optional[float]
    :param data_directory: The data directory , defaults to pathlib.Path()
    :type data_directory: Union[pathlib.Path, str], optional
    :return: The sampling and filtering properties
    :rtype: Dict[str, Union[List[int], Dict[str, float]]]
    """
    data_directory = pathlib.Path(data_directory)
    dict1 = sampling(tensor_info, dt_cgps)
    dict2 = filter_tele(tensor_info)
    dict3 = filter_surf()
    dict4 = filter_strong(tensor_info)
    dict5 = filter_strong(tensor_info, cgps=True)
    scales = wavelet_scales()
    dict6: Dict[str, Union[List[int], Dict[str, float]]] = {
        "sampling": dict1,
        "tele_filter": dict2,
        "surf_filter": dict3,
        "strong_filter": dict4,
        "cgps_filter": dict5,
        "wavelet_scales": scales,
    }
    with open(data_directory / "sampling_filter.json", "w") as f:
        json.dump(
            dict6,
            f,
            sort_keys=True,
            indent=4,
            separators=(",", ": "),
            ensure_ascii=False,
        )
    return dict6


def wavelet_scales() -> List[int]:
    """Set scales for the wavelet transformation

    :return: The scales
    :rtype: List[int]
    """
    n_begin = 1
    n_end = 8
    return [n_begin, n_end]
