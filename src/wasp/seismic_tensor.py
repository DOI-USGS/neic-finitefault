# -*- coding: utf-8 -*-
"""
Module for retrieving data from moment tensor files in CMTSOLUTION format, and
retrieve the nodal planes from such a moment tensor.
"""


import json
import os
import pathlib
import xml.etree.ElementTree as ET
from copy import deepcopy
from datetime import datetime
from typing import Optional, Tuple, Union

import numpy as np
from obspy import read_events  # type:ignore
from obspy.core.utcdatetime import UTCDateTime  # type:ignore


def get_tensor(
    cmt_file: Optional[Union[pathlib.Path, str]] = None,
    quake_file: Optional[Union[pathlib.Path, str]] = None,
    directory: Union[pathlib.Path, str] = pathlib.Path(),
) -> dict:
    """Parse tensor information from a cmt file into a dictionary

    :param cmt_file: The path to the cmt file, defaults to None
    :type cmt_file: Optional[Union[pathlib.Path, str]], optional
    :param quake_file: The path to the quake file file, defaults to None
    :type quake_file: Optional[Union[pathlib.Path, str]], optional
    :param directory: Where the tensor file should be read/written to,
                        defaults to pathlib.Path()
    :type directory: Union[pathlib.Path, str], optional
    :raises RuntimeError: If both the cmt file and quake file are undefined
    :return: The parsed tensor
    :rtype: dict
    """
    directory = pathlib.Path(directory)
    if cmt_file:
        tensor_info = read_gcmt_file(cmt_file)
        tensor_info = modify_tensor(tensor_info)
        delta = UTCDateTime.utcnow() - tensor_info["date_origin"]
        tensor_info["timedelta"] = delta  # delta.total_seconds()
    if quake_file:
        tensor_info = read_quake_file(quake_file)
        tensor_info = modify_tensor(tensor_info)
        delta = UTCDateTime.utcnow() - tensor_info["date_origin"]
        tensor_info["timedelta"] = delta  # .total_seconds()
    if not cmt_file and not quake_file:
        if not os.path.isfile(directory / "tensor_info.json"):
            raise RuntimeError(
                "No file named tensor_info.json located in "
                "folder {}".format(os.getcwd())
            )
        with open(directory / "tensor_info.json") as t:
            tensor_info = json.load(t)
        tensor_info = modify_tensor(tensor_info)
    return tensor_info


def write_tensor(
    tensor_info: dict, directory: Union[pathlib.Path, str] = pathlib.Path()
):
    """Write moment tensor properties to tensor_info.json

    :param tensor_info: The moment tensor properties
    :type tensor_info: dict
    :param directory: Where to write the tensor file, defaults to pathlib.Path()
    :type directory: Union[pathlib.Path, str], optional
    """
    directory = pathlib.Path(directory)
    tensor_info2 = deepcopy(tensor_info)
    with open(directory / "tensor_info.json", "w") as f:
        try:
            del tensor_info2["date_origin"]
        except KeyError:
            pass
        json.dump(
            tensor_info2,
            f,
            sort_keys=True,
            indent=4,
            separators=(",", ": "),
            ensure_ascii=False,
        )
    return


def read_gcmt_file(cmt_file: Union[pathlib.Path, str]) -> dict:
    """Read a GCMT moment tensor and parse the information into a dictionary

    :param cmt_file: The full path to the cmt file
    :type cmt_file: Union[pathlib.Path, str]
    :return: The parsed moment tensor properties
    :rtype: dict
    """
    with open(cmt_file, "r") as input_file:
        lines = [line.split() for line in input_file]

    if "Q" in lines[0][0]:
        new_line = lines[0]
        code, year = lines[0][0].split("Q")
        lines[0] = [code, year] + new_line[1:]
    if "W" in lines[0][0]:
        new_line = lines[0]
        code, year = lines[0][0].split("W")
        lines[0] = [code, year] + new_line[1:]

    year, month, day, hour, minute, second, lat, lon, hyp_depth = lines[0][1:10]
    origin_time = UTCDateTime(
        int(year), int(month), int(day), int(hour), int(minute), int(float(second))
    )
    date_time = origin_time.isoformat()
    time_shift = max(float(lines[2][2]), 5)
    half_duration = float(lines[3][2])

    centroid_lat = float(lines[4][1])
    centroid_lon = float(lines[5][1])
    centroid_depth = float(lines[6][1])

    mrr = float(lines[7][1])
    mtt = float(lines[8][1])
    mpp = float(lines[9][1])
    mrt = float(lines[10][1])
    mrp = float(lines[11][1])
    mtp = float(lines[12][1])

    input_dict = {
        "mrr": mrr,
        "mtt": mtt,
        "mpp": mpp,
        "mrt": mrt,
        "mrp": mrp,
        "mtp": mtp,
        "datetime": date_time,
        "date_origin": origin_time,
        "lat": lat,
        "lon": lon,
        "depth": hyp_depth,
        "time_shift": time_shift,
        "half_duration": half_duration,
        "centroid_lat": centroid_lat,
        "centroid_lon": centroid_lon,
        "centroid_depth": centroid_depth,
    }

    return input_dict


def read_quake_file(quake_file: Union[pathlib.Path, str]) -> dict:
    """Read a QuakeML file and parse the information into a dictionary

    :param quake_file: The full path to the quakeml file
    :type quake_file: Union[pathlib.Path, str]
    :return: The parsed moment tensor properties
    :rtype: dict
    """
    event = read_events(quake_file, format="QUAKEML")[0]
    moment_tensor = event.focal_mechanisms[0].moment_tensor
    tensor = moment_tensor.tensor
    mrr = tensor.m_rr * 10**7
    mtt = tensor.m_tt * 10**7
    mpp = tensor.m_pp * 10**7
    mrt = tensor.m_rt * 10**7
    mrp = tensor.m_rp * 10**7
    mtp = tensor.m_tp * 10**7
    m0 = moment_tensor.scalar_moment * 10**7
    is_rise_time = False
    if (
        moment_tensor.source_time_function is not None
        and moment_tensor.source_time_function.rise_time is not None
    ):
        is_rise_time = True
    if is_rise_time:
        time_shift = moment_tensor.source_time_function.rise_time
    else:
        time_shift = 5
    half_duration = 1.2 * 10**-8 * m0 ** (1 / 3)

    origins = event.origins
    # TODO: investigate why this assumes there are only two origins (old script got first and last).
    first_origin = origins[0]
    second_origin = origins[1]
    event_lat = first_origin.latitude
    event_lon = first_origin.longitude
    depth = first_origin.depth / 1000
    centroid_lat = second_origin.latitude
    centroid_lon = second_origin.longitude
    centroid_depth = second_origin.depth / 1000
    time = first_origin.time
    origin_time = UTCDateTime(time)
    date_time = origin_time.isoformat()

    input_dict = {
        "mrr": mrr,
        "mtt": mtt,
        "mpp": mpp,
        "mrt": mrt,
        "mrp": mrp,
        "mtp": mtp,
        "datetime": date_time,
        "date_origin": origin_time,
        "lat": event_lat,
        "lon": event_lon,
        "depth": depth,
        "time_shift": time_shift,
        "half_duration": half_duration,
        "centroid_lat": centroid_lat,
        "centroid_lon": centroid_lon,
        "centroid_depth": centroid_depth,
    }

    return input_dict


def modify_tensor(tensor_info: dict) -> dict:
    """Add more information to the tensor dictionary

    :param tensor_info: The tensor dictionary
    :type tensor_info: dict
    :return: The updated tensor information
    :rtype: dict
    """
    date_time = tensor_info["datetime"]
    tensor_info = {
        key: float(value) for (key, value) in tensor_info.items() if __is_number(value)
    }

    tensor_info["datetime"] = date_time
    tensor_info["date_origin"] = UTCDateTime(date_time)

    mzz = tensor_info["mrr"]
    mxx = tensor_info["mtt"]
    myy = tensor_info["mpp"]
    mxz = tensor_info["mrt"]
    myz = -tensor_info["mrp"]
    mxy = -tensor_info["mtp"]
    moment_tensor = np.zeros((3, 3))
    moment_tensor[0, :] = np.array([mxx, mxy, mxz])
    moment_tensor[1, :] = np.array([mxy, myy, myz])
    moment_tensor[2, :] = np.array([mxz, myz, mzz])

    tensor_info["moment_mag"] = __get_moment_mag(moment_tensor)
    return tensor_info


def planes_from_tensor(
    tensor_info: dict, dip_sens: float = 0, stk_sens: float = 0
) -> Tuple[dict, dict]:
    """Get nodal planes from tensor information

    :param tensor_info: The moment tensor properties
    :type tensor_info: dict
    :param dip_sens: Modify dip angle by this amount, defaults to 0
    :type dip_sens: float, optional
    :param stk_sens: Modify strike angle by this amount, defaults to 0
    :type stk_sens: float, optional
    :return: The nodal planes
    :rtype: Tuple[dict, dict]
    """
    mzz = tensor_info["mrr"]
    mxx = tensor_info["mtt"]
    myy = tensor_info["mpp"]
    mxz = tensor_info["mrt"]
    myz = -tensor_info["mrp"]
    mxy = -tensor_info["mtp"]
    moment_tensor = np.zeros((3, 3))
    moment_tensor[0, :] = np.array([mxx, mxy, mxz])
    moment_tensor[1, :] = np.array([mxy, myy, myz])
    moment_tensor[2, :] = np.array([mxz, myz, mzz])

    mt_eigenvalue, mt_eigenvectors = np.linalg.eigh(moment_tensor)
    imax = np.argmax(mt_eigenvalue)
    imin = np.argmin(mt_eigenvalue)
    slip_vector = mt_eigenvectors[:, imax] + mt_eigenvectors[:, imin]
    fault_normal = mt_eigenvectors[:, imax] - mt_eigenvectors[:, imin]
    if fault_normal[2] > 0:
        fault_normal = -fault_normal
        slip_vector = -slip_vector
    strike1, dip1, rake1 = __strike_dip_rake_from_ln(slip_vector, fault_normal)
    if slip_vector[2] > 0:
        fault_normal = -fault_normal
        slip_vector = -slip_vector
    strike2, dip2, rake2 = __strike_dip_rake_from_ln(fault_normal, slip_vector)
    strike1 = np.remainder(strike1, 360) + stk_sens
    strike2 = np.remainder(strike2, 360) + stk_sens
    dip1 = dip1 + dip_sens
    dip2 = dip2 + dip_sens
    rake1 = np.remainder(rake1, 360)
    rake2 = np.remainder(rake2, 360)
    if rake1 > 180:
        rake1 = rake1 - 360
    if rake2 > 180:
        rake2 = rake2 - 360
    np1_plane_info = {"strike": strike1, "dip": dip1, "rake": rake1}
    ffm_np1_data = {"plane_info": np1_plane_info}
    np2_plane_info = {"strike": strike2, "dip": dip2, "rake": rake2}
    ffm_np2_data = {"plane_info": np2_plane_info}
    return ffm_np1_data, ffm_np2_data


def __strike_dip_rake_from_ln(
    slip_vector: np.ndarray, fault_normal: np.ndarray
) -> Tuple[float, float, float]:
    """Calculate strike, dip, and rake from the slip vector and fault plane's
       normal vector (Following Udias Fig 16.19)

    :param slip_vector: The slip vectory
    :type slip_vector: np.ndarray
    :param fault_normal: The vector normal to the fault plane
    :type fault_normal: np.ndarray
    :return: The strike, dip, and rake
    :rtype: Tuple[float,float,float]

    .. warning::
        This routine is copied from the same routine in ``instaseis`` library.
        The only modification is here we work in xyz coordinates, as opposed
        to mrt as in said routine.
        https://github.com/krischer/instaseis/blob/master/instaseis/source.py#L85
    """
    l_norm = slip_vector / np.linalg.norm(slip_vector, 2)
    n_norm = fault_normal / np.linalg.norm(fault_normal, 2)
    delta = np.arccos(-n_norm[2])
    phi = np.arctan2(-n_norm[0], n_norm[1])

    # needs two different formulas, because the first is unstable for dip = 0
    # and the second for dip = 90
    if delta > 0.1:
        lambd = np.arctan2(
            -l_norm[2],
            np.sin(delta) * (l_norm[0] * np.cos(phi) + l_norm[1] * np.sin(phi)),
        )
    else:
        lambd = np.arctan2(
            (l_norm[0] * np.sin(phi) - l_norm[1] * np.cos(phi)),
            np.cos(delta) * (l_norm[0] * np.cos(phi) + l_norm[1] * np.sin(phi)),
        )

    strike = np.rad2deg(phi)
    dip = np.rad2deg(delta)
    rake = np.rad2deg(lambd)

    return strike, dip, rake


def __get_moment_mag(moment_tensor: np.ndarray) -> float:
    """Calculate the moment magnitude from the moment tensor following Herrmann, 1989

    :param moment_tensor: The moment tensor array
    :type moment_tensor: np.ndarray
    :return: The moment magnitude
    :rtype: float
    """
    mt_eigenvalue = np.linalg.eigh(moment_tensor)[0]
    imax = np.argmax(mt_eigenvalue)
    amax = mt_eigenvalue[imax]
    imin = np.argmin(mt_eigenvalue)
    amin = mt_eigenvalue[imin]
    moment_mag = (np.abs(amax) + np.abs(amin)) / 2
    return moment_mag


def __is_number(value: Union[str, float, int]) -> bool:
    """Determine if the provided value is a number

    :param value: The provided value
    :type value: Union[str,float,int]
    :return: Whether the value is a number
    :rtype: bool
    """
    try:
        float(value)
        return True
    except ValueError:
        return False
    except TypeError:
        return False


if __name__ == "__main__":
    import argparse

    parser = argparse.ArgumentParser()
    parser.add_argument(
        "-gcmt", "--gcmt_tensor", help="location of GCMT moment tensor file"
    )
    args = parser.parse_args()

    tensor_info = read_gcmt_file(args.gcmt_tensor)
    tensor_info = modify_tensor(tensor_info)
    delta = datetime.utcnow() - tensor_info["date_origin"]
    tensor_info["timedelta"] = delta.total_second
    moment_mag = tensor_info["moment_mag"]
    moment_mag = 2 * np.log10(moment_mag) / 3 - 10.7
    print(tensor_info)
    write_tensor(tensor_info)
