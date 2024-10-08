#!/usr/bin/python
# -*- coding: utf-8 -*-
"""Script for automatic creation of a fault plane, and for defining properties
of the fault plane.
"""


import errno
import json
import os
import pathlib
from typing import Any, Dict, List, Optional, Tuple, Union

import numpy as np
from obspy import UTCDateTime  # type: ignore

import wasp.management as mng
import wasp.seismic_tensor as tensor


def create_finite_fault(
    tensor_info: dict,
    np_plane_info: dict,
    data_type: List[str],
    water_level: float = 0,
    rupture_vel: Optional[float] = None,
    directory: Union[pathlib.Path, str] = pathlib.Path(),
) -> dict:
    """Create the fault plane and rise time function for the given plane and tensor

    :param tensor_info: The tensor information
    :type tensor_info: dict
    :param np_plane_info: The nodal plane properties
    :type np_plane_info: dict
    :param data_type: The data types used
    :type data_type: List[str]
    :param water_level: The water level depth, defaults to 0
    :type water_level: float, optional
    :param rupture_vel: The rupture velocity, defaults to None
    :type rupture_vel: Optional[float], optional
    :param directory: Where the file(s) should be read/written, defaults to pathlib.Path()
    :type directory: Union[pathlib.Path, str], optional
    :return: The segments data
    :rtype: dict

    .. rubric:: Example:

    >>> tensor_info = {
            'time_shift': 10.0,
            'depth': 25.0,
            'moment_mag': 10 ** 28,
            'lat': -10,
            'lon': -70,
            'centroid_lat': -11,
            'centroid_lon': -69
            }
    >>> np_plane_info = {'strike': 350, 'dip': 15, 'rake': 90}
    >>> data_type = ['strong_motion']
    >>> create_finite_fault(tensor_info, np_plane_info, data_type)
    """
    directory = pathlib.Path(directory)
    print("Create fault segments from input parameters")
    time_shift = tensor_info["time_shift"]
    strike = np_plane_info["strike"]
    dip = np_plane_info["dip"]
    rake = np_plane_info["rake"]
    rupture_vel = __default_vel_of_eq(tensor_info) if not rupture_vel else rupture_vel
    plane_info = __plane_tensor_def(strike, dip, rake, rupture_vel)
    plane_info2 = plane_info.copy()

    eq_time = 2 * time_shift + 0.75 * time_shift
    subfaults = __fault_plane_properties(eq_time, tensor_info, plane_info2, water_level)
    rise_time = __rise_time_parameters(tensor_info, subfaults, data_type)
    hyp_location = __hypocenter_location2(
        plane_info2,
        subfaults,
        tensor_info,
        water_level,
        rise_time,
    )

    plane_info2.update(subfaults)
    plane_info2.update(hyp_location)
    segments_data = __save_plane_data(
        plane_info, subfaults, hyp_location, rise_time, output_directory=directory
    )
    return segments_data


def point_sources_param(
    segments: dict,
    tensor_info: dict,
    rise_time: dict,
    connections: Optional[List[dict]] = None,
) -> list:
    """Define the point sources for fault segments for the given segments,
       hypocenter location, and moment tensor

    :param segments: The segments information (from segments_data.json "segments")
    :type segments: dict
    :param tensor_info: The tensor information
    :type tensor_info: dict
    :param rise_time: The rise time properties(from segments_data.json "rise_time")
    :type rise_time: dict
    :param connections: Connections properties, defaults to None
    :type connections: Optional[List[dict]], optional
    :raises Exception: If the point source is above the ground
    :return: The point sources
    :rtype: list
    """
    event_lat = tensor_info["lat"]
    event_lon = tensor_info["lon"]
    depth = tensor_info["depth"]
    delta_strike = segments[0]["delta_strike"]
    delta_dip = segments[0]["delta_dip"]
    rupture_vel = segments[0]["rupture_vel"]
    subfaults = {"delta_strike": delta_strike, "delta_dip": delta_dip}
    subfaults2 = _point_sources_def(rise_time, rupture_vel, subfaults)
    strike_ps = subfaults2["strike_ps"]
    dip_ps: int = subfaults2["dip_ps"]
    dx = subfaults2["dx"]
    dy = subfaults2["dy"]
    nx = int(strike_ps / 2.0 + 0.51)
    ny = int(dip_ps / 2.0 + 0.51)
    deg2rad = np.pi / 180.0

    point_sources: list = [[]] * len(segments)
    ref_coords: list = [[]] * len(segments)
    hypocenters: list = [[]] * len(segments)
    #
    # first we define reference coordinates!
    #
    segment = segments[0]
    strike = segment["strike"]
    dip = segment["dip"]
    hyp_stk = segment["hyp_stk"]
    hyp_dip = segment["hyp_dip"]
    x_ref = hyp_stk * delta_strike + nx * dx
    y_ref = hyp_dip * delta_dip + ny * dy
    ref_coords[0] = [event_lat, event_lon, depth, x_ref, y_ref]
    for i, segment in enumerate(segments):
        if not "hypocenter" in segment:
            hypocenters[i] = [event_lat, event_lon, depth]
        else:
            new_hypocenter = segment["hypocenter"]
            hyp_stk = segment["hyp_stk"]
            hyp_dip = segment["hyp_dip"]
            x1_ref = hyp_stk * delta_strike + nx * dx
            y1_ref = hyp_dip * delta_dip + ny * dy
            new_lat = new_hypocenter["lat"]
            new_lon = new_hypocenter["lon"]
            new_depth = new_hypocenter["depth"]
            ref_coords[i] = [new_lat, new_lon, new_depth, x1_ref, y1_ref]
            hypocenters[i] = [new_lat, new_lon, new_depth]
    if connections:
        for connection in connections:
            segment1 = connection["segment1"]
            segment2 = connection["segment2"]
            subfault1 = connection["segment1_subfault"]
            subfault2 = connection["segment2_subfault"]
            lat0, lon0, depth0, x_ref0, y_ref0 = ref_coords[segment1 - 1]
            first_segment = segments[segment1 - 1]
            strike = first_segment["strike"]
            dip = first_segment["dip"]
            n1_stk, n1_dip = subfault1
            x = n1_stk * delta_strike - x_ref0 + dx * 0
            y = n1_dip * delta_dip - y_ref0 + dy * 0
            dep_ref = y * np.sin(dip * deg2rad) + depth0
            lat_ref, lon_ref = __lat_lon(strike, dip, x, y, lat0, lon0)
            second_segment = segments[segment2 - 1]
            n2_stk, n2_dip = subfault2
            x2_ref = n2_stk * delta_strike
            y2_ref = n2_dip * delta_dip
            ref_coords[segment2 - 1] = [lat_ref, lon_ref, dep_ref, x2_ref, y2_ref]
    #
    # now we define the point sources for the segments
    #
    point_sources = [[]] * len(segments)
    zipped = zip(segments, ref_coords, hypocenters)
    for i, (segment, ref_coord, hypocenter) in enumerate(zipped):
        strike = segment["strike"]
        dip = segment["dip"]
        delay_segment = 0
        if "delay_segment" in segment:
            delay_segment = segment["delay_segment"]
        stk_subfaults: int = segment["stk_subfaults"]
        dip_subfaults: int = segment["dip_subfaults"]
        hyp_stk = segment["hyp_stk"]
        hyp_dip = segment["hyp_dip"]
        hypo_lat, hypo_lon, hypo_depth = hypocenter
        matrix = np.zeros((dip_subfaults, stk_subfaults, dip_ps, strike_ps, 7))  # type: ignore
        #
        # we give location of hypocenter relative to the fault segment
        #
        x_center = hyp_stk * delta_strike + nx * dx
        y_center = hyp_dip * delta_dip + ny * dy
        lat0, lon0, depth0, x0_ref, y0_ref = ref_coord
        for k2 in range(dip_subfaults):
            for j2 in range(stk_subfaults):
                for k1 in range(dip_ps):
                    for j1 in range(strike_ps):
                        #
                        # distance from the point source to the hypocenter over rupture surface
                        #
                        x = (j2 + 1) * delta_strike + (j1 + 1) * dx - x_center
                        y = (k2 + 1) * delta_dip + (k1 + 1) * dy - y_center
                        distance = np.sqrt(x**2 + y**2)
                        t1 = distance / rupture_vel
                        #
                        # depth of point source
                        #
                        x = (j2 + 1) * delta_strike + (j1 + 1) * dx - x0_ref
                        y = (k2 + 1) * delta_dip + (k1 + 1) * dy - y0_ref
                        dep = y * np.sin(dip * deg2rad) + depth0
                        if dep < 0.1:
                            raise Exception("Point source is above the ground!")
                        lat, lon = __lat_lon(strike, dip, x, y, lat0, lon0)
                        #
                        # distance over earth surface
                        #
                        dist, az, baz = mng._distazbaz(lat, lon, hypo_lat, hypo_lon)
                        matrix[k2, j2, k1, j1, :] = (
                            lat,
                            lon,
                            dep,
                            distance,
                            t1,
                            dist,
                            az,
                        )
        point_sources[i] = matrix
    return point_sources


def __lat_lon(
    strike: float, dip: float, x: float, y: float, lat0: float, lon0: float
) -> Tuple[float, float]:
    """Get the reverence latitude/longitude

    :param strike: The fault segment strike
    :type strike: float
    :param dip: The fault segment dip
    :type dip: float
    :param x: x distance from the point source to the hypocenter over rupture surface
    :type x: float
    :param y: y distance from the point source to the hypocenter over rupture surface
    :type y: float
    :param lat0: Reference latitude
    :type lat0: float
    :param lon0: Reference longitude
    :type lon0: float
    :return: The new coordinates
    :rtype: Tuple[float, float]
    """
    deg2rad = np.pi / 180.0
    cos_stk = np.cos(strike * deg2rad)
    sin_stk = np.sin(strike * deg2rad)
    cos_dip = np.cos(dip * deg2rad)
    degree = 111.12
    lat_ref = lat0 + (x * cos_stk - y * cos_dip * sin_stk) / degree
    lon_ref = lon0 + (x * sin_stk + y * cos_dip * cos_stk) / degree / np.cos(
        lat0 * deg2rad
    )
    return lat_ref, lon_ref


def shear_modulous(point_sources: list, velmodel: Optional[dict] = None) -> list:
    """Calculate the shear modulus for each subfault

    :param point_sources: Point source array for fault segments
    :type point_sources: list
    :param velmodel: The velocity model properties, defaults to None
    :type velmodel: Optional[dict], optional
    :return: The calculated shear moduli
    :rtype: list
    """
    if not velmodel:
        p_vel = np.array(
            [
                5.800,
                6.800,
                8.080,
                8.594,
                8.732,
                8.871,
                9.219,
                9.561,
                9.902,
                10.073,
                10.212,
                10.791,
                10.869,
            ]
        )
        s_vel = np.array(
            [
                3.200,
                3.900,
                4.473,
                4.657,
                4.707,
                4.757,
                4.981,
                5.176,
                5.370,
                5.467,
                5.543,
                5.982,
                6.056,
            ]
        )
        dens = np.array(
            [
                2.600,
                2.900,
                3.3754,
                3.4465,
                3.4895,
                3.5325,
                3.7448,
                3.8288,
                3.9128,
                3.9548,
                3.9840,
                4.3886,
                4.4043,
            ]
        )
        thick = np.array(
            [
                12.000,
                9.400,
                196.000,
                36.000,
                108.00,
                36.000,
                33.333,
                100.00,
                33.333,
                33.333,
                70.000,
                25.250,
                0.0,
            ]
        )
        velmodel = {"p_vel": p_vel, "s_vel": s_vel, "dens": dens, "thick": thick}
    vel_s = velmodel["s_vel"]
    dens = velmodel["dens"]
    thick = velmodel["thick"]
    #
    # now we compute the shear modulous at every subfault
    #
    shear: list = [[]] * len(point_sources)
    for segment, point_sources_seg in enumerate(point_sources):
        n_dip, n_stk, dip_ps, strike_ps, etc = point_sources_seg.shape
        depth_sources = point_sources_seg[:, :, :, :, 2]
        matrix = np.zeros((n_dip, n_stk))
        for i in range(n_dip):
            for j in range(n_stk):
                dep_p = depth_sources[i, j, dip_ps // 2, strike_ps // 2]
                source_layer = __source_layer(thick, dep_p)
                niu = (
                    float(vel_s[source_layer]) ** 2 * float(dens[source_layer]) * 10**10
                )
                matrix[i, j] = niu
        shear[segment] = matrix
    return shear


def __default_vel_of_eq(tensor_info: dict) -> float:
    """Get a default rupture velocity
       Rupture velocity is roughly 0.7*shear wave velocity, Vs. Vs increases with depth.
       Rupture velocity can be modified manually in segments_data.json.

    :param tensor_info: Dictionary with plane tensor information
    :type tensor_info: dict
    :return: Rupture velocity
    :rtype: float
    """
    #
    #  Default rupture velocity for shallow crustal events is 2.5 km/sec.
    #
    time_shift = tensor_info["time_shift"]
    moment_mag = tensor_info["moment_mag"]
    depth = tensor_info["depth"]
    default_vel = 2.5
    #
    # Default rupture velocity for intermediate depth earthquakes (100-300 km) is 3.0 km/sec.
    #
    if depth > 100:
        default_vel = 3.0
    #
    # Default rupture velocity for deep earthquakes (>300 km) is 3.6 km/sec.
    #
    if depth > 300:
        default_vel = 3.6
    #
    # Default rupture velocity for "slow" earthquakes follows Duputel et al. (2013),
    # Earth and Planetary Science Letters, in which earthquake is classified as slow if the
    # centroid time delay is much larger than a magnitude-based estimate of the half-duration
    #
    if time_shift / (1.2 * 10**-8 * moment_mag ** (1 / 3)) > 3:
        default_vel = 1.5

    return default_vel


def __fault_plane_properties(
    eq_time: int, tensor_info: dict, plane_info: dict, water_level: float
) -> dict:
    """Define dimensions of fault plane and subfaults

    :param eq_time: The earthquake time
    :type eq_time: int
    :param tensor_info: Dictionary with moment tensor information
    :type tensor_info: Dict[str, float]
    :param plane_info: Dictionary with plane tensor information
    :type plane_info: Dict[str, float]
    :param water_level: The water level
    :type water_level: float
    :return: The fault/subfault dimensions
    :rtype: dict
    """
    #
    # Fault dimensions
    #  The fault plane is constrained by following conditions (from Goldberg et al., 2022 SRL):
    #  1. Length = 3.3*Vr*Tc  where Vr is rupture velocity and Tc is the centroid time
    #  2. Width < 2z/sin(d)  where z is hypocentral depth and d is fault dip
    #  3. Width < Length
    #  4. Width < 300 km
    #  5. Fault plane is divided into ~225 subfaults of equal size
    #
    dip = plane_info["dip"]
    default_vel = plane_info["rupture_vel"]
    depth = tensor_info["depth"]
    dist_hypo_surface = max(depth - water_level, 0.8 * depth) / np.sin(
        dip * np.pi / 180.0
    )
    max_length = default_vel * eq_time
    max_width = min(300.0, max_length, 2 * dist_hypo_surface)
    max_width = max(max_width, 30)
    size0 = np.sqrt(max_width * max_length / 225.0)
    delta_strike = max(size0, 1.0)
    delta_dip = max(size0, 1.0)
    stk_subfaults = int(min(int(max_length / delta_strike), 45))
    if stk_subfaults % 2 == 0:
        stk_subfaults = stk_subfaults + 1
    max_length = 1.2 * max_length
    delta_strike = max_length / stk_subfaults
    dip_subfaults = int(min(max(int(max_width / delta_dip), 3), 15))
    if dip_subfaults % 2 == 0:
        dip_subfaults = dip_subfaults + 1
    if dist_hypo_surface < delta_dip:
        delta_dip = max(max_width / dip_subfaults, 1.9 * dist_hypo_surface)
    else:
        delta_dip = 0.99 * max_width / dip_subfaults

    fault_dimensions = __subfaults_properties(
        delta_strike, delta_dip, stk_subfaults, dip_subfaults
    )
    return fault_dimensions


def __rise_time_parameters(
    tensor_info: Dict[str, float],
    fault_dimensions: Dict[str, float],
    data_type: List[str],
) -> Dict[str, float]:
    """Automatically generate the rise time information

    :param tensor_info: Dictionary with plane tensor information
    :type tensor_info: Dict[str, float]
    :param fault_dimensions:  Dictionary with subfault information
    :type fault_dimensions: Dict[str, float]
    :param data_type: They data type
    :type data_type: List[str]
    :return: The rise time information
    :rtype: Dict[str, float]
    """
    delta_strike = fault_dimensions["delta_strike"]
    delta_dip = fault_dimensions["delta_dip"]

    #
    # finite fault
    #
    if tensor_info["time_shift"] <= 10:
        windows = int(1.5 * max(delta_strike, delta_dip) / 2)
        delta_rise = 1.0
    elif tensor_info["time_shift"] <= 24:
        windows = int(1.5 * max(delta_strike, delta_dip) / 4)
        delta_rise = 2.0
    elif tensor_info["time_shift"] >= 48:
        windows = int(1.5 * max(delta_strike, delta_dip) / 8)
        delta_rise = 4.0
    else:
        # TODO: Check whether this is ever evaluated and if it should be
        windows = int(
            1.5 * max(delta_strike, delta_dip) * 6 / tensor_info["time_shift"]
        )
        delta_rise = tensor_info["time_shift"] / 12
    if "body" in data_type:
        windows = max(int(1.5 * max(delta_strike, delta_dip) / 3), windows)
        delta_rise = min(1.5, delta_rise)
    if tensor_info["depth"] > 200:
        windows = int(1.5 * max(delta_strike, delta_dip) / 2)
        delta_rise = 1.0

    min_rise = delta_rise
    windows = windows + 2  # 3

    rise_time_param = {
        "min_rise": min_rise,
        "delta_rise": delta_rise,
        "windows": windows,
    }
    return rise_time_param


def _point_sources_def(
    rise_time_param: Dict[str, float],
    rupture_vel: float,
    fault_dimensions: Dict[str, float],
) -> dict:
    """Determine the amount of point sources per subfault from the subfault
       dimensions and rise time information

    :param rise_time_param: Dictionary with rise time information
    :type rise_time_param: Dict[str, float]
    :param rupture_vel: The rupture velocity
    :type rupture_vel: float
    :param fault_dimensions: Dictionary with subfault information
    :type fault_dimensions: Dict[str, float]
    :return: The general point source information
    :rtype: dict
    """
    delta_strike = fault_dimensions["delta_strike"]
    delta_dip = fault_dimensions["delta_dip"]
    t1 = rise_time_param["delta_rise"]
    delta = t1 * rupture_vel

    strike_ps = int(delta_strike / delta) + 1
    dip_ps = min(int(delta_dip / delta) + 1, 17)

    strike_ps = strike_ps + 1 if strike_ps % 2 == 0 else strike_ps
    dip_ps = dip_ps + 1 if dip_ps % 2 == 0 else dip_ps

    dx = delta_strike / strike_ps
    dy = delta_dip / dip_ps
    extra_info = __point_sources_general(strike_ps, dip_ps, dx, dy)
    return extra_info


def __hypocenter_location2(
    plane_info: Dict[str, float],
    fault_dimensions: Dict[str, float],
    tensor_info: Dict[str, float],
    water_level: float,
    rise_time: Dict[str, float],
) -> Dict[str, float]:
    """Determine in which subfault the hypocenter is located.
       Currently, we center the plane at the centroid in strike direction,
       and at the hypocenter, in dip direction.

    :param plane_info: Dictionary with plane tensor information
    :type plane_info: Dict[str, float]
    :param fault_dimensions: Dictionary with subfault information
    :type fault_dimensions: Dict[str, float]
    :param tensor_info: Dictionary with moment tensor information
    :type tensor_info: Dict[str, float]
    :param water_level: The water level
    :type water_level: float
    :param rise_time: Dictionary with rise time information
    :type rise_time: Dict[str, float]
    :return: The hypocenter location
    :rtype: Dict[str,float]
    """
    dip = plane_info["dip"]
    strike = plane_info["strike"]
    deg2rad = np.pi / 180
    degree = 111.19
    event_lat = tensor_info["lat"]
    event_lon = tensor_info["lon"]
    depth = tensor_info["depth"]
    centroid_lat = tensor_info["centroid_lat"]
    centroid_lon = tensor_info["centroid_lon"]
    stk_subfaults = fault_dimensions["stk_subfaults"]
    dip_subfaults = fault_dimensions["dip_subfaults"]
    delta_strike = fault_dimensions["delta_strike"]
    delta_dip = fault_dimensions["delta_dip"]
    rupture_vel = plane_info["rupture_vel"]
    subfaults = {"delta_strike": delta_strike, "delta_dip": delta_dip}
    subfaults2 = _point_sources_def(rise_time, rupture_vel, subfaults)
    strike_ps = subfaults2["strike_ps"]
    dip_ps = subfaults2["dip_ps"]
    cos_stk = np.cos(strike * deg2rad)
    sin_stk = np.sin(strike * deg2rad)
    cos_dip = np.cos(dip * deg2rad)
    matrix = np.array(
        [
            [cos_stk / degree, -cos_dip * sin_stk / degree],
            [
                sin_stk / (degree * np.cos(event_lat * deg2rad)),
                cos_dip * cos_stk / (degree * np.cos(event_lat * deg2rad)),
            ],
        ]
    )
    matrix = np.linalg.inv(matrix)
    vector = np.array([[centroid_lat - event_lat], [centroid_lon - event_lon]])
    solution = np.dot(matrix, vector)
    x, y = solution.flatten()

    hyp_stk: Union[float, int] = int(-x // delta_strike) + int(stk_subfaults / 2.0) + 1
    hyp_stk = max(1, min(stk_subfaults, hyp_stk))
    surface_dist = max(depth - water_level, 0.8 * depth) / np.sin(dip * np.pi / 180.0)
    hyp_dip = int(dip_subfaults / 2.0) + 1
    if delta_dip * dip_subfaults / 2.0 > surface_dist:
        for j in range(hyp_dip):
            if delta_dip * (j + 0.5) > 1.01 * surface_dist:
                break
        hyp_dip = j
    hyp_stk = hyp_stk if stk_subfaults > 1 else 1
    hyp_dip = hyp_dip if dip_subfaults > 1 else 1
    nx_hyp = int(strike_ps / 2.0 + 0.51)
    ny_hyp = int(dip_ps / 2.0 + 0.51)
    hyp_location = __epicenter_location(hyp_stk, hyp_dip)
    return hyp_location


def __source_layer(thick: np.ndarray, source_depth: float) -> int:
    """Get the index of the layer that the source is in

    :param thick: List of layer thicknesses
    :type thick: np.ndarray
    :param source_depth: The depth of the source
    :type source_depth: float
    :return: The index of the layer where the source exists
    :rtype: int
    """
    n_layers = len(thick)
    cumul_depth = np.zeros(n_layers + 1)

    for j in range(n_layers):
        cumul_depth[j + 1] = cumul_depth[j] + float(thick[j])

    for j in range(n_layers):
        if (source_depth >= cumul_depth[j]) and (source_depth <= cumul_depth[j + 1]):
            source_layer = j
            break
    return source_layer


def __plane_tensor_def(
    strike: float, dip: float, rake: float, rupture_vel: float
) -> Dict[str, float]:
    """Format tensor properties

    :param strike: The strike
    :type strike: float
    :param dip: The dip
    :type dip: float
    :param rake: The rake
    :type rake: float
    :param rupture_vel: The rupture velocity
    :type rupture_vel: float
    :return: The formatted properties
    :rtype: Dict[str, float]
    """
    values = {
        "strike": strike,
        "dip": dip,
        "rake": rake,
        "rupture_vel": rupture_vel,
        "min_vel": 0.4 * rupture_vel,
        "max_vel": 1.25 * rupture_vel,
    }
    return values


def __subfaults_properties(
    delta_strike: float, delta_dip: float, stk_subfaults: int, dip_subfaults: int
) -> Dict[str, Union[float, int]]:
    """Format subfault properties

    :param delta_strike: The strike delta
    :type delta_strike: float
    :param delta_dip: The dip delta
    :type delta_dip: float
    :param stk_subfaults: The number of strike subfaults
    :type stk_subfaults: int
    :param dip_subfaults: The number of dip subfaults
    :type dip_subfaults: int
    :return: The formatted properties
    :rtype: Dict[str, Union[float, int]]
    """
    values = {
        "delay_segment": 0,
        "delta_strike": delta_strike,
        "delta_dip": delta_dip,
        "stk_subfaults": stk_subfaults,
        "dip_subfaults": dip_subfaults,
    }
    return values


def __point_sources_general(
    strike_ps: float, dip_ps: float, dx: float, dy: float
) -> Dict[str, float]:
    """Format properties of point sources for each subfault

    :param strike_ps: The point source strike
    :type strike_ps: float
    :param dip_ps: The point source dip
    :type dip_ps: float
    :param dx: Delta x
    :type dx: float
    :param dy: Delta y
    :type dy: float
    :return: The formatted properties
    :rtype: Dict[str,float]
    """
    values = {"strike_ps": strike_ps, "dip_ps": dip_ps, "dx": dx, "dy": dy}
    return values


def __epicenter_location(hyp_stk: float, hyp_dip: float) -> Dict[str, float]:
    """Format epicenter information

    :param hyp_stk: Hypocenter strike
    :type hyp_stk: float
    :param hyp_dip: Hypocenter dip
    :type hyp_dip: float
    :return: Epicenter location
    :rtype: Dict[str, float]
    """
    values = {
        "hyp_stk": hyp_stk,  # hypocenter subfault
        "hyp_dip": hyp_dip,  #
    }
    return values


def __save_plane_data(
    plane_tensor: Dict[str, float],
    subfaults: Dict[str, float],
    epicenter_loc: Dict[str, float],
    rise_time: Dict[str, float],
    output_directory: Union[pathlib.Path, str] = pathlib.Path(),
) -> dict:
    """Save fault plane properties to json file

    :param plane_tensor: Dictionary with plane tensor information
    :type plane_tensor: Dict[str, float]
    :param subfaults: Dictionary with subfault information
    :type subfaults: Dict[str, float]
    :param epicenter_loc: Dictionary with epicenter information
    :type epicenter_loc: Dict[str, float]
    :param rise_time: Dictionary with rise time information
    :type rise_time: Dict[str, float]
    :param output_directory: The directory where modelling outputs exist,
                            defaults to pathlib.Path()
    :type output_directory: Union[pathlib.Path, str], optional
    :return: The plane properties
    :rtype: dict
    """
    output_directory = pathlib.Path(output_directory)
    segment_info: dict = {"neighbours": []}
    segment_info.update(plane_tensor)
    segment_info.update(subfaults)
    segment_info.update(epicenter_loc)
    segments_info = [segment_info]
    dictionary = {"segments": segments_info, "rise_time": rise_time}
    with open(output_directory / "segments_data.json", "w") as f:
        json.dump(
            dictionary,
            f,
            sort_keys=True,
            indent=4,
            separators=(",", ": "),
            ensure_ascii=False,
        )
    return dictionary


def is_fault_correct(
    tensor_info: Dict[str, Union[float, str, UTCDateTime]], segment: Dict[str, float]
) -> bool:
    """Check if fault will not touch surface

    :param tensor_info: Dictionary with moment tensor information
    :type tensor_info: Dict[str, Union[float, str, UTCDateTime]]
    :param segment: The segment properties
    :type segment: Dict[str, float]
    :return: True if the fault will not touch the surface
    :rtype: bool
    """
    depth = tensor_info["depth"]
    dip = segment["dip"]
    delta_dip = segment["delta_dip"]
    hyp_dip = segment["hyp_dip"]
    length = (hyp_dip - 0.2) * delta_dip
    height = np.sin(dip * np.pi / 180) * length
    return depth > height
