# -*- coding: utf-8 -*-
"""
The routines here calcaulte the shakemap polygon from the FFM 
"""
from typing import Tuple, Union

import numpy as np
from pyproj import Geod

geod = Geod(ellps="WGS84")


def equivalent_slip_length(slip_array: np.ndarray, subfault_len: float) -> int:
    """Get the length of the equivalent slip

    :param slip_array: The array of slip
    :type slip_array: np.ndarray
    :param subfault_len: The length of subfaults
    :type subfault_len: float
    :return: The length of equivalent slip
    :rtype: int
    """
    ### ALONG STRIKE AUTOCORRELATION ###
    autocorr = np.zeros(len(slip_array))
    for i in range(len(slip_array)):
        sum_val = 0.0
        for j in range(len(slip_array)):
            v1 = j - i
            if v1 >= 0:
                v2 = slip_array[v1]
            else:
                v2 = 0.0
            sum_val = sum_val + slip_array[j] * v2
        autocorr[i] = sum_val

    ### MIRROR FUNCTION ###
    autocorr_2 = np.zeros(2 * len(autocorr))
    autocorr_2[0 : len(autocorr)] = autocorr[::-1]
    autocorr_2[len(autocorr) : len(autocorr_2)] = autocorr
    ### INTEGRATE ###
    area = 0.0
    for k in range(len(autocorr_2) - 1):
        sq = autocorr_2[k] * subfault_len
        tri = 0.5 * subfault_len * (autocorr_2[k + 1] - autocorr_2[k])
        area = area + sq + tri
    ### DIVIDE BY ZERO LAG CORRELATION ###
    equivalent_len = area / autocorr[0]
    equivalent_len = round(equivalent_len, 2)

    return equivalent_len


def locate_equivalent_slip(
    slip_array: np.ndarray, subfault_len: float, eq_len: int
) -> float:
    """Locate the left edge of the equivalent slip

    :param slip_array: The array of slip
    :type slip_array: np.ndarray
    :param subfault_len: The length of subfaults
    :type subfault_len: float
    :param eq_len: The length of the equivalent slip
    :type eq_len: int
    :return: the location of the equivalent slip
    :rtype: float
    """
    num_subfaults = eq_len / subfault_len
    num_subfaults = int(round(num_subfaults, 0))
    excess_len = eq_len - (num_subfaults * subfault_len)
    left_edge_ind = 0
    max_slip = 0
    for edge_location in range(len(slip_array) - num_subfaults):
        total_slip = sum(slip_array[edge_location : edge_location + num_subfaults])
        if total_slip > max_slip:
            max_slip = total_slip
            left_edge_ind = edge_location
    left_edge: Union[float, int] = left_edge_ind * subfault_len
    left_edge = left_edge - (excess_len / 2)
    left_edge = left_edge - (subfault_len / 2)

    return left_edge


def translate_xy_to_latlondep(
    segment: dict,
    hyp_lon: float,
    hyp_lat: float,
    hyp_dep: float,
    eq_len_AS: int,
    eq_len_AD: int,
    left_edge_AS: float,
    left_edge_AD: float,
) -> Tuple[str, str, str, str]:
    strk = segment["strike"]
    dip = segment["dip"]
    ### calculate upper and lower latitudes at same along-strike location as hypocenter###
    dip_azimuth = strk + 90.0
    if dip_azimuth > 360.0:
        dip_azimuth = dip_azimuth - 360.0
    lon0, lat0, _ = geod.fwd(
        hyp_lon, hyp_lat, dip_azimuth, np.cos(np.radians(dip)) * (left_edge_AD * 1000.0)
    )
    lon1, lat1, _ = geod.fwd(
        hyp_lon,
        hyp_lat,
        dip_azimuth,
        np.cos(np.radians(dip)) * (left_edge_AD * 1000.0 + eq_len_AD * 1000.0),
    )
    ### move along-dip to get corner locations ###
    lon_c1, lat_c1, _ = geod.fwd(lon0, lat0, strk, left_edge_AS * 1000.0)
    lon_c4, lat_c4, _ = geod.fwd(lon1, lat1, strk, left_edge_AS * 1000.0)
    lon_c3, lat_c3, _ = geod.fwd(
        lon1, lat1, strk, left_edge_AS * 1000.0 + eq_len_AS * 1000.0
    )
    lon_c2, lat_c2, _ = geod.fwd(
        lon0, lat0, strk, left_edge_AS * 1000.0 + eq_len_AS * 1000.0
    )
    lon_c1 = "{:.2f}".format(lon_c1)
    lat_c1 = "{:.2f}".format(lat_c1)
    lon_c2 = "{:.2f}".format(lon_c2)
    lat_c2 = "{:.2f}".format(lat_c2)
    lon_c3 = "{:.2f}".format(lon_c3)
    lat_c3 = "{:.2f}".format(lat_c3)
    lon_c4 = "{:.2f}".format(lon_c4)
    lat_c4 = "{:.2f}".format(lat_c4)
    ### calculate top and bottom depths ###
    dist_dep0 = left_edge_AD
    dist_dep1 = left_edge_AD + eq_len_AD
    dep0 = hyp_dep + dist_dep0 * np.sin(np.radians(dip))
    dep1 = hyp_dep + dist_dep1 * np.sin(np.radians(dip))
    dep0 = round(dep0, 2)
    dep1 = round(dep1, 2)
    ### put corners together ###
    corner_1 = str(lon_c1) + " " + str(lat_c1) + " " + str(dep0)
    corner_2 = str(lon_c2) + " " + str(lat_c2) + " " + str(dep0)
    corner_3 = str(lon_c3) + " " + str(lat_c3) + " " + str(dep1)
    corner_4 = str(lon_c4) + " " + str(lat_c4) + " " + str(dep1)

    return corner_1, corner_2, corner_3, corner_4
