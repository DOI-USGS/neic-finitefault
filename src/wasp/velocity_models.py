#!/usr/bin/python
# -*- coding: utf-8 -*-
"""Module for finding the required velocity model, and creating json and
dictionaries with this velocity model.
"""

import json
import pathlib
from typing import List, Union

import numpy as np
from netCDF4 import Dataset  # type:ignore


def select_velmodel(
    tensor_info: dict,
    default_dirs: dict,
    directory: Union[pathlib.Path, str] = pathlib.Path(),
) -> dict:
    """Select a velocity model given the hypocenter location

    :param tensor_info: The moment tensor information
    :type tensor_info: dict
    :param default_dirs: The paths to default directories
    :type default_dirs: dict
    :param directory: Where the file(s) should be read/written, defaults to pathlib.Path()
    :type directory: Union[pathlib.Path, str], optional
    :return: The velocity model
    :rtype: dict

    .. rubric:: Example:

    >>> tensor_info = {
            'moment_mag': 7 * 10 ** 27,
            'date_origin': UTCDateTime(2014, 4, 1, 23, 46, 47),
            'lat': -19.5,
            'lon': -70.5,
            'depth': 25,
            'time_shift': 44,
            'half_duration': 26,
            'centroid_lat': -21,
            'centroid_lon': -70,
            'centroid_depth': 35,
            'timedelta': 21 * 60
        }
    >>> velmodel = select_velmodel(tensor_info)
    >>> print(velmodel)
        {'p_vel': array([5.21, 5.37, 5.55, 5.72, 5.89, 5.98, 6.8, 7.01, 7.55,
                         8.05, 8.08]),
         's_vel': array([2.99, 3.09, 3.19, 3.29, 3.39, 3.44, 3.81, 3.95, 4.24,
                         4.39, 4.473]),
         'qb': array([150., 150., 150., 150., 150., 150., 300., 300., 300.,
                      300., 500.]),
         'dens': array([2.5, 2.5, 2.6, 2.7, 2.7, 2.8, 2.8, 2.9, 3., 3.4,
                        3.3754]),
         'water_level': 0,
         'qa': array([300., 300., 300., 300., 300., 300., 600., 600., 600.,
                      600., 1200.]),
         'thick': array([2.5, 2., 2., 2., 2. , 4.5, 10., 15., 10. , 20., 196.])
         }

    .. note::

        The locations of the velocity models can be modified in the routine
        ``default_dirs()`` of the module ``management.py``.

    """
    directory = pathlib.Path(directory)
    crust_model = __litho_crust_velmodel(tensor_info, default_dirs)
    velmodel = __crust_mantle_model(crust_model, tensor_info["depth"])
    velmodel2json(velmodel, directory=directory)
    return velmodel


def velmodel2json(velmodel: dict, directory: Union[pathlib.Path, str] = pathlib.Path()):
    """Write the velocity model dictionary to velmodel_data.json

    :param velmodel: The velocity model
    :type velmodel: dict
    :param directory: The directory to write to, defaults to pathlib.Path()
    :type directory: Union[pathlib.Path, str], optional
    """
    directory = pathlib.Path(directory)
    with open(directory / "velmodel_data.json", "w") as f:
        json.dump(
            velmodel,
            f,
            sort_keys=True,
            indent=4,
            separators=(",", ": "),
            ensure_ascii=False,
        )
    return


def __crust_crust_velmodel(tensor_info: dict, default_dirs: dict) -> dict:
    """Get the velocity model interpolated from crust2.0, for the location of
        the hypocenter of the earthquake

    :param tensor_info: The moment tensor information
    :type tensor_info: dict
    :param default_dirs: The paths to default directories
    :type default_dirs: dict
    :return: The velocity model
    :rtype: dict
    """
    #
    # code for crust type
    #
    lat = tensor_info["lat"]
    lon = tensor_info["lon"]
    crust_codes = default_dirs["crust_codes"]
    models_codes = default_dirs["models_codes"]
    with open(crust_codes, "r") as input_file:
        data = input_file.readlines()
    codes = [__process_line(line)[1:] for line in data[1:]]

    dx = 2
    cola = 90.0 - lat
    if lon >= 180:
        lon = lon - 360.0
    ilat = int(cola / dx)
    ilon = int((lon + 180.0) / dx)

    p_vel_crust: Union[np.ndarray, list] = np.zeros(7)
    s_vel_crust: Union[np.ndarray, list] = np.zeros(7)
    dens_crust: Union[np.ndarray, list] = np.zeros(7)
    thick_crust: Union[np.ndarray, list] = np.zeros(7)
    #
    # crust velocity model
    #
    with open(models_codes, "r") as input_file:
        data = input_file.readlines()
    j = 0
    for i in range(360):
        j = j + 5
        if data[j][:2] == codes[ilat][ilon]:
            p_vel_crust = __process_line(data[j + 1])
            s_vel_crust = __process_line(data[j + 2])
            dens_crust = __process_line(data[j + 3])
            thick_crust = __process_line(data[j + 4])[:7]
            #
            # we flip the first two layers.
            #
            aux = p_vel_crust[0]
            p_vel_crust[0] = p_vel_crust[1]
            p_vel_crust[1] = aux
            aux = s_vel_crust[0]
            s_vel_crust[0] = s_vel_crust[1]
            s_vel_crust[1] = aux
            aux = dens_crust[0]
            dens_crust[0] = dens_crust[1]
            dens_crust[1] = aux
            aux = thick_crust[0]
            thick_crust[0] = thick_crust[1]
            thick_crust[1] = aux
            break
    #
    # remove water layer
    #
    if s_vel_crust[0] <= 0.5:
        p_vel_crust = p_vel_crust[1:]
        s_vel_crust = s_vel_crust[1:]
        dens_crust = dens_crust[1:]
        thick_crust = thick_crust[1:]

    s_vel_crust = [max(s_vel, 0.01) for s_vel in s_vel_crust]
    indexes = [i for i, thick in enumerate(thick_crust) if thick > 0.0001]

    p_vel_crust = np.array([p_vel_crust[i] for i in indexes])
    s_vel_crust = np.array([s_vel_crust[i] for i in indexes])
    dens_crust = np.array([dens_crust[i] for i in indexes])
    thick_crust = np.array([thick_crust[i] for i in indexes])
    qa_crust = 1000 * np.ones(len(p_vel_crust))
    qb_crust = 500 * np.ones(len(p_vel_crust))
    crust_velmodel = __dict_velmodel(
        p_vel_crust, s_vel_crust, dens_crust, thick_crust, qa_crust, qb_crust
    )

    return crust_velmodel


def __crust_mantle_model(crust_model: dict, depth: float) -> dict:
    """Add the PREM model for the mantle with a given crust velocity model

    :param crust_model: The crust velocity model
    :type crust_model: dict
    :param depth: The depth
    :type depth: float
    :return: The velocity model
    :rtype: dict
    """
    max_depth = depth + 60
    #
    # PREM velocity model
    #
    p_vel_mantle = np.array(
        [
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
    s_vel_mantle = np.array(
        [4.473, 4.657, 4.707, 4.757, 4.981, 5.176, 5.370, 5.467, 5.543, 5.982, 6.056]
    )
    dens_mantle = np.array(
        [
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
    thick_mantle = np.array(
        [
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
    qa_mantle = np.array(
        [
            1.2e03,
            3.6e02,
            3.7e02,
            3.7e02,
            3.7e02,
            3.6e02,
            3.6e02,
            3.6e02,
            3.6e02,
            7.6e02,
            7.5e02,
        ]
    )
    qb_mantle = np.array(
        [
            5.0e02,
            1.4e02,
            1.4e02,
            1.4e02,
            1.4e02,
            1.4e02,
            1.4e02,
            1.4e02,
            1.4e02,
            3.1e02,
            3.1e02,
        ]
    )
    mantle_model = __dict_velmodel(
        p_vel_mantle, s_vel_mantle, dens_mantle, thick_mantle, qa_mantle, qb_mantle
    )

    p_vel = np.concatenate([crust_model["p_vel"], mantle_model["p_vel"]])
    s_vel = np.concatenate([crust_model["s_vel"], mantle_model["s_vel"]])
    dens = np.concatenate([crust_model["dens"], mantle_model["dens"]])
    thick = np.concatenate([crust_model["thick"], mantle_model["thick"]])
    qa = np.concatenate([crust_model["qa"], mantle_model["qa"]])
    qb = np.concatenate([crust_model["qb"], mantle_model["qb"]])

    depth = 0

    j = 0
    for i, thick_layer in enumerate(thick):
        depth = depth + float(thick_layer)
        if depth > max_depth:
            j = i + 2
            break
    j = len(thick) if j == 0 else j

    velmodel = __dict_velmodel(
        p_vel[:j], s_vel[:j], dens[:j], thick[:j], qa[:j], qb[:j]
    )
    return velmodel


def model2dict(model_file: Union[pathlib.Path, str]) -> dict:
    """Read the custom crust velocity model into a dictionary

    :param model_file: Path to the model file
    :type model_file: Union[pathlib.Path, str]
    :return: The parsed crust velocity model
    :rtype:
    """
    with open(model_file, "r") as input_file:
        lines = input_file.readlines()

    if len(__process_line(lines[0])) == 1:
        del lines[0]

    p_vel_crust = np.array([__process_line(line)[0] for line in lines])
    s_vel_crust = np.array([__process_line(line)[1] for line in lines])
    dens_crust = np.array([__process_line(line)[2] for line in lines])
    thick_crust = np.array([__process_line(line)[3] for line in lines])
    qa_crust = np.array([__process_line(line)[4] for line in lines])
    qb_crust = np.array([__process_line(line)[5] for line in lines])

    crust_model = __dict_velmodel(
        p_vel_crust, s_vel_crust, dens_crust, thick_crust, qa_crust, qb_crust
    )
    return crust_model


def __dict_velmodel(
    p_vel: np.ndarray,
    s_vel: np.ndarray,
    dens: np.ndarray,
    thick: np.ndarray,
    qa: np.ndarray,
    qb: np.ndarray,
) -> dict:
    """Format the velocity model into a dictionary

    :param p_vel: The p velocity
    :type p_vel: np.ndarray
    :param s_vel: The s velocity
    :type s_vel: np.ndarray
    :param dens: The density of layers
    :type dens: np.ndarray
    :param thick: The thickness of layer
    :type thick: np.ndarray
    :param qa: qa
    :type qa: np.ndarray
    :param qb: qb
    :type qb: np.ndarray
    :return: The velocity model
    :rtype: dict
    """
    velmodel_dict = {
        "p_vel": [str(v) for v in p_vel],
        "s_vel": [str(v) for v in s_vel],
        "dens": [str(v) for v in dens],
        "thick": [str(v) for v in thick],
        "qa": [str(v) for v in qa],
        "qb": [str(v) for v in qb],
    }
    return velmodel_dict


def __litho_crust_velmodel(tensor_info: dict, default_dirs: dict) -> dict:
    """Get the velocity model interpolated from litho1.0 for the location of
       the hypocenter of the earthquake

    :param tensor_info: The moment tensor information
    :type tensor_info: dict
    :param default_dirs: The paths to default directories
    :type default_dirs: dict
    :return: The velocity model
    :rtype: dict
    """
    #
    # code for crust type
    #
    lat = tensor_info["lat"]
    lon = tensor_info["lon"]
    litho_model = default_dirs["litho_model"]
    rootgrp = Dataset(litho_model, "r", format="NETCDF4")

    vars = rootgrp.variables
    latitudes = vars["latitude"]
    latitudes = np.array([val for val in latitudes])
    longitudes = vars["longitude"]
    longitudes = np.array([val for val in longitudes])

    latitudes2 = (latitudes - lat) ** 2
    longitudes2 = (longitudes - lon) ** 2
    index_lat = np.argmin(latitudes2)
    index_lon = np.argmin(longitudes2)

    layers = [
        "water_top",
        "ice_top",
        "upper_sediments_top",
        "middle_sediments_top",
        "lower_sediments_top",
        "upper_crust_top",
        "middle_crust_top",
        "lower_crust_top",
        "lid_top",
        "asthenospheric_mantle_top",
    ]

    p_vel_crust: Union[np.ndarray, list] = [
        vars[layer + "_vp"][index_lat, index_lon] for layer in layers
    ]
    p_vel_crust = np.array([val for val in p_vel_crust if not np.isnan(val)])

    s_vel_crust: Union[np.ndarray, list] = [
        vars[layer + "_vs"][index_lat, index_lon] for layer in layers
    ]
    s_vel_crust = np.array([val for val in s_vel_crust if not np.isnan(val)])

    dens_crust: Union[np.ndarray, list] = [
        vars[layer + "_density"][index_lat, index_lon] for layer in layers
    ]
    dens_crust = np.array([val for val in dens_crust if not np.isnan(val)]) / 1000

    depth_crust: Union[np.ndarray, list] = [
        vars[layer + "_depth"][index_lat, index_lon] for layer in layers
    ]
    depth_crust = np.array([val for val in depth_crust if not np.isnan(val)])

    qb_crust: Union[np.ndarray, list] = [
        vars[layer + "_qmu"][index_lat, index_lon] for layer in layers
    ]
    qb_crust = np.array([val for val in qb_crust if not np.isnan(val)])
    qa_crust = 2 * qb_crust
    #
    # remove water layer
    #
    if s_vel_crust[0] <= 0.1:
        p_vel_crust = p_vel_crust[1:]
        s_vel_crust = s_vel_crust[1:]
        dens_crust = dens_crust[1:]
        depth_crust = depth_crust[1:]

    model = __depth2thick(
        p_vel_crust, s_vel_crust, dens_crust, depth_crust, qa_crust, qb_crust
    )

    p_vel_crust = model["p_vel"][:-1]
    s_vel_crust = model["s_vel"][:-1]
    dens_crust = model["dens"][:-1]
    thick_crust = model["thick"][:-1]
    qa_crust = model["qa"][:-1]
    qb_crust = model["qb"][:-1]

    indexes = [i for i, thick in enumerate(thick_crust) if float(thick) > 0.0001]

    p_vel_crust = np.array([p_vel_crust[i] for i in indexes])
    s_vel_crust = np.array([s_vel_crust[i] for i in indexes])
    dens_crust = np.array([dens_crust[i] for i in indexes])
    thick_crust = np.array([thick_crust[i] for i in indexes])
    qa_crust = np.array([qa_crust[i] for i in indexes])
    qb_crust = np.array([qb_crust[i] for i in indexes])
    crust_velmodel = __dict_velmodel(
        p_vel_crust, s_vel_crust, dens_crust, thick_crust, qa_crust, qb_crust
    )

    return crust_velmodel


def model_depth2thick(model_file: Union[pathlib.Path, str]) -> dict:
    """Read the custom crust velocity model with depth into a dictionary

    :param model_file: Path to the model file
    :type model_file: Union[pathlib.Path, str]
    :return: The parsed crust velocity model
    :rtype:
    """
    with open(model_file, "r") as input_file:
        lines = input_file.readlines()

    if len(__process_line(lines[0])) == 1:
        del lines[0]

    p_vel_crust = np.array([__process_line(line)[0] for line in lines])
    s_vel_crust = np.array([__process_line(line)[1] for line in lines])
    dens_crust = np.array([__process_line(line)[2] for line in lines])
    depth_crust = np.array([__process_line(line)[3] for line in lines])
    qa_crust = np.array([__process_line(line)[4] for line in lines])
    qb_crust = np.array([__process_line(line)[5] for line in lines])

    model = __depth2thick(
        p_vel_crust, s_vel_crust, dens_crust, depth_crust, qa_crust, qb_crust
    )

    return model


def __depth2thick(
    p_vel: np.ndarray,
    s_vel: np.ndarray,
    dens: np.ndarray,
    depth: np.ndarray,
    qa: np.ndarray,
    qb: np.ndarray,
) -> dict:
    """Convert depths to thicknesses

    :param p_vel: The p velocity
    :type p_vel: np.ndarray
    :param s_vel: The s velocity
    :type s_vel: np.ndarray
    :param dens: The density of layers
    :type dens: np.ndarray
    :param depth: The depth of layer
    :type thick: np.ndarray
    :param qa: qa
    :type qa: np.ndarray
    :param qb: qb
    :type qb: np.ndarray
    :return: The velocity model
    :rtype: dict
    """
    thick = np.diff(depth)
    p_vel2 = 2 / (1 / p_vel[:-1] + 1 / p_vel[1:])
    p_vel2 = np.round(p_vel2, 2)
    s_vel2 = 2 / (1 / s_vel[:-1] + 1 / s_vel[1:])
    s_vel2 = np.round(s_vel2, 2)

    model = __dict_velmodel(p_vel2, s_vel2, dens[:-1], thick, qa[:-1], qb[:-1])
    return model


def __process_line(line: str) -> List[Union[float, int]]:
    """Parse numbers in a line

    :param line: The line of a file
    :type line: str
    :return: The line with numbers parsed
    :rtype: Union[float,int]
    """
    line = line.replace("\n", "")
    line = line.replace("\t", " ")
    line_split: list = line.split()
    line_split = [string for string in line_split if string]
    for i in range(len(line_split)):
        try:
            line_split[i] = float(line_split[i])
            if line_split[i] == int(line_split[i]):
                line_split[i] = int(line_split[i])
        except:
            pass
    return line_split
