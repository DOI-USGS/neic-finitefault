#!/usr/bin/python
# -*- coding: utf-8 -*-
"""Module for finding the required velocity model, and creating json and
dictionaries with this velocity model.
"""

import numpy as np
import json
from netCDF4 import Dataset

    
def select_velmodel(tensor_info, default_dirs):
    """Script to select velocity model given hypocenter location

    :param tensor_info: dictionary with moment tensor information
    :type tensor_info: dict
    :returns: dictionary with properties of velocity model

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
    # crust_model = __crust_crust_velmodel(tensor_info, default_dirs)
    crust_model = __litho_crust_velmodel(tensor_info, default_dirs)
    velmodel = __crust_mantle_model(crust_model, tensor_info['depth'])
    velmodel2json(velmodel)
    return velmodel


def velmodel2json(velmodel):
    """Write a dictionary with a velocity model to a JSON file
    
    :param default_dirs: dictionary with velocity model
    :type default_dirs: dict
    """
    with open('velmodel_data.json','w') as f:
         json.dump(
                 velmodel, f, sort_keys=True, indent=4,
                 separators=(',', ': '), ensure_ascii=False)
    return


def water_depth(tensor_info):
    """We guess the water depth at a location from the crust velocity model.
    """
    crust_velmodel = __crust_crust_velmodel(tensor_info)
    s_vel = crust_velmodel['s_vel']
    thick = crust_velmodel['thick']
    water_depth = float(thick[0])
    return float(water_depth) if float(s_vel[0]) < 0.1 else 0
    
    
def __crust_crust_velmodel(tensor_info, default_dirs):
    """We find the velocity model interpolated from crust2.0, for the location
    of the hypocenter of the earthquake
    
    :param tensor_info: dictionary with moment tensor information
    :param default_dirs: dictionary with default directories to be used
    :type default_dirs: dict
    :type tensor_info: dict
    """
#
# code for crust type
#
    lat = tensor_info['lat']
    lon = tensor_info['lon']
    crust_codes = default_dirs['crust_codes']
    models_codes = default_dirs['models_codes']
    with open(crust_codes, 'r') as input_file:
        data = input_file.readlines()
    codes = [__process_line(line)[1:] for line in data[1:]]

    dx = 2
    cola = 90.0 - lat
    if lon >= 180:
        lon = lon - 360.0
    ilat = int(cola / dx)
    ilon = int((lon + 180.0) / dx)

    p_vel_crust = np.zeros(7)
    s_vel_crust = np.zeros(7)
    dens_crust = np.zeros(7)
    thick_crust = np.zeros(7)
#
# cruts velocity model
#
    with open(models_codes, 'r') as input_file:
        data = input_file.readlines()
    j = 0
    for i in range(360):
        j = j + 5
        if data[j][:2] == codes[ilat][ilon]:
            print(codes[ilat][ilon])
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
            p_vel_crust, s_vel_crust, dens_crust, thick_crust,
            qa_crust, qb_crust)
                
    return crust_velmodel
    
    
def __crust_mantle_model(crust_model, depth):
    """We add the prem model for the mantle, with a given crust velocity model.
    """
    max_depth = depth + 60
#
# PREM velocity model
#        
    p_vel_mantle = np.array(
        [8.080, 8.594, 8.732, 8.871, 9.219, 9.561, 9.902, 10.073, 10.212,
         10.791, 10.869])
    s_vel_mantle = np.array(
        [4.473, 4.657, 4.707, 4.757, 4.981, 5.176, 5.370, 5.467, 5.543, 5.982,
         6.056])
    dens_mantle = np.array(
        [3.3754, 3.4465, 3.4895, 3.5325, 3.7448, 3.8288, 3.9128, 3.9548,
         3.9840, 4.3886, 4.4043])
    thick_mantle = np.array(
        [196.000, 36.000, 108.00, 36.000, 33.333, 100.00, 33.333, 33.333,
         70.000, 25.250, 0.0])
    qa_mantle = np.array(
        [1.2e+03, 3.6e+02, 3.7e+02, 3.7e+02, 3.7e+02, 3.6e+02, 3.6e+02,
         3.6e+02, 3.6e+02, 7.6e+02, 7.5e+02])
    qb_mantle = np.array(
        [5.0e+02, 1.4e+02, 1.4e+02, 1.4e+02, 1.4e+02, 1.4e+02, 1.4e+02,
         1.4e+02, 1.4e+02, 3.1e+02, 3.1e+02])
    mantle_model = __dict_velmodel(
        p_vel_mantle, s_vel_mantle, dens_mantle,
        thick_mantle, qa_mantle, qb_mantle)

    p_vel = np.concatenate([crust_model['p_vel'], mantle_model['p_vel']])
    s_vel = np.concatenate([crust_model['s_vel'], mantle_model['s_vel']])
    dens = np.concatenate([crust_model['dens'], mantle_model['dens']])
    thick = np.concatenate([crust_model['thick'], mantle_model['thick']])
    qa = np.concatenate([crust_model['qa'], mantle_model['qa']])
    qb = np.concatenate([crust_model['qb'], mantle_model['qb']])

    depth = 0
    
    j = 0
    for i, thick_layer in enumerate(thick):
        depth = depth + float(thick_layer)
        if depth > max_depth:
            j = i + 2
            break
    j = len(thick) if j == 0 else j
    
    velmodel = __dict_velmodel(
        p_vel[:j], s_vel[:j], dens[:j], thick[:j], qa[:j], qb[:j])
    return velmodel
    
    
def model2dict(model_file):
    """We get a dictionary with a custom crust velocity model.
    """
    with open(model_file, 'r') as input_file:
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
        p_vel_crust, s_vel_crust, dens_crust, thick_crust,
        qa_crust, qb_crust)
    return crust_model
    
    
def __prem_model(depth, water_depth):
    """Prem velocity model for the crust
    """
    p_vel_crust = np.array([5.800, 6.800])
    s_vel_crust = np.array([3.200, 3.900])
    dens_crust = np.array([2.600, 2.900])
    thick_crust = np.array([12.000, 9.400])
    qa_crust = np.array([1.5e+03, 1.4e+03])
    qb_crust = np.array([6.0e+02, 6.0e+02])
    if water_depth > 0.5:
        p_vel_crust = np.concatenate([[1.500], p_vel_crust])
        s_vel_crust = np.concatenate([[0.01], s_vel_crust])
        dens_crust = np.concatenate([[1.000], dens_crust])
        thick_crust = np.concatenate([[3.000], thick_crust])
        qa_crust = np.concatenate([[1.5e+03], qa_crust])
        qb_crust = np.concatenate([[1.5e+03], qa_crust])
        
    crust_prem = __dict_velmodel(
        p_vel_crust, s_vel_crust, dens_crust, thick_crust, qa_crust, qb_crust)
    
    prem_velmodel = __crust_mantle_model(crust_prem, depth)
    return prem_velmodel
    
    
def __dict_velmodel(p_vel, s_vel, dens, thick, qa, qb):
    """Helper routine. We create a dictionary for a given velocity model.
    """
    velmodel_dict = {
        'p_vel': [str(v) for v in p_vel],
        's_vel': [str(v) for v in s_vel],
        'dens': [str(v) for v in dens],
        'thick': [str(v) for v in thick],
        'qa': [str(v) for v in qa],
        'qb': [str(v) for v in qb]
    }
    return velmodel_dict


def __litho_crust_velmodel(tensor_info, default_dirs):
    """We find the velocity model interpolated from litho1.0, for the location
    of the hypocenter of the earthquake
    
    :param tensor_info: dictionary with moment tensor information
    :param default_dirs: dictionary with default directories to be used
    :type default_dirs: dict
    :type tensor_info: dict
    """
#
# code for crust type
#
    lat = tensor_info['lat']
    lon = tensor_info['lon']
    litho_model = default_dirs['litho_model']
    rootgrp = Dataset(litho_model, "r", format="NETCDF4")
    
    vars = rootgrp.variables
    latitudes = vars['latitude']
    latitudes = np.array([val for val in latitudes])
    longitudes = vars['longitude']
    longitudes = np.array([val for val in longitudes])
    
    latitudes2 = (latitudes - lat) ** 2
    longitudes2 =  (longitudes - lon) ** 2
    index_lat = np.argmin(latitudes2)
    index_lon = np.argmin(longitudes2)
    
    layers = [
        'water_top',
        'ice_top',
        'upper_sediments_top',
        'middle_sediments_top',
        'lower_sediments_top',
        'upper_crust_top',
        'middle_crust_top',
        'lower_crust_top',
        'lid_top',
        'asthenospheric_mantle_top'
    ]
    
    p_vel_crust = [vars[layer + '_vp'][index_lat, index_lon] for layer in layers]
    p_vel_crust = np.array([val for val in p_vel_crust if not np.isnan(val)])
    
    s_vel_crust = [vars[layer + '_vs'][index_lat, index_lon] for layer in layers]
    s_vel_crust = np.array([val for val in s_vel_crust if not np.isnan(val)])
    
    dens_crust = [vars[layer + '_density'][index_lat, index_lon] for layer in layers]
    dens_crust = np.array([val for val in dens_crust if not np.isnan(val)])/1000
    
    depth_crust = [vars[layer + '_depth'][index_lat, index_lon] for layer in layers]
    depth_crust = np.array([val for val in depth_crust if not np.isnan(val)])
    
    qb_crust = [vars[layer + '_qmu'][index_lat, index_lon] for layer in layers]
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
        p_vel_crust, s_vel_crust, dens_crust, depth_crust, qa_crust, qb_crust)
    
    p_vel_crust = model['p_vel'][:-1]
    s_vel_crust = model['s_vel'][:-1]
    dens_crust = model['dens'][:-1]
    thick_crust = model['thick'][:-1]
    qa_crust = model['qa'][:-1]
    qb_crust = model['qb'][:-1]

    indexes = [i for i, thick in enumerate(thick_crust) if float(thick) > 0.0001]

    p_vel_crust = np.array([p_vel_crust[i] for i in indexes])
    s_vel_crust = np.array([s_vel_crust[i] for i in indexes])
    dens_crust = np.array([dens_crust[i] for i in indexes])
    thick_crust = np.array([thick_crust[i] for i in indexes])
    qa_crust = np.array([qa_crust[i] for i in indexes])
    qb_crust = np.array([qb_crust[i] for i in indexes])
    crust_velmodel = __dict_velmodel(
        p_vel_crust, s_vel_crust, dens_crust, thick_crust,
        qa_crust, qb_crust)
                
    return crust_velmodel


def model_depth2thick(model_file):
    """
    """
    with open(model_file, 'r') as input_file:
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
        p_vel_crust, s_vel_crust, dens_crust, depth_crust, qa_crust, qb_crust)
    
    return model


def __depth2thick(p_vel, s_vel, dens, depth, qa, qb):
    """
    """
    thick = np.diff(depth)
    p_vel2 = 2 / (1/p_vel[:-1] + 1/p_vel[1:])
    p_vel2 = np.round(p_vel2, 2)
    s_vel2 = 2 / (1/s_vel[:-1] + 1/s_vel[1:])
    s_vel2 = np.round(s_vel2, 2)
    
    model = __dict_velmodel(p_vel2, s_vel2, dens[:-1], thick, qa[:-1], qb[:-1])
    return model


def __process_line(line):
    line = line.replace('\n', '')
    line = line.replace('\t', ' ')
    line = line.split()
    line = [string for string in line if string]
    for i in range(len(line)):
        try:
            line[i] = float(line[i])
            if line[i] == int(line[i]):
                line[i] = int(line[i])
        except:
            pass
    return line
    
    
if __name__ == '__main__':
    velmodel = model2dict('vel_model')
    velmodel2json(velmodel)
    import argparse
    import management as mng
    parser = argparse.ArgumentParser()
    parser.add_argument("-o", "--output", help="name of output file")
    parser.add_argument("--lat", help="latitude of hypocenter location")
    parser.add_argument("--lon", help="longitude of hypocenter location")
    parser.add_argument("--depth", help="depth of hypocenter location")
    args = parser.parse_args()
    default_dirs = mng.default_dirs()
    tensor_info = {
        'lat': float(args.lat),
        'lon': float(args.lon),
        'depth': float(args.depth)
    }
    velmodel = select_velmodel(tensor_info, default_dirs)
    with open(args.output, 'w') as outf:
        nlen = len(velmodel['p_vel'])
        outf.write('{}\n'.format(nlen))
        zipped = zip(velmodel['p_vel'], velmodel['s_vel'], velmodel['thick'],
                     velmodel['dens'], velmodel['qa'], velmodel['qb'])
        for p, s, thick, dens, qa, qb in zipped:
            outf.write('{} {} {} {} {} {}\n'.format(p, s, dens, thick, qa, qb))
