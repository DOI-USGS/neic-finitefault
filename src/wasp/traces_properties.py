# -*- coding: utf-8 -*-
"""Module with automatic settings for sampling, length and frequency filtering
for teleseismic and near field data
"""


import numpy as np
import json


###################################
# automatic settings
###################################
    
    
def nyquist_frequency(dt):
    """
    """
    return 1 / 2 / dt
    

def sampling(tensor_info, dt_cgps=None):
    """Automatic setting for data sampling
    
    :param tensor_info: dictionary with moment tensor information
    :type tensor_info: dict
    """
    moment_mag = tensor_info['moment_mag']
    time_shift = tensor_info['time_shift'] 
    depth = tensor_info['depth']
    dt_tele = 0.2    
    if moment_mag > 10 ** 29:
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
    return {'dt_tele': dt_tele, 'dt_strong': dt_strong, 'dt_cgps': dt_cgps}


def filtro_tele(tensor_info):
    """Automatic settings for filter of teleseismic body waves
    
    :param tensor_info: dictionary with moment tensor information
    :type tensor_info: dict
    """
    time_shift = tensor_info['time_shift']
    depth = tensor_info['depth']
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
    filtro = {
            'freq0': freq0, 
            'low_freq': freq1, 
            'high_freq': freq2, 
            'freq3': freq3        
            }
    return filtro


def filtro_surf(tensor_info):
    """Automatic settings for filter of teleseismic surface waves

    :param tensor_info: dictionary with moment tensor information
    :type tensor_info: dict
    """
    filtro = {
        'freq1': 0.003,
        'freq2': 0.004,
        'freq3': 0.006,
        'freq4': 0.007
    }
    return filtro


def filtro_strong(tensor_info, cgps=False):
    """Automatic settings for filter of strong motion data.
    
    :param tensor_info: dictionary with moment tensor information
    :type tensor_info: dict
    """
    time_shift = tensor_info['time_shift']

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
    
    max_freq = 0.125# if tensor_info['time_shift'] > 10 else 0.25
    if cgps:
        max_freq = 0.3
    max_freq = max_freq if tensor_info['depth'] < 300 else 0.05
    filtro = {
        'low_freq': min_freq,
        'high_freq': max_freq
    }
    return filtro


def properties_json(tensor_info, dt_cgps):
    """We set automatic properties for waveform data to be used in modelling
    and write those to a JSON file.
    
    :type tensor_info: dict
    :param tensor_info: dictionary with moment tensor information
    """
    dict1 = sampling(tensor_info, dt_cgps)
    dict2 = filtro_tele(tensor_info)
    dict3 = filtro_surf(tensor_info)
    dict4 = filtro_strong(tensor_info)
    dict5 = filtro_strong(tensor_info, cgps=True)
    scales = wavelet_scales(tensor_info)
    dict6 = {
        'sampling': dict1,
        'tele_filter': dict2,
        'surf_filter': dict3,
        'strong_filter': dict4,
        'cgps_filter': dict5,
        'wavelet_scales': scales
    }
    with open('sampling_filter.json', 'w') as f:
         json.dump(
             dict6, f, sort_keys=True, indent=4, separators=(',', ': '),
             ensure_ascii=False)
    return dict6
    
    
def wavelet_scales(tensor_info):
    """Scales of wavelet transform to be used in modelling
    
    :param tensor_info: dictionary with moment tensor information
    :type tensor_info: dict
    """
    n_begin = 1
    n_end = 8
    return n_begin, n_end


if __name__ == '__main__':
    import argparse
    import seismic_tensor as tensor
    import os
    
    parser = argparse.ArgumentParser()
    parser.add_argument("-f", "--folder", default=os.getcwd(),
                        help="folder where there are input files")
    parser.add_argument("-gcmt", "--gcmt_tensor",
                        help="location of GCMT moment tensor file")
    args = parser.parse_args()
    os.chdir(args.folder)
    if args.gcmt_tensor:
        cmt_file = args.gcmt_tensor
        tensor_info = tensor.get_tensor(cmt_file=cmt_file)
    else:
        tensor_info = tensor.get_tensor() 
    properties_json(tensor_info)
