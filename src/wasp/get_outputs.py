#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""Script with routines for retrieving both the kinematic model which solves
the inverse problem, and the synthetic waveforms produced by such model,
"""


import numpy as np
import json
import os
import errno
import plane_management as pl_mng
from obspy import read


##########################
# Get FFM model
##########################


def read_solution_static_format(segments):
    """We read a solution file in static format.

    :param segments: dictionary with properties of the fault segments
    :type segments: dict
    """
    lat = []
    lon = []
    depth = []
    slip = []
    rake = []
    trup = []
    trise = []
    tfall = []
    moment = []

    if not os.path.isfile('Solucion.txt'):
        raise FileNotFoundError(
            errno.ENOENT, os.strerror(errno.ENOENT), 'Solucion.txt')
    with open('Solucion.txt', 'r') as input_file:
        jk = [line.split() for line in input_file]

    faults_data = [index + 1 for index, line in enumerate(jk)\
                   if set(['#Lat.', 'Lon.', 'depth', 'slip']) <= set(line)]
    headers = [index for index, line in enumerate(jk)\
               if set(['#Fault_segment', 'nx(Along-strike)=', 'ny(downdip)='])\
               <= set(line)]
    headers = headers[1:] + [len(jk)]
    if not len(headers) == len(segments):
        raise RuntimeError(
            'Inconsistency between Fault.time and Solucion.txt.'
            ' Different amount of fault segments')
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
        stk_subfaults, dip_subfaults, delta_strike, delta_dip, hyp_stk, hyp_dip\
            = pl_mng.__unpack_plane_data(segment)
#
# Reshape the rupture process
#
        if not slip_fault.size == stk_subfaults * dip_subfaults:
            raise RuntimeError(
                'Inconsistency between Fault.time and Solucion.txt.'\
                ' Different size of fault segment')
        lat_fault.shape = dip_subfaults, stk_subfaults
        lon_fault.shape = dip_subfaults, stk_subfaults
        depth_fault.shape = dip_subfaults, stk_subfaults
        slip_fault.shape = dip_subfaults, stk_subfaults
        rake_fault.shape = dip_subfaults, stk_subfaults
        trup_fault.shape = dip_subfaults, stk_subfaults
        trise_fault.shape = dip_subfaults, stk_subfaults
        tfall_fault.shape = dip_subfaults, stk_subfaults
        moment_fault.shape = dip_subfaults, stk_subfaults
        lat = lat + [lat_fault]
        lon = lon + [lon_fault]
        depth = depth + [depth_fault]
        slip = slip + [slip_fault]
        rake = rake + [rake_fault]
        trup = trup + [trup_fault]
        trise = trise + [trise_fault]
        tfall = tfall + [tfall_fault]
        moment = moment + [moment_fault]

    solution = {
        'slip': slip,
        'rake': rake,
        'rupture_time': trup,
        'trise': trise,
        'tfall': tfall,
        'lat': lat,
        'lon': lon,
        'depth': depth,
        'moment': moment
    }
    return solution


def read_solution_fsp_format(file_name, custom=False):
    """We read a solution file in fsp format

    :param file_name: string with location of solution file in fsp format
    :param custom:
    :type file_name: string
    :type custom: bool, optional
    """
    lat = []
    lon = []
    depth = []
    slip = []
    rake = []
    trup = []
    trise = []
    width = []

    if not os.path.isfile(file_name):
        raise FileNotFoundError(
            errno.ENOENT, os.strerror(errno.ENOENT), file_name)
    with open(file_name, 'r') as input_file:
        jk = [line.split() for line in input_file]

    tensor_info = {
        'lat': float(jk[5][5]),
        'lon': float(jk[5][8]),
        'depth': float(jk[5][11])
    }
    n_segments = int(jk[14][8])
    subfaults_data = {
        'stk_subfaults': int(jk[12][5]),
        'dip_subfaults': int(jk[12][8]),
        'delta_strike': float(jk[13][5]),
        'delta_dip': float(jk[13][9]),
        'strike': float(jk[7][5])
    }

    line0 = [i for i, line in enumerate(jk)\
             if {'SOURCE', 'MODEL', 'PARAMETERS'} < set(line)]
    line0 = line0[0] + 9
    if n_segments == 1:
        for line in jk[line0:]:
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
            width0 = 0 if not custom else float(jk[line0 + 2][7])
            subfaults_seg = int(jk[line0 + 6][3])\
                if not custom else int(jk[line0 + 7][3])
            line0 = line0 + 10 if not custom else line0 + 11
            for line in jk[line0:line0 + subfaults_seg]:
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
            line0 = line0 + subfaults_seg + 1

    solution = {
        'slip': slip,
        'rake': rake,
        'trup': trup,
        'trise': trise,
        'lat': lat,
        'lon': lon,
        'depth': depth,
        'new_width': width
    }
    return tensor_info, solution, subfaults_data


##########################
# Get data
##########################


def get_data_dict(traces_info, syn_file=None, observed=True, checker=False,
                  obs_file=None, margin=10):
    """Fills dictionary with synthetic data at station and channel

    :param traces_info: list of dictionaries with stations and channels metadata
    :param syn_file: string with location of file with synthetic data
    :param observed: whether to add observed waveform or not
    :param margin: shifts observed waveforms to the left by `margin` seconds
    :type traces_info: list
    :type syn_file: string, optional
    :type observed: bool, optional
    :type margin: float, optional
    """
    used_channels = []
    if syn_file:
        with open(syn_file, 'r') as infile:
            lines = [line.split() for line in infile]

        lines0 = [i for i, line in enumerate(lines) if len(line) > 2]
        for file in traces_info:
            name = file['name']
            channel = file['component']
            channel = __get_channel(channel)
            indexes = [i for i in lines0 if lines[i][2]==name\
                and lines[i][3] in channel]
            this_channel = [name, channel]
            if not this_channel in used_channels:
                used_channels = used_channels + [this_channel]
                index = indexes[0]
            else:
                index = indexes[1]
            npts = int(lines[index][0])
            synthetic = [float(real) for real, imag in lines[index + 1:index + npts]]
            file['synthetic'] = np.array(synthetic)
    if observed:
        if obs_file:
            for file in traces_info:
                file['observed'] = _get_observed_from_chen(file, obs_file)
        else:
            for file in traces_info:
                file['observed'] = _get_observed_from_chen2(file, margin=margin)
                derivative = False if not 'derivative' in file else file['derivative']
                dt = file['dt']
                file['observed'] = np.gradient(
                    file['observed'], dt, edge_order=2)\
                    if derivative else file['observed']

    else:
        for file in traces_info:
            file['observed'] = [0 for i in range(1024)]
    return traces_info


def _get_synthetic_from_chen(file, syn_file):
    """Gets synthetic waveform from a file with synthetic data, for a given
    station
    """
    name = file['name']
    channel = file['component']
    channel = __get_channel(channel)
    with open(syn_file, 'r') as infile:
        lines = [line.split() for line in infile]

    lines0 = [i for i, line in enumerate(lines) if len(line) > 2]
    index = next(
        i for i in lines0 if lines[i][2]==name and lines[i][3] in channel
    )
    npts = int(lines[index][0])
    synthetic = [float(real) for real, imag in lines[index + 1:index + npts]]
    return np.array(synthetic)


def _get_observed_from_chen(file, obse_file):
    """Gets observed waveform from a file with observed data, for a given
    station
    """
    name = file['name']
    channel = file['component']
    channel = __get_channel(channel)
    with open(obse_file, 'r') as infile:
        lines = [line.split() for line in infile]

    lines0 = [i for i, line in enumerate(lines)
              if not __is_number(line[0])]
    indexes = [i for i in lines0
               if (len(lines[i]) > 1 and lines[i][1] == name)]
    index = next(i for i in indexes if lines[i + 1][1] in channel)
    npts = int(lines[index + 3][1])
    npts = min(file['duration'], npts)
    observed = [float(real[0]) for real in lines[index + 6:index + 5 + npts]]
    return np.array(observed)


def _get_observed_from_chen2(file, margin=10):
    """Fills dictionary with observed data at station and channel
    """
    file2 = file['file']
    st = read(file2)
    data = st[0].data
    dt = file['dt']
    index0 = file['start_signal'] - int(margin // dt)
    # index0 = max(index0, 0)
    index1 = file['start_signal'] + file['duration']
    index1 = max(index1, index0 + 1)
    if index0 >= 0:
        data = data[index0:index1]
    else:
        data = data[index0:]
    return data


def _old_get_observed_from_chen(file, obse_file):
    """Fills dictionary with observed data at station and channel

    :param file: dictionary with properties of the fault segments
    :param syn_file: dictionary with moment tensor information
    :type file: dict
    :type syn_file: dict
    """
    name = file['name']
    component = file['component']
    component = __get_component(component)
    start = file['start_signal']
    with open(obse_file, 'r') as infile:
        lines = [line.split() for line in infile]

    lines0 = [i for i, line in enumerate(lines) if len(line) > 2]
    indexes = [i for i in lines0 if lines[i][2] == name]
    index = next(i for i in indexes if lines[i][3] in component)
    npts = int(lines[index][0])
    observed = [float(real[0]) for real in lines[index + 1:index + npts]]
    return np.array(observed[start:])


def __get_channel(channel):
    """Auxiliary routine
    """
    if channel in ['P', 'BHZ']:
        channel = ['P', 'BHZ']
    if channel in ['SH', None, 'BH1', 'BH2', 'BHE', 'BHN']:
        channel = ['SH']
    if channel in ['HNZ', 'HLZ', 'BNZ']:
        channel = ['HNZ', 'HLZ', 'BNZ']
    if channel in ['HNE', 'HLE', 'BNE']:
        channel = ['HNE', 'HLE', 'BNE']
    if channel in ['HNN', 'HLN', 'BNN']:
        channel = ['HNN', 'HLN', 'BNN']
    if channel in ['LXZ', 'LHZ', 'LYZ']:
        channel = ['LXZ', 'LHZ', 'LYZ']
    if channel in ['LXE', 'LHE', 'LYE']:
        channel = ['LXE', 'LHE', 'LYE']
    if channel in ['LXN', 'LHN', 'LYN']:
        channel = ['LXN', 'LHN', 'LYN']
    if channel in ['dart']:
        channel = ['dart']
    return channel


def retrieve_gps(syn_name='static_synthetics.txt'):
    """Get inverted and observed GPS data with station location

    :param syn_name: name of data with synthetic GPS data
    :type syn_name: strong, optional
    """
    obs_data = json.load(open('static_data.json'))

    lats = [data['location'][0] for data in obs_data]
    lons = [data['location'][1] for data in obs_data]
    names = [data['name'] for data in obs_data]
    observed = [data['observed'] for data in obs_data]
    error = [data['data_error'] for data in obs_data]

    with open(syn_name, 'r') as inf:
        lines = [line.split() for line in inf]

    name_synthetic = [[line[1], line[4:]] for line in lines[1:]]
    synthetic = []
    for index, name in enumerate(names):
        synthetic_gps = [
            gps_syn for gps_name, gps_syn in name_synthetic if gps_name == name]
        synthetic = synthetic + synthetic_gps
    return names, lats, lons, observed, synthetic, error


def get_insar():
    """
    """
    insar_data = json.load(open('insar_data.json'))

    with open('insar_synthetics.txt', 'r') as syn_file:
        lines_syn = [line.split() for line in syn_file]

    lines_ramp = []
    ramps = []
    lines0 = 0
    lines1 = 1
    if 'ascending' in insar_data:
        asc_properties = insar_data['ascending']
        ramps = [asc_property['ramp'] for asc_property in asc_properties]
    if 'descending' in insar_data:
        desc_properties = insar_data['descending']
        ramps = ramps + [desc_property['ramp'] for desc_property in desc_properties]
    lines_ramp = []
    if any(ramps):
        with open('insar_ramp.txt', 'r') as ramp_file:
            lines_ramp = [line.split() for line in ramp_file]
    if 'ascending' in insar_data:
        asc_properties = insar_data['ascending']
        for asc_property in asc_properties:
            insar_points = []
            insar_asc = asc_property['name']
            ramp_asc = asc_property['ramp']
            with open(insar_asc, 'r') as asc_file:
                lines_asc = [line.split() for line in asc_file]
            lines_asc = [line for line in lines_asc if not '#' in ''.join(line)]
            lines1 = lines0 + len(lines_asc)
            ramp_track = [[0] * 5] * len(lines_asc)
            if ramp_asc:
                ramp_track = lines_ramp[lines0 + 1:lines1]
            zipped = zip(
                lines_asc[:],
                lines_syn[lines0 + 1:lines1],
                ramp_track[:])
            for line1, line2, line3 in zipped:
                lat = float(line1[1])
                lon = float(line1[0])
                observed = 100*float(line1[2])
                synthetic = float(line2[4])
                ramp = float(line3[4])
                new_dict = {
                    'lat': lat,
                    'lon': lon,
                    'observed': observed,
                    'synthetic': synthetic,
                    'ramp': ramp
                }
                insar_points = insar_points + [new_dict]
            asc_property['points'] = insar_points
            lines0 = lines0 + len(lines_asc)
    if 'descending' in insar_data:
        desc_properties = insar_data['descending']
        for desc_property in desc_properties:
            insar_points = []
            insar_desc = desc_property['name']
            ramp_desc = desc_property['ramp']
            with open(insar_desc, 'r') as desc_file:
                lines_desc = [line.split() for line in desc_file]
            lines_desc = [line for line in lines_desc if not '#' in ''.join(line)]
            lines1 = lines0 + len(lines_desc)
            ramp_track = [[0] * 5] * len(lines_desc)
            if ramp_desc:
                ramp_track = lines_ramp[lines0 + 1:lines1]
            zipped = zip(
                lines_desc[:],
                lines_syn[lines0 + 1:lines1],
                ramp_track[:])
            for line1, line2, line3 in zipped:
                lat = float(line1[1])
                lon = float(line1[0])
                observed = 100*float(line1[2])
                synthetic = float(line2[4])
                ramp = float(line3[4])
                new_dict = {
                    'lat': lat,
                    'lon': lon,
                    'observed': observed,
                    'synthetic': synthetic,
                    'ramp': ramp
                }
                insar_points = insar_points + [new_dict]
            desc_property['points'] = insar_points
            lines0 = lines0 + len(lines_desc)
    # else:
    #     with open('insar_ramp.txt', 'r') as ramp_file:
    #         lines_ramp = [line.split() for line in ramp_file]
    #     if 'ascending' in insar_data:
    #         asc_properties = insar_data['ascending']
    #         for asc_property in asc_properties:
    #             insar_points = []
    #             insar_asc = asc_property['name']
    #             with open(insar_asc, 'r') as asc_file:
    #                 lines_asc = [line.split() for line in asc_file]
    #             lines_asc = [line for line in lines_asc if not '#' in ''.join(line)]
    #             lines1 = lines0 + len(lines_asc)
    #             zipped = zip(
    #                 lines_asc[:],
    #                 lines_syn[lines0 + 1:lines1],
    #                 lines_ramp[lines0 + 1:lines1])
    #             for line1, line2, line3 in zipped:
    #                 lat = float(line1[1])
    #                 lon = float(line1[0])
    #                 observed = 100*float(line1[2])
    #                 synthetic = float(line2[4])
    #                 ramp = float(line3[4])
    #                 new_dict = {
    #                     'lat': lat,
    #                     'lon': lon,
    #                     'observed': observed,
    #                     'synthetic': synthetic,
    #                     'ramp': ramp
    #                 }
    #                 insar_points = insar_points + [new_dict]
    #             asc_property['points'] = insar_points
    #             lines0 = lines0 + len(lines_asc)
    #     if 'descending' in insar_data:
    #         desc_properties = insar_data['descending']
    #         for desc_property in desc_properties:
    #             insar_points = []
    #             insar_desc = desc_property['name']
    #             with open(insar_desc, 'r') as desc_file:
    #                 lines_desc = [line.split() for line in desc_file]
    #             lines_desc = [line for line in lines_desc if not '#' in ''.join(line)]
    #             lines1 = lines0 + len(lines_desc)
    #             zipped = zip(
    #                 lines_desc[:],
    #                 lines_syn[lines0 + 1:lines1],
    #                 lines_ramp[lines0 + 1:lines1])
    #             for line1, line2, line3 in zipped:
    #                 lat = float(line1[1])
    #                 lon = float(line1[0])
    #                 observed = 100*float(line1[2])
    #                 synthetic = float(line2[4])
    #                 ramp = float(line3[4])
    #                 new_dict = {
    #                     'lat': lat,
    #                     'lon': lon,
    #                     'observed': observed,
    #                     'synthetic': synthetic,
    #                     'ramp': ramp
    #                 }
    #                 insar_points = insar_points + [new_dict]
    #             desc_property['points'] = insar_points
    #             lines0 = lines0 + len(lines_desc)
    return insar_data


def __is_number(string):
    try:
        float(string)
        return True
    except ValueError:
        return False
