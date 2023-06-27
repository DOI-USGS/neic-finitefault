# -*- coding: utf-8 -*-
"""Routine for performing administrative tasks, such as changing folders, or
moving to a different folder.
"""


import os
import numpy as np
import json
from obspy.geodetics import locations2degrees, degrees2kilometers
from obspy.core.utcdatetime import UTCDateTime
import subprocess
from read_config import read_config


def theoretic_arrivals(model, dist, depth):
    """Store theoretic arrivals of some phases.
    
    :param model: traveltime model.
    :param dist: distance of event to station in degrees
    :param depth: depth of event in km
    :type model: TauPyModel
    :type dist: float
    :type depth: float
    """
    p_arrival = model.get_travel_times(
            source_depth_in_km=depth, distance_in_degree=dist,
            phase_list=['P'])
    pp_arrival = model.get_travel_times(
            source_depth_in_km=depth, distance_in_degree=dist,
            phase_list=['PP'])
    s_arrival = model.get_travel_times(
            source_depth_in_km=depth, distance_in_degree=dist,
            phase_list=['S'])
    ss_arrival = model.get_travel_times(
            source_depth_in_km=depth, distance_in_degree=dist,
            phase_list=['SS'])
    p_arrivals = model.get_travel_times(
            source_depth_in_km=depth, distance_in_degree=dist + 0.1,
            phase_list=['P'])
    p_slowness = (p_arrivals[0].time - p_arrival[0].time) / 0.1
    s_arrivals = model.get_travel_times(
            source_depth_in_km=depth, distance_in_degree=dist + 0.1,
            phase_list=['S'])
    s_slowness = (s_arrivals[0].time - s_arrival[0].time) / 0.1
    arrivals = {
        'p_arrival': p_arrival,
        'pp_arrival': pp_arrival,
        's_arrival': s_arrival,
        'ss_arrival': ss_arrival,
        'p_slowness': p_slowness,
        's_slowness': s_slowness
    }
    return arrivals
    
    
def default_dirs():
    """Environment variables.
    """
    config = read_config()
    paths = config['PATHS']
    info = paths['info']
    compute_near_gf = paths['compute_near_gf']
    get_near_gf = paths['get_near_gf']
    modelling = paths['modelling']
    default_dirs = {
        'root_dir': paths['code_path'],
        'long_gf_bank': paths['surf_gf_bank'],
        'crust_codes': os.path.join(info, 'CNtype2.txt'),
        'models_codes': os.path.join(info, 'CNtype2_key.txt'),
        'litho_model': os.path.join(info, 'LITHO1.0.nc'),
        'gf_bank': paths['surf_gf_bank'],
        'strong_motion_gf_bank2': os.path.join(
            compute_near_gf, 'green_bank_openmp_f95'),
        'strong_motion_gf': os.path.join(get_near_gf, 'get_strong_motion'),
        'cgps_gf_bank': os.path.join(get_near_gf, 'cgps'),
        'gps_gf': os.path.join(compute_near_gf, 'gf_static_f95'),
        'tele_gf': os.path.join(modelling, 'green_tele'),
        'finite_fault': os.path.join(modelling, 'run_modelling'),
        'forward': os.path.join(modelling, 'run_forward'),
        'trench_graphics': os.path.join(
            paths['cartopy_files'], 'PB2002_plates'),
        'sac_exec': paths['sac_exec']
    }
    return default_dirs


def start_time_id(tensor_info):
    """We get a string, by concatenating the data for the origin time of
    an earthquake
    
    :param tensor_info: dictionary with moment tensor information
    :type tensor_info: dict
    """
    datetime = tensor_info['datetime']
    chunk1, chunk2 = datetime.split('T')
    year, month, day = chunk1.split('-')[:3]
    hour, minute, second = chunk2.split(':')

    time_origin = '{}{}{}{}{}{}'.format(year, month, day, hour, minute, second)
    return time_origin


def _distazbaz(station_lat, station_lon, event_lat, event_lon):
    """We compute the distance, azimuth and back_azimuth, between two pairs
    lat lon.
    """
    if np.abs(event_lat - station_lat) < 10 ** -7\
    and np.abs(event_lon - station_lon) < 10 ** -7:
        return 0, 0, 0
    degrees2rad = np.pi / 180.0
    flattening = (298.257223563)**-1
    event_lat2 = (1 - flattening) * np.tan(event_lat * degrees2rad)
    event_lat2 = np.arctan(event_lat2)
    station_lat2 = (1 - flattening) * np.tan(station_lat * degrees2rad)
    station_lat2 = np.arctan(station_lat2)
    sigma = np.sin(event_lat2) * np.sin(station_lat2)
    sigma2 = np.cos(event_lat2) * np.cos(station_lat2)
    sigma2 = sigma2 * np.cos(np.abs(event_lon - station_lon) * degrees2rad)
    sigma = np.arccos(sigma + sigma2)
    p = (event_lat2 + station_lat2) / 2
    q = (station_lat2 - event_lat2) / 2
    x = (sigma - np.sin(sigma)) * np.sin(p)**2 * np.cos(q)**2 / (np.cos(sigma/2))**2
    y = (sigma + np.sin(sigma)) * np.cos(p)**2 * np.sin(q)**2 / (np.sin(sigma/2))**2
    distance = 6378 * (sigma - 0.5*flattening*(x + y))

    delta = np.abs(station_lon - event_lon)
    if event_lon <= station_lon:
        azimuth = np.arctan2(
            np.cos(station_lat * degrees2rad) * np.cos(event_lat * degrees2rad)\
            * np.sin(delta * degrees2rad),
            np.sin(station_lat * degrees2rad) - np.cos(sigma)\
            * np.sin(event_lat * degrees2rad))
        azimuth = azimuth / degrees2rad
    
        back_azimuth = np.arctan2(
            np.cos(station_lat * degrees2rad) * np.cos(event_lat * degrees2rad)\
            * np.sin(delta * degrees2rad),
            np.sin(event_lat * degrees2rad) - np.cos(sigma)\
            * np.sin(station_lat * degrees2rad))
        back_azimuth = 360 - back_azimuth / degrees2rad
    else:
        back_azimuth = np.arctan2(
            np.cos(station_lat * degrees2rad) * np.cos(event_lat * degrees2rad)\
            * np.sin(delta * degrees2rad),
            np.sin(event_lat * degrees2rad) - np.cos(sigma)\
            * np.sin(station_lat * degrees2rad))
        back_azimuth = back_azimuth / degrees2rad
    
        azimuth = np.arctan2(
            np.cos(station_lat * degrees2rad) * np.cos(event_lat * degrees2rad)\
            * np.sin(delta * degrees2rad),
            np.sin(station_lat * degrees2rad) - np.cos(sigma)\
            * np.sin(event_lat * degrees2rad))
        azimuth = 360 - azimuth / degrees2rad
    return distance, azimuth, back_azimuth


def coords2utm(lat, lon, ref_lon):
    """Transform (lat, lon) pair to utm coordinates. See also 
    """
    deg2rad = np.pi / 180
    lat2 = deg2rad * lat
    lon2 = deg2rad * lon
    ref_lon2 = ref_lon * deg2rad
    radius = 6378.137
    flattening = 1 / 298.257223563
    k0 = 0.9996
    easting0 = 500
    northing0 = 0 if lat >= 0 else 10000
    n = flattening / (2 - flattening)
    A = radius * (1 + (n**2) / 4 + (n**4) / 64 + (n**6) / 256) / (1 + n)
    alfa1 = 0.5*n - (2/3)*n**2 + (5/16)*n**3 + (41/180)*n**4\
        - (127/288)*n**5 + (7891/37800)*n**6
    alfa2 = (13/48)*n**2 - (3/5)*n**3 + (557/1440)*n**4\
        + (281/630)*n**5 - (1983433/1935360)*n**6
    alfa3 = (61/240)*n**3 - (103/140)*n**4 + (15061/26880)*n**5\
        + (167603/181440)*n**6
    alfa4 = (49561/161280)*n**4 - (179/168)*n**5 + (6601661/7257600)*n**6
    alfa5 = (34729/80640)*n**5 - (3418889/1995840)*n**6
    alfa6 = (212378941/319334400)*n**6
    t = np.arctanh(np.sin(lat2))\
        - 2*np.sqrt(n)*np.arctanh(2*np.sqrt(n)*np.sin(lat2)/(n + 1))/(n + 1)
    t = np.sinh(t)
    delta = lon2 - ref_lon2
    psi = np.arctan(t / np.cos(delta))
    eta = np.arcsinh(np.sin(delta) / np.sqrt(t**2 + np.cos(delta)**2))
    easting = eta\
        + alfa1 * np.cos(2*psi) * np.sinh(2*eta)\
        + alfa2 * np.cos(4*psi) * np.sinh(4*eta)\
        + alfa3 * np.cos(6*psi) * np.sinh(6*eta)\
        + alfa4 * np.cos(8*psi) * np.sinh(8*eta)\
        + alfa5 * np.cos(10*psi) * np.sinh(10*eta)\
        + alfa6 * np.cos(12*psi) * np.sinh(12*eta)
    easting = easting0 + k0*A*easting
    northing = psi\
        + alfa1 * np.sin(2*psi) * np.cosh(2*eta)\
        + alfa2 * np.sin(4*psi) * np.cosh(4*eta)\
        + alfa3 * np.sin(6*psi) * np.cosh(6*eta)\
        + alfa4 * np.sin(8*psi) * np.cosh(8*eta)\
        + alfa5 * np.sin(10*psi) * np.cosh(10*eta)\
        + alfa6 * np.sin(12*psi) * np.cosh(12*eta)
    northing = northing0 + k0*A*northing
    return easting, northing
    


def correct_response_file(tensor_info, pzfile):
    """Routine for selecting instrumental response only in period of earthquake
    """
    date_origin = tensor_info['date_origin']
    with open(pzfile, 'r') as infile:
        lines = [line.split() for line in infile]
    start_times = [line[-1] for line in lines if 'START' in line]
    end_times = [line[-1] for line in lines if 'END' in line]
    end_times[0] = '2030-12-31T23:59:59'
    print(start_times)
    print(end_times)
    print('\n\n\n\n')
    network_lines = [index for index, line in enumerate(lines) if 'NETWORK' in line]
    constant_lines = [index for index, line in enumerate(lines) if 'CONSTANT' in line]
    
    start_times2 = [UTCDateTime(time) for time in start_times]
    end_times2 = [UTCDateTime(time) for time in end_times]
    
    zipped1 = zip(start_times2, end_times2)
    
    index = next(i for i, (start, end) in enumerate(zipped1)\
                 if start <= date_origin <= end)
    index = max(0, index - 1)
    index1 = network_lines[index] - 1
    index2 = constant_lines[index] + 1
    
    with open(pzfile, 'w') as outfile:
        for line in lines[index1:index2]:
            line2 = ' '.join(line)
            outfile.write('{}\n'.format(line2))


def use_waveforms(data_type):
    """
    """
    use_waveforms = False
    if 'strong_motion' in data_type:
        use_waveforms = True
    if 'cgps' in data_type:
        use_waveforms = True
    if 'tele_body' in data_type:
        use_waveforms = True
    if 'surf_tele' in data_type:
        use_waveforms = True
    return use_waveforms


# def run_sac(command):
#     """Routine to execute sac with a given input command from python
#     """
#     dirs = default_dirs()
#     sac_exec = dirs['sac_exec']
#     p = subprocess.Popen(
#         sac_exec, stdin=subprocess.PIPE, stdout=subprocess.PIPE,
#         stderr=subprocess.PIPE)
#     out, err = p.communicate(input=command.encode('utf-8'))
#     p.terminate()
#     return out, err


if __name__ == '__main__':
    import glob
    import argparse
    import seismic_tensor as tensor

    parser = argparse.ArgumentParser()
    parser.add_argument("-f", "--folder", default=os.getcwd(),
                        help="folder where there are input files")
    parser.add_argument("-gcmt", "--gcmt_tensor",
                        help="location of GCMT moment tensor file")
    
    args = parser.parse_args()
    os.chdir(args.folder)
    files = glob.glob('SACPZ*')
    tensor_info = tensor.get_tensor(cmt_file=args.gcmt_tensor)
    for pzfile in files:
        correct_response_file(tensor_info, pzfile)
        
    

 
