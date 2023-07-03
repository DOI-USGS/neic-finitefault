# -*- coding: utf-8 -*-
"""
"""


import os
import subprocess
import json
import logging
import modulo_logs as ml


def gf_retrieve(used_data_type, default_dirs):
    """Compute and store in binary files Green functions for each station, both
    for teleseismic body waves, as for strong motion, cGPS and static data.
    
    :param used_data_type: list with data types used in modelling
    :param default_dirs: dictionary with default directories to be used
    :type used_data_type: list
    :type default_dirs: dict
    """
    green_fun_tele = default_dirs['tele_gf']
    green_fun_str = default_dirs['strong_motion_gf']
    green_fun_gps = default_dirs['gps_gf']
    
    processes = []
    loggers = []
    data_types = []
    ch = logging.StreamHandler()
    formatter = logging.Formatter('%(levelname)s - %(message)s')
    ch.setFormatter(formatter)
    ch.setLevel(logging.ERROR)
    
    if 'tele_body' in used_data_type:
        print('Computing teleseismic GFs')
        logger1 = ml.create_log(
            'body_wave_GF', os.path.join('logs', 'green_tele_log'))
        logger1.addHandler(ch)
        p1 = subprocess.Popen(
            [green_fun_tele], stdout=subprocess.PIPE, stderr=subprocess.PIPE)
        # with open(os.path.join('logs', 'green_tele_log'), "w") as out_tele:
        #     p1 = subprocess.Popen([green_fun_tele], stdout=out_tele)
        processes = processes + [p1]
        loggers = loggers + [logger1]
        data_types = data_types + ['body waves']
    if 'strong_motion' in used_data_type:
        print('Computing strong motion GFs')
        logger2 = ml.create_log(
            'get_strong_motion_GF', os.path.join('logs', 'green_str_log'))
        logger2.addHandler(ch)
        p2 = subprocess.Popen(
            [green_fun_str], stdout=subprocess.PIPE, stderr=subprocess.PIPE)
        # with open(os.path.join('logs', 'green_np1_str_log'), "w") as out_strong:
        #     p2 = subprocess.Popen([green_fun_str], stdout=out_strong)
        processes = processes + [p2]
        loggers = loggers + [logger2]
        data_types = data_types + ['strong motion']
    if 'cgps' in used_data_type:
        print('Computing cGPS GFs')
        logger3 = ml.create_log(
            'get_cgps_GF', os.path.join('logs', 'green_cgps_log'))
        logger3.addHandler(ch)
        p3 = subprocess.Popen(
            [green_fun_str, 'cgps'], stdout=subprocess.PIPE,
            stderr=subprocess.PIPE)
        # with open(os.path.join('logs', 'green_np1_cgps_log'), "w") as out_cgps:
            # p3 = subprocess.Popen([green_fun_str, 'cgps'], stdout=out_cgps)
        processes = processes + [p3]
        loggers = loggers + [logger3]
        data_types = data_types + ['cgps']
    if 'gps' in used_data_type:
        print('Computing static GPS GFs')
        logger4 = ml.create_log(
            'GPS_GF', os.path.join('logs', 'green_gps_log'))
        logger4.addHandler(ch)
        p4 = subprocess.Popen(
            [green_fun_gps, 'gps'], stdout=subprocess.PIPE,
            stderr=subprocess.PIPE)
        processes = processes + [p4]
        loggers = loggers + [logger4]
        data_types = data_types + ['gps']
    if 'insar' in used_data_type:
        print('Computing InSAR GFs')
        logger5 = ml.create_log(
            'INSAR_GF', os.path.join('logs', 'green_insar_log'))
        logger5.addHandler(ch)
        p5 = subprocess.Popen(
            [green_fun_gps, 'insar'], stdout=subprocess.PIPE,
            stderr=subprocess.PIPE)
        processes = processes + [p5]
        loggers = loggers + [logger5]
        data_types = data_types + ['insar']

    # [p.wait() for p in processes]
    for p, log, data_type in zip(processes, loggers, data_types):
        out, err = p.communicate(timeout=50 * 60)
        log.info(out.decode('utf-8'))
        if err: log.error(err.decode('utf-8', 'ignore'))
        ml.close_log(log)
        if err:
            raise RuntimeError(
                'Got following error while retrieving GF '\
                'for {}:\n{}'.format(data_type, err)
            )
    
    
def fk_green_fun0(dt, tensor_info, default_dirs, gf_bank=None):
    """We write a file with important data for computing or retrieving strong
    motion Green functions.
    
    :param tensor_info: dictionary with moment tensor information
    :param dt: sampling rate of data
    :param default_dirs: dictionary with default directories to be used
    :param gf_bank: location of GF bank to be used.
    :type tensor_info: dict
    :type tensor_info: float
    :type default_dirs: dict
    :type gf_bank: string, optional
    """
#    dirs = mng.default_dirs()
    min_depth = 3
    if not gf_bank:
#
# strong motion
#
        min_depth = 3
        lat = tensor_info['lat']
        depth = tensor_info['depth']
        print(dt, lat, depth)
        if lat >= -28 and depth <= 200 and dt <= 0.2:
            gf_bank = 'Hussen2_0.2_600_3_300'
        elif lat >= -28 and depth <= 200 and dt > 0.2:
            gf_bank = 'Hussen2_0.4_1000_3_300'
        elif lat >= -28 and depth > 500:
            gf_bank = 'FAPE_0.8_1000_500_800'
            dt = 0.8
            min_depth = 500
#        if lat >= -32 and depth > 80 and dt <= 0.2:
#            gf_bank = 'FAPE_0.2_600_3_300'
        if lat < -40 and dt <= 0.2:
            gf_bank = 'crust_T6_0.2_600_3_300'
        elif lat < -40 and dt > 0.2:
            gf_bank = 'crust_T6_0.4_1000_3_300'
        if -28 > lat >= -32 and depth <= 140 and dt <= 0.2:
            gf_bank = 'Hussen_0.2_600_3_300'#zona_central
        elif -28 > lat >= -32 and depth <= 140 and dt > 0.2:
            gf_bank = 'Hussen_0.4_1000_3_300'
        if -32 > lat >= -36 and depth <= 140 and dt <= 0.2:
            gf_bank = 'Hicks_0.2_600_3_300'#zona_central
        elif -32 > lat >= -36 and depth <= 140 and dt > 0.2:
            gf_bank = 'Hicks_0.4_1000_3_300'#zona_central
        if -36 > lat >= -40 and depth <= 200 and dt <= 0.2:
            gf_bank = 'Hicks_0.2_600_3_300'
        elif -36 > lat >= -40 and depth <= 200 and dt > 0.2:
            gf_bank = 'Hicks_0.4_1000_3_300'#Maule
        if 30 < lat < 40 and dt <= 0.2:
            min_depth = 1
            gf_bank = 'california1_0.2_500_1_60'
#
# cGPS
#
        if lat >= -32 and depth <= 200 and abs(dt - 1.0) < 0.01:
            gf_bank = 'Hussen_1.0_1000_3_300'
        if lat < -40 and abs(dt - 1.0) < 0.01:
            gf_bank = 'crust_T6_1.0_1000_3_300'
        if -32 > lat >= -36 and depth <= 140 and abs(dt - 1.0) < 0.01:
            gf_bank = 'zona_central_1.0_1000_3_300'
        if -36 > lat >= -40 and depth <= 200 and abs(dt - 1.0) < 0.01:
            gf_bank = 'Hicks_1.0_1000_3_300'#Maule_crust
        location = os.path.join(default_dirs['strong_motion_gf_bank'], gf_bank)

    if dt < 0.9:
        with open(os.path.join(default_dirs['strong_motion_gf_bank'],
                               'gf_banks_data.json'),'r') as f:
            dicts = json.load(f)
    else:
        with open(os.path.join(default_dirs['cgps_gf_bank'],
                               'gf_banks_data.json'),'r') as f:
            dicts = json.load(f)
    
    max_depth = next(d['max_depth'] for d in dicts if d['file']==gf_bank)
    max_dist = next(d['max_dist'] for d in dicts if d['file']==gf_bank)

    if abs(dt - 1.0) < 0.01:
        location = os.path.join(default_dirs['cgps_gf_bank'], gf_bank)
    
    green_dict = {
        'location': location,
        'min_depth': min_depth,
        'max_depth': max_depth,
        'min_dist': 0,
        'max_dist': max_dist,
        'dt': dt
    }
    
    name = 'strong_motion_gf.json' if dt < 0.9 else 'cgps_gf.json'
    with open(name, 'w') as f:
        json.dump(
            green_dict, f, sort_keys=True, indent=4,
            separators=(',', ': '), ensure_ascii=False)
    return green_dict


def fk_green_fun1(data_prop, tensor_info, location, cgps=False, max_depth=None):
    """We write a file with important data for computing or retrieving strong
    motion Green functions.
    
    :param tensor_info: dictionary with moment tensor information
    :param dt: sampling rate of data
    :param location: location of GF bank to be used.
    :param cgps: whether GF are of cgps data or not.
    :param max_depth: maximum depth of point sources.
    :type tensor_info: dict
    :type tensor_info: float
    :type location: string
    :type cgps: bool, optional
    :type max_depth: float, optional
    """
    sampling = data_prop['sampling']
    dt = sampling['dt_strong'] if not cgps else sampling['dt_cgps']
    depth = tensor_info['depth']
    time_shift = tensor_info['time_shift']
    min_depth = max(1, depth - 100) 
    if not max_depth:
        max_depth = max(30, 2 * depth)
        max_depth = min(max_depth, depth + 60)
    max_depth = max_depth + 5
    min_dist = 0
    max_dist = 600 if time_shift < 40 else 1000
    time_corr = 10 if not cgps else 25
    
    green_dict = {
        'location': location,
        'min_depth': min_depth,
        'max_depth': max_depth,
        'min_dist': min_dist,
        'max_dist': max_dist,
        'dt': dt,
        'time_corr': time_corr
    }
    
    name = 'strong_motion_gf.json' if not cgps else 'cgps_gf.json'
    with open(name, 'w') as f:
        json.dump(
            green_dict, f, sort_keys=True, indent=4,
            separators=(',', ': '), ensure_ascii=False)
    return green_dict


if __name__ == '__main__':
    import argparse
    import seismic_tensor as tensor
    import management as mng
    import manage_parser as mp

    parser = argparse.ArgumentParser()
    parser.add_argument(
        "-f", "--folder", default=os.getcwd(),
        help="folder where there are input files")
    parser = mp.parser_add_tensor(parser)
    parser = mp.parser_add_gf(parser)
    parser.add_argument(
        "-dt", type=float, default=0.2,
        help="sampling step of strong motion data")
    args = parser.parse_args()
    os.chdir(args.folder)
    if args.gcmt_tensor:
        cmt_file = args.gcmt_tensor
        tensor_info = tensor.get_tensor(cmt_file=cmt_file)
    else:
        tensor_info = tensor.get_tensor()
    tensor_info['timedelta'] = 81 * 90
    used_data = mp.get_used_data(args)
    default_dirs = mng.default_dirs()
    if 'strong_motion' in used_data and not os.path.isfile('strong_motion_gf.json'):
        green_dict = fk_green_fun0(args.dt, tensor_info, default_dirs)
#        write_green_file(green_dict)
    if 'cgps' in used_data and not os.path.isfile('cgps_gf.json'):
        green_dict = fk_green_fun0(1.0, tensor_info)
#        write_green_file(green_dict)
    gf_retrieve(used_data, default_dirs)
