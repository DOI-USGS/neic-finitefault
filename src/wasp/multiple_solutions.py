import inversion_chen_new as inv
import os
import numpy as np
from multiprocessing import Pool, cpu_count
import pandas as pd
import seismic_tensor as tensor
import management as mng
import json
import shutil
from itertools import product
import glob
import fault_plane as pf


def multiple_solutions(tensor_info, data_type, default_dirs, folders,
                       strike=None, dip=None, rupt_vel=None):
    """Run multiple solutions, each with specified parameters

    :param tensor_info: dictionary with moment tensor properties
    :param data_type: list with data types to be used in modelling
    :param default_dirs: dictionary with default directories to be used
    :param folders: name of folders where to store parallel solutions
    :param strike: list of strike values to use
    :param dip: list of dip values to use
    :param rupt_vel: list of rupture velocity values to use
    :type tensor_info: dict
    :type data_type: list
    :type default_dirs: dict
    :type folders: string
    :type strike: list, optional
    :type dip: list, optional
    :type rupt_vel: list, optional
    """
    this_folder = os.path.abspath(os.getcwd())
    event_folder = os.path.dirname(this_folder)
    segments_data = json.load(open('segments_data.json'))
    segments = segments_data['segments']
    strike = strike if len(strike)>0 else [segments[0]['strike']]
    dip = dip if len(dip)>0 else [segments[0]['dip']]
    rupt_vel = rupt_vel if len(rupt_vel)>0 else [segments[0]['rupture_vel']]
    os.chdir(event_folder)

    new_iter = product(strike, dip, rupt_vel)
    subfolders = []
    for i, (strike1, dip1, rupt_vel1) in enumerate(new_iter):
        shutil.copytree(this_folder, '{}.{}'.format(folders, i))
        os.chdir('{}.{}'.format(folders, i))
        subfolders = subfolders + [os.path.abspath(os.getcwd())]
        new_segments_data = segments_data.copy()
        new_segments_data['segments'][0]['strike'] = strike1
        new_segments_data['segments'][0]['dip'] = dip1
        new_segments_data['segments'][0]['rupture_vel'] = rupt_vel1
        segments2 = new_segments_data['segments']
        rise_time = new_segments_data['rise_time']
        force_plane_above_ground(tensor_info, segments2)
        with open('segments_data.json', 'w') as f:
            json.dump(
                new_segments_data, f, sort_keys=True, indent=4,
                separators=(',', ': '), ensure_ascii=False)
        os.chdir(event_folder)

    cpus = cpu_count()
    processes = int(cpus / 3)
    print(processes)
    with Pool(processes=processes) as pool:
        results = pool.starmap(
            __worker,
            [
                [tensor_info, data_type, default_dirs, subfolder]\
                for subfolder in subfolders
            ]
        )


def __worker(tensor_info, data_type, default_dirs, subfolder):
    """
    """
    os.chdir(subfolder)
    segments_data = json.load(open('segments_data.json'))
    inv.manual_modelling(tensor_info, data_type, default_dirs, segments_data)
    return


def force_plane_above_ground(tensor_info, segments):
    """As name says

    :param tensor_info: dictionary with moment tensor properties
    :param segments: list of fault segments used
    :type tensor_info: dict
    :type segments: list
    """
    segment = segments[0]
    correct = pf.is_fault_correct(tensor_info, segment)
    if correct:
        return segments
    else:
        depth = tensor_info['depth']
        dip = segment['dip']
        hyp_dip = segment['hyp_dip']
        delta_dip = 0.99*depth/np.sin(dip*np.pi/180)/hyp_dip
        segments[0]['delta_dip'] = delta_dip
        return segments


def get_summary_all_models(subfolders):
    """Get a table with misfits and parameters in all subfolders

    :param subfolders: list with directions where results are located
    :type subfolders: list
    """
    summary_table = {}
    strike = []
    dip = []
    rupt_vel = []
    objective_error = []
    misfit_error = []
    for subfolder in subfolders:
        os.chdir(subfolder)
        segments_data = json.load(open('segments_data.json'))
        segments = segments_data['segments']
        strike = strike + [segments[0]['strike']]
        dip = dip + [segments[0]['dip']]
        rupt_vel = rupt_vel + [segments[0]['rupture_vel']]
        errors = get_summary()
        objective_error = objective_error + [errors['objective_error']]
        misfit_error = misfit_error + [errors['misfit_error']]
        os.chdir('..')
    df = pd.DataFrame(
        {
            'subfolders': subfolders,
            'strike': strike,
            'dip': dip,
            'rupt_vel': rupt_vel,
            'objective_error': objective_error,
            'misfit_error': misfit_error
        }
    )
    print(df)


def get_summary():
    """Get summary from FFM model
    """
    with open('modelling_summary.txt', 'r') as infile:
        lines = [line.split() for line in infile]

    error1 = 0
    error2 = 0
    for line in lines:
        if 'averaged' in line:
            error1 = float(line[-1])
        if 'objective' in line:
            error2 = float(line[-1])
    errors = {'misfit_error': error1, 'objective_error': error2}
    return errors


if __name__ == '__main__':
    import argparse
    import manage_parser as mp

    parser = argparse.ArgumentParser()
    parser.add_argument(
        "-f", "--folder", default=os.getcwd(),
        help="folder where there are input files")
    parser = mp.parser_add_tensor(parser)
    parser = mp.parser_ffm_data(parser)
    parser.add_argument(
        "--strike", default=[], nargs='*', type=float,
        help="strike values to model")
    parser.add_argument(
        "--dip", default=[], nargs='*', type=float,
        help="dip values to model")
    parser.add_argument(
        "--rupt_vel", default=[], nargs='*', type=float,
        help="rupture velocity values to model")
    parser.add_argument(
        "-r", "--run", action="store_true",
        help="run multiple models")
    parser.add_argument(
        "-v", "--view_results", action="store_true",
        help="view results from multiple models")
    parser.add_argument(
        "--folders", default='NP3',
        help="where to store results")
    args = parser.parse_args()

    this_folder = os.path.abspath(args.folder)
    os.chdir(this_folder)
    data_type = mp.get_used_data(args)

    default_dirs = mng.default_dirs()
    segments_data = json.load(open('segments_data.json'))
    if args.gcmt_tensor:
        cmt_file = args.gcmt_tensor
        tensor_info = tensor.get_tensor(cmt_file=cmt_file)
    else:
        tensor_info = tensor.get_tensor()
    strike = args.strike
    dip = args.dip
    rupt_vel = args.rupt_vel

    if args.run:
        multiple_solutions(
            tensor_info, data_type, default_dirs, args.folders, strike=strike,
            dip=dip, rupt_vel=rupt_vel)
    if args.view_results:
        os.chdir(this_folder)
        os.chdir('..')
        subfolders = glob.glob('{}.*'.format(args.folders))
        get_summary_all_models(subfolders)

