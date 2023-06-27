#!/usr/bin/python
# -*- coding: utf-8 -*-
"""Routine for doing a time shift of observed data.
"""


import numpy as np
import os
import json
import glob
import get_outputs
from obspy import read
import matplotlib.pyplot as plt
import seismic_tensor as tensor
from waveform_plots import plot_waveforms
from many_events import select_waveforms_event


def shift_match(data_type, plot=False, method='full', zero_start=True):
    """We shift synthetic data to maximize cross-correlation with observed data.

    :param data_type: list of data types to be used in modelling.
    :param plot: whether to plot results of waveform shift
    :param method: method for cross-correlation
    :type data_type: list
    :type tensor_info: bool, optional
    :type method: string, optional
    """
    if data_type == 'tele':
        json_file = 'tele_waves.json'
    if data_type == 'strong':
        json_file = 'strong_motion_waves.json'
    if data_type == 'cgps':
        json_file = 'cgps_waves.json'
    if data_type == 'surf':
        json_file = 'surf_waves.json'
    files = json.load(open(json_file))
    plot_folder = 'tele_shift' if data_type == 'tele' else 'strong_shift'\
        if data_type == 'strong' else 'surf_shift' if data_type == 'surf'\
        else 'cgps_shift'
    if not os.path.isdir(plot_folder):
        os.mkdir(plot_folder)
    synthetics_file = 'synthetics_body.txt' if data_type == 'tele' else 'synthetics_strong.txt'\
        if data_type == 'strong' else 'synthetics_surf.txt' if data_type == 'surf'\
        else 'synm.cgps'

    dt = float(files[0]['dt'])
    files = get_outputs.get_data_dict(files, syn_file=synthetics_file)
    for file in files:
        derivative = False if not 'derivative' in file else file['derivative']
        synt_tr = file['synthetic']
        stream = read(file['file'])
        nshift = int(5 / dt) if data_type == 'tele' else int(5 / dt)\
            if data_type == 'strong' else int(12 / dt) if data_type == 'surf'\
            else int(4 / dt)
        obser_tr = stream[0].data
        if derivative:
            obser_tr = np.diff(obser_tr)
        length = int(float(file['duration'])) if method == 'full' else int(25 / dt)
        start = int(file['start_signal'])
        tr_shift = _shift(obser_tr, synt_tr, nshift, length, start)
        file['start_signal'] = start + tr_shift
        new_baseline = stream[0].data[start + tr_shift]

        file['synthetic'] = []
        file['observed'] = []
        if plot:
            length2 = int(10 / dt)
            name = file['name']
            component = file['component']
            synthetic = synt_tr[:length]
            syn_waveforms = [synthetic, synthetic]
            time0 = np.arange(len(synthetic)) * dt
            fig, axes = plt.subplots(2, 1)
            axes2 = axes.ravel()
            fig.suptitle('{} {}'.format(name, component))
            derivative = False if not 'derivative' in file else file['derivative']
            time1, observed0 = get_observed(file, start, length, margin=10)
            observed0 = np.gradient(observed0, dt, edge_order=2) if derivative\
            else observed0
            axes[0].axvline(0)
            axes[0].set_title('Before Shift')
            start4 = start + tr_shift
            time2, observed1 = get_observed(
                file, start4, length, margin=10, zero_start=zero_start)
            observed1 = np.gradient(observed1, dt) if derivative\
            else observed1
            obs_times = [time1, time2]
            syn_times = [time0, time0]
            obs_waveforms = [observed0, observed1]
            axes2 = plot_waveforms(axes2, obs_times, obs_waveforms, color='black')
            axes2 = plot_waveforms(axes2, syn_times, syn_waveforms, color='red',
                                   custom='fill')
            axes[1].axvline(0)
            axes[1].set_title('After Shift')
            name_file = os.path.join(
                plot_folder, '{}_{}.png'.format(name, component))
            plt.savefig(name_file)
            plt.close(fig)

        if zero_start:
            stream[0].data = stream[0].data - new_baseline
            stream.write(file['file'], format='SAC', byteorder=0)
    return files


def shift_match2(data_type, plot=False, method='full', zero_start=True, event=None):
    """We shift synthetic data to maximize cross-correlation with observed data.
    This can hadle the case when one needs to shift wavefomrs both in displacement
    and velocity at the same time.

    :param data_type: list of data types to be used in modelling.
    :param plot: whether to plot results of waveform shift
    :param method: method for cross-correlation
    :type data_type: list
    :type tensor_info: bool, optional
    :type method: string, optional
    """
    if data_type == 'tele':
        json_file = 'tele_waves.json'
    if data_type == 'strong':
        json_file = 'strong_motion_waves.json'
    if data_type == 'cgps':
        json_file = 'cgps_waves.json'
    if data_type == 'surf':
        json_file = 'surf_waves.json'
    files = json.load(open(json_file))
    synthetics_file = 'synthetics_body.txt' if data_type == 'tele' else 'synthetics_strong.txt'\
        if data_type == 'strong' else 'synthetics_surf.txt' if data_type == 'surf'\
        else 'synm.cgps'

    dt = float(files[0]['dt'])
    plot_folder = 'tele_shift' if data_type == 'tele' else 'strong_shift'\
        if data_type == 'strong' else 'surf_shift' if data_type == 'surf'\
        else 'cgps_shift'
    files = get_outputs.get_data_dict(files, syn_file=synthetics_file)
    if event is not None:
        files = select_waveforms_event(files, event)
        plot_folder = '{}_event{}'.format(plot_folder, event)
    if not os.path.isdir(plot_folder):
        os.mkdir(plot_folder)
    used_channels = []
    for file in files:
        name = file['name']
        channel = file['component']
        this_channel = [name, channel]
        if not this_channel in used_channels:
            used_channels = used_channels + [this_channel]
    for station, channel in used_channels:
        files2 = [file for file in files if file['name'] == station]
        files2 = [file for file in files2 if file['component'] == channel]
        synthetics = [file['synthetic'] for file in files2]
        derivatives = [False if not file['derivative'] else file['derivative']\
                       for file in files2]
        titles = ['Disp' if not derivative else 'Vel' for derivative in derivatives]
        streams = [read(file['file']) for file in files2]
        streams = [stream.differentiate() if derivative else stream\
                   for stream, derivative in zip(streams, derivatives)]
        waveforms = [stream[0].data for stream in streams]
        nshift = int(5 / dt) if data_type == 'tele' else int(5 / dt)\
            if data_type == 'strong' else int(12 / dt) if data_type == 'surf'\
            else int(4 / dt)
        lengths = [int(float(file['duration'])) for file in files2]
        length = np.min(np.array(lengths))
        start = int(files2[0]['start_signal'])
        tr_shift = _shift2(waveforms, synthetics, nshift, length, start)
        for file in files2:
            file['start_signal'] = start + tr_shift

        for file in files2:
            file['synthetic'] = []
            file['observed'] = []
        if plot:
            length2 = int(10 / dt)
            start0 = 0
            start00 = 0
            synthetics = [synthetic[:length] for synthetic in synthetics]
            fig, axes = plt.subplots(2, len(synthetics), figsize=(10, 5))
            fig.text(0.04, 0.6, 'Before Shift', va='center', rotation='vertical')
            fig.text(0.04, 0.3, 'After Shift', va='center', rotation='vertical')
            fig.suptitle('{} {}'.format(station, channel))
            obs_times = []
            syn_times = []
            obs_waveforms = []
            syn_waveforms = []
            if len(files2) > 1:
                zipped = zip(files2, synthetics, axes[0, :], titles)
                for file, synthetic, ax, title in zipped:
                    time1, observed0 = get_observed(
                        file, start, length, margin=10)
                    time0 = np.arange(len(synthetic)) * dt
                    obs_times = obs_times + [time1]
                    syn_times = syn_times + [time0]
                    obs_waveforms = obs_waveforms + [observed0]
                    syn_waveforms = syn_waveforms + [synthetic]
                    ax.axvline(0)
                    ax.set_title(title)
                zipped = zip(files2, synthetics, axes[1, :])
                for file, synthetic, ax in zipped:
                    start4 = start + tr_shift
                    time2, observed1 = get_observed(
                        file, start4, length, margin=10, zero_start=zero_start)
                    time0 = np.arange(len(synthetic)) * dt
                    obs_times = obs_times + [time2]
                    syn_times = syn_times + [time0]
                    obs_waveforms = obs_waveforms + [observed1]
                    syn_waveforms = syn_waveforms + [synthetic]
                    ax.axvline(0)
            else:
                file = files2[0]
                synthetic = synthetics[0]
                time1, observed0 = get_observed(file, start, length, margin=10)
                time0 = np.arange(len(synthetic)) * dt
                min_val = np.minimum(np.min(observed0), np.min(synthetic))
                max_val = np.maximum(np.max(observed0), np.max(synthetic))
                axes[0].set_title(file['component'])
                start4 = start + tr_shift
                time2, observed1 = get_observed(
                    file, start4, length, margin=10, zero_start=zero_start)
                time0 = np.arange(len(synthetic)) * dt
                obs_times = [time1, time2]
                syn_times = [time0, time0]
                obs_waveforms = [observed0, observed1]
                syn_waveforms = [synthetic, synthetic]
                axes[1].axvline(0)
                axes[1].set_title(file['component'])
            axes2 = axes.ravel()
            axes2 = plot_waveforms(axes2, obs_times, obs_waveforms, color='black')
            axes2 = plot_waveforms(axes2, syn_times, syn_waveforms, color='red',
                                   custom='fill')
            name_file = os.path.join(
                plot_folder, '{}_{}.png'.format(station, channel))
            plt.savefig(name_file)
            plt.close(fig)

        if zero_start:
            stream = read(files2[0]['file'])
            new_baseline = stream[0].data[start + tr_shift]
            stream[0].data = stream[0].data - new_baseline
            stream.write(file['file'], format='SAC', byteorder=0)
    return files


def shift_match_regional(data_type, plot=False, method='full', zero_start=True, event=None):
    """We shift regional synthetic data which satisfies pareto-optimal
    cross-correlation with observed data, for a given station.

    :param data_type: list of data types to be used in modelling.
    :param plot: whether to plot results of waveform shift
    :param method: method for cross-correlation
    :type data_type: list
    :type tensor_info: bool, optional
    :type method: string, optional
    """
    if data_type == 'strong':
        json_file = 'strong_motion_waves.json'
    if data_type == 'cgps':
        json_file = 'cgps_waves.json'
    files = json.load(open(json_file))
    synthetics_file = 'synthetics_strong.txt'\
    if data_type == 'strong' else 'synthetics_cgps.txt'

    dt = float(files[0]['dt'])
    plot_folder = 'strong_shift' if data_type == 'strong' else 'cgps_shift'
    files = get_outputs.get_data_dict(files, syn_file=synthetics_file)
    if event is not None:
        files = select_waveforms_event(files, event)
        plot_folder = '{}_event{}'.format(plot_folder, event)
    if not os.path.isdir(plot_folder):
        os.mkdir(plot_folder)
    stations = [file['name'] for file in files]
    stations = list(set(stations))
    for station in stations:
        files2 = [file for file in files if file['name'] == station]
        components = [file['component'] for file in files2]
        synthetics = [file['synthetic'] for file in files2]
        streams = [read(file['file']) for file in files2]
        waveforms = [stream[0].data for stream in streams]
        nshift = int(5 / dt) if data_type == 'strong' else int(4 / dt)
        lengths = [int(float(file['duration'])) for file in files2]
        length = np.min(np.array(lengths))
        start = int(files2[0]['start_signal'])
        tr_shift = _shift2(waveforms, synthetics, nshift, length, start)
        for file in files2:
            file['start_signal'] = start + tr_shift
        # new_baselines = [st[0].data[start + tr_shift] for st in streams]

        for file in files2:
            file['synthetic'] = []
            file['observed'] = []
        if plot:
            length2 = int(10 / dt)
            start0 = 0
            start00 = 0
            fig, axes = plt.subplots(2, len(synthetics), figsize=(30, 10))
            fig.text(0.04, 0.6, 'Before Shift', va='center', rotation='vertical')
            fig.text(0.04, 0.3, 'After Shift', va='center', rotation='vertical')
            fig.suptitle(station)
            if len(files2) > 1:
                zipped = zip(files2, synthetics, axes[0, :])
                for file, synthetic, ax in zipped:
                    time1, observed0 = get_observed(
                        file, start, length, margin=10)
                    time0 = np.arange(len(synthetic)) * dt
                    ax.plot(time1, observed0)
                    ax.plot(time0, synthetic, 'r')
                    ax.axvline(0)
                    ax.set_title(file['component'])
                zipped = zip(files2, synthetics, axes[1, :])
                for file, synthetic, ax in zipped:
                    start4 = start + tr_shift
                    time2, observed1 = get_observed(
                        file, start4, length, margin=10, zero_start=zero_start)
                    time0 = np.arange(len(synthetic)) * dt
                    ax.plot(time2, observed1)
                    ax.plot(time0, synthetic, 'r')
                    ax.axvline(0)
            else:
                file = files2[0]
                synthetic = synthetics[0]
                time1, observed0 = get_observed(file, start, length, margin=10)
                time0 = np.arange(len(synthetic)) * dt
                min_val = np.minimum(np.min(observed0), np.min(synthetic))
                max_val = np.maximum(np.max(observed0), np.max(synthetic))
                axes[0].plot(time1, observed0)
                axes[0].plot(time0, synthetic, 'r')
                axes[0].vlines(0, min_val, max_val)
                axes[0].set_title(file['component'])
                start4 = start + tr_shift
                time2, observed1 = get_observed(
                    file, start4, length, margin=10, zero_start=zero_start)
                time0 = np.arange(len(synthetic)) * dt
                axes[1].plot(time2, observed1)
                axes[1].plot(time0, synthetic, 'r')
                axes[1].axvline(0)
                axes[1].set_title(file['component'])
            name_file = os.path.join(plot_folder, '{}.png'.format(station))
            plt.savefig(name_file)
            plt.close(fig)

        if zero_start:
            for stream in streams:
                new_baseline = stream[0].data[start + tr_shift]
                stream[0].data = stream[0].data - new_baseline
                stream.write(file['file'], format='SAC', byteorder=0)
    return files


def save_waveforms(data_type, files):
    """
    """
    if data_type == 'tele':
        json_file = 'tele_waves.json'
    if data_type == 'strong':
        json_file = 'strong_motion_waves.json'
    if data_type == 'cgps':
        json_file = 'cgps_waves.json'
    if data_type == 'surf':
        json_file = 'surf_waves.json'
    with open(json_file,'w') as f:
         json.dump(
             files, f, sort_keys=True, indent=4,
             separators=(',', ': '), ensure_ascii=False)


def get_observed(file_dict, start, length, margin=10, zero_start=False):
    """
    """
    stream = read(file_dict['file'])
    dt = stream[0].stats.delta
    length2 = int(margin / dt)
    waveform = stream[0].data
    new_baseline = waveform[start]
    if start >= 0:
        start2 = max(0, start - length2)
        start3 = start - start2
        waveform1 = np.array(
            [val for i, val in enumerate(waveform[start2:]) if i < length])
    else:
        waveform1 = np.array(
            [val for i, val in enumerate(waveform) if i < length])
        waveform1 = np.concatenate(([0] * -start, waveform1))
    time = np.arange(-start3, len(waveform1) - start3) * dt
    if zero_start:
        waveform1 = waveform1 - new_baseline
    derivative = False if not 'derivative' in file_dict else file_dict['derivative']
    if derivative:
        waveform1 = np.gradient(waveform1, dt)
    return time, waveform1


def _shift(obser_tr, syn_tr, nshift, length, start_pos):
    """Routine for shifting an observed waveform to maximize the
    cross-correlation to a synthetic waveform.
    """
    err_max = 0
    j_min = 0
    synthetic = syn_tr[:length]
    for j in range(-nshift, nshift + 1):
        start2 = start_pos + j
        exy = 0
        if start2 < 0:
            continue
        else:
            observed = np.array(
                [val for i, val in enumerate(obser_tr[start_pos + j:])\
                if i < length])
            synthetic2 = synthetic[:len(observed)]
            exy = np.sum(observed * synthetic2)
        err = 2 * exy

        if err_max <= err:
            err_max = err
            j_min = j
    return j_min


def _shift2(waveforms, synthetics, nshift, length, start_pos):
    """Routine for finding the shift with pareto-optimal cross-correlation for
    all channels of a certain station.
    """
    synthetics = [synthetic[:length] for synthetic in synthetics]

    j_min = 0
    err_max = 0
    for j in range(-nshift, nshift + 1):
        start2 = start_pos + j
        err = 0
        if start2 < 0:
            continue
        else:
            zipped = zip(waveforms, synthetics)
            for i, (observed, synthetic) in enumerate(zipped):
                observed2 = np.array(
                    [val for i, val in enumerate(observed[start_pos + j:])\
                    if i < length])
                synthetic2 = synthetic[:len(observed2)]
                err = err + 2*np.sum(observed2 * synthetic2)

        if err_max <= err:
            err_max = err
            j_min = j

    return j_min


def _shift3(waveforms, synthetics, nshift, length, start_pos):
    """Routine for finding the shift with pareto-optimal cross-correlation for
    all channels of a certain station.
    """
    synthetics = [synthetic[:length] for synthetic in synthetics]

    cross_corr = [[]] * (2 * nshift + 1)
    for j in range(-nshift, nshift + 1):
        start2 = start_pos + j
        err = np.zeros(len(synthetics))
        if start2 < 0:
            continue
        else:
            zipped = zip(waveforms, synthetics)
            for i, (observed, synthetic) in enumerate(zipped):
                observed2 = np.array(
                    [val for i, val in enumerate(observed[start_pos + j:])\
                    if i < length])
                synthetic2 = synthetic[:len(observed2)]
                err[i] = 2 * np.sum(observed2 * synthetic2)
        cross_corr[j + nshift] = err

    pareto_set = []
    j_min = 0
    for j in range(-nshift, nshift + 1):
        candidate = cross_corr[j + nshift]
        pareto_max = True
        for i in range(-nshift, nshift + 1):
            err = cross_corr[i + nshift]
            print(err)
            print(candidate)
            if np.any(err > candidate):
                pareto_max = False
                break
        if pareto_max:
            pareto_set = pareto_set + [j]
    print(pareto_set)
    if pareto_set:
        j_min = pareto_set[0]
    return j_min


def _print_arrival(tensor_info):
    """
    """
    date_origin = tensor_info['date_origin']
    other_files = glob.glob('../data/*BHZ*sac')
    files = json.load(open('tele_waves.json'))
    plot_folder = 'tele_arrival'
    if not os.path.isdir(plot_folder):
        os.mkdir(plot_folder)
    for file in files:
        start = file['start_signal']
        sac = file['file']
        stream = read(sac)
        trace = stream[0]
        delta = trace.stats.delta
        stat = trace.stats.station
        chan = trace.stats.channel
        if not chan == 'BHZ':
            continue
        other_file = [v for v in other_files if stat in v]
        other_stream = read(other_file[0])
        trace2 = other_stream[0]
        if not 't1' in trace2.stats.sac:
            continue
        arrival = trace2.stats.sac['t1']
        start = date_origin + arrival - 15
        end = date_origin + arrival + 15
        trace2.trim(starttime=start, endtime=end)
        fig = plt.figure()
        data = trace2.data
        time = np.linspace(-15, 15, len(data))
        plt.plot(time, data)
        plt.title(stat)
        plt.axvline(x=0, color='r')
#        plt.xticks(time)
        fig.savefig(os.path.join(plot_folder, '{}_pick'.format(stat)))
        del fig
    return


if __name__ == '__main__':
    import argparse

    parser = argparse.ArgumentParser()
    parser.add_argument("-f", "--folder", default=os.getcwd(),
                        help="folder where there are input files")
    parser.add_argument('-o', '--option', choices=['match', 'manual'],
                        help='choose whether to shift by cross-correlation\
                        or plot to pick manually')
    parser.add_argument("-t", "--type", choices=['tele', 'strong', 'cgps', 'surf'],
                        help="type of data to apply shift")
    parser.add_argument("-p", "--plot", action="store_true",
                        help="plot or not results of shift")
    parser.add_argument("-gcmt", "--gcmt_tensor",
                        help="location of GCMT moment tensor file")
    parser.add_argument(
        "-me", "--many_events", action="store_true", help="plots for many events")
    args = parser.parse_args()
    os.chdir(args.folder)
    if args.gcmt_tensor:
        cmt_file = args.gcmt_tensor
        tensor_info = tensor.get_tensor(cmt_file=cmt_file)
    if args.option == 'match':
        plot = False if not args.plot else True
        if not args.many_events:
            if args.type in ['tele', 'surf']:
                files = shift_match2(args.type, plot=plot)
            else:
                files = shift_match_regional(args.type, plot=plot)
        else:
            if args.type in ['tele', 'surf']:
                files1 = shift_match2(args.type, plot=plot, event=1)
                files2 = shift_match2(args.type, plot=plot, event=2)
                files = files1 + files2
            else:
                files1 = shift_match_regional(args.type, plot=plot, event=1)
                files2 = shift_match_regional(args.type, plot=plot, event=2)
                files = files1 + files2
        save_waveforms(args.type, files)
    else:
        _print_arrival(tensor_info)

