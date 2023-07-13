import sys
import os
import numpy as np
import json
import time
#import matplotlib
#matplotlib.use('Qt5Agg')
from matplotlib import pyplot as plt
from matplotlib import gridspec, ticker, patches
from scipy.signal import butter, filtfilt
from obspy import read


def plot_waveforms(axes, times, waveforms, color='blue', custom=None):
    """Method to plot some (waveform) time series.

    :param axes: where to plot the waveform time series
    :param times: times of data to be plotted
    :param waveforms: waveforms to plot
    :param color: color of the plotted waveform
    :param custom: fancier plotting methods
    :type axes: list
    :type times: list
    :type waveforms: list
    :type color: string
    :type custom: string
    """
    for ax, time, waveform in zip(axes, times, waveforms):
        if waveform is None:
            continue
        if len(waveform) == 0:
            continue
        ax.plot(time, waveform, color=color, linewidth=0.6)
        min_time, max_time = ax.get_xlim()
        min_time = np.minimum(np.min(time), min_time)
        max_time = np.maximum(np.max(time), max_time)
        ax.set_xlim([min_time, max_time])
        if custom == 'fill':
            min_val, max_val = ax.get_ylim()
            min_val = np.minimum(np.min(waveform), min_val)
            max_val = np.maximum(np.max(waveform), max_val)
            ax.set_ylim([min_val, max_val])
            ax.vlines(0, min_val, max_val)
        ax.xaxis.set_major_locator(
            ticker.MaxNLocator(nbins=3, min_n_ticks=3))
        ax.yaxis.set_major_locator(
            ticker.MaxNLocator(nbins=3, min_n_ticks=3))
    return axes


def add_metadata(axes, **kwargs):
    """Add metadata to some axes.
    """
    if 'names' in kwargs:
        for ax, name in zip(axes, kwargs['names']):
            if name is None:
                continue
            ax.text(
                0.9, 0.9, name, ha='center', va='center', transform=ax.transAxes)
    if 'distances' in kwargs:
        for ax, dist in zip(axes, kwargs['distances']):
            if dist is None:
                continue
            ax.text(
                0.1, 0.1, '{:0.1f}'.format(dist), ha='center',
                va='center', transform=ax.transAxes)
    if 'azimuths' in kwargs:
        for ax, az in zip(axes, kwargs['azimuths']):
            if az is None:
                continue
            ax.text(
                0.1, 0.9, '{:0.1f}'.format(az), ha='center',
                va='center', transform=ax.transAxes)
    if 'weights' in kwargs:
        for ax, weight in zip(axes, kwargs['weights']):
            if weight is None:
                continue
            alpha = 1 if weight > 0 else 0.1
            lines = ax.get_lines()
            for line in lines:
                line.set_alpha(alpha)
    return axes


def plot_waveform_fits(files, components, type_str, start_margin=10,
                       test=False, event=None):
    """Plot fit of observed to synthetic data for selected channels.

    :param files: waveform files to plot
    :param components: components (channels) of data selected for plotting
    :param type_str: data type of given waveform files
    :param start_margin: start margin of data for plotting
    :type files: list
    :type components: list
    :type type_str: string
    :type start_margin: float, optinoal
    """
    files = [file for file in files if file['component'] in components]
    files = sorted(files, key=lambda k: k['azimuth'])
    sampling = [file['dt'] for file in files]
    names = [file['name'] for file in files]
    azimuths = [file['azimuth'] for file in files]
    distances = [file['distance'] for file in files]
    weights = [file['trace_weight'] for file in files]
    obs_waveforms = [file['observed'] for file in files]
    syn_waveforms = [file['synthetic'] for file in files]
    zipped = zip(obs_waveforms, syn_waveforms)
    syn_waveforms = [syn_waveform[:len(obs_waveform)]\
                     for obs_waveform, syn_waveform in zipped]
    zipped = zip(sampling, syn_waveforms)
    syn_times = [dt * np.arange(0, len(synthetic)) for dt, synthetic in zipped]
    start_waveform = []
    for file in files:
        dt = file['dt']
        nstart = file['start_signal']
        margin = int(start_margin / dt)
        margin = min(nstart, margin)
        start_waveform = start_waveform + [margin]
    zipped = zip(sampling, start_waveform, obs_waveforms)
    obs_times = [dt * np.arange(-start, len(observed) - start)\
                 for dt, start, observed in zipped]
    numrows_phase = len(files) // 4 + 1
    fig, axes = plt.subplots(max(4, numrows_phase), 4, figsize=(13, 9))
    axes2 = axes.ravel()
    for ax in axes2[len(files):]:
        ax.axis('off')
    obs_waveforms = (waveform for waveform in obs_waveforms)
    syn_waveforms = (waveform for waveform in syn_waveforms)

    axes2 = plot_waveforms(axes2, obs_times, obs_waveforms, color='black')
    axes2 = plot_waveforms(axes2, syn_times, syn_waveforms, color='red',
                           custom='fill')
    dict = {
        'weights': weights,
        'azimuths': azimuths,
        'names': names,
        'distances': distances
    }
    axes2 = add_metadata(axes2, **dict)

    if type_str == 'cgps':
        if 'LXZ' in components: plot_name = 'LXZ_cgps_waves'
        if 'LXN' in components: plot_name = 'LXN_cgps_waves'
        if 'LXE' in components: plot_name = 'LXE_cgps_waves'

    if type_str == 'strong_motion':
        if 'HNZ' in components: plot_name = 'HNZ_strong_motion_waves'
        if 'HNN' in components: plot_name = 'HNN_strong_motion_waves'
        if 'HNE' in components: plot_name = 'HNE_strong_motion_waves'

    if type_str == 'tele_body':
        if 'BHZ' in components: plot_name = 'P_body_waves'
        if 'SH' in components: plot_name = 'SH_body_waves'

    if type_str == 'surf_tele':
        if 'BHZ' in components: plot_name = 'Rayleigh_surf_waves'
        if 'SH' in components: plot_name = 'Love_surf_waves'

    if type_str == 'dart':
        plot_name = 'DART_waves.png'

    if event is not None:
        plot_name = '{}_event{}'.format(plot_name, event)
    plt.savefig(plot_name, bbox_inches='tight')
    plt.close()
    return


# def filt_waveform(file, high_freq):
#     """
#     """
#     stream = read(file['file'])
#     data = stream[0][1000:]
#     high_freq = high_freq / 25
#     b, a = butter(2, high_freq, btype='lowpass')
#     filt_data = filtfilt(b, a, data)
#     file['synthetic'] = filt_data
#     return file


# def plot_spectra(files, dt):
#     """
#     """
#     for file in files:
#         waveform = file['synthetic']
#         name = file['name']
#         comp = file['component']
#         fft = np.fft.fft(waveform)
#         n = len(waveform)
#         freq = np.fft.fftfreq(n, d=dt)
#         plt.loglog(freq[:n//2], np.abs(fft[:n//2]))
#         plt.title('{} {}'.format(name, comp))
#         plt.savefig('spectra_{}_{}'.format(name, comp))
#         plt.close()


if __name__ == '__main__':
    files = [
        {
            'file': '/home/pkoch/folder_plot16/STR.MT07.HNE.C1.ACC',
            'name': 'MT07',
            'component': 'HNE',
            'synthetic': []
        },
        {
            'file': '/home/pkoch/folder_plot16/STR.MT07.HNN.C1.ACC',
            'name': 'MT07',
            'component': 'HNN',
            'synthetic': []
        },
        {
            'file': '/home/pkoch/folder_plot16/STR.MT07.HNZ.C1.ACC',
            'name': 'MT07',
            'component': 'HNZ',
            'synthetic': []
        },
    ]
    high_freq = 15
    for file in files:
        file = filt_waveform(file, high_freq)
    plot_spectra(files, 0.01)
    fig, axes = plt.subplots(1, 3, figsize=(15, 5))
    waveforms = [file['synthetic'] for file in files]
    times = [np.arange(len(waveform)) * 0.01 for waveform in waveforms]
    axes = plot_waveforms(axes, times, waveforms)
    axes[0].set_title('N')
    axes[1].set_title('E')
    axes[2].set_title('Z')
    plot_name = 'MT07_lowpass_{}'.format(high_freq)
    plt.savefig(plot_name, bbox_inches='tight')
    print(1)
