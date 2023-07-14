import sys
import os
import numpy as np
import json
#import matplotlib
#matplotlib.use('Qt5Agg')
from matplotlib import pyplot as plt
from matplotlib import gridspec, ticker, patches
from scipy.signal import butter, filtfilt
from obspy import read

plt.rc('axes', titlesize=14)
plt.rc('axes', labelsize=12)
plt.rc('xtick', labelsize=12)
plt.rc('ytick',labelsize=12)
plt.rc('font', size=12)

def plot_waveforms(axes, times, waveforms, weights, type_str=None, comp=None, color='blue', custom=None):
    """
    """
    nPlot=0
    for ax, time, waveform, weight in zip(axes, times, waveforms, weights):
        nPlot+=1
        if weight == 0.0:
            ax.plot(time, waveform, color=color, linewidth=2, linestyle='dashed')
        else:
            ax.plot(time, waveform, color=color, linewidth=2*weight)
        min_time, max_time = ax.get_xlim()
        min_time = np.minimum(np.min(time), min_time)
        max_time = np.maximum(np.max(time), max_time)
        if custom == 'fill':
            min_val, max_val = ax.get_ylim()
            min_val = np.minimum(np.min(waveform), min_val)
            max_val = np.maximum(np.max(waveform), max_val)
            #if type_str == 'cgps' or type_str == 'strong_motion':
            #    if max_val < 1 and min_val > -1:
            #        max_val = 1
            #        min_val = -1
            ax.set_ylim([-(max(abs(min_val), max_val)),max(abs(min_val), max_val)])
            min_val, max_val = ax.get_ylim()
            ax.vlines(0, min_val, max_val,'k', lw=1)
            if type_str == 'tele_body':
                ax.text(np.max(time), 0.6*max_val, 
                   '{:0.1f}'.format(max(abs(min_val),max_val)),ha='right',va='center')
                ax.hlines(0, -20, np.max(time), 'k', lw=1)
                ax.set_xlim([-20, np.max(time)])
            elif type_str == 'surf_tele':
                ax.text(np.max(time), 0.6*max_val, 
                   '{:0.3f}'.format(max(abs(min_val),max_val)),ha='right',va='center')
                ax.hlines(0, -350, np.max(time), 'k', lw=1)
                #ax.set_xlim([-350,3000])
                ax.set_xlim([-350, np.max(time)])
            elif type_str == 'cgps' or type_str == 'strong_motion':
                min_wval = np.min(waveform)
                max_wval = np.max(waveform)
                if max_wval > abs(min_wval):
                    ax.text(np.max(time), 0.6*max_val,
                       '{:0.2f}'.format(max_wval),ha='right',va='center')
                else:
                    ax.text(np.max(time), 0.6*max_val,
                       '{:0.2f}'.format(min_wval),ha='right',va='center')
                ax.hlines(0, -15, np.max(time), 'k', lw=1)
                ax.set_xlim([-15,np.max(time)])
            min_time, max_time = ax.get_xlim()
            if type_str == 'tele_body' and comp == 'BHZ':
                ax.text(1.1*min_time,0.2*max(abs(min_val),max_val),'P',ha='right',va='bottom')
            if type_str == 'tele_body' and comp == 'SH':
                ax.text(1.1*min_time,0.2*max(abs(min_val),max_val),'SH',ha='right',va='bottom')
            if type_str == 'surf_tele' and comp == 'BHZ':
                ax.text(1.2*min_time,0.2*max(abs(min_val),max_val),'Z',ha='right',va='bottom')
            if type_str == 'surf_tele' and comp == 'SH':
                ax.text(1.2*min_time,0.2*max(abs(min_val),max_val),'T',ha='right',va='bottom')
            #if type_str == 'cgps' or type_str == 'strong_motion':
            #    ax.text(1.4*min_time, 0.2*max(abs(min_val),max_val),None,ha='right',va='bottom')
        if custom == 'syn':
            max_val = np.maximum(abs(min(waveform)),max(waveform))
            tmin, tmax = ax.get_xlim()
            ymin, ymax = ax.get_ylim()
            if type_str == 'tele_body':
                ax.text(tmax, 0.6*ymin,
                    '{:0.1f}'.format(max_val),ha='right',va='center',color='red')
            elif type_str == 'surf_tele':
                ax.text(tmax, 0.6*ymin,
                    '{:0.3f}'.format(max_val),ha='right',va='center',color='red')
            elif type_str == 'cgps' or type_str == 'strong_motion':
                ax.text(tmax, 0.6*ymin,
                    '{:0.2f}'.format(max_val),ha='right',va='center',color='red')
        #ax.xaxis.set_major_locator(
        #    ticker.MaxNLocator(nbins=3, min_n_ticks=3))
        #ax.yaxis.set_major_locator(
        #    ticker.MaxNLocator(nbins=3, min_n_ticks=3))
        if type_str == 'tele_body':
            ax.xaxis.set_major_locator(
                 ticker.MaxNLocator(nbins=5, min_n_ticks=5))
            ax.yaxis.set_major_locator(
                 ticker.NullLocator())
            ax.xaxis.set_minor_locator(
                 ticker.MultipleLocator(5))
        elif type_str == 'surf_tele':
            ax.xaxis.set_major_locator(
                 ticker.MultipleLocator(500))
            ax.yaxis.set_major_locator(
                 ticker.NullLocator())
        elif type_str == 'cgps' or type_str == 'strong_motion':
            ax.xaxis.set_major_locator(
                 ticker.MultipleLocator(40))
            ax.yaxis.get_major_locator().set_params(integer=True)
            #ax.yaxis.tick_right()
            #ax.yaxis.set_major_locator(
            #     ticker.NullLocator())
            #ax.xaxis.set_minor_locator(
            #     ticker.MultipleLocator(500))
        if nPlot > len(weights) - 3:
            ax.set_xlabel('Time After OT (s)')
        ax.grid(axis='x',which='both',linestyle='dotted', color='0.5')
    return axes


def add_metadata(axes, **kwargs):
    """
    """
    if 'type_str' in kwargs:
        if kwargs['type_str'] == 'cgps' or kwargs['type_str'] == 'strong_motion':
#            if 'names' in kwargs:
#                for ax, name in zip(axes, kwargs['names']):
#                    ax.text(-0.07, 0.50, name, ha='right', va='center', transform=ax.transAxes)
            if 'distances' in kwargs:
               for ax, dist in zip(axes, kwargs['distances']):
                   ax.text(0.01, 0.46, '{:0.2f}'.format(dist), ha='left',
                           va='top', transform=ax.transAxes)
            if 'comps' in kwargs:
               if 'names' in kwargs:
                   for ax, comp, name in zip(axes, kwargs['comps'], kwargs['names']):
                       if comp[-1] == 'E':
                           ax.text(-0.18, 0.5, name, ha='right', va='center', transform=ax.transAxes, rotation=90, fontweight='bold')
               for ax, comp in zip(axes, kwargs['comps']):
                   ax.text(-0.13, 0.5, comp, ha='right', va='center', transform=ax.transAxes, rotation=90)
        else:
            if 'names' in kwargs:
                for ax, name in zip(axes, kwargs['names']):
                    ax.text(-0.02, 0.50, name, ha='right', va='center', transform=ax.transAxes)
            if 'distances' in kwargs:
                for ax, dist in zip(axes, kwargs['distances']):
                    ax.text(
                        0.01, 0.46, '{:0.0f}'.format(dist), ha='left',
                        va='top', transform=ax.transAxes)
            #if 'comps' in kwargs:
            #    for ax, comp in zip(axes, kwargs['comps']):
            #        ax.text(-0.02, 0.40, comp, ha='right', va='center', transform=ax.transAxes)
    if 'azimuths' in kwargs:
        for ax, az in zip(axes, kwargs['azimuths']):
            ax.text(
                0.01, 0.5, '{:0.0f}'.format(az), ha='left',
                va='bottom', transform=ax.transAxes)
    if 'weights' in kwargs:
        for ax, weight in zip(axes, kwargs['weights']):
            alpha = 1 if weight > 0 else 0.25
            lines = ax.get_lines()
            for line in lines:
                line.set_alpha(alpha)
    return axes


def plot_waveform_fits(files, components, type_str, start_margin=10,
                       test=False, forward=False):
    """
    """
    if type_str == 'tele_body' or type_str == 'surf_tele':
        files = [file for file in files if file['component'] in components]
        print('Creating Waveform Fit Plot: ' + str(type_str) + ' ' + str(components[0]))
    if type_str == 'cgps' or type_str == 'strong_motion':
        files = [file for file in files]
        print('Creating Waveform Fit Plot: ' + str(type_str))
    files = sorted(files, key=lambda k: (k['azimuth'], k['component']))
    sampling = [file['dt'] for file in files]
    names = [file['name'] for file in files]
    azimuths = [file['azimuth'] for file in files]
    distances = [file['distance'] for file in files]
    weights = [file['trace_weight'] for file in files]
    obs_waveforms = [file['observed'] for file in files]
    syn_waveforms = [file['synthetic'] for file in files]
    comp = [file['component'] for file in files]
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
    numrows_phase = len(files) // 3 + 1
    fig, axes = plt.subplots(max(3, numrows_phase), 3, figsize=(20, int(2.2*numrows_phase)))
    axes2 = axes.ravel()
    for ax in axes2[len(files):]:
        ax.axis('off')

    if type_str == 'tele_body' or type_str == 'surf_tele':
        axes2 = plot_waveforms(axes2, obs_times, obs_waveforms, weights, type_str=type_str, 
                comp=comp[0], color='black', custom='fill')
        axes2 = plot_waveforms(axes2, syn_times, syn_waveforms, weights, type_str=type_str,
                comp=comp[0], color='red', custom='syn')
    if type_str == 'cgps' or type_str == 'strong_motion':
        axes2 = plot_waveforms(axes2, obs_times, obs_waveforms, weights, type_str=type_str,
                comp=comp, color='black', custom='fill')
        axes2 = plot_waveforms(axes2, syn_times, syn_waveforms, weights, type_str=type_str,
                comp=comp, color='red', custom='syn')
   
    
    dict = {
        'weights': weights,
        'azimuths': azimuths,
        'names': names,
        'distances': distances,
        'type_str': type_str,
        'comps':comp
    }
    axes2 = add_metadata(axes2, **dict)

    if type_str == 'tele_body':
        if 'BHZ' in components: plot_name = 'P_body_waves'
        if 'SH' in components: plot_name = 'SH_body_waves'

    if type_str == 'surf_tele':
        if 'BHZ' in components: plot_name = 'Rayleigh_surf_waves'
        if 'SH' in components: plot_name = 'Love_surf_waves'

    if type_str == 'cgps':
        plot_name = 'cGPS_waves'

    if type_str == 'strong_motion':
        plot_name = 'strong_motion_waves'

    plt.savefig(plot_name+'.png', dpi=300)#bbox_inches='tight')
    plt.savefig(plot_name+'.ps')
    plt.close()
    return


def filt_waveform(file, high_freq):
    """
    """
    stream = read(file['file'])
    data = stream[0][1000:]
    high_freq = high_freq / 25
    b, a = butter(2, high_freq, btype='lowpass')
    filt_data = filtfilt(b, a, data)
    file['synthetic'] = filt_data
    return file


def plot_spectra(files, dt):
    """
    """
    for file in files:
        waveform = file['synthetic']
        name = file['name']
        comp = file['component']
        fft = np.fft.fft(waveform)
        n = len(waveform)
        freq = np.fft.fftfreq(n, d=dt)
        plt.loglog(freq[:n//2], np.abs(fft[:n//2]))
        plt.title('{} {}'.format(name, comp))
        plt.savefig('spectra_{}_{}'.format(name, comp))
        plt.close()


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
