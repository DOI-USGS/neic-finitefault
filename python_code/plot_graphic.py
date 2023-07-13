# -*- coding: utf-8 -*-
"""The routines here allow to plot the solution model FFM modelling, as well
as moment rate function, and waveform fits.
"""


import argparse
from matplotlib import pyplot as plt
from matplotlib import gridspec, ticker, patches, colors
import cartopy.crs as ccrs
import cartopy.io.shapereader as shpreader
import cartopy.feature as cf
from obspy.imaging.beachball import beach
import numpy as np
import os
import get_outputs
import load_ffm_model
import json
import errno
from shutil import copy2, move
import glob
#from clawpack.geoclaw import dtopotools
#
# local modules
#
import fault_plane as pf
import velocity_models as mv
import plane_management as pl_mng
import seismic_tensor as tensor
from waveform_plots import plot_waveform_fits
from plot_maps import plot_map, set_map_cartopy, plot_borders


def plot_ffm_sol(tensor_info, segments_data, point_sources, shear, solution,
                 vel_model, default_dirs, use_waveforms=True, event=None):
    """Main routine. Allows to coordinate execution of different plotting
    routines.

    :param tensor_info: dictionary with moment tensor information
    :param segments: list of dictionaries with properties of fault segments
    :param point_sources: properties of point sources of the fault plane
    :param shear: values of shear modulous for each subfault
    :param solution: dictionary with output kinematic model properties
    :param vel_model: dictionary with velocity model properties
    :param default_dirs: dictionary with default directories to be used
    :type default_dirs: dict
    :type tensor_info: dict
    :type segments: list
    :type point_sources: array
    :type shear: array
    :type solution: dict
    :type vel_model: dict

    .. rubric:: Example:

    First, we load necessary modules.

    >>> import json
    >>> import get_outputs # Allows us to get properties of inverted model
    >>> import management as mng # Allows us to load location of plotting files
    >>> import fault_plane as pf
    >>> import plane_management as pl_mng
    Next, we load necessary data for plots.
    >>> vel_model = json.load(open('velmodel_data.json')) # Assume vel_model stored in file 'velmodel_data.json'
    >>> segments, rise_time, point_sources = pl_mng.__read_planes_info() # Loads point sources and segments information
    >>> solution = get_outputs.read_solution_static_format(segments, point_sources)
    >>> shear = pf.shear_modulous(point_sources, velmodel=vel_model)
    >>> tensor_info = {
            'moment_mag': 7 * 10 ** 27,
            'date_origin': UTCDateTime(2014, 04, 01, 23, 46, 47)
            'lat': -19.5,
            'lon': -70.5,
            'depth': 25,
            'time_shift': 44,
            'half_duration': 40,
            'centroid_lat': -21,
            'centroid_lon': -70,
            'centroid_depth': 35
        }
    Next, we plot solution
    >>> default_dirs = mng.default_dirs()
    >>> plot_ffm_sol(tensor_info, segments, point_sources, shear, solution,
    >>>              vel_model, default_dirs)

    .. note::

        To plot the results of the FFM modelling, we need to run this code
        in a folder whih contains files Solucion.txt, Fault.time, Fault.pos,
        Event_mult.in, and some among the files synm.tele, synm.str_low,
        synm.str and synm.cgps.

    .. note::

        When running this code manually, it is good idea to check if
        the information among files Solucion.txt, Fault.pos, Fault.time,
        and Event_mult.in is consistent.

    """
    segments = segments_data['segments']
    _plot_vel_model(vel_model, point_sources)
    if use_waveforms:
        _plot_moment_rate_function(segments_data, shear, point_sources, event=event)
        _PlotRiseTime(segments, point_sources, solution)
        _PlotRuptTime(segments, point_sources, solution)
    _PlotSlipDistribution(segments, point_sources, solution)
    _PlotMap(tensor_info, segments, point_sources, solution, default_dirs, event=event)


def plot_beachballs(segments, data_type):
    """Main routine. Allows to coordinate execution of different plotting
    routines.

    :param segments: list of dictionaries with properties of fault segments
    :param data_type: list with data types used in modelling
    :type data_type: list
    :type segments: list
    """
    if 'tele_body' in data_type:
        if not os.path.isfile('tele_waves.json'):
            raise FileNotFoundError(
                errno.ENOENT, os.strerror(errno.ENOENT), 'tele_waves.json')
        traces_info = json.load(open('tele_waves.json'))
        plot_beachball(segments, files=traces_info, phase='P')
        plot_beachball(segments, files=traces_info, phase='SH')


def plot_misfit(used_data_type, forward=False, event=None):
    """Plot misfit of observed and synthetic data

    :param used_data_type: list with data types used in modelling
    :param forward: whether model is result of kinematic modelling or not
    :type used_data_type: list
    :type forward: bool, optional
    """
    from many_events import select_waveforms_event

    if 'dart' in used_data_type:
        if not os.path.isfile('dart_waves.json'):
            raise FileNotFoundError(
                errno.ENOENT, os.strerror(errno.ENOENT), 'dart_waves.json')
        traces_info = json.load(open('dart_waves.json'))
        traces_info = get_outputs.get_data_dict(
            traces_info, obs_file='waveforms_dart.txt',
            syn_file='synthetics_dart.txt')
        values = ['dart']
        for components in values:
            plot_waveform_fits(
                traces_info, components, 'dart', start_margin=0)
    if 'tele_body' in used_data_type:
        if not os.path.isfile('tele_waves.json'):
            raise FileNotFoundError(
                errno.ENOENT, os.strerror(errno.ENOENT), 'tele_waves.json')
        traces_info = json.load(open('tele_waves.json'))
        traces_info = get_outputs.get_data_dict(
            traces_info, syn_file='synthetics_body.txt')
        if event:
            traces_info = select_waveforms_event(traces_info, event)
        values = [['BHZ'], ['SH']]
        for components in values:
            plot_waveform_fits(
                traces_info, components, 'tele_body', event=event)
    if 'surf_tele' in used_data_type:
        if not os.path.isfile('surf_waves.json'):
            raise FileNotFoundError(
                errno.ENOENT, os.strerror(errno.ENOENT), 'surf_waves.json')
        traces_info = json.load(open('surf_waves.json'))
        traces_info = get_outputs.get_data_dict(
            traces_info, syn_file='synthetics_surf.txt', margin=0)
        if event:
            traces_info = select_waveforms_event(traces_info, event)
        values = [['BHZ'], ['SH']]
        for components in values:
            plot_waveform_fits(
                traces_info, components, 'surf_tele', event=event)
    if 'strong_motion' in used_data_type:
        if not os.path.isfile('strong_motion_waves.json'):
            raise FileNotFoundError(
                errno.ENOENT, os.strerror(errno.ENOENT),
                'strong_motion_waves.json')
        traces_info = json.load(open('strong_motion_waves.json'))
        traces_info = get_outputs.get_data_dict(
            traces_info, syn_file='synthetics_strong.txt')
        if event:
            traces_info = select_waveforms_event(traces_info, event)
        values = [['HLZ', 'HNZ'], ['HLE', 'HNE'], ['HLN', 'HNN']]
        for components in values:
            plot_waveform_fits(
                traces_info, components, 'strong_motion', event=event)
    if 'cgps' in used_data_type:
        if not os.path.isfile('cgps_waves.json'):
            raise FileNotFoundError(
                errno.ENOENT, os.strerror(errno.ENOENT), 'cgps_waves.json')
        traces_info = json.load(open('cgps_waves.json'))
        traces_info = get_outputs.get_data_dict(
            traces_info, syn_file='synthetics_cgps.txt')
        if event:
            traces_info = select_waveforms_event(traces_info, event)
        values = [
            ['LXZ', 'LHZ', 'LYZ'],
            ['LXE', 'LHE', 'LYE'],
            ['LXN', 'LHN', 'LYN']
        ]
        for components in values:
            plot_waveform_fits(
                traces_info, components, 'cgps', event=event)
    return


def _plot_vel_model(velmodel, point_sources):
    """We plot the seismic velocity model as a function of depth

    :param point_sources: properties of point sources of the fault plane
    :param vel_model: dictionary with velocity model properties
    :type point_sources: array
    :type vel_model: dict
    """
    max_depth = [max(ps_segment[:, :, :, :, 2].flatten())\
        for ps_segment in point_sources]
    max_depth = max(max_depth)
    p_vel = np.array(velmodel['p_vel']).astype(float)
    sh_vel = np.array(velmodel['s_vel']).astype(float)
    thick = np.array(velmodel['thick']).astype(float)

    depths = np.zeros(len(thick) + 1)

    depths[1:] = np.cumsum(thick)
    depths = np.array([depth for depth in depths if depth < 70])
    depths = np.append([depths], [70])#[max_depth])
    plt.plot((p_vel[0], p_vel[0]), (depths[0], depths[1]), 'b-', label='P')
    plt.plot((sh_vel[0], sh_vel[0]), (depths[0], depths[1]), 'r-', label='SH')
    j = len(depths) - 3#2
    for i in range(j):
        plt.plot((p_vel[i], p_vel[i]), (depths[i], depths[i + 1]), 'b-')
        plt.plot(
            (p_vel[i], p_vel[i + 1]), (depths[i + 1], depths[i + 1]), 'b-')
        plt.plot((sh_vel[i], sh_vel[i]), (depths[i], depths[i + 1]), 'r-')
        plt.plot(
            (sh_vel[i], sh_vel[i + 1]), (depths[i + 1], depths[i + 1]), 'r-')

    plt.plot((p_vel[j], p_vel[j]), (depths[j], depths[j + 1]), 'b-')
    plt.plot((sh_vel[j], sh_vel[j]), (depths[j], depths[j + 1]), 'r-')

    plt.title('Crust model for north of Chile')#'Body wave velocity model')
    plt.xlabel('Body wave velocity $(km/s)$')
    plt.ylabel('Depth $(km)$')
    plt.legend(loc='upper right')

    ax = plt.gca()
    ax.invert_yaxis()
    plt.savefig('crust_body_wave_vel_model.png', bbox_inches='tight')
    plt.close()


def _PlotRuptTime(segments, point_sources, solution):
    """We plot time distribution based on the FFM solution model

    :param segments: list of dictionaries with properties of fault segments
    :param point_sources: properties of point sources of the fault plane
    :param solution: dictionary with output kinematic model properties
    :type segments: list
    :type point_sources: array
    :type solution: dict
    """
    rupt_time = solution['rupture_time']
    max_rupt_time = [np.max(rupt_time_seg.flatten()) for rupt_time_seg in rupt_time]
    max_rupt_time = np.max(max_rupt_time)
    slip = solution['slip']
    max_slip = [np.max(slip_seg.flatten()) for slip_seg in slip]
    max_slip = np.max(max_slip)
    x_label = 'Distance along strike $(km)$'
    y_label = 'Distance along dip $(km)$'
    for i_segment, (segment, rupt_time_seg, slip_seg, ps_seg)\
    in enumerate(zip(segments, rupt_time, slip, point_sources)):
#
# Plot the slip distribution
#
        indexes = np.where(slip_seg < 0.1 * max_slip)
        rupt_time_seg[indexes] = 0
        fig = plt.figure(figsize=(15, 8))
        ax = fig.add_subplot(111)
        ax.set_ylabel(y_label)
        ax.set_xlabel(x_label)
        fig.subplots_adjust(right=0.75)
        if i_segment == 0:
            ax.plot(0, 0, 'w*', ms=20)
        ax, im = __several_axes(
            rupt_time_seg, segment, ps_seg, ax, max_val=max_rupt_time)
        cbar_ax = fig.add_axes([0.85, 0.15, 0.05, 0.7])
        cb = fig.colorbar(im, cax=cbar_ax)
        cb.set_label('Rupt_time (s)')
        plt.savefig(
            'RuptTime_plane{}.png'.format(i_segment), bbox_inches='tight')
        plt.close()
    return


def _PlotRiseTime(segments, point_sources, solution):
    """We plot rise time distribution based on the FFM solution model

    :param segments: list of dictionaries with properties of fault segments
    :param point_sources: properties of point sources of the fault plane
    :param solution: dictionary with output kinematic model properties
    :type segments: list
    :type point_sources: array
    :type solution: dict
    """
    rise_time = solution['trise']
    max_trise = [np.max(trise_seg.flatten()) for trise_seg in rise_time]
    max_trise = np.max(max_trise)
    fall_time = solution['tfall']
    max_tfall = [np.max(tfall_seg.flatten()) for tfall_seg in fall_time]
    max_tfall = np.max(max_tfall)
    slip = solution['slip']
    max_slip = [np.max(slip_seg.flatten()) for slip_seg in slip]
    max_slip = np.max(max_slip)
    x_label = 'Distance along strike $(km)$'
    y_label = 'Distance along dip $(km)$'
    for i_segment, (segment, trise_seg, tfall_seg, slip_seg, ps_seg)\
    in enumerate(zip(segments, rise_time, fall_time, slip, point_sources)):
#
# Plot the slip distribution
#
        indexes = np.where(slip_seg < 0.1 * max_slip)
        trise_seg[indexes] = 0
        tfall_seg[indexes] = 0
        fig, axes = plt.subplots(1, 2, figsize=(20, 10), sharex=True, sharey=True)
        fig.subplots_adjust(bottom=0.15)
        axes[0].set_ylabel(y_label)
        axes[0].set_xlabel(x_label)
        if i_segment == 0:
            axes[0].plot(0, 0, 'w*', ms=20)
        axes[0], im = __several_axes(
            trise_seg, segment, ps_seg, axes[0],
            max_val=max_trise, autosize=False)
        axes[1].set_ylabel(y_label)
        axes[1].set_xlabel(x_label)
        if i_segment == 0:
            axes[1].plot(0, 0, 'w*', ms=20)
        axes[1], im = __several_axes(
            tfall_seg, segment, ps_seg, axes[1],
            max_val=max_tfall, autosize=False)
        cbar_ax = fig.add_axes([0.1, 0.05, 0.8, 0.05])
        cb = fig.colorbar(im, cax=cbar_ax, orientation='horizontal')
        cb.set_label('Rise_time (s)')
        plt.savefig(
            'RiseTime_plane{}.png'.format(i_segment), bbox_inches='tight')
        plt.close()
    return


def _PlotSlipDistribution(segments, point_sources, solution):
    """We plot slip distribution based on the FFM solution model

    :param segments: list of dictionaries with properties of fault segments
    :param point_sources: properties of point sources of the fault plane
    :param solution: dictionary with output kinematic model properties
    :type segments: list
    :type point_sources: array
    :type solution: dict
    """
    slip = solution['slip']
    rake = solution['rake']
    max_slip = [np.max(slip_seg.flatten()) for slip_seg in slip]
    max_slip = np.max(max_slip)
    x_label = 'Distance along strike $(km)$'
    y_label = 'Distance along dip $(km)$'
    for i_segment, (segment, slip_seg, rake_seg, ps_seg)\
    in enumerate(zip(segments, slip, rake, point_sources)):
        max_slip_seg = np.max(slip_seg.flatten())
        u = slip_seg * np.cos(rake_seg * np.pi / 180.0) / max_slip_seg
        v = slip_seg * np.sin(rake_seg * np.pi / 180.0) / max_slip_seg
#
# Plot the slip distribution
#
        fig = plt.figure(figsize=(15, 8))
        ax = fig.add_subplot(111)
        ax.set_ylabel(y_label)
        ax.set_xlabel(x_label)
        fig.subplots_adjust(right=0.75)
        stk_subfaults, dip_subfaults, delta_strike, delta_dip, hyp_stk, hyp_dip\
            = pl_mng.__unpack_plane_data(segment)
        x = np.arange(stk_subfaults) * delta_strike - hyp_stk * delta_strike
        y = np.arange(dip_subfaults) * delta_dip - hyp_dip * delta_dip
        ax.quiver(x, y, u, v, scale=15.0, width=0.003)
        if i_segment == 0:
            ax.plot(0, 0, 'w*', ms=20)
        ax, im = __several_axes(
            slip_seg, segment, ps_seg, ax, max_val=max_slip)
        cbar_ax = fig.add_axes([0.85, 0.15, 0.05, 0.7])
        cb = fig.colorbar(im, cax=cbar_ax)
        cb.set_label('Slip (cm)')
        plt.savefig(
            'SlipDist_plane{}.png'.format(i_segment), bbox_inches='tight')
        plt.close()
    return


def _PlotSlipDist_Compare(segments, point_sources, input_model,
                          solution):
    """We plot slip distribution based on the FFM solution model and compare
    inverted model to an input model

    :param segments: list of dictionaries with properties of fault segments
    :param point_sources: properties of point sources of the fault plane
    :param solution: dictionary with output kinematic model properties
    :param input_model: input kinematic model
    :type segments: list
    :type point_sources: array
    :type solution: dict
    :type input_model: dict
    """
    slip = solution['slip']
    rake = solution['rake']
    slip2 = input_model['slip']
    rake2 = input_model['rake']
    max_slip = [np.max(slip_seg.flatten()) for slip_seg in slip]
    max_slip = np.max(max_slip)
    max_slip2 = [np.max(slip_seg2.flatten()) for slip_seg2 in slip2]
    max_slip2 = np.max(max_slip2)
    max_slip = np.maximum(max_slip, max_slip2)
    x_label = 'Distance along strike $(km)$'
    y_label = 'Distance along dip $(km)$'
    zipped = zip(segments, slip, rake, slip2, rake2, point_sources)
    for i_segment, (segment, slip_seg, rake_seg, slip_seg2, rake_seg2, ps_seg)\
    in enumerate(zipped):
        max_slip_seg = np.max(slip_seg.flatten())
        max_slip_seg2 = np.max(slip_seg2.flatten())
        max_slip_seg = np.maximum(max_slip_seg, max_slip_seg2)
        u = slip_seg * np.cos(rake_seg * np.pi / 180.0) / max_slip_seg
        v = slip_seg * np.sin(rake_seg * np.pi / 180.0) / max_slip_seg
        u2 = slip_seg2 * np.cos(rake_seg2 * np.pi / 180.0) / max_slip_seg
        v2 = slip_seg2 * np.sin(rake_seg2 * np.pi / 180.0) / max_slip_seg
#
# Plot the slip distribution
#
        fig, (ax0, ax1) = plt.subplots(1, 2, figsize=(30, 8))
        ax0.set_ylabel(y_label)
        ax0.set_xlabel(x_label)
        ax1.set_ylabel(y_label)
        ax1.set_xlabel(x_label)
        fig.subplots_adjust(right=0.75)
        stk_subfaults, dip_subfaults, delta_strike, delta_dip, hyp_stk, hyp_dip\
            = pl_mng.__unpack_plane_data(segment)
        x = np.arange(stk_subfaults) * delta_strike - hyp_stk * delta_strike
        y = np.arange(dip_subfaults) * delta_dip - hyp_dip * delta_dip
        ax0.quiver(x, y, u, v, scale=15.0, width=0.003)
        ax1.quiver(x, y, u2, v2, scale=15.0, width=0.003)
        if i_segment == 0:
            ax0.plot(0, 0, 'w*', ms=20)
            ax1.plot(0, 0, 'w*', ms=20)
        ax0, im = __several_axes(
            slip_seg, segment, ps_seg, ax0,
            max_val=max_slip, autosize=False)
        ax1, im = __several_axes(
            slip_seg2, segment, ps_seg, ax1,
            max_val=max_slip, autosize=False)
        ax0.set_title('Inverted model', fontsize=20)
        ax1.set_title('Original model', fontsize=20)
        cbar_ax = fig.add_axes([0.85, 0.15, 0.05, 0.7])
        cb = fig.colorbar(im, cax=cbar_ax)
        cb.set_label('Slip (cm)')
        plt.savefig('SlipDist_Compare_plane{}.png'.format(i_segment),
                    bbox_inches='tight')
        plt.close()
    return


def _PlotMap(tensor_info, segments, point_sources, solution, default_dirs,
             files_str=None, stations_gps=None, max_slip=None, event=None):
    """We plot slip map.

    :param tensor_info: dictionary with moment tensor information
    :param segments: list of dictionaries with properties of fault segments
    :param point_sources: properties of point sources of the fault plane
    :param files_str: location of regional data to add
    :param station_gps: static data to add
    :param solution: dictionary with output kinematic model properties
    :param max_slip: maximum slip to appear in colorbar
    :param default_dirs: dictionary with default directories to be used
    :type default_dirs: dict
    :type tensor_info: dict
    :type segments: list
    :type point_sources: array
    :type files_str: list
    :type station_gps: list
    :type solution: dict
    :type max_slip: float, optional
    """
    plane_info = segments[0]
    stk_subfaults, dip_subfaults, delta_strike, delta_dip, hyp_stk, hyp_dip\
        = pl_mng.__unpack_plane_data(plane_info)
    slip = solution['slip']
    
    point_sources2 = point_sources.copy()
    segments2 = segments.copy()
    if event is not None:
        zipped = zip(segments, slip)
        slip = [slip_seg for segment, slip_seg in zipped if segment['event']==event]
        zipped = zip(segments, point_sources)
        point_sources2 = [ps_seg for segment, ps_seg in zipped if segment['event']==event]
        segments2 = [segment for segment in segments if segment['event']==event]
#
# accurate plot coordinates
#
    segments_lats, segments_lons = __redefine_lat_lon(segments2, point_sources2)
    min_lats = [np.min(segment_lat.flatten()) for segment_lat in segments_lats]
    max_lats = [np.max(segment_lat.flatten()) for segment_lat in segments_lats]
    min_lons = [np.min(segment_lon.flatten()) for segment_lon in segments_lons]
    max_lons = [np.max(segment_lon.flatten()) for segment_lon in segments_lons]
    min_lat = np.min(min_lats)# - 0.5
    max_lat = np.max(max_lats)# + 0.5
    min_lon = np.min(min_lons)# - 0.5
    max_lon = np.max(max_lons)# + 0.5

    margin = 1.3 * (stk_subfaults * delta_strike) / 111.19
    lat0 = tensor_info['lat']
    lon0 = tensor_info['lon']
    tectonic = '{}.shp'.format(default_dirs['trench_graphics'])
    dictn = {
        'projection': ccrs.PlateCarree(),
        'facecolor': '#eafff5'
    }

    fig, ax = plt.subplots(1, 1, figsize=(15, 15), subplot_kw=dictn)
    fig.subplots_adjust(hspace=0, wspace=0, top=0.9, bottom=0.1, right=0.8)
    tectonic = cf.ShapelyFeature(
        shpreader.Reader(tectonic).geometries(), ccrs.PlateCarree(),
        edgecolor='red', facecolor=(198/255., 236/255., 253/255.))
    shpfilename = shpreader.natural_earth(
        resolution='10m', category='cultural', name='admin_0_countries')
    countries = cf.ShapelyFeature(
        shpreader.Reader(shpfilename).geometries(), ccrs.PlateCarree(),
        edgecolor='black', facecolor='lightgray')

    if files_str is not None:
        for file in files_str:
            name = file['name']
            latp, lonp = file['location']
            min_lat = min(min_lat, latp)
            max_lat = max(max_lat, latp)
            min_lon = min(min_lon, lonp)
            max_lon = max(max_lon, lonp)
            distance = max(np.abs(latp - lat0), np.abs(lonp - lon0))
            margin = max(margin, 1.2 * distance)
            ax.plot(
                lonp, latp, 'wo', markersize=10,
                transform=ccrs.PlateCarree(), zorder=4)
            ax.text(
                lonp + 0.1, latp + 0.1, '{}'.format(name),
                transform=ccrs.PlateCarree(), zorder=4)
    if stations_gps is not None:
        max_obs = np.zeros(3)
        stations_gps2 = []
        for name, sta_lat, sta_lon, obs, syn, error in stations_gps:
            min_lat = min(min_lat, sta_lat)
            max_lat = max(max_lat, sta_lat)
            min_lon = min(min_lon, sta_lon)
            max_lon = max(max_lon, sta_lon)
            stations_gps2 = stations_gps2\
                + [[name, sta_lat, sta_lon, obs, syn, error]]
            max_obs = np.maximum([abs(float(v)) for v in obs], max_obs)
            distance = max(np.abs(sta_lat - lat0), np.abs(sta_lon - lon0))
            margin = max(margin, 1.2 * distance)
        max_obs = np.max(max_obs)
        plt.text(
            lon0 + margin - 2, lat0 + margin - 0.25,
            '{:.2f} cm'.format(max_obs), transform=ccrs.PlateCarree())
        plt.text(
            lon0 + margin - 2, lat0 + margin - 0.45,
            '{:.2f} cm'.format(max_obs), transform=ccrs.PlateCarree())
        for name, sta_lat, sta_lon, obs, syn, error in stations_gps2:
            plt.plot(
                sta_lon, sta_lat, 'ks', transform=ccrs.PlateCarree(),
                markersize=14)
            gps_z, gps_n, gps_e = syn
            east_west = float(gps_e) / max_obs
            north_south = float(gps_n) / max_obs
            plt.arrow(
                sta_lon, sta_lat, east_west, north_south, color='r',
                zorder=3, linewidth=2, head_width=0.05, head_length=0.05,
                transform=ccrs.PlateCarree())
            up_down = float(gps_z) / max_obs
            plt.arrow(
                sta_lon, sta_lat, 0.0, up_down, color='r', zorder=3,
                linewidth=2, head_width=0.05, head_length=0.05,
                transform=ccrs.PlateCarree())
            gps_z, gps_n, gps_e = obs
            east_west = float(gps_e) / max_obs
            north_south = float(gps_n) / max_obs
            plt.arrow(
                sta_lon, sta_lat, east_west, north_south, zorder=3,
                linewidth=2, head_width=0.05, head_length=0.05,
                transform=ccrs.PlateCarree())
            up_down = float(gps_z) / max_obs
            plt.arrow(
                sta_lon, sta_lat, 0.0, up_down, zorder=3,
                linewidth=2, head_width=0.05, head_length=0.05,
                transform=ccrs.PlateCarree())
            plt.text(
                sta_lon + 0.1, sta_lat + 0.1, '{}'.format(name),
                transform=ccrs.PlateCarree())
            err_z, err_n, err_e = error
            width = float(err_e) / max_obs#/ 100
            height = float(err_n) / max_obs#/ 100
            ellipse = patches.Ellipse(
                (sta_lon + east_west, sta_lat + north_south), width,
                height, zorder=4, color='k', linewidth=10,
                transform=ccrs.PlateCarree())
            plt.gca().add_patch(ellipse)
        plt.arrow(
            lon0 + margin - 0.2, lat0 + margin - 0.2, -1, 0, color='r',
            zorder=3, linewidth=2, head_width=0.05, head_length=0.05,
            transform=ccrs.PlateCarree())
        plt.arrow(
            lon0 + margin - 0.2, lat0 + margin - 0.4, -1, 0, color='k',
            zorder=3, linewidth=2, head_width=0.05, head_length=0.05,
            transform=ccrs.PlateCarree())
    max_slip = max([np.amax(slip_fault) for slip_fault in slip])\
        if not max_slip else max_slip
    margins = [min_lon - 0.5, max_lon + 0.5, min_lat - 0.5, max_lat + 0.5]
    ax = set_map_cartopy(ax, margins, tectonic=tectonic, countries=countries)
    ax.plot(
        lon0, lat0, 'w*', markersize=15, transform=ccrs.PlateCarree(), zorder=4)
#
# plot slip map
#
    ax, cs = plot_map(
        ax, segments_lats, segments_lons, slip,
        max_val=max_slip, transform=dictn['projection'])
    cbar_ax = fig.add_axes([0.85, 0.1, 0.05, 0.8])
    cbar = plt.colorbar(cs, cax=cbar_ax)
    cbar.set_label('Slip (cm)', size=15)
    cbar.ax.yaxis.set_ticks_position('left')
    ax.set_aspect('auto', adjustable=None)
    plot_name = 'Map'
    if event is not None:
        plot_name = '{}_event{}'.format(plot_name, event)
    plt.savefig(plot_name, bbox_inches='tight')
    plt.close()
    return


def _PlotInsar(tensor_info, segments, point_sources, default_dirs,
               insar_points, los='ascending'):
    """We compare insar data with insar produced by the inverted earthquake
    model.

    :param tensor_info: dictionary with moment tensor information
    :param segments: list of dictionaries with properties of fault segments
    :param point_sources: properties of point sources of the fault plane
    :param insar_points: location, observed and inverted data for insar
    :param solution: dictionary with output kinematic model properties
    :param default_dirs: dictionary with default directories to be used
    :param los: whether insar track is ascending or descending
    :type default_dirs: dict
    :type tensor_info: dict
    :type segments: list
    :type point_sources: array
    :type insar_points: list
    :type solution: dict
    :type los: string, optional
    """
    plane_info = segments[0]
    stk_subfaults, dip_subfaults, delta_strike, delta_dip, hyp_stk, hyp_dip\
        = pl_mng.__unpack_plane_data(plane_info)
    # slip = solution['slip']
#
# accurate plot coordinates
#
    segments_lats, segments_lons = __redefine_lat_lon(segments, point_sources)
    min_lats = [np.min(segment_lat.flatten()) for segment_lat in segments_lats]
    max_lats = [np.max(segment_lat.flatten()) for segment_lat in segments_lats]
    min_lons = [np.min(segment_lon.flatten()) for segment_lon in segments_lons]
    max_lons = [np.max(segment_lon.flatten()) for segment_lon in segments_lons]
    min_lat = np.min(min_lats)# - 0.5
    max_lat = np.max(max_lats)# + 0.5
    min_lon = np.min(min_lons)# - 0.5
    max_lon = np.max(max_lons)# + 0.5

    margin = 1.3 * (stk_subfaults * delta_strike) / 111.19
    lat0 = tensor_info['lat']
    lon0 = tensor_info['lon']
    tectonic = '{}.shp'.format(default_dirs['trench_graphics'])
    dictn = {
        'projection': ccrs.PlateCarree(),
        'facecolor': '#eafff5'
    }

    fig, axes = plt.subplots(2, 3, figsize=(30, 15), subplot_kw=dictn)
    tectonic = cf.ShapelyFeature(
        shpreader.Reader(tectonic).geometries(), ccrs.PlateCarree(),
        edgecolor='red', facecolor=(198/255., 236/255., 253/255.))
    shpfilename = shpreader.natural_earth(
        resolution='10m', category='cultural', name='admin_0_countries')
    countries = cf.ShapelyFeature(
        shpreader.Reader(shpfilename).geometries(), ccrs.PlateCarree(),
        edgecolor='black', facecolor='lightgray')

    max_diff = -1
    min_diff = 1
    lats = [point['lat'] for point in insar_points]
    lons = [point['lon'] for point in insar_points]
    min_lat = min(min_lat, np.min(lats))
    max_lat = max(max_lat, np.max(lats))
    min_lon = min(min_lon, np.min(lons))
    max_lon = max(max_lon, np.max(lons))
    observed = [point['observed'] for point in insar_points]
    synthetic = [point['synthetic'] for point in insar_points]
    ramp = [point['ramp'] for point in insar_points]
    diffs = [obs - syn for obs, syn in zip(observed, synthetic)]
    obs_no_ramp = [obs - ramp for obs, ramp in zip(observed, ramp)]
    syn_no_ramp = [syn - ramp for syn, ramp in zip(synthetic, ramp)]

    margins = [min_lon - 0.5, max_lon + 0.5, min_lat - 0.5, max_lat + 0.5]

    values = [observed, synthetic, diffs, obs_no_ramp, syn_no_ramp, ramp]
    titles = ['Observed', 'Synthetic', 'Misfit','Observed - Ramp', 'Synthetic - Ramp', 'Ramp']
    labels = ['Observed LOS (m)', 'Modeled LOS (m)', 'Residual (m)', 'Observed LOS (m)', 'Modeled LOS (m)', 'Modeled Ramp (m)']
    rows = [0, 0, 0, 1, 1, 1]
    cols = [0, 1, 2, 0, 1, 2]
    zipped = zip(values, titles, labels, rows, cols)
    for value, title, label, row, col in zipped:
        axes[row][col].set_title(title, fontdict={'fontsize': 20})
        max_abs = np.max(np.abs(value))
        vmin = -max_abs - 0.001# if not title == 'Misfit' else -40
        vmax = max_abs + 0.001# if not title == 'Misfit' else 40
        norm = colors.TwoSlopeNorm(vmin=vmin, vcenter=0, vmax=vmax)
        cs = axes[row][col].scatter(
            lons, lats, zorder=4, c=value, cmap='bwr',
            norm=norm, transform=ccrs.PlateCarree())
        ax = set_map_cartopy(
            axes[row][col], margins, tectonic=tectonic, countries=countries)
        ax.plot(
            lon0, lat0, 'y*', markersize=15,
            transform=ccrs.PlateCarree(), zorder=4)
        ax = plot_borders(
            ax, segments_lats, segments_lons, transform=dictn['projection'])
        fig.colorbar(cs, ax=ax, orientation="horizontal")
        ax.set_aspect('auto', adjustable=None)

    fig.tight_layout()
    plt.savefig('Insar_{}_fit.png'.format(los), bbox_inches='tight')
    plt.close()
    return


def _PlotComparisonMap(tensor_info, segments, point_sources, input_model,
                       solution):
    """We plot slip map and compare with map from input model.

    :param tensor_info: dictionary with moment tensor information
    :param segments: list of dictionaries with properties of fault segments
    :param point_sources: properties of point sources of the fault plane
    :param solution: dictionary with output kinematic model properties
    :param input_model: input kinematic model
    :type tensor_info: dict
    :type segments: list
    :type point_sources: array
    :type solution: dict
    :type input_model: dict
    """
    input_slip = input_model['slip']
    plane_info = segments[0]
    stk_subfaults, dip_subfaults, delta_strike, delta_dip, hyp_stk, hyp_dip\
        = pl_mng.__unpack_plane_data(plane_info)
    if stk_subfaults * dip_subfaults == 1:
        return
    slip = solution['slip']
#
# accurate plot coordinates
#
    segments_lats, segments_lons = __redefine_lat_lon(segments, point_sources)

    margin = min(1.5 * (stk_subfaults * delta_strike) / 111.19, 10)#3
    lat0 = tensor_info['lat']
    lon0 = tensor_info['lon']
    margins = [lon0 - margin, lon0 + margin, lat0 - margin, lat0 + margin]
    dictn = {'projection': ccrs.PlateCarree(), 'facecolor': '#eafff5'}
    shpfilename = shpreader.natural_earth(
        resolution='10m', category='cultural', name='admin_0_countries')
    countries = cf.ShapelyFeature(
        shpreader.Reader(shpfilename).geometries(), ccrs.PlateCarree(),
        edgecolor='black', facecolor='lightgray')

    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(30, 15), subplot_kw=dictn)
    ax1.set_title('Inverted model', fontsize=22)
    ax2.set_title('Original model', fontsize=22)
    fig.subplots_adjust(hspace=0, wspace=0.1, top=0.9, bottom=0.3)
    for ax in [ax1, ax2]:
        ax = set_map_cartopy(ax, margins, countries=countries)
        ax.plot(
            lon0, lat0, 'w*', markersize=15,
            transform=ccrs.PlateCarree(), zorder=4)
    max_slip = max([np.amax(slip_fault) for slip_fault in slip])
    max_slip2 = max([np.amax(input_slip2) for input_slip2 in input_slip])
    max_slip = max(max_slip, max_slip2)
#
# plot slip map
#
    ax1, cs1 = plot_map(
        ax1, segments_lats, segments_lons, slip,
        max_val=max_slip, transform=dictn['projection'])
    ax2, cs2 = plot_map(
        ax2, segments_lats, segments_lons, input_slip,
        max_val=max_slip, transform=dictn['projection'])
    fig.subplots_adjust(bottom=0.15)
    cbar_ax = fig.add_axes([0.1, 0.05, 0.8, 0.05])
    cb = fig.colorbar(cs2, cax=cbar_ax, orientation='horizontal')
    cb.set_label('Slip (cm)')
    plt.savefig('Comparison.png', bbox_inches='tight')
    plt.close()
    return


def __redefine_lat_lon(segments, point_sources):
    """
    """
    segments_lats = [[]] * len(segments)
    segments_lons = [[]] * len(segments)
    for i, point_sources_seg in enumerate(point_sources):
        lat = point_sources_seg[:, :, :, :, 0]
        lon = point_sources_seg[:, :, :, :, 1]
        ny, nx, a, b = lat.shape
        new_lat = np.zeros((ny + 1, nx + 1))
        new_lon = np.zeros((ny + 1, nx + 1))
        for j in range(ny):
            for k in range(nx):
                new_lat[j, k] = lat[j, k, 0, 0]
                new_lon[j, k] = lon[j, k, 0, 0]
        for k in range(nx):
            new_lat[-1, k] = lat[-1, k, -1, 0]
            new_lon[-1, k] = lon[-1, k, -1, 0]
        for j in range(ny):
            new_lat[j, -1] = lat[j, -1, 0, -1]
            new_lon[j, -1] = lon[j, -1, 0, -1]
        new_lat[-1, -1] = lat[-1, -1, -1, -1]
        new_lon[-1, -1] = lon[-1, -1, -1, -1]
        segments_lats[i] = new_lat
        segments_lons[i] = new_lon
    return segments_lats, segments_lons


def _plot_moment_rate_function(segments_data, shear, point_sources, event=None):
    r"""We plot moment rate function
    """
    print('Creating Moment Rate Plot...')
    segments = segments_data['segments']
    plane_info = segments[0]
    dt = 0.01
    model = load_ffm_model.load_ffm_model(segments_data, point_sources)
    slip = model['slip']
    trup = model['trup']
    tl = model['trise']
    tr = model['tfall']
    properties = pl_mng.__unpack_plane_data(plane_info)
    delta_strike, delta_dip = [properties[2], properties[3]]
    stk_subfaults, dip_subfaults = [properties[0], properties[1]]
    tmax = 1.5 * np.max([np.amax(trup_seg + tl_seg + tr_seg)\
                         for trup_seg, tl_seg, tr_seg in zip(trup, tl, tr)])
    tmax = tmax if stk_subfaults * dip_subfaults > 1\
        else (tl[0][0, 0] + tr[0][0, 0]) * 8.0
    nmax = int((tmax/dt + 1))
    mr = np.zeros(nmax)
    seismic_moment = 0

    shear2 = shear.copy()
    point_sources2 = point_sources.copy()
    segments2 = segments.copy()
    if event is not None:
        zipped = zip(segments, slip)
        slip = [slip_seg for segment, slip_seg in zipped if segment['event']==event]
        zipped = zip(segments, trup)
        trup = [trup_seg for segment, trup_seg in zipped if segment['event']==event]
        zipped = zip(segments, tl)
        tl = [trise_seg for segment, trise_seg in zipped if segment['event']==event]
        zipped = zip(segments, tr)
        tr = [tfall_seg for segment, tfall_seg in zipped if segment['event']==event]
        zipped = zip(segments, shear)
        shear2 = [shear_seg for segment, shear_seg in zipped if segment['event']==event]
        zipped = zip(segments, point_sources)
        point_sources2 = [ps_seg for segment, ps_seg in zipped if segment['event']==event]
        segments2 = [segment for segment in segments if segment['event']==event]

    zipped = zip(segments2, slip, trup, tl, tr, shear2, point_sources2)

    for segment, slip_seg, trup_seg, trise_seg, tfall_seg, shear_seg,\
    point_sources_seg in zipped:
        dip_subfaults, stk_subfaults = np.shape(slip_seg)
        moment_rate = np.zeros(nmax)
        for iy in range(dip_subfaults):
            for ix in range(stk_subfaults):
                rupt_vel = segment['rupture_vel']
                rise_time = np.zeros(nmax)
                tfall = tfall_seg[iy, ix]
                trise = trise_seg[iy, ix]
                array1 = np.arange(0, trise, dt)
                tmid = len(array1)
                rise_time[:tmid]\
                    = (1 - np.cos(np.pi * array1 / trise)) / (trise + tfall)
                array2 = np.arange(0, tfall, dt)
                tend = tmid + len(array2)
                rise_time[tmid:tend]\
                    = (1 + np.cos(np.pi * array2 / tfall)) / (tfall + trise)
                duration = int(max(delta_strike, delta_dip) / dt / rupt_vel)
                source_dur = np.ones(duration) / duration
                start_index = max(0, int(trup_seg[iy, ix] / dt))
                product = slip_seg[iy, ix] * shear_seg[iy, ix] / 100 / 10
                sub_rise_time = rise_time[:tend]
                convolve = np.convolve(source_dur, sub_rise_time)
                moment_rate[start_index:start_index + len(convolve)]\
                = moment_rate[start_index:start_index + len(convolve)]\
                + convolve * product

        seismic_moment = seismic_moment\
            + np.sum((slip_seg / 100) * (shear_seg / 10)\
                     * (delta_strike * 1000) * (delta_dip * 1000))

#
# find moment rate function
#
        for i in range(nmax):
            time = i * dt
            mr[i] = mr[i]\
                + moment_rate[i] * (delta_strike * 1000) * (delta_dip * 1000)

    time = np.arange(nmax) * dt
    with open('STF.txt', 'w') as outf:
        outf.write('dt: {}\n'.format(dt))
        outf.write('Time[s]     Moment_Rate [Nm]\n')
        for t, val in zip(time, mr):
            outf.write('{:8.2f}:   {:8.4e}\n'.format(t, val))

    seismic_moment = np.trapz(mr, dx=0.01)
    magnitude = 2.0 * (np.log10(seismic_moment * 10 ** 7) - 16.1) / 3.0
    plt.xlabel('Time $(s)$')
    plt.ylabel('Moment rate $(Nm/s)$')
    plt.text(
        0.5 * max(time), 0.95 * max(mr),
        '$M_0$: {:.2E} $Nm$'.format(seismic_moment))
    plt.text(
        0.5 * max(time), 0.85 * max(mr),
        '$M_w$: {:.2f}'.format(magnitude))
    plt.grid('on')
    plt.fill_between(time, mr)
    plot_name = 'MomentRate'
    if event is not None:
        plot_name = '{}_event{}'.format(plot_name, event)
    plt.savefig(plot_name, bbox_inches='tight')
    plt.close()
    return


def _PlotSnapshotSlip(tensor_info, segments, point_sources, solution):
    """we plot snapshots of the rupture process.

    :param tensor_info: dictionary with moment tensor information
    :param segments: list of dictionaries with properties of fault segments
    :param point_sources: properties of point sources of the fault plane
    :param solution: dictionary with output kinematic model properties
    :type tensor_info: dict
    :type segments: list
    :type point_sources: array
    :type solution: dict
    """
    plane_info = segments[0]
    dt = 0.01
    slip = solution['slip']
    trup = solution['rupture_time']
    tl = solution['t_rise']
    tr = solution['t_fall']
    stk_subfaults, dip_subfaults, delta_strike, delta_dip, hyp_stk, hyp_dip\
        = pl_mng.__unpack_plane_data(plane_info)
    if [stk_subfaults, dip_subfaults] == [1, 1]:
        return
    tmid = trup + tl
    tstop = trup + tl + tr
    srmax = slip / (tr + tl) * 2.
#
# Define the vector field of the slip and plotting parameters.
#
    x = np.arange(stk_subfaults) * delta_strike - hyp_stk * delta_strike
    y = np.arange(dip_subfaults) * delta_dip - hyp_dip * delta_dip
    ratio = max(1,
        (9 / 16) * ((stk_subfaults * delta_strike) / (dip_subfaults * delta_dip)))
    vmax = np.amax(slip)
    tmax = np.amax(trup)
    step = int((tmax/dt + 1) / 9.)
#
# snapshots of rupture process
#
    fig, axes = plt.subplots(nrows=3, ncols=3, figsize=(9, 9))
    i = 1

    for ax in axes.flat:
        time = i * step * dt
        srate, cslip, broken = __rupture_process(
            time, slip, srmax, trup, tl, tr, tmid, tstop)
        ax.set_yticklabels([])
        ax.set_xticklabels([])
        if np.max(broken) > np.min(broken):
            ax.contour(
                -broken, 1, colors='k', extent=__extent_plot(plane_info))
        ax.contourf(
            cslip, cmap='jet', vmin=0, vmax=vmax,
            extent=__extent_plot(plane_info))
        ax.plot(0, 0, 'r*', ms=20)
        ax.invert_yaxis()
        ax.set_aspect(ratio)
        ax.set_title('Time: {0:.2f} s'.format(time))
        if i == 9: im = ax.contourf(x, y, cslip, cmap='jet')
        i = i + 1

    fig.subplots_adjust(right=0.8)
    cbar_ax = fig.add_axes([0.85, 0.15, 0.05, 0.7])
    cb = fig.colorbar(im, cax=cbar_ax)
    cb.set_label('Slip (cm)')
    fig.text(
        0.1, 0.1, 'CSN Automatic \nSolution', fontsize=50,
        color='gray', ha='left', va='bottom', alpha=0.5, wrap=True)
    plt.savefig('SlipSnapshot.png', bbox_inches='tight')
    plt.close()
    return


def plot_beachball(segments, files=None, phase=None):
    """Here we plot the beachball for the event. Optionally, we add the
    location of teleseismic data used in the FFM modelling.
    """
    segment = segments[0]
#
# Get the focal mechanism
#
    strike = segment['strike']
    dip = segment['dip']
    rake = segment['rake']
    fm = [strike, dip, rake]
#
# Plot the beach ball
#
    fig = plt.figure(figsize=(7, 7))
    bb = beach(fm, width=3.0, zorder=1)
    ax = plt.gca()
    ax.add_collection(bb)
    plt.plot(0, 0, 'b*', markersize=20)
    ax.set_aspect('equal')
#
# Load the station informations
#
    if files:
        plt.plot(0, 0, 'r*', markersize=20)
        for file in files:
            comp = file['component']
            if comp == 'BHZ': comp = 'P'
            dist = file['distance']
            az = file['azimuth']
            nam = file['name']
            if comp == phase:
                r = dist * (3.0 / 2.) / 90.
                x1 = np.sin(az * np.pi / 180.) * r
                y1 = np.cos(az * np.pi / 180.) * r
                plt.plot(x1, y1, 'ro')
                plt.text(x1 + 0.01, y1 + 0.01, nam)

    fig.patch.set_visible(False)
    plt.gca().axis('off')
    name_plot = '{}_azimuthcover.png'.format(phase) if phase else 'Tensor.png'
    plt.savefig(name_plot)
    plt.close()
    return


def __extent_plot(plane_info):
    """
    """
    stk_subfaults, dip_subfaults, delta_strike, delta_dip, hyp_stk, hyp_dip\
        = pl_mng.__unpack_plane_data(plane_info)
    return [-hyp_stk * delta_strike, (stk_subfaults - hyp_stk) * delta_strike,
            -hyp_dip * delta_dip, (dip_subfaults - hyp_dip) * delta_dip]


def __several_axes(data, segment, point_source_seg, ax, max_val=None,
                   autosize=True):
    """
    """
    stk_subfaults, dip_subfaults, delta_strike, delta_dip, hyp_stk, hyp_dip\
        = pl_mng.__unpack_plane_data(segment)
    min_dist = - (hyp_dip + 0.5) * delta_dip
    max_dist = (dip_subfaults - hyp_dip - 0.5) * delta_dip
    min_strike = - (hyp_stk + 0.5) * delta_strike
    max_strike = (stk_subfaults - hyp_stk - 0.5) * delta_strike
    dep = point_source_seg[:, :, :, :, 2]
    min_depth = dep[0, 0, 0, 0]
    max_depth = dep[-1, 0, -1, 0]
    if stk_subfaults * dip_subfaults == 1:
        dip = segment['dip']
        delta_z = delta_dip * np.sin(dip * np.pi / 180.0) / 2
        min_depth = min_depth - delta_z
        max_depth = max_depth + delta_z
    max_val = np.max(data.flatten()) if not max_val else max_val
    im = ax.imshow(
        data, cmap='jet', origin='lower', vmax=max_val, aspect='auto',
        extent=[min_strike, max_strike, min_dist, max_dist])
    ax2 = ax.twinx()
    ax2.set_xlim([min_strike, max_strike])
    ax2.set_ylim([min_depth, max_depth])
    ax.set(adjustable='datalim')
    ax2.set(adjustable='datalim')
    if autosize:
        ax.figure.set_size_inches(
            4 * stk_subfaults * delta_strike / dip_subfaults / delta_dip, 4)
        ax2.figure.set_size_inches(
            4 * stk_subfaults * delta_strike / dip_subfaults / delta_dip, 4)
    ax2.set_ylabel('Depth $(km)$')
    ax.invert_yaxis()
    ax2.invert_yaxis()
    return ax, im


def __rupture_process(time, slip, trup, tmid, tstop, rise_time,
                      point_sources, dt):
    """We give slip rate, rupture front, and accumulated slip at a certain
    time ``time``.
    """
    dip_subfaults, stk_subfaults = np.shape(slip)
    srate = np.zeros((dip_subfaults, stk_subfaults))
    cslip = np.zeros((dip_subfaults, stk_subfaults))
    broken = np.ones((dip_subfaults, stk_subfaults))

    for i in range(dip_subfaults):
        for j in range(stk_subfaults):
            convolve = rise_time[i, j, :]
            index = int(time / dt - trup[i, j] / dt)
            if (time < trup[i, j]):
                broken[i, j] = 0.
            elif index < len(convolve):
                srate[i, j] = convolve[index] * slip[i, j]#srmax[iy, ix]
            if trup[i, j] < time <= tmid[i, j]:
                cslip[i, j] = (time - trup[i, j]) * srate[i, j] / 2.
            if tmid[i, j] < time <= tstop[i, j]:
                cslip[i, j]\
                    = slip[i, j] - (tstop[i, j] - time) * srate[i, j] / 2.
            if (time > tstop[i, j]):
                cslip[i, j] = slip[i, j]
    return srate, cslip, broken


def __add_watermark(fig):
    """
    """
    fig.text(0.1, 0.1, 'CSN Automatic \nSolution', fontsize=50, color='gray',
             ha='left', va='bottom', alpha=0.5, wrap=True)
    return fig


if __name__ == '__main__':
    """
    """
    import management as mng
    import manage_parser as mp
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "-f", "--folder", default=os.getcwd(),
        help="folder where there are input files")
    parser = mp.parser_add_tensor(parser)
    parser = mp.parser_data_plot(parser)
    parser.add_argument(
        "-ffms", "--ffm_solution", action="store_true",
        help="plot FFM solution slip maps, rise time")
    parser.add_argument(
        "-bb", "--beachballs", action="store_true", help="plot beachballs")
    parser.add_argument(
        "-me", "--many_events", action="store_true", help="plots for many events")
    args = parser.parse_args()
    os.chdir(args.folder)
    if args.gcmt_tensor:
        args.gcmt_tensor = os.path.abspath(args.gcmt_tensor)
    used_data = mp.get_used_data(args)
    use_waveforms = False
    use_waveforms = use_waveforms if not 'strong_motion' in used_data else True
    use_waveforms = use_waveforms if not 'cgps' in used_data else True
    use_waveforms = use_waveforms if not 'tele_body' in used_data else True
    use_waveforms = use_waveforms if not 'surf_tele' in used_data else True
    default_dirs = mng.default_dirs()
    if args.gcmt_tensor:
        cmt_file = args.gcmt_tensor
        tensor_info = tensor.get_tensor(cmt_file=cmt_file)
    else:
        tensor_info = tensor.get_tensor()
    segments_data = json.load(open('segments_data.json'))
    segments = segments_data['segments']
    rise_time = segments_data['rise_time']
    connections = None
    if 'connections' in segments_data:
        connections = segments_data['connections']
    point_sources = pf.point_sources_param(
        segments, tensor_info, rise_time, connections=connections)
    if args.ffm_solution:
        solution = get_outputs.read_solution_static_format(segments)
        if not os.path.isfile('velmodel_data.json'):
            vel_model = mv.select_velmodel(tensor_info, default_dirs)
        else:
            vel_model = json.load(open('velmodel_data.json'))
        shear = pf.shear_modulous(point_sources, velmodel=vel_model)
        plot_ffm_sol(
            tensor_info, segments_data, point_sources, shear,
            solution, vel_model, default_dirs, use_waveforms=use_waveforms)
        if args.many_events:
            plot_ffm_sol(
                tensor_info, segments_data, point_sources, shear,
                solution, vel_model, default_dirs,
                event=1, use_waveforms=use_waveforms)
            plot_ffm_sol(
                tensor_info, segments_data, point_sources, shear,
                solution, vel_model, default_dirs,
                event=2, use_waveforms=use_waveforms)

    traces_info, stations_gps = [None, None]
    if args.gps:
        names, lats, lons, observed, synthetic, error\
            = get_outputs.retrieve_gps()
        stations_gps = zip(names, lats, lons, observed, synthetic, error)
    if args.strong:
        traces_info = json.load(open('strong_motion_waves.json'))
    if args.strong or args.gps:
        solution = get_outputs.read_solution_static_format(segments)
        _PlotMap(
            tensor_info, segments, point_sources, solution, default_dirs,
            files_str=traces_info, stations_gps=stations_gps)
        input_model = load_ffm_model.load_ffm_model(
            segments_data, point_sources, option='fault&rise_time.txt')
        _PlotSlipDist_Compare(segments, point_sources, input_model, solution)
        _PlotComparisonMap(
            tensor_info, segments, point_sources, input_model, solution)
    if args.insar:
        solution = get_outputs.read_solution_static_format(segments)
        _PlotMap(
            tensor_info, segments, point_sources, solution, default_dirs)
        insar_data = get_outputs.get_insar()
        if 'ascending' in insar_data:
            asc_properties = insar_data['ascending']
            for i, asc_property in enumerate(asc_properties):
                insar_points = asc_property['points']
                _PlotInsar(
                    tensor_info, segments, point_sources,
                    default_dirs, insar_points, los='ascending{}'.format(i))
        if 'descending' in insar_data:
            desc_properties = insar_data['descending']
            for i, desc_property in enumerate(desc_properties):
                insar_points = desc_property['points']
                _PlotInsar(
                    tensor_info, segments, point_sources,
                    default_dirs, insar_points, los='descending{}'.format(i))

    if args.beachballs:
        plot_beachballs(segments, used_data)
    plot_misfit(used_data)#, forward=True)
    if args.many_events:
        plot_misfit(used_data, event=1)
        plot_misfit(used_data, event=2)
    plot_files = glob.glob(os.path.join('plots', '*png'))
    for plot_file in plot_files:
        os.remove(plot_file)
    plot_files = glob.glob('*png')
    for plot_file in plot_files:
        move(plot_file, 'plots')
