#!/usr/bin/python
# -*- coding: utf-8 -*-
"""Module with routines to perform data processing (bandpass filtering,
downsampling, removal of instrumental response) to use such data in finite
fault modelling
"""

import subprocess
from shutil import copy2
from multiprocessing import Pool#, cpu_count
from obspy import read
from obspy.io.sac import SACTrace
from obspy.taup import TauPyModel
from obspy.core.utcdatetime import UTCDateTime
from obspy.geodetics import locations2degrees, kilometers2degrees, gps2dist_azimuth
from obspy.signal.invsim import simulate_seismometer
from obspy.signal import rotate
import numpy as np
import os
import glob
import logging
import management as mng
import modulo_logs as ml
import wang_baseline_removal_v1 as wang1
import time
import seismic_tensor as tensor
from scipy.integrate import cumtrapz
from scipy.signal import butter, filtfilt
import errno
import json
import remove_response as response


###################################
# teleseismic body waves
###################################


def select_process_tele_body(tele_files0, tensor_info, data_prop):
    """Module for automatic selection and processing of teleseismic body waves.

    :param tensor_info: dictionary with moment tensor information
    :param data_prop: dictionary with waveform properties
    :param tele_files0: files with body waves to be selected and processed
    :type tensor_info: dict
    :type data_prop: dict
    :type tele_files0: list

    .. note::

        Currently, this code allows only to process files in sac format, where the
        instrumental response is in a SACPZ file. Other data formats are not
        supported.
    """
    time1 = time.time()
    depth = tensor_info['depth']
    timedelta = tensor_info['timedelta']
    tele_files1 = __pre_select_tele(tele_files0, tensor_info, timedelta)
    if not tele_files1:
        return
    if not os.path.isdir('logs'):
        os.mkdir('logs')
    if not os.path.isdir('P'):
        os.mkdir('P')
    if not os.path.isdir('SH'):
        os.mkdir('SH')
    logger1 = ml.create_log('body_tele_processing',
                            os.path.join('logs', 'body_tele_processing.log'))
    logger1 = ml.add_console_handler(logger1)
    logger1.info('Start selection of teleseismic body data')

    data_folder = os.path.abspath(os.getcwd())

    logger1.info('Get theoretical arrivals of P and S waves')
    model = TauPyModel(model="ak135f_no_mud")
    __picker(tensor_info, tele_files1, depth, model)
    logger1.info('time spent getting theoretic arrivals: {}'.format(
        time.time() - time1))
    logger1.info('Remove response for selected teleseismic traces')
    response_files = glob.glob('SACPZ*') + glob.glob('SAC_PZs*')
    __remove_response_body(tele_files1, response_files, tensor_info, data_prop,
                           logger=logger1)
    logger1.info('time until remove response: {}'.format(time.time() - time1))
    logger1.info('Rotate horizontal components body waves')
#
    os.chdir(os.path.join(data_folder, 'SH'))
    horizontal_sacs = glob.glob('*sac')
    __rotation(tensor_info, horizontal_sacs, logger=logger1)

    logger1.info('Process body waves')
    os.chdir(data_folder)
    os.chdir(os.path.join(data_folder, 'P'))
    p_files = glob.glob('*_BHZ*sac')
    process_body_waves(p_files, 'P', data_prop, tensor_info)
    os.chdir(os.path.join(data_folder, 'SH'))
    s_files = glob.glob('*_SH*sac')
    process_body_waves(s_files, 'SH', data_prop, tensor_info)

    logger1.info('Body waves have been succesfully processed')
    os.chdir(data_folder)
    logger1.info('time until create files: {}'.format(time.time() - time1))
    logger1.info('total time spent: {}'.format(time.time() - time1))
    ml.close_log(logger1)
    return


def __pre_select_tele(tele_files0, tensor_info, timedelta):
    """Module which pre-selects teleseismic waveforms for processing
    """
    date_origin = tensor_info['date_origin']
    sacheaders = (SACTrace.read(sac) for sac in tele_files0)
    zipped = zip(tele_files0, sacheaders)
    fun1 = lambda header: header.stla and header.stlo
    zipped2 = ((sac, header) for sac, header in zipped if fun1(header))

    used_tele_sacs = []
    for sac, sacheader in zipped2:
        if not __select_distance(tensor_info, sacheader.stla, sacheader.stlo,
                                 max_distance=89, min_distance=31):
            continue
        if not __select_time(sac, tensor_info):
            continue
        if not sacheader.kcmpnm in sac:
            continue
        if not sacheader.khole in ['00', '', None]:
            continue
        used_tele_sacs = used_tele_sacs + [sac]
        sacheader.reftime = UTCDateTime(date_origin)
        sacheader.write(sac, byteorder='little')
    return used_tele_sacs


def __picker(tensor_info, used_tele_sacs, depth, model):
    """Obspy taupy picker.
    """
    event_lat = tensor_info['lat']
    event_lon = tensor_info['lon']
    for sac in used_tele_sacs:
        sacheader = SACTrace.read(sac)
        dist = locations2degrees(
                event_lat, event_lon, sacheader.stla, sacheader.stlo)
        arrivals = mng.theoretic_arrivals(model, dist, depth)
        sacheader.t1 = arrivals['p_arrival'][0].time
        sacheader.t2 = arrivals['pp_arrival'][0].time
        sacheader.t3 = arrivals['p_slowness']
        sacheader.t5 = arrivals['s_arrival'][0].time
        sacheader.t6 = arrivals['ss_arrival'][0].time
        sacheader.t7 = arrivals['s_slowness']
        sacheader.evla = event_lat
        sacheader.evlo = event_lon
        sacheader.write(sac, byteorder = 'little')
    return


def __remove_response_body(used_tele_sacs, response_files, tensor_info,
                           data_prop, logger=None):
    """We remove instrumental response for data, and replace it for a common
    instrumental response
    """
    string = [
        'ZEROS             2\n',
        '0.               0.\n',
        '0.               0.\n',
        'POLES             4\n',
        '-6.17E-03  6.17E-03\n',
        '-6.17E-03 -6.17E-03\n',
        '-39.18        49.12\n',
        '-39.18       -49.12\n',
        'CONSTANT         3948'
    ]
    string = ''.join(string)
    with open('response_file.txt', 'w') as resp_file:
        resp_file.write(string)
    paz_dict2, is_paz = __read_paz('response_file.txt')
    filtro = data_prop['tele_filter']
    depth = tensor_info['depth']
    time_shift = tensor_info['time_shift']
    t_min = 120.0
    t_max = 400.0
    if depth < 300:
        t_max = 240.0 if time_shift < 80 else 400.0
    freq0 = filtro['freq0']
    freq1 = filtro['low_freq']
    freq2 = filtro['high_freq']
    freq3 = filtro['freq3']
    freqs = [freq0, freq1, freq2, freq3]

    full_signal1 = lambda arrival, begin, end:\
        arrival - t_min >= begin and arrival + t_max <= end
    new_name = lambda network, station, channel:\
        os.path.join('{}_{}_{}.sac'.format(network, station, channel))

    for sac in used_tele_sacs:
        sacheader = SACTrace.read(sac)
        stream = read(sac)
        channel = stream[0].stats.channel
        phase = 'P' if channel == 'BHZ' else 'SH'
        arrival = sacheader.t1 if channel == 'BHZ' else sacheader.t5
        begin = sacheader.b
        end = sacheader.e
        full_signal = full_signal1(arrival, begin, end)
        name = stream[0].stats.station
        network = stream[0].stats.network
        str_name = os.path.join(phase, new_name(network, name, channel))
        loc_code = stream[0].stats.location
        pz_files0 = [resp for resp in response_files if name in resp]
        pz_files0 = [resp for resp in pz_files0 if channel in resp]
        if loc_code == '00':
            pz_files = [response for response in pz_files0 if '00' in response]
        if loc_code in [None, '']:
            pz_files = [resp for resp in pz_files0 if '--' in resp or '__' in resp]
        if not pz_files:
            continue
        pzfile_val = next(iter(pz_files))#
        if not os.path.isfile(pzfile_val):
            continue
        paz_dict, is_paz = __read_paz(pzfile_val)
        if not is_paz:
            if logger:
                logger.warning(
                    'PZ response unavailable at station {}, channel {}'.format(
                        name, channel))
            continue

        if full_signal:
            stream[0].detrend('linear')
            stream[0].taper(max_percentage=0.05)
            nyq = 0.5 * stream[0].stats.sampling_rate
            low_freq2 = freq1 / nyq
            b, a = butter(2, low_freq2, btype='highpass')
            filt_data = filtfilt(b, a, stream[0].data)
            stream[0].data = filt_data
            stream[0].data = simulate_seismometer(
                stream[0].data,
                stream[0].stats.sampling_rate,
                paz_remove=paz_dict,
                paz_simulate=paz_dict2,
                pre_filt=freqs,
                taper=False
            )
            stream.write(str_name, format='SAC', byteorder = 0)
    return


def __rotation(tensor_info, used_sacs, rotate_RT=True, logger=None):
    """We rotate horizontal body waves to match transverse and radial
    directions to the event.
    """
    event_lat = tensor_info['lat']
    event_lon = tensor_info['lon']
    string = lambda x, y, z: '{}_{}_{}.sac'.format(x, y, z)

    streams = [read(sac) for sac in used_sacs]
    streams = tuple(streams)
    stations = [st[0].stats.station for st in streams]
    stations = list(set(stations))

    for sta in stations:
        streams2 = [st for st in streams if st[0].stats.station == sta]

        streams_n1 = [st for st in streams2 if st[0].stats.channel[-1] in ['N', '1']]
        streams_e2 = [st for st in streams2 if st[0].stats.channel[-1] in ['E', '2']]
        if len(streams_n1) + len(streams_e2) < 2:
            continue
        stream = streams_n1[0] + streams_e2[0]
        starttime1 = stream[0].stats.starttime
        starttime2 = stream[1].stats.starttime
        endtime1 = stream[0].stats.endtime
        endtime2 = stream[1].stats.endtime
        diff1 = starttime1 - starttime2
        diff2 = endtime1 - endtime2
        starttime = starttime1 if diff1 >= 0 else starttime2
        endtime = endtime1 if diff2 <= 0 else endtime2
        stream.trim(starttime=starttime, endtime=endtime)
        station_lat = stream[0].stats.sac.stla
        station_lon = stream[0].stats.sac.stlo
#
# previous rotation for not NE components
#
        if stream[0].stats.channel[-1] == '1':
            cmpaz1 = stream[0].stats.sac.cmpaz
            cmpaz2 = stream[1].stats.sac.cmpaz
            dip1 = stream[0].stats.sac.cmpinc
            dip2 = stream[1].stats.sac.cmpinc
            data1 = stream[0].data
            data2 = stream[1].data
            data0 = np.ones(len(data1))
            data0, data1, data2 = rotate.rotate2zne(
                data0,
                0,
                -90,
                data1,
                cmpaz1,
                0,
                data2,
                cmpaz2,
                0
            )
            stream[0].data = data1
            stream[1].data = data2
            channel1 = stream[0].stats.channel
            channel2 = stream[1].stats.channel
            stream[0].stats.channel = channel1[:-1] + 'N'
            stream[1].stats.channel = channel2[:-1] + 'E'
        net = stream[0].stats.network
        if rotate_RT:
            transverse = string(net, sta, 'SH')
            radial = string(net, sta, 'RD')
            dist, az, baz = gps2dist_azimuth(
                event_lat, event_lon, station_lat, station_lon)
            stream.rotate(method='NE->RT', back_azimuth=baz)
            stream_radial = stream.select(channel="*R")
            stream_transverse = stream.select(channel="*T")
            stream_transverse.write(transverse, format='SAC', byteorder = 0)
            stream_radial.write(radial, format='SAC', byteorder = 0)
        else:
            channel1 = stream[0].stats.channel
            channel2 = stream[1].stats.channel
            north = string('acc', sta, channel1)
            east = string('acc', sta, channel2)
            stream_north = stream.select(channel="*N")
            stream_east = stream.select(channel="*E")
            stream_east.write(east, format='SAC', byteorder = 0)
            stream_north.write(north, format='SAC', byteorder = 0)
    return


def process_body_waves(tele_files, phase, data_prop, tensor_info):
    r"""We high pass filter and downsample teleseismic data.

    :param phase: string indicating whether phase is P or SH
    :param tensor_info: dictionary with moment tensor information
    :param data_prop: dictionary with waveform properties
    :param tele_files: body wave files to be resampled and filtered
    :type phase: string
    :type tensor_info: dict
    :type data_prop: dict
    :type tele_files: list
    """
    dt = data_prop['sampling']['dt_tele']
    date_origin = tensor_info['date_origin']
    filtro = data_prop['tele_filter']
    low_freq = filtro['low_freq']
    high_freq = filtro['high_freq']
    t_min = 120.0
    t_max = 360.0

    for sac in tele_files:
        sacheader = SACTrace.read(sac)
        sacheader.t8 = low_freq
        sacheader.t9 = high_freq
        if phase == 'SH':
            sacheader.kcmpnm = 'SH'
        sacheader.write(sac, byteorder='little')
    if phase == 'SH':
        for sac in tele_files:
            sacheader = SACTrace.read(sac)
            sacheader.kcmpnm = 'SH'
            sacheader.write(sac, byteorder='little')

    high_freq = 1 if phase == 'P' else 0.5
    high_freq = min(high_freq, 1 / 2.5 / dt)
    for sac in tele_files:
        sacheader = SACTrace.read(sac)
        stream = read(sac)
        arrival = sacheader.t1 if phase == 'P' else sacheader.t5
        start = date_origin + arrival - t_min
        end = date_origin + arrival + t_max
        tr = stream[0]
        tr.trim(starttime=start, endtime=end)
        nyq = 0.5 * tr.stats.sampling_rate
        high_freq2 = high_freq / nyq
        b, a = butter(2, high_freq2, btype='lowpass')
        filt_data = filtfilt(b, a, tr.data)
        tr.data = filt_data
        factor = dt * tr.stats.sampling_rate
        if abs(int(round(factor)) - factor) > 10 ** - 2: continue
        tr.decimate(int(round(factor)), no_filter=True)
        tr.data = tr.data - tr.data[0]
        if phase == 'SH': tr.detrend(type='linear')
        tr.data = 10 ** 6 * tr.data
        if np.max(np.abs(tr.data)) > 5 * 10 ** 5: continue
        stream[0] = tr
        stream.write('final_{}'.format(sac), format='SAC', byteorder = 0)
    return


##################################
# teleseismic surface waves
##################################


def select_process_surf_tele(tele_files0, tensor_info, data_prop):
    """Module for selecting and processing surface wave data.

    :param tensor_info: dictionary with moment tensor information
    :param tele_files0: files with surface wave to be selected and processed
    :type tensor_info: dict
    :type tele_files0: list

    .. rubric:: Example:

    >>> from obspy.core.utcdatetime import UTCDateTime
    >>> import glob
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
            'centroid_depth': 35,
            'timedelta': 30 * 60
        }
    >>> tele_files0 = glob.glob('*sac')
    >>> select_process_surf_tele(tele_files0, tensor_info)
    .. note::

        Currently, this code allows only to process files in sac format, where the
        instrumental response is in a SACPZ file. Other data formats are not
        supported.
    """
    timedelta = tensor_info['timedelta']
    depth = tensor_info['depth']
    tele_files1 = __pre_select_tele(tele_files0, tensor_info, timedelta)
    if not tele_files1:
        return
    if not os.path.isdir('logs'):
        os.mkdir('logs')
    if not os.path.isdir('LONG'):
        os.mkdir('LONG')
    logger1 = ml.create_log('surf_tele_processing',
                            os.path.join('logs', 'surf_tele_processing.log'))
    logger1 = ml.add_console_handler(logger1)
    logger1.info('Start processing of long period surface waves')

    data_folder = os.path.abspath(os.getcwd())

    os.chdir(data_folder)
    logger1.info('Get theoretical arrivals of P and S waves')
    model = TauPyModel(model="ak135f_no_mud")
    __picker(tensor_info, tele_files1, depth, model)
    logger1.info('Remove instrumental response for surface waves')
    response_files = glob.glob('SACPZ*') + glob.glob('SAC_PZs*')
    __remove_response_surf(tele_files1, response_files, data_prop, logger=logger1)
    logger1.info('Rotate horizontal components for surface waves')
    os.chdir(os.path.join(data_folder, 'LONG'))
    horizontal_sacs = glob.glob('*_BH[EN12]*sac')
    __rotation(tensor_info, horizontal_sacs)
    logger1.info('Get final surface waves')
    os.chdir(os.path.join(data_folder, 'LONG'))
    surf_files = glob.glob('*_BHZ*.sac') + glob.glob('*_SH*.sac')
    __create_long_period(surf_files, tensor_info)

    logger1.info('Long period surface waves have been succesfully processed')
    os.chdir(data_folder)
    ml.close_log(logger1)
    return


def __remove_response_surf(used_tele_sacs, response_files, data_prop, logger=None):
    """We remove instrumental response for surface wave data.
    """
    new_name = lambda network, station, channel:\
        os.path.join('{}_{}_{}.sac'.format(network, station, channel))

    filtro = data_prop['surf_filter']
    freq1 = filtro['freq1']
    freq2 = filtro['freq2']
    freq3 = filtro['freq3']
    freq4 = filtro['freq4']
    freqs = [freq1, freq2, freq3, freq4]

    for sac in used_tele_sacs:
        stream = read(sac)
        name = stream[0].stats.station
        channel = stream[0].stats.channel
        loc_code = stream[0].stats.location
        starttime1 = stream[0].stats.starttime
        endtime1 = stream[0].stats.endtime
        diff1 = endtime1 - starttime1
        if diff1 < 45 * 60:
            continue
        pz_files0 = [resp for resp in response_files if name in resp]
        pz_files0 = [resp for resp in pz_files0 if channel in resp]
        if loc_code == '00':
            pz_files = [response for response in pz_files0 if '00' in response]
        if loc_code in [None, '']:
            pz_files = [resp for resp in pz_files0 if '--' in resp or '__' in resp]
        if not pz_files:
            continue
        pzfile_val = next(iter(pz_files))#
        if not os.path.isfile(pzfile_val):
            continue
        paz_dict, is_paz = __read_paz(pzfile_val)
        if not is_paz:
            if logger:
                logger.warning(
                    'PZ response unavailable at station {}, channel {}'.format(
                        name, channel))
            continue
        netwk = stream[0].stats.network
        str_name = os.path.join('LONG', new_name(netwk, name, channel))

        stream[0].differentiate()
        stream[0].taper(max_percentage=0.05)
        stream[0].integrate()
        stream[0].data = simulate_seismometer(
            stream[0].data,
            stream[0].stats.sampling_rate,
            paz_remove=paz_dict,
            pre_filt=freqs,
            taper=False
        )
        stream.write(str_name, format='SAC', byteorder = 0)
    return


def __create_long_period(surf_files, tensor_info):
    """We downsample long period data.
    """
    for sac in surf_files:
        sacheader = SACTrace.read(sac)
        sacheader.kcmpnm = 'SH'\
        if not sacheader.kcmpnm == 'BHZ' else sacheader.kcmpnm
        sacheader.write(sac, byteorder = 'little')

    for sac in surf_files:
        st = read(sac)
        tr = st[0]
        factor = 4 * tr.stats.sampling_rate
        tr.decimate(int(round(factor)), no_filter=True)
        tr.data = 1000 * tr.data
        if np.max(np.abs(tr.data)) < 10 ** -3: continue
        st[0] = tr
        st.write('final_{}'.format(sac), format='SAC', byteorder = 0)
    return


###################################
# cgps
###################################


def select_process_cgps(cgps_files, tensor_info, data_prop):
    """Routine for selecting and processing cgps data

    :param tensor_info: dictionary with moment tensor information
    :param data_prop: dictionary with waveform properties
    :param cgps_files: files with cGPS data to be selected and processed
    :type tensor_info: dict
    :type data_prop: dict
    :type cgps_files: list

    .. note::

        Currently, this code allows only to process files in sac format.
        Other data formats are not supported.
    """
    if not os.path.isdir('logs'):
        os.mkdir('logs')
    logger1 = ml.create_log('cgps_processing',
                            os.path.join('logs', 'cgps_processing.log'))
    logger1 = ml.add_console_handler(logger1)
    logger1.info('Process cGPS data')
    data_folder = os.path.abspath(os.getcwd())
    separator = '\n*************************************************\n'
    logger1.info('Select cGPS traces')
    __select_cgps_files(cgps_files, tensor_info)
    logger1.info('Process selected cGPS traces')
    os.chdir(os.path.join(os.getcwd(), 'cGPS'))
    final_files = glob.glob('final*')
    for final_file in final_files:
        if os.path.isfile(final_file):
            os.remove(final_file)
    cgps_files1 = glob.glob('*L[HXY]*.sac') + glob.glob('*L[HXY]*.SAC')
    __change_start(cgps_files1, tensor_info, cgps=True)
    new_process_cgps(tensor_info, cgps_files1, data_prop, logger=logger1)
    os.chdir(data_folder)
    logger1.info('cGPS waves have been succesfully processed')
    strong_files3 = glob.glob(os.path.join('cGPS', 'final*'))
    select_strong_stations(tensor_info, strong_files3)
    ml.close_log(logger1)
    return


def __select_distance(tensor_info, station_lat, station_lon,
                      max_distance=180, min_distance=0, use_centroid=False):
    """
    """
    event_lat = tensor_info['lat']
    event_lon = tensor_info['lon']
    distance = locations2degrees(event_lat, event_lon, station_lat, station_lon)
    if use_centroid:
        centroid_lat = tensor_info['centroid_lat']
        centroid_lon = tensor_info['centroid_lon']
        distance2 = locations2degrees(
            centroid_lat, centroid_lon, station_lat, station_lon)
        distance = min(distance, distance2)
    return min_distance <= distance < max_distance


def __select_time(sac, tensor_info):
    """
    """
    datetime = tensor_info['datetime']
    origin_time = UTCDateTime(datetime)
    st = read(sac)
    return st[0].stats.starttime - origin_time < 20 * 60


def __select_cgps_files(cgps_files, tensor_info):
    """We select cgps data. We must make sure our cgps data is not too bad
    """
    event_lat = tensor_info['lat']
    event_lon = tensor_info['lon']
    depth = tensor_info['depth']
    time_shift = tensor_info['time_shift']
    centroid_lat = tensor_info['centroid_lat']
    event_lon = tensor_info['lon']
    date_origin = tensor_info['date_origin']
    depth = tensor_info['depth']
    distance = 2 if time_shift < 50 else 4

    for sac in cgps_files:
        stream = read(sac)
        station_lat = stream[0].stats.sac.stla
        station_lon = stream[0].stats.sac.stlo
        if station_lat == None or station_lon == None:
            continue
        syn_len = 2 * time_shift + 7 * depth / 50
        start = date_origin
        starttime = stream[0].stats.starttime
        delta = min(start - starttime, 20)
        stream.trim(starttime=start - delta)
        start, end = [stream[0].stats.starttime, stream[0].stats.endtime]
        if end - start < syn_len:
            continue
        if not __select_distance(tensor_info, station_lat, station_lon,
                                 max_distance=distance, use_centroid=True):
            continue
        indexes = np.isfinite(stream[0].data)
        if np.max(stream[0].data[indexes]) == np.min(stream[0].data[indexes]):
            continue

        stream[0].data = stream[0].data\
            - np.sum(stream[0].data[np.isfinite(stream[0].data)])\
            / np.sum(np.isfinite(stream[0].data))

        if np.max(np.abs(stream[0].data)) > 10 ** 4:
            continue

        if not np.sum(np.isnan(stream[0].data)) == 0:
            station_lat = stream[0].stats.sac.stla
            station_lon = stream[0].stats.sac.stlo
            delta = stream[0].stats.delta
            distance = locations2degrees(
                event_lat, event_lon, station_lat, station_lon)
            distance = 111.11 * distance
            distance2 = np.sqrt(distance ** 2 + depth ** 2)

            index0 = distance2 / 5 / delta
            indexes2 = [i for i, value in enumerate(stream[0].data)\
                if not np.isfinite(value)]
            indexes3 = [i for i in indexes2 if index0 <= i < index0 + 100]
            gaps = __get_gaps(stream[0].data)
            if not max(indexes2) < index0\
            and len(indexes3) >= 5 and max([len(gap) for gap in gaps]) >= 4:
                continue
            else:
                for gap in gaps:
                    begin = gap[0] - 1
                    end = gap[-1] + 1
                    end = min(end, len(stream[0].data) - 1)
                    stream[0].data = __linear_fill(stream[0].data, begin, end)

        station = stream[0].stats.station
        channel = stream[0].stats.channel
        name = os.path.join('cGPS', '{}_{}.sac'.format(station, channel))
        stream.write(name, format='SAC', byteorder=0)


def __change_start(stations_str, tensor_info, cgps=False):
    """Routine for modifying start time of cGPS data
    """
    event_lat = tensor_info['lat']
    event_lon = tensor_info['lon']
    depth = tensor_info['depth']
    for sac in stations_str:
        st = read(sac)
        sacheader = SACTrace.read(sac)
        station_lat = st[0].stats.sac.stla
        station_lon = st[0].stats.sac.stlo
        dist = locations2degrees(event_lat, event_lon, station_lat, station_lon)
        dist = dist * 111.12
        dist = np.sqrt(dist**2 + depth**2)
        est_arrival = dist / 5
        begin = sacheader.b
        reftime = sacheader.reftime
        st[0].stats.starttime = UTCDateTime(reftime) + begin
        diff = st[0].stats.starttime - tensor_info['date_origin']
        if diff > est_arrival:
            continue
        if cgps:
            st[0].trim(endtime=tensor_info['date_origin'] + 1100)  # - 60 * 10)
            st[0].data = st[0].data - np.mean(st[0].data[:10])
        st[0].trim(starttime=tensor_info['date_origin'] - 20, pad=True,
                   fill_value=0)
        st.write(sac, format='SAC', byteorder=0)


def new_process_cgps(tensor_info, stations_str, data_prop, logger=None):
    """Routine for processing (filtering and filling header) selected cGPS data.

    :param tensor_info: dictionary with moment tensor information
    :param stations_str: list with waveforms (in sac format) to use
    :param data_prop: dictionary with waveform properties
    :param logger: where to log results of processing data
    :type tensor_info: dict
    :type stations_str: list
    :type data_prop: dict
    :type logger: logging, optional
    """
    start = tensor_info['date_origin']
    filtro_cgps = data_prop['strong_filter']
    if 'cgps_filter' in data_prop:
        filtro_cgps = data_prop['cgps_filter']
    lat_ep = tensor_info['lat']
    lon_ep = tensor_info['lon']
    sampling = data_prop['sampling']
    dt_cgps = sampling['dt_cgps']

    high_freq = min(filtro_cgps['high_freq'], 0.5)

    for sac in stations_str:
        st = read(sac)
        if dt_cgps > st[0].stats.delta:
            ratio = int(dt_cgps / st[0].stats.delta)
            st[0].decimate(ratio)
        if st[0].stats.delta > dt_cgps:
            st[0].interpolate(1/dt_cgps)
            # st[0].resample(sampling_rate=1/dt_cgps)
        sacheader = SACTrace.from_obspy_trace(st[0])
        sacheader.t9 = high_freq
        sacheader.write(sac, byteorder='little')

    for sac in stations_str:
        sacheader = SACTrace.read(sac)
        begin = sacheader.b
        start_time = UTCDateTime(sacheader.reftime) + begin
        sacheader.o = start - start_time
        sacheader.write(sac, byteorder='little')

    _filter_decimate(stations_str, filtro_cgps, dt_cgps, corners=4, passes=2,
                     decimate=False, filter0='lowpass', logger=logger)
    return


def __linear_fill(data, begin, end):
    """
    """
    slope = (data[begin] - data[end]) / (begin - end) if begin >= 0 else 0
    intercept = data[begin] - slope * begin if begin >= 0 else data[end]
    for index in range(begin + 1, end):
        data[index] = slope * index + intercept
    return data


def __get_gaps(data):
    """
    """
    indexes = [i for i, value in enumerate(data) if not np.isfinite(value)]
    if not indexes:
        return None
    start_gaps = [index for index in indexes if index - 1 not in indexes]
    end_gaps = [index for index in indexes if index + 1 not in indexes]
    gaps = []
    for start, end in zip(start_gaps, end_gaps):
        gaps = gaps + [[index for index in range(start, end + 1)]]
    return gaps


###################################
# custom strong motion
###################################


def _filter_decimate(sac_files, filtro, dt, corners=4, passes=2,
                     decimate=True, logger=None, filter0='bandpass'):
    """
    """
    low_freq = filtro['low_freq']
    high_freq = filtro['high_freq']
    power = lambda m, n: max([d for d in range(10) if m % (n**d) == 0])
    for sac in sac_files:
        st = read(sac)
        if decimate:
            st[0].interpolate(1/dt)
        delta = st[0].stats.delta
        nyq = 0.5 / delta
        high_freq2 = high_freq / nyq
        low_freq2 = low_freq / nyq
        if filter0 == 'bandpass':
            b, a = butter(corners, [low_freq2, high_freq2], btype='bandpass')
        else:
            b, a = butter(corners, high_freq2, btype='lowpass')
        filt_data = filtfilt(b, a, st[0].data)
        st[0].data = 100*filt_data
        st.write('final_{}'.format(sac), format='SAC', byteorder=0)
    return


#####################
# automatic strong motion
#####################



def select_process_strong(strong_files0, tensor_info, data_prop,
                          remove_response=True):
    """Module for selecting and processing strong motion data

    :param tensor_info: dictionary with moment tensor information
    :param data_prop: dictionary with waveform properties
    :param strong_files0: files with strong motions to be selected and processed
    :param remove_response: whether to remove paz response in processing
    :type tensor_info: dict
    :type data_prop: dict
    :type strong_files0: list
    :type remove_repsonse: bool, optional

    .. note::
        Currently, this code allows only to process files in sac format, where the
        instrumental response is in a SACPZ file. Other data formats are not
        supported.
    """
    data_folder = os.path.abspath(os.getcwd())
    if not os.path.isdir('logs'):
        os.mkdir('logs')
    if not os.path.isdir('STR'):
        os.mkdir('STR')
    logger1 = ml.create_log('strong_motion_processing',
                            os.path.join('logs', 'strong_motion_processing.log'))
    logger1 = ml.add_console_handler(logger1)
    logger1.info('Process strong motion data')

    move_str_files(strong_files0)
    response_files = glob.glob('SACPZ*') + glob.glob('SAC_PZs*')
    response_files = [os.path.abspath(response) for response in response_files]
    os.chdir(os.path.join(data_folder, 'STR'))
    final_files = glob.glob('final*')
    for final_file in final_files:
        if os.path.isfile(final_file):
            os.remove(final_file)
    strong_files1 = glob.glob('acc*')
    if remove_response:
        logger1.info('Remove response for selected strong motion traces')
        __remove_response_str(strong_files1, response_files, logger=logger1)

    logger1.info('Select strong motion traces')
    os.chdir(os.path.join(data_folder, 'STR'))
    horizontal_sacs = glob.glob('*HN1*sac') + glob.glob('*HN2*sac')
    __rotation(tensor_info, horizontal_sacs, rotate_RT=False, logger=logger1)
    for file in horizontal_sacs:
        if os.path.isfile(file):
            os.remove(file)
    strong_files1 = glob.glob('acc*')
    strong_files2 = __select_str_files(strong_files1, tensor_info)
    logger1.info('Update duration of traces')
    __change_start(strong_files2, tensor_info)
    logger1.info('Process selected strong motion traces')
    process_strong_motion(strong_files2, tensor_info, data_prop, logger=logger1)
    os.chdir(data_folder)
    strong_files3 = glob.glob(os.path.join('STR', 'final*'))
    select_strong_stations(tensor_info, strong_files3)
    logger1.info('Strong motion waves have been succesfully processed')
    ml.close_log(logger1)
    return


def __convert_response_acc(resp_file):
    """Modify SACPZ response to ensure output is in units of acceleration after
    response removal.
    """
    with open(resp_file, 'r') as infile:
        lines = [line for line in infile]

    indexes0 = [i for i, line in enumerate(lines) if 'CONSTANT' in line.split()]
    index0 = 0
    temp_resp_file = 'temp.txt'
    lines3 = []
    for index in indexes0:
        lines1 = lines[index0:index + 1]
        index0 = index
        index2 = next(i for i, line in enumerate(lines1) if 'ZEROS' in line.split())
        index3 = next(i for i, line in enumerate(lines1) if 'POLES' in line.split())
        lines2 = lines1[:index2]
        with open(temp_resp_file, 'w') as outfile:
            for line in lines1:
                outfile.write(line)
        paz_dict, is_paz = __read_paz(temp_resp_file)
        zeros = paz_dict['zeros']
        zeros1 = [zero for zero in zeros if abs(zero - 0j) < 1e-10]
        input_unit = [line for line in lines1 if 'INPUT UNIT' in line]
        if len(input_unit) > 0:
            input_unit = input_unit[0]
            input_unit = input_unit.split()[-1]
            input_unit = input_unit.lower()
            if input_unit in ['m']:
                nzeros = 2
            elif input_unit in ['m/s']:
                nzeros = 1
            elif input_unit in ['m/s**2', 'm/s/s']:
                nzeros = 0
            else:
                nzeros = len(zeros1)
            zeros1 = zeros1[:nzeros]
        for zero in zeros1:
            zeros.remove(zero)
        paz_dict['zeros'] = zeros
        nzeros1 = len(zeros)
        lines2 = lines2 + ['ZEROS\t {}\n'.format(nzeros1)]
        for zero in zeros:
            real = np.real(zero)
            imag = np.imag(zero)
            lines2 = lines2 + ['\t{:.4e}\t{:.4e}\n'.format(real, imag)]
        lines2 = lines2 + lines1[index3:]
        lines3 = lines3 + lines2

    with open(resp_file, 'w') as outfile:
        for line in lines3:
            outfile.write(line)
    return resp_file


def __select_str_files(strong_files, tensor_info):
    r"""We pre-select strong motion data
    """
    time_shift = tensor_info['time_shift']
    centroid_lat = tensor_info['centroid_lat']
    centroid_lon = tensor_info['centroid_lon']
    depth = tensor_info['depth']
    strong_files2 = []
    distance = 2 if time_shift < 50 else 4
    syn_len = 20 + 4*time_shift + 7*depth/50

    for sac in strong_files:
        select_stat = True
        sacfile = SACTrace.read(sac)
        station_lat = sacfile.stla
        station_lon = sacfile.stlo

        st = read(sac)
        start, end = [st[0].stats.starttime, st[0].stats.endtime]
        if end - start < syn_len:
            continue
        if not __select_distance(tensor_info, station_lat, station_lon,
                                 max_distance=distance, use_centroid=True):
            continue
        strong_files2 = strong_files2 + [sac]
    if len(strong_files2) == 0:
        raise RuntimeError(
            'No strong motion waveforms selected. Either strong motion '\
            + 'waveforms are too distant, or they are shorter than '\
            + '{} '.format(syn_len)\
            + 'seconds, which is the required length for strong motions '\
            + 'for this event'
        )

    return strong_files2


def move_str_files(stations_str):
    """
    """
    new_name = lambda network, station, channel:\
        'acc_{}_{}_{}.sac'.format(network, station, channel)

    for sac in stations_str:
        st = read(sac)
        name = st[0].stats.station
        channel = st[0].stats.channel
        network = st[0].stats.network
        filename = os.path.join('STR', new_name(name, channel, network))
        st.write(filename, format='SAC', byteorder=0)
    return


def __read_paz(paz_file):
    """
    """
    is_paz = True
    with open(paz_file, 'r') as infile:
        lines = [line.split() for line in infile]
    indexes = [i for i, line in enumerate(lines) if 'ZEROS' in line]
    if len(indexes) == 0:
        return None, False
    index = indexes[-1]
    lines2 = lines[index:]
    n_zeros = int(lines2[0][1])
    n_poles = int(lines2[n_zeros + 1][1])
    gain_line = n_zeros + n_poles + 2
    if len(lines2) <= gain_line:
        gain = 1
    elif len(lines2[gain_line]) <= 1:
        gain = 1
    else:
        gain = float(lines2[gain_line][1])
    zeros = [float(real)+1j*float(imag) for real, imag in lines2[1:n_zeros + 1]]
    poles = [float(real)+1j*float(imag) for real, imag in lines2[n_zeros+2:gain_line]]
    paz_dict = {
        'zeros': zeros,
        'poles': poles,
        'gain': gain,
        'sensitivity': 1
    }
    return paz_dict, is_paz


def __remove_response_str(stations_str, response_files, logger=None):
    """Remove instrumental response.
    """
    for sac in stations_str:
        st = read(sac)
        name = st[0].stats.station
        channel = st[0].stats.channel
        network = st[0].stats.network
        pz_files = [resp for resp in response_files if name in resp]
        pz_files = [resp for resp in pz_files if channel in resp]

        if not pz_files:
            os.remove(sac)
            continue
        pzfile2 = next(iter(pz_files))
        if not os.path.isfile(pzfile2):
            os.remove(sac)
            continue

        paz_dict, is_paz = __read_paz(pzfile2)
        if not is_paz:
            if logger:
                logger.warning(
                    'PZ response unavailable at station {}, channel {}'.format(
                        name, channel))
            continue
        __convert_response_acc(pzfile2)
        paz_dict, is_paz = __read_paz(pzfile2)
        st[0].simulate(paz_remove=paz_dict)
        st.write(sac, format='SAC', byteorder=0)
    return


def process_strong_motion(strong_files, tensor_info, data_prop, logger=None):
    r"""We integrate strong motion to velocity and proceed to downsample it.

    :param tensor_info: dictionary with moment tensor information
    :param data_prop: dictionary with waveform properties
    :param strong_files: strong motion to be resampled and integrated to velocity
    :param logger: where to log results of strong motion resampling and integration
    :type tensor_info: dict
    :type data_prop: dict
    :type strong_files: list
    :type logger: logging, optional
    """
    muestreo = data_prop['sampling']
    dt_strong = muestreo['dt_strong']
    filtro_strong = data_prop['strong_filter']

    low_freq = filtro_strong['low_freq']
    high_freq = filtro_strong['high_freq']

    for sac in strong_files:
        sacheader = SACTrace.read(sac)
        sacheader.t8 = low_freq
        sacheader.t9 = high_freq
        sacheader.write(sac, byteorder='little')

    print('Baseline removal procedure')
    with Pool(processes=6) as pool:
        results = pool.map(__worker, [[str_file] for str_file in strong_files])
    strong_files2 = glob.glob('vel*')

    _filter_decimate(strong_files2, filtro_strong, dt_strong)
    return


def small_events(select_str_data):
    """A simple method for removing strong motion baselines. Based on Boore et
    al. [2002]_

    :param select_st_data: list with waveforms (in sac format) to use
    :type select_st_data: list

    This method consists on padding the record with zeros on the extremes and
    applying a low pass filter.

    .. [2002] Boore DM: On pads and filters: processing strong-motion data.
        Bulletin of the Seismological Society of America 2005,95(2):745-750.
        10.1785/0120040160
    """
    order = 4
    for sac in select_str_data:
        stream = read(sac)
        delta = stream[0].stats.delta
        t_pad = 1.5 * order / 0.02 / 2
        t_pad = int(t_pad / delta)
        fs = stream[0].stats['sampling_rate']
        if _delete_criteria(stream[0].data):
            continue
        data = stream[0].data
        data = np.pad(data, t_pad, 'constant')
        data = cumtrapz(data, dx=delta, initial=0)
        st_vel = stream.copy()
        nyq = 0.5 * fs
        low = 0.02 / nyq
        b, a = butter(4, low, btype='highpass')
#
# apply filter
#
        data = filtfilt(b, a, data)
        st_vel[0].data = data[t_pad:-t_pad]
        st_vel.write('vel_{}'.format(sac[4:]), format='SAC')
    return


def __worker(select_str_data):
    """Helper routine. Method for baseline removal based on Wang et al.
    """
    for sac in select_str_data:
        stream = read(sac)
        if _delete_criteria(stream[0].data):
            continue
        st_vel = wang1.wang_process(sac)
        if not st_vel:
            continue
        st_vel.write('vel_{}'.format(sac[4:]), format='SAC')
    return


def _delete_criteria(data):
    """
    """
    data = data - np.mean(data)
    if np.max(np.abs(data)) < 10 ** -6:
        return True
    return False


##########################
# We proceed to select strong motion stations
##########################


def select_strong_stations(tensor_info, files):
    """We select strong motion data to use in finite fault modelling

    :param files: list of waveform files (in sac format) to select
    :type files: list
    """
    if len(files) < 150:
        return files
    lat = tensor_info['lat']
    lon = tensor_info['lon']
    sacheaders = tuple([SACTrace.read(sac) for sac in files])
    streams = tuple([read(sac) for sac in files])
    names = [header.kstnm for header in sacheaders]
    azimuth0 = []
    for header in sacheaders:
        station_lat = header.stla
        station_lon = header.stlo
        dist, az, baz = mng._distazbaz(station_lat, station_lon, lat, lon)
        azimuth0 = azimuth0 + [az]
    azimuth0 = np.array(azimuth0)#[header.az for header in sacheaders])
    az0 = np.amin(azimuth0)
    az1 = np.amax(azimuth0)
    jump = int((az1 - az0) / 60) + 1
    names = list(set(names))
    select_stations = []

# Limiting the number of data. Choose at most one station (3 channels)
# for every degree in azimuth. Preferably use data with higher PGV

    for az0 in range(0, 360, jump):
        az1 = az0 + jump
        best_pgv = 0.0
        for name in names:
            fun1 = lambda header: name in header.kstnm
            header0 = next(filter(fun1, sacheaders))
            station_lat = header0.stla
            station_lon = header0.stlo
            dist, azimuth, baz = mng._distazbaz(
                station_lat, station_lon, lat, lon)
            if not az0 <= azimuth < az1:
                continue
            zipped = zip(streams, sacheaders)
            streams2 = (st for st, header in zipped if name in header.kstnm)
            pgvs = np.array([np.max(np.abs(st[0].data)) for st in streams2])
            pgv = np.sqrt(np.sum(pgvs ** 2))
            if pgv > best_pgv:
                best_pgv = pgv
                zipped = zip(sacheaders, files)
                chosen_sta = [sac for header, sac in zipped if name in header.kstnm]

        if best_pgv > 0.0:
            select_stations = select_stations + chosen_sta

    for file in files:
        if not file in select_stations:
            os.remove(file)
    return select_stations


if __name__ == '__main__':
    import argparse
    import manage_parser as mp

    parser = argparse.ArgumentParser()
    parser.add_argument("-f", "--folder", default=os.getcwd(),
                        help="folder where there are input files")
    parser = mp.parser_add_tensor(parser)
    parser = mp.parser_data_process(parser)
    parser.add_argument("-dt", "--sampling_delta",
                        help="process cGPS data")
    parser.add_argument("--pick", action="store_true",
                        help="pick of teleseismic waves")
    args = parser.parse_args()
    if not os.path.isfile('sampling_filter.json'):
        raise FileNotFoundError(
                errno.ENOENT, os.strerror(errno.ENOENT),
                'sampling_filter.json')
    data_prop = json.load(open('sampling_filter.json'))
    os.chdir(args.folder)
    if args.gcmt_tensor:
        cmt_file = args.gcmt_tensor
        tensor_info = tensor.get_tensor(cmt_file=cmt_file)
    else:
        tensor_info = tensor.get_tensor()
    if args.tele:
        tele_files = glob.glob('*BH*SAC') + glob.glob('*BH*sac')
        select_process_tele_body(tele_files, tensor_info, data_prop)
    if args.surface:
        tele_files = glob.glob('*BH*SAC') + glob.glob('*BH*sac')
        select_process_surf_tele(tele_files, tensor_info, data_prop)
    if args.strong:
        strong_files = glob.glob('*.HN*SAC') + glob.glob('*.HL*SAC')\
                    + glob.glob('*.HN*sac') + glob.glob('*.HL*sac')\
                    + glob.glob('*.AH?.*') + glob.glob('*.AH?.*')\
                    + glob.glob('*_HN*sac') + glob.glob('*_HL*sac')
        select_process_strong(strong_files, tensor_info, data_prop)
    if args.cgps:
        cgps_files = glob.glob('*.L[HXY]*SAC') + glob.glob('*.L[HXY]*sac')
        select_process_cgps(cgps_files, tensor_info, data_prop)

