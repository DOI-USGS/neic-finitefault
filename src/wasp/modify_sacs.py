import json
import argparse
import os
from obspy import read
import numpy as np
from matplotlib import pyplot as plt
import get_outputs


def correct_waveforms(json_file, plot=True):
    """
    """
    channels = json.load(open(json_file))
    message1 = '\nSelect station for modification. To exit, type exit: '
    message2 = \
            '\nSelect channel from station {} to modify. To exit, '\
            + 'write exit: '
    message3 = '\nSpecify time correction. To shift waveform backward, '\
                + 'enter a negative number: '
    message4 = '\nNo time correction selected. '\
                + 'Proceed to baseline correction.'
    message5 = '\nSpecify baseline correction. To shift waveform down, '\
                + 'enter a negative number: '
    message6 = '\nNo baseline correction selected. '\
                + 'Going back to station selection.'
    message7 = 'No channel was selected. Going back to station selection.'
    print('List of available stations\n')
    stations = [channel['name'] for channel in channels]
    stations = list(set(stations))
    print(*stations)
    while True:
        stations = [channel['name'] for channel in channels]
        stations = list(set(stations))
        station = input(message1)
        if station.lower() == 'exit':
            print('Exit\n')
            break
        if station not in stations:
            raise RuntimeError(
                'Selected station does not belong to list of available stations')
        channels2 = [channel for channel in channels if channel['name'] == station]
        channels_station = [channel['component'] for channel in channels2]
        print('List of available channels for station {}\n'.format(station))
        print(*channels_station)
        channel0 = input(message2.format(station))
        if channel0.lower() == 'exit':
            print('Exit\n')
            break
        if channel0 == '':
            print(message7)
            continue
        if channel0 not in channels_station:
            raise RuntimeError(
                'Selected channel does not belong to list of available channels')
        channel1 = next(channel for channel in channels\
                        if channel['name'] == station\
                        and channel['component'] == channel0)
        file = channel1['file']
        st = read(file)
        delta = st[0].stats.delta
        start = channel1['start_signal']
        time_shift = input(message3)
        shift = False
        if not __is_number(time_shift):
            if time_shift.lower() == 'exit':
                print('Exit\n')
                break
            elif time_shift == '':
                print(message4)
            else:
                raise RuntimeError('Invalid time correction value.')
        else:
            time_shift = float(time_shift)
            shift = True
        baseline_shift = input(message5)
        shift2 = False
        if not __is_number(baseline_shift):
            if baseline_shift.lower() == 'exit':
                print('Exit\n')
                break
            elif baseline_shift == '':
                print(message6)
            else:
                raise RuntimeError('Invalid baseline correction value.')
        else:
            baseline_shift = float(baseline_shift)
            shift2 = True
        if shift:
            time_shift2 = int(time_shift / delta)
            channel1['start_signal'] = channel1['start_signal'] - time_shift2
        if shift2:
            st[0].data = st[0].data + baseline_shift
        st.write(file, format='SAC', byteorder = 0)
        channels.remove(channel1)
        channels.append(channel1)

    with open(json_file,'w') as f:
         json.dump(
                 channels, f, sort_keys=True, indent=4, separators=(',', ': '),
                 ensure_ascii=False)
    return


def plot_channels(data_type):
    """
    """
    files = json.load(open(json_file))
    plot_folder = 'review_tele' if json_file == 'tele_waves.json'\
        else 'review_strong' if json_file == 'strong_motion_waves.json'\
        else 'review_surf' if json_file == 'surf_waves.json'\
        else 'review_manual'
    if not os.path.isdir(plot_folder):
        os.mkdir(plot_folder)

    dt = float(files[0]['dt'])
    for file in files:
        length = int(float(file['duration']))
        length2 = int(10 / dt)
        start0 = 0
        start00 = 0
        name = file['name']
        component = file['component']
        stream = read(file['file'])
        obser_tr = stream[0].data
        start = int(file['start_signal'])
        fig = plt.figure(figsize=(10, 10))
        ax = fig.gca()
        fig.suptitle('{} {}'.format(name, component))
        start2 = max(0, start - length2)
        start3 = start - start2
        observed0 = obser_tr[start2:start2 + length]
        time1 = np.arange(-start3, len(observed0) - start3) * dt
        min_val = np.min(observed0)
        max_val = np.max(observed0)
        ax.plot(time1, observed0)
        ax.axhline(color='k')
        ax.axvline(color='k')
        ax.set_title('Comparison of observed and synthetic')
        name_file = os.path.join(
            plot_folder, '{}_{}.png'.format(name, component))
        plt.savefig(name_file)
        plt.close(fig)


def main(json_file, plot=False, modify=False):
    """
    """
    if plot:
        plot_channels(json_file)
    if modify:
        correct_waveforms(json_file)
    return


def __is_number(string):
    try:
        float(string)
        return True
    except ValueError:
        return False


if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument("-f", "--folder", default=os.getcwd(),
                        help="folder where there are input files")
    parser.add_argument("-t", "--tele", action="store_true",
                        help="select teleseismic data")
    parser.add_argument("-su", "--surface", action="store_true",
                        help="select surface waves data")
    parser.add_argument("-st", "--strong", action="store_true",
                        help="select strong motion data")
    parser.add_argument("--cgps", action="store_true",
                        help="select cGPS data")
    parser.add_argument("-p", "--plot", action="store_true",
                        help="create more detailed plots of waveforms")
    parser.add_argument("-m", "--modify", action="store_true",
                        help="manual modification of waveforms")
    args = parser.parse_args()
    os.chdir(args.folder)
    if args.tele:
        if not os.path.isfile('tele_waves.json'):
            raise FileNotFoundError(
                    errno.ENOENT, os.strerror(errno.ENOENT), 'tele_waves.json')
        json_file = 'tele_waves.json'
        main(json_file, plot=args.plot, modify=args.modify)
    if args.surface:
        if not os.path.isfile('surf_waves.json'):
            raise FileNotFoundError(
                    errno.ENOENT, os.strerror(errno.ENOENT), 'surf_waves.json')
        json_file = 'surf_waves.json'
        main(json_file, plot=args.plot, modify=args.modify)
    if args.strong:
        if not os.path.isfile('strong_motion_waves.json'):
            raise FileNotFoundError(
                    errno.ENOENT, os.strerror(errno.ENOENT),
                    'strong_motion_waves.json')
        json_file = 'strong_motion_waves.json'
        main(json_file, plot=args.plot, modify=args.modify)
    if args.cgps:
        if not os.path.isfile('cgps_waves.json'):
            raise FileNotFoundError(
                    errno.ENOENT, os.strerror(errno.ENOENT), 'cgps_waves.json')
        json_file = 'cgps_waves.json'
        main(json_file, plot=args.plot, modify=args.modify)