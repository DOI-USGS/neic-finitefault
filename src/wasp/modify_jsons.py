import json
import argparse
import os


def modify_channels(json_file, method='delete'):
    """Method to select interactively channels to remove or downweight in
    modelling.

    :param json_file: file with different channels
    :param method: whether to downweight or remove data
    :type json_file: list
    :type method: string, optional
    """
    channels = json.load(open(json_file))
    first_message = '\nWrite station for channel removal. To exit, type exit: '
    second_message = \
            '\nWrite channel from station {} to remove. To exit, '\
            + 'write exit: '
    third_message = 'No channel was selected. Going back to station selection.'
    if method == 'downweight':
        first_message = \
                '\nWrite station for channel downweight to 0. To exit, type exit: '
        second_message = \
                '\nWrite channel from station {} to downweight to 0. To exit, '\
                + 'type exit: '
    print('List of available stations\n')
    stations = [channel['name'] for channel in channels]
    stations = list(set(stations))
    print(*stations)
    while True:
        stations = [channel['name'] for channel in channels]
        stations = list(set(stations))
        station = input(first_message)
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
        channel0 = input(second_message.format(station))
        if channel0.lower() == 'exit':
            print('Exit\n')
            break
        if channel0 == '':
            print(third_message)
            continue
        if channel0 not in channels_station:
            raise RuntimeError(
                'Selected channel does not belong to list of available channels')
        channel1 = next(channel for channel in channels\
                        if channel['name'] == station\
                        and channel['component'] == channel0)
        if method == 'downweight':
            channel1['trace_weight'] = 0.0
            channels.remove(channel1)
            channels.append(channel1)
        elif method == 'delete':
            channels.remove(channel1)
    with open(json_file,'w') as f:
         json.dump(
                 channels, f, sort_keys=True, indent=4, separators=(',', ': '),
                 ensure_ascii=False)
    return


if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument("-f", "--folder", default=os.getcwd(),
                        help="folder where there are input files")
    parser.add_argument("-del", "--delete", action="store_true",
                        help="delete selected channels")
    parser.add_argument("-dn", "--downweight", action="store_true",
                        help="give selected channels zero weight")
    parser.add_argument("-t", "--tele", action="store_true",
                        help="compute files with teleseismic data")
    parser.add_argument("-su", "--surface", action="store_true",
                        help="compute files with surface waves data")
    parser.add_argument("-st", "--strong", action="store_true",
                        help="compute files with strong motion data")
    parser.add_argument("--cgps", action="store_true",
                        help="compute files with cGPS data")
    parser.add_argument("--gps", action="store_true",
                        help="compute files with static GPS data")
    args = parser.parse_args()
    os.chdir(args.folder)
    if args.tele:
        if not os.path.isfile('tele_waves.json'):
            raise FileNotFoundError(
                    errno.ENOENT, os.strerror(errno.ENOENT), 'tele_waves.json')
        json_file = 'tele_waves.json'
    if args.surface:
        if not os.path.isfile('surf_waves.json'):
            raise FileNotFoundError(
                    errno.ENOENT, os.strerror(errno.ENOENT), 'surf_waves.json')
        json_file = 'surf_waves.json'
    if args.strong:
        if not os.path.isfile('strong_motion_waves.json'):
            raise FileNotFoundError(
                    errno.ENOENT, os.strerror(errno.ENOENT),
                    'strong_motion_waves.json')
        json_file = 'strong_motion_waves.json'
    if args.cgps:
        if not os.path.isfile('cgps_waves.json'):
            raise FileNotFoundError(
                    errno.ENOENT, os.strerror(errno.ENOENT), 'cgps_waves.json')
        json_file = 'cgps_waves.json'
    if args.delete:
        modify_channels(json_file, method='delete')
    if args.downweight:
        modify_channels(json_file, method='downweight')
