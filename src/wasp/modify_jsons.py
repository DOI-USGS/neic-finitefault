import argparse
import errno
import json
import os
import pathlib
from copy import copy
from typing import Dict, List, Literal, Union


def _modify_by_dict(
    json_file: Union[pathlib.Path, str],
    method: Literal["downweight", "delete"],
    station_dict: Dict[str, List[str]],
):
    """Modify the channels using realtime a station dictionary

    :param json_file: File with channels
    :type json_file: Union[pathlib.Path, str]
    :param method: Whether to downweight or delete the channel
    :type method: Literal[&quot;downweight&quot;, &quot;delete&quot;]
    :param station_dict: The station dictionary to use instead of user input,
                            dict format: <station>: [<channel(s)>]
    :type station_dict: Dict[str, List[str]]
    :raises ValueError: If a specified station does not exist
    :raises ValueError: If a specified channel does not exist
    """
    with open(json_file, "r") as f:
        channels = json.load(f)
    for station, chosen_channels in station_dict.items():
        channels_to_search = copy(chosen_channels)
        for idx in range(len(channels) - 1, -1, -1):
            channel = channels[idx]
            if channel["name"] == station and channel["component"] in chosen_channels:
                if method == "downweight":
                    channels[idx]["trace_weight"] = 0
                else:
                    del channels[idx]
                channels_to_search = [
                    c for c in channels_to_search if c != channel["component"]
                ]
            if not len(channels_to_search):
                break
    with open(json_file, "w") as f:
        json.dump(channels, f)


def _modify_by_input(
    json_file: Union[pathlib.Path, str],
    method: Literal["downweight", "delete"],
):
    """Modify the channels using realtime user input

    :param json_file: File with channels
    :type json_file: Union[pathlib.Path, str]
    :param method: Whether to downweight or delete the channel
    :type method: Literal[&quot;downweight&quot;, &quot;delete&quot;]
    :raises ValueError: If a specified station does not exist
    :raises ValueError: If a specified channel does not exist
    """
    channels = json.load(open(json_file))
    first_message = "\nWrite station for channel removal. To exit, type exit: "
    second_message = (
        "\nWrite channel from station {} to remove. To exit, " + "write exit: "
    )
    third_message = "No channel was selected. Going back to station selection."
    if method == "downweight":
        first_message = (
            "\nWrite station for channel downweight to 0. To exit, type exit: "
        )
        second_message = (
            "\nWrite channel from station {} to downweight to 0. To exit, "
            + "type exit: "
        )
    print("List of available stations\n")
    stations = [channel["name"] for channel in channels]
    stations = list(set(stations))
    print(*stations)
    while True:
        stations = [channel["name"] for channel in channels]
        stations = list(set(stations))
        station = input(first_message)
        if station.lower() == "exit":
            print("Exit\n")
            break
        if station not in stations:
            raise ValueError(
                "Selected station does not belong to list of available stations"
            )
        channels2 = [channel for channel in channels if channel["name"] == station]
        channels_station = [channel["component"] for channel in channels2]
        print("List of available channels for station {}\n".format(station))
        print(*channels_station)
        channel0 = input(second_message.format(station))
        if channel0.lower() == "exit":
            print("Exit\n")
            break
        if channel0 == "":
            print(third_message)
            continue
        if channel0 not in channels_station:
            raise ValueError(
                "Selected channel does not belong to list of available channels"
            )
        channel1 = next(
            channel
            for channel in channels
            if channel["name"] == station and channel["component"] == channel0
        )
        if method == "downweight":
            channel1["trace_weight"] = 0.0
            channels.remove(channel1)
            channels.append(channel1)
        elif method == "delete":
            channels.remove(channel1)
    with open(json_file, "w") as f:
        json.dump(
            channels,
            f,
            sort_keys=True,
            indent=4,
            separators=(",", ": "),
            ensure_ascii=False,
        )
    return


def modify_channels(
    json_file: Union[pathlib.Path, str],
    method: Literal["downweight", "delete"],
    input: bool = True,
    station_dict: Dict[str, List[str]] = {},
):
    """Method to select interactively channels to remove or downweight in
       modelling

    :param json_file: File with channels
    :type json_file: Union[pathlib.Path, str]
    :param method: Whether to downweight or delete the channel
    :type method: Literal[&quot;downweight&quot;, &quot;delete&quot;]
    :param input: Whether or not to use user input or the station_dict, defaults to True
    :type input: bool, optional
    :param station_dict: The station dictionary to use instead of user input,
                            dict format: <station>: [<channel(s)>], defaults to {}
    :type station_dict: Dict[str, List[str]], optional
    :raises ValueError: If input is False but no station dict is provided
    """
    if input:
        _modify_by_input(json_file=json_file, method=method)
    elif station_dict == {}:
        raise ValueError(f"If input is not used, a station_dict must be provided.")
    else:
        _modify_by_dict(json_file=json_file, method=method, station_dict=station_dict)


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "-f", "--folder", default=os.getcwd(), help="folder where there are input files"
    )
    parser.add_argument(
        "-del", "--delete", action="store_true", help="delete selected channels"
    )
    parser.add_argument(
        "-dn",
        "--downweight",
        action="store_true",
        help="give selected channels zero weight",
    )
    parser.add_argument(
        "-t", "--tele", action="store_true", help="compute files with teleseismic data"
    )
    parser.add_argument(
        "-su",
        "--surface",
        action="store_true",
        help="compute files with surface waves data",
    )
    parser.add_argument(
        "-st",
        "--strong",
        action="store_true",
        help="compute files with strong motion data",
    )
    parser.add_argument(
        "--cgps", action="store_true", help="compute files with cGPS data"
    )
    parser.add_argument(
        "--gps", action="store_true", help="compute files with static GPS data"
    )
    args = parser.parse_args()
    os.chdir(args.folder)
    if args.tele:
        if not os.path.isfile("tele_waves.json"):
            raise FileNotFoundError(
                errno.ENOENT, os.strerror(errno.ENOENT), "tele_waves.json"
            )
        json_file = "tele_waves.json"
    if args.surface:
        if not os.path.isfile("surf_waves.json"):
            raise FileNotFoundError(
                errno.ENOENT, os.strerror(errno.ENOENT), "surf_waves.json"
            )
        json_file = "surf_waves.json"
    if args.strong:
        if not os.path.isfile("strong_motion_waves.json"):
            raise FileNotFoundError(
                errno.ENOENT, os.strerror(errno.ENOENT), "strong_motion_waves.json"
            )
        json_file = "strong_motion_waves.json"
    if args.cgps:
        if not os.path.isfile("cgps_waves.json"):
            raise FileNotFoundError(
                errno.ENOENT, os.strerror(errno.ENOENT), "cgps_waves.json"
            )
        json_file = "cgps_waves.json"
    if args.delete:
        modify_channels(json_file, method="delete")
    if args.downweight:
        modify_channels(json_file, method="downweight")
