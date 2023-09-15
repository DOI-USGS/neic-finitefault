import json
import os
import pathlib
from typing import Dict, List, Literal, Union

import numpy as np
from matplotlib import pyplot as plt  # type: ignore
from obspy import read  # type: ignore


def _modify_by_dict(
    json_file: Union[pathlib.Path, str],
    station_dict: Dict[str, List[Dict[str, Union[float, int, str]]]],
):
    """Modify the channels using a station dictionary

    :param json_file: File with channels
    :type json_file: Union[pathlib.Path, str]
    :param station_dict: The station dictionary to use instead of user input,
                            dict format:
                            {
                                <station>: {
                                    <channel>: {
                                        "baseline_correction": float or None,
                                        "time_correction": float or None,
                                    },
                                    ...
                                },
                                ...
                            }
    :type station_dict: Dict[str, List[Dict[str, Union[float, int, str]]]]
    """
    with open(json_file, "r") as f:
        channels = json.load(f)
    for idx, channel in enumerate(channels):
        name = channel["name"]
        component = channel["component"]
        if name in station_dict and component in station_dict[name]:
            modifications = station_dict[name][component]
            file_name = channel["file"]
            st = read(file_name)
            delta = st[0].stats.delta
            if modifications.get("baseline_correction") is not None:
                st[0].data = st[0].data + modifications["baseline_correction"]
                st.write(file_name, format="SAC", byteorder=0)
            if modifications.get("time_correction") is not None:
                time_shift = int(modifications["time_correction"] / delta)
                channels[idx]["start_signal"] = channel["start_signal"] - time_shift
    with open(json_file, "w") as f:
        json.dump(channels, f)


def _modify_by_input(
    json_file: Union[pathlib.Path, str],
):
    """Modify the channels using realtime user input

    :param json_file: File with channels
    :type json_file: Union[pathlib.Path, str]
    :raises ValueError: If a specified station does not exist
    :raises ValueError: If a specified channel does not exist
    :raises ValueError: If a correction provided is not a number
    """
    with open(json_file, "r") as f:
        channels = json.load(f)
    message1 = "\nSelect station for modification. To exit, type exit: "
    message2 = "\nSelect channel from station {} to modify. To exit, " + "write exit: "
    message3 = (
        "\nSpecify time correction. To shift waveform backward, "
        + "enter a negative number: "
    )
    message4 = "\nNo time correction selected. " + "Proceed to baseline correction."
    message5 = (
        "\nSpecify baseline correction. To shift waveform down, "
        + "enter a negative number: "
    )
    message6 = (
        "\nNo baseline correction selected. " + "Going back to station selection."
    )
    message7 = "No channel was selected. Going back to station selection."
    print("List of available stations\n")
    stations = [channel["name"] for channel in channels]
    stations = list(set(stations))
    print(*stations)
    while True:
        stations = [channel["name"] for channel in channels]
        stations = list(set(stations))
        station = input(message1)
        if station.lower() == "exit":
            print("Exit\n")
            break
        if station not in stations:
            raise RuntimeError(
                "Selected station does not belong to list of available stations"
            )
        channels2 = [channel for channel in channels if channel["name"] == station]
        channels_station = [channel["component"] for channel in channels2]
        print("List of available channels for station {}\n".format(station))
        print(*channels_station)
        channel0 = input(message2.format(station))
        if channel0.lower() == "exit":
            print("Exit\n")
            break
        if channel0 == "":
            print(message7)
            continue
        if channel0 not in channels_station:
            raise RuntimeError(
                "Selected channel does not belong to list of available channels"
            )
        channel1 = next(
            channel
            for channel in channels
            if channel["name"] == station and channel["component"] == channel0
        )
        file = channel1["file"]
        st = read(file)
        delta = st[0].stats.delta
        start = channel1["start_signal"]
        time_shift: Union[float, str] = input(message3)
        shift = False
        if not __is_number(time_shift):
            if str(time_shift).lower() == "exit":
                print("Exit\n")
                break
            elif time_shift == "":
                print(message4)
            else:
                raise RuntimeError("Invalid time correction value.")
        else:
            time_shift = float(time_shift)
            shift = True
        baseline_shift: Union[float, str] = input(message5)
        shift2 = False
        if not __is_number(baseline_shift):
            if str(baseline_shift).lower() == "exit":
                print("Exit\n")
                break
            elif baseline_shift == "":
                print(message6)
            else:
                raise RuntimeError("Invalid baseline correction value.")
        else:
            baseline_shift = float(baseline_shift)
            shift2 = True
        if shift:
            time_shift2 = int(time_shift / delta)
            channel1["start_signal"] = channel1["start_signal"] - time_shift2
        if shift2:
            st[0].data = st[0].data + baseline_shift
        st.write(file, format="SAC", byteorder=0)
        channels.remove(channel1)
        channels.append(channel1)

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


def correct_waveforms(
    json_file: Union[pathlib.Path, str],
    input: bool = True,
    station_dict: Dict[str, List[Dict[str, Union[float, int, str]]]] = {},
):
    """_summary_

    :param json_file: File with channels
    :type json_file: Union[pathlib.Path, str]
    :param input: Whether or not to use realtime user input, defaults to True
    :type input: bool, optional
    :param station_dict: The station dictionary to use instead of user input,
                            dict format:
                            {
                                <station>: {
                                    <channel>: {
                                        "baseline_correction": float or None,
                                        "time_correction": float or None,
                                    },
                                    ...
                                },
                                ...
                            }
                            defaults to {}
    :type station_dict: Dict[str, List[Dict[str, Union[float, int, str]]]], optional
    :raises ValueError: _description_
    """
    if input:
        _modify_by_input(json_file=json_file)
    elif station_dict == {}:
        raise ValueError(f"If input is not used, a station_dict must be provided.")
    else:
        _modify_by_dict(json_file=json_file, station_dict=station_dict)


def plot_channels(
    json_file: pathlib.Path,
    plot_directory: Union[pathlib.Path, str],
):
    """Plot channels

    :param json_file: The full path to the json file with the available files
    :type json_file: str
    :param plot_directory: Where plots (and plot folders) should be saved
    :type plot_directory: Union[pathlib.Path, str]
    """
    with open(json_file, "r") as f:
        files = json.load(f)

    plot_folder: Union[pathlib.Path, str] = (
        "review_tele"
        if os.path.basename(json_file) == "tele_waves.json"
        else "review_strong"
        if os.path.basename(json_file) == "strong_motion_waves.json"
        else "review_surf"
        if os.path.basename(json_file) == "surf_waves.json"
        else "review_manual"
    )
    plot_folder = pathlib.Path(plot_directory) / plot_folder
    if not os.path.isdir(plot_folder):
        os.mkdir(plot_folder)

    dt = float(files[0]["dt"])
    for file in files:
        length = int(float(file["duration"]))
        length2 = int(10 / dt)
        start0 = 0
        start00 = 0
        name = file["name"]
        component = file["component"]
        stream = read(file["file"])
        obser_tr = stream[0].data
        start = int(file["start_signal"])
        fig = plt.figure(figsize=(10, 10))
        ax = fig.gca()
        fig.suptitle("{} {}".format(name, component))
        start2 = max(0, start - length2)
        start3 = start - start2
        observed0 = obser_tr[start2 : start2 + length]
        time1 = np.arange(-start3, len(observed0) - start3) * dt
        min_val = np.min(observed0)
        max_val = np.max(observed0)
        ax.plot(time1, observed0)
        ax.axhline(color="k")
        ax.axvline(color="k")
        ax.set_title("Comparison of observed and synthetic")
        name_file = os.path.join(plot_folder, "{}_{}.png".format(name, component))
        plt.savefig(name_file)
        plt.close(fig)


def __is_number(value: Union[str, float]) -> bool:
    """Determine if the provided value is a number

    :param value: The value
    :type value: str
    :return: Whether the provided value is a number
    :rtype: bool
    """
    try:
        float(value)
        return True
    except ValueError:
        return False
