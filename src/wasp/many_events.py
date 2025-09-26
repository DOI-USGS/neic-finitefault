import json
import os
import pathlib
from typing import List, Union

import wasp.seismic_tensor as tensor


def get_waveforms_events(
    waveforms_events: List[List[dict]],
    data_type: str,
    directory: Union[pathlib.Path, str] = pathlib.Path(),
):
    """Write files that label waveforms with their event number

    :param waveforms_events: The events' waveform properties
    :type waveforms_events: List[List[dict]]
    :param data_type: The data type being associated with events
    :type data_type: str
    :param directory: The directory to write to, defaults to pathlib.Path()
    :type directory: Union[pathlib.Path, str], optional
    """
    directory = pathlib.Path(directory)
    if "cgnss" in data_type:
        file_name = "cgnss_events.txt"
        dict_name = "cgnss_waves.json"
    if "strong" in data_type:
        file_name = "strong_motion_events.txt"
        dict_name = "strong_motion_waves.json"
    if "body" in data_type:
        file_name = "tele_events.txt"
        dict_name = "tele_waves.json"
    if "surf" in data_type:
        file_name = "surf_events.txt"
        dict_name = "surf_waves.json"
    if "gnss" in data_type:
        file_name = "static_events.txt"
        dict_name = "static_data.json"
    get_waveforms_events2(waveforms_events, directory / file_name)
    save_dict(waveforms_events, directory / dict_name)


def get_waveforms_events2(
    waveforms_events: List[List[dict]], file_name: Union[pathlib.Path, str]
):
    """Write a text file listing the waveforms with their component and event

    :param waveforms_events: The list of events' waveform properties
    :type waveforms_events: List[List[dict]]
    :param file_name: The full path to file to write
    :type file_name: Union[pathlib.Path, str]
    """
    with open(file_name, "w") as outf:
        for i, traces in enumerate(waveforms_events):
            for trace in traces:
                name = trace["name"]
                component = trace["component"]
                if len(component) == 0:
                    component = "None"
                outf.write("{} {} {}\n".format(name, component, i + 1))


def save_dict(
    elements: List[List[dict]], file_name: Union[pathlib.Path, str]
) -> List[dict]:
    """Save a file that associates waveforms with the event number

    :param elements: The list of events' waveform lists
    :type elements: List[List[dict]]
    :param file_name: The full path to the file to be written
    :type file_name: Union[pathlib.Path,str]
    :return: The updated list of waveforms associated to their events
    :rtype: List[dict]
    """
    new_list: List[dict] = []
    for i, element in enumerate(elements):
        for value in element:
            value["event"] = i + 1
        new_list = new_list + element
    with open(file_name, "w") as write_file:
        json.dump(
            new_list, write_file, indent=4, separators=(",", ": "), sort_keys=True
        )
    return new_list


def get_segments_events(
    segments_events: List[dict],
    tensors: List[dict],
    directory: Union[pathlib.Path, str] = pathlib.Path(),
):
    """Write a list of segments and their corresponding events and a
    full segments_data dictionary

    :param segments_events: The list of segments data for each event
    :type segments_events: List[dict]
    :param tensors: The list of tensors for each event
    :type tensors: List[dict]
    :param directory: The directory to write the files to, defaults to pathlib.Path()
    :type directory: Union[pathlib.Path, str], optional
    """
    directory = pathlib.Path(directory)
    new_connections: list = []
    rise_time = segments_events[0]["rise_time"]
    new_segments_data = {
        "rise_time": rise_time,
        "segments": [],
    }
    new_segments: List[dict] = []
    segments_events2: List[dict] = []
    zipped = zip(segments_events, tensors)
    segment0 = 0
    for i, (segments_data, tensor) in enumerate(zipped):
        lat = tensor["lat"]
        lon = tensor["lon"]
        depth = tensor["depth"]
        hypocenter = {"lat": lat, "lon": lon, "depth": depth}
        segments = segments_data["segments"]
        for j, segment in enumerate(segments):
            new_dict = {"segment": j + segment0 + 1, "event": i + 1}
            segments_events2 = segments_events2 + [new_dict]
            segment["event"] = i + 1
        segments[0]["hypocenter"] = hypocenter
        new_segments = new_segments + segments
        if not "connections" in segments_data:
            segment0 = segment0 + len(segments)
        else:
            connections = segments_data["connections"]
            for connection in connections:
                segment1 = connection["segment1"]
                segment2 = connection["segment2"]
                connection["segment1"] = segment1 + segment0
                connection["segment2"] = segment2 + segment0
            new_connections = new_connections + connections
            segment0 = segment0 + len(segments)
    new_segments_data["segments"] = new_segments
    if len(new_connections) > 0:
        new_segments_data["connections"] = new_connections

    with open(directory / "segments_events.txt", "w") as outf:
        outf.write("segment event\n")
        for segment_event in segments_events2:
            segment = segment_event["segment"]
            event = segment_event["event"]
            outf.write("{} {}\n".format(segment, event))
    with open(directory / "segments_data.json", "w") as write_file:
        json.dump(
            new_segments_data,
            write_file,
            indent=4,
            separators=(",", ": "),
            sort_keys=True,
        )


def get_model_space_events(
    model_events: List[List[dict]], directory: Union[pathlib.Path, str] = pathlib.Path()
):
    """Get model information for events

    :param model_events: The models corresponding to events
    :type model_events: List[List[dict]]
    :param directory: The directory where to write the file, defaults to pathlib.Path()
    :type directory: Union[pathlib.Path, str], optional
    """
    directory = pathlib.Path(directory)
    new_model_spaces: List[dict] = []
    segment0 = 0
    for i, model_space in enumerate(model_events):
        new_model_space: List[dict] = []
        for model_segment in model_space:
            regularization = model_segment["regularization"]
            neighbour_down = regularization["neighbour_down"]
            if neighbour_down is not None:
                segment = neighbour_down["segment"]
                segment = segment + segment0
                neighbour_down["segment"] = segment
            neighbour_left = regularization["neighbour_left"]
            if neighbour_left is not None:
                segment = neighbour_left["segment"]
                segment = segment + segment0
                neighbour_left["segment"] = segment
            neighbour_right = regularization["neighbour_right"]
            if neighbour_right is not None:
                segment = neighbour_right["segment"]
                segment = segment + segment0
                neighbour_right["segment"] = segment
            neighbour_up = regularization["neighbour_up"]
            if neighbour_up is not None:
                segment = neighbour_up["segment"]
                segment = segment + segment0
                neighbour_up["segment"] = segment
            regularization = {
                "neighbour_down": neighbour_down,
                "neighbour_left": neighbour_left,
                "neighbour_right": neighbour_right,
                "neighbour_up": neighbour_up,
            }
            model_segment["regularization"] = regularization
            new_model_space = new_model_space + [model_segment]
        segment0 = segment0 + len(model_space)
        new_model_spaces = new_model_spaces + new_model_space
    with open(directory / "model_space.json", "w") as write_file:
        json.dump(
            new_model_spaces,
            write_file,
            indent=4,
            separators=(",", ": "),
            sort_keys=True,
        )


def get_moment_events(
    tensors: List[dict], directory: Union[pathlib.Path, str] = pathlib.Path()
):
    """Write a list of events' seismic moments

    :param tensors: The list of events' tensor properties
    :type tensors: List[dict]
    :param directory: The directory to write to, defaults to pathlib.Path()
    :type directory: Union[pathlib.Path, str], optional
    """
    directory = pathlib.Path(directory)
    seismic_moments: List[dict] = []
    for tensor_info in tensors:
        moment = tensor_info["seismic_moment"]
        seismic_moments = seismic_moments + [moment]

    with open(directory / "moment_events.txt", "w") as outf:
        for moment in seismic_moments:
            outf.write("{}\n".format(moment))


def select_segments_event(segments_data: dict, event: int) -> dict:
    """Select segments corresponding to the specified event

    :param segments_data: The segments dictionary
    :type segments_data: dict
    :param event: The event number
    :type event: int
    :return: The Updated segments data
    :rtype: dict
    """
    rise_time = segments_data["rise_time"]
    new_segments_data = {
        "rise_time": rise_time,
        "segments": [],
    }
    new_segments: List[dict] = []
    segments = segments_data["segments"]
    for segment in segments:
        if segment["event"] == event:
            new_segments = new_segments + [segment]
    new_segments_data["segments"] = new_segments
    return new_segments_data


def select_waveforms_event(traces_info: List[dict], event: int) -> List[dict]:
    """Select waveforms corresponding to the specified event

    :param traces_info: The trace properties
    :type traces_info: List[dict]
    :param event: The event number
    :type event: int
    :return: The updated list of waveforms
    :rtype: List[dict]
    """
    new_traces: List[dict] = []
    for trace in traces_info:
        if trace["event"] == event:
            new_traces = new_traces + [trace]
    return new_traces
