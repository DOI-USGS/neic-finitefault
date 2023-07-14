import json
import os
import seismic_tensor as tensor


def get_waveforms_events(waveforms_events, data_type):
    """
    """
    if 'cgps' in data_type:
        file_name = 'cgps_events.txt'
        dict_name = 'cgps_waves.json'
    if 'strong_motion' in data_type:
        file_name = 'strong_motion_events.txt'
        dict_name = 'strong_motion_waves.json'
    if 'tele_body' in data_type:
        file_name = 'tele_events.txt'
        dict_name = 'tele_waves.json'
    if 'surf_tele' in data_type:
        file_name = 'surf_events.txt'
        dict_name = 'surf_waves.json'
    if 'gps' in data_type:
        file_name = 'static_events.txt'
        dict_name = 'static_data.json'
    get_waveforms_events2(waveforms_events, file_name)
    save_dict(waveforms_events, dict_name)


def get_waveforms_events2(waveforms_events, file_name):
    """
    """
    with open(file_name, 'w') as outf:
        for i, traces in enumerate(waveforms_events):
            for trace in traces:
                name = trace['name']
                component = trace['component']
                if len(component) == 0:
                    component = 'None'
                outf.write('{} {} {}\n'.format(name, component, i + 1))


def save_dict(elements, file_name):
    """
    """
    new_list = []
    for i, element in enumerate(elements):
        for value in element:
            value['event'] = i + 1
        new_list = new_list + element
    with open(file_name, "w") as write_file:
        json.dump(new_list, write_file, indent=4,
                  separators=(',', ': '), sort_keys=True)
    return new_list


def get_segments_events(segments_events, tensors):
    """
    """
    print(json.dumps(segments_events, indent=4,
                     separators=(',', ': '), sort_keys=True))
    new_connections = []
    rise_time = segments_events[0]['rise_time']
    new_segments_data = {
        'rise_time': rise_time,
        'segments': [],
    }
    new_segments = []
    segments_events2 = []
    zipped = zip(segments_events, tensors)
    segment0 = 0
    for i, (segments_data, tensor) in enumerate(zipped):
        lat = tensor['lat']
        lon = tensor['lon']
        depth = tensor['depth']
        hypocenter = {'lat': lat, 'lon': lon, 'depth': depth}
        segments = segments_data['segments']
        for j, segment in enumerate(segments):
            new_dict = {
                'segment': j + segment0 + 1,
                'event': i + 1
            }
            segments_events2 = segments_events2 + [new_dict]
            segment['event'] = i + 1
        segments[0]['hypocenter'] = hypocenter
        new_segments = new_segments + segments
        if not 'connections' in segments_data:
            segment0 = segment0 + len(segments)
        else:
            connections = segments_data['connections']
            for connection in connections:
                segment1 = connection['segment1']
                segment2 = connection['segment2']
                connection['segment1'] = segment1 + segment0
                connection['segment2'] = segment2 + segment0
            new_connections = new_connections + connections
            segment0 = segment0 + len(segments)
    new_segments_data['segments'] = new_segments
    if len(new_connections) > 0:
        new_segments_data['connections'] = new_connections

    with open('segments_events.txt', 'w') as outf:
        outf.write('segment event\n')
        for segment_event in segments_events2:
            segment = segment_event['segment']
            event = segment_event['event']
            outf.write('{} {}\n'.format(segment, event))
    with open('segments_data.json', "w") as write_file:
        json.dump(new_segments_data, write_file, indent=4,
                  separators=(',', ': '), sort_keys=True)


def get_model_space_events(model_events):
    """
    """
    new_model_spaces = []
    segment0 = 0
    for i, model_space in enumerate(model_events):
        new_model_space = []
        for model_segment in model_space:
            regularization = model_segment['regularization']
            neighbour_down = regularization["neighbour_down"]
            if neighbour_down is not None:
                segment = neighbour_down['segment']
                segment = segment + segment0
                neighbour_down['segment'] = segment
            neighbour_left = regularization["neighbour_left"]
            if neighbour_left is not None:
                segment = neighbour_left['segment']
                segment = segment + segment0
                neighbour_left['segment'] = segment
            neighbour_right = regularization["neighbour_right"]
            if neighbour_right is not None:
                segment = neighbour_right['segment']
                segment = segment + segment0
                neighbour_right['segment'] = segment
            neighbour_up = regularization["neighbour_up"]
            if neighbour_up is not None:
                segment = neighbour_up['segment']
                segment = segment + segment0
                neighbour_up['segment'] = segment
            regularization = {
                "neighbour_down" : neighbour_down,
                "neighbour_left" : neighbour_left,
                "neighbour_right" : neighbour_right,
                "neighbour_up" : neighbour_up
            }
            model_segment['regularization'] = regularization
            new_model_space = new_model_space + [model_segment]
        segment0 = segment0 + len(model_space)
        new_model_spaces = new_model_spaces + new_model_space
    with open('model_space.json', "w") as write_file:
        json.dump(new_model_spaces, write_file, indent=4,
                  separators=(',', ': '), sort_keys=True)


def get_moment_events(tensors):
    """
    """
    seismic_moments = []
    for tensor_info in tensors:
        moment = tensor_info['seismic_moment']
        seismic_moments = seismic_moments + [moment]

    with open('moment_events.txt', 'w') as outf:
        for moment in seismic_moments:
            outf.write('{}\n'.format(moment))


def select_segments_event(segments_data, event):
    """
    """
    rise_time = segments_data['rise_time']
    new_segments_data = {
        'rise_time': rise_time,
        'segments': [],
    }
    new_segments = []
    segments = segments_data['segments']
    for segment in segments:
        if segment['event'] == event:
            new_segments = new_segments + segment
    new_segments_data['segments'] = new_segments
    return new_segments_data


def select_waveforms_event(traces_info, event):
    """
    """
    new_traces = []
    for trace in traces_info:
        if trace['event'] == event:
            new_traces = new_traces + [trace]
    return new_traces


if __name__ == '__main__':
    import argparse
    import manage_parser as mp

    parser = argparse.ArgumentParser()
    parser.add_argument(
        "-f", "--folder", default=os.getcwd(),
        help="folder where there are input files")
    parser = mp.parser_add_tensor(parser)
    parser = mp.parser_fill_data_files(parser)
    parser.add_argument(
        "-pevs", "--path_events", action="append",
        help="add inputs for modelling many events")
    parser.add_argument(
        "-mgcmt", "--many_gcmt_tensor", action='append',
        help="location of GCMT moment tensor file")
    args = parser.parse_args()
    this_folder = os.path.abspath(args.folder)
    os.chdir(args.folder)
    data_type = mp.get_used_data(args)
    if args.gcmt_tensor:
        cmt_file = args.gcmt_tensor
        tensor_info = tensor.get_tensor(cmt_file=cmt_file)
    else:
        tensor_info = tensor.get_tensor()
    if args.path_events:
        tensors = []
        segments_events = []
        model_events = []
        annealing_events = []
        strong_motion_events = []
        tele_events = []
        surf_events = []
        cgps_events = []
        static_events = []
        for path_event in args.path_events:
            os.chdir(path_event)
            tensor = json.load(open('tensor_info.json'))
            segments_data = json.load(open('segments_data.json'))
            annealing_data = json.load(open('annealing_prop.json'))
            model_space = json.load(open('model_space.json'))
            if 'strong_motion' in data_type:
                traces_strong = json.load(open('strong_motion_waves.json'))
                strong_motion_events = strong_motion_events + [traces_strong]
            if 'tele_body' in data_type:
                traces_tele = json.load(open('tele_waves.json'))
                tele_events = tele_events + [traces_tele]
            if 'surf_tele' in data_type:
                traces_surf = json.load(open('surf_waves.json'))
                surf_events = surf_events + [traces_surf]
            if 'cgps' in data_type:
                traces_cgps = json.load(open('cgps_waves.json'))
                cgps_events = cgps_events + [traces_cgps]
            if 'gps' in data_type:
                static_data = json.load(open('static_data.json'))
                static_events = static_events + [static_data]
            segments_events = segments_events + [segments_data]
            annealing_events = annealing_events + [annealing_data]
            model_events = model_events + [model_space]
            tensors = tensors + [tensor]
            os.chdir(this_folder)
        get_moment_events(annealing_events)
        get_model_space_events(model_events)
        get_segments_events(segments_events, tensors)
        if 'strong_motion' in data_type:
            get_waveforms_events(strong_motion_events, 'strong_motion')
        if 'cgps' in data_type:
            get_waveforms_events(cgps_events, 'cgps')
        if 'tele_body' in data_type:
            get_waveforms_events(tele_events, 'tele_body')
        if 'surf_tele' in data_type:
            get_waveforms_events(surf_events, 'surf_tele')
        if 'gps' in data_type:
            get_waveforms_events(static_events, 'gps')
