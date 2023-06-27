# -*- coding: utf-8 -*-
"""Management of fault planes files and solution model
"""


import numpy as np
import json
import errno
import os


def __get_planes_json():
    """
    """
    if not os.path.isfile('segments_data.json'):
        raise FileNotFoundError(
                errno.ENOENT, os.strerror(errno.ENOENT), 'segments_data.json')
    data = json.load(open('segments_data.json'))
    rise_time = data['rise_time']
    segments = data['segments']
    return segments, rise_time


def __read_planes_info():
    """
    """

    segments, rise_time = __get_planes_json()
    n_segments = len(segments)

    with open('fault&rise_time.txt', 'r') as subfaults_file:
        lines = [line.split() for line in subfaults_file]

    strike_ps = int(lines[1][3])
    dip_ps = int(lines[1][4])

    with open('fault&rise_time.txt', 'r') as point_sources_file:
        lines2 = [line.split() for line in point_sources_file]
#
# read Fault.pos later
#
    point_sources = [[]] * n_segments
    line = 0
    for i_segment, segment in enumerate(segments):
        dip_subfaults = segment['dip_subfaults']
        stk_subfaults = segment['stk_subfaults']
        matrix = np.zeros((dip_subfaults, stk_subfaults, dip_ps, strike_ps, 7))
        line = line + 1
        for iys in range(dip_subfaults):
            for ixs in range(stk_subfaults):
                for j in range(dip_ps):
                    for i in range(strike_ps):
                        array = [float(val) for val in lines2[line]]
                        matrix[iys, ixs, j, i, :] = np.array(array)
                        line = line + 1
        point_sources[i_segment] = matrix
    return segments, rise_time, point_sources


def __unpack_plane_data(plane_info):
    """Auxiliary routine. We extract information from a fault segment.
    """
    stk_subfaults = plane_info['stk_subfaults']
    dip_subfaults = plane_info['dip_subfaults']
    delta_strike = plane_info['delta_strike']
    delta_dip = plane_info['delta_dip']
    hyp_stk = plane_info['hyp_stk'] - 1
    hyp_dip = plane_info['hyp_dip'] - 1
    return stk_subfaults, dip_subfaults, delta_strike, delta_dip, hyp_stk, hyp_dip


if __name__ == '__main__':
    import os
#    segments, rise_time, point_sources = __read_planes_info()
#    print(segments)
    os.chdir('/home/pk/Inversion_Chen_Ji/Inversiones/20140401234647/paper_challa/pk.5/NP2/srcmod_solutions/')
    tensor_info, solution = read_solution_fsp_format('hayes')
    print(solution)
