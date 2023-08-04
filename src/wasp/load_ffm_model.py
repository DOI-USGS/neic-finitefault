#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
"""


import numpy as np

import wasp.plane_management as pl_mng
from wasp import get_outputs


def load_ffm_model(
    segments_data,
    point_sources,
    option="Solucion.txt",
    max_slip=1000,
    len_stk=4,
    len_dip=4,
):
    """Load FFM model from some input file.

    :param segments_data: properties of fault segments and rise time
    :param point_sources: properties of point sources
    :param option: file from where to load the input kinematic model
    :param max_slip: largest slip of kinematic model, when creating custom
     kinematic model
    :param len_stk: length of patches when creating checkerboard model
    :param len_dip: width of patches when creating checkerboard model
    :type point_sources: list
    :type segments_data: dict
    :type option: string, optional
    :type max_slip: float, optional
    :type len_stk: int, optional
    :type len_dip: int, optional
    """
    from random import randint

    segments = segments_data["segments"]
    rise_time = segments_data["rise_time"]

    slip = []
    rake = []
    trup = []
    trise = []
    tfall = []

    if option == "Solucion.txt":
        solution = get_outputs.read_solution_static_format(segments)
        slip = solution["slip"]
        rake = solution["rake"]
        trup = solution["rupture_time"]
        trise = solution["trise"]
        tfall = solution["tfall"]

    with open("fault&rise_time.txt", "r") as input_file:
        jk = [line.split() for line in input_file]

    faults_data = [
        index + 4
        for index, (line0, line1) in enumerate(zip(jk[3:-1], jk[4:]))
        if len(line0) <= 3 and len(line1) >= 5
    ]
    headers = [
        index + 4
        for index, (line0, line1) in enumerate(zip(jk[3:-1], jk[4:]))
        if len(line0) >= 5 and len(line1) <= 3
    ]
    headers = headers[:] + [len(jk)]
    if option == "fault&rise_time.txt":
        for segment, point_source_seg, start, end in zip(
            segments, point_sources, faults_data, headers
        ):
            #
            # load FFM input model
            #
            slip_fault = np.array([float(line[0]) for line in jk[start:end]])
            rake_fault = np.array([float(line[1]) for line in jk[start:end]])
            trup_fault = np.array([float(line[2]) for line in jk[start:end]])
            trise_fault = np.array([float(line[3]) for line in jk[start:end]])
            tfall_fault = np.array([float(line[4]) for line in jk[start:end]])
            (
                n_sub_stk,
                n_sub_dip,
                delta_x,
                delta_y,
                hyp_stk,
                hyp_dip,
            ) = pl_mng.__unpack_plane_data(segment)
            #
            # Reshape the rupture process
            #
            slip_fault.shape = n_sub_dip, n_sub_stk
            rake_fault.shape = n_sub_dip, n_sub_stk
            trup_fault.shape = n_sub_dip, n_sub_stk
            trise_fault.shape = n_sub_dip, n_sub_stk
            tfall_fault.shape = n_sub_dip, n_sub_stk
            slip = slip + [slip_fault]
            rake = rake + [rake_fault]
            trup = trup + [trup_fault]
            trise = trise + [trise_fault]
            tfall = tfall + [tfall_fault]

    if option == "point_source":
        min_rise = rise_time["min_rise"]
        for i_segment, (segment, point_source_seg, start, end) in enumerate(
            zip(segments, point_sources, faults_data, headers)
        ):
            rake_value = segment["rake"]
            a, b, c, d, e = point_sources[0].shape
            ny = int(c / 2)
            nx = int(d / 2)
            #
            # load FFM input model
            #
            slip_fault = np.array([float(line[0]) for line in jk[start:end]])
            rake_fault = np.array([rake_value for line in jk[start:end]])
            trup_fault = point_source_seg[:, :, ny, nx, 4]
            trise_fault = np.array([min_rise for line in jk[start:end]])
            tfall_fault = np.array([min_rise for line in jk[start:end]])
            (
                n_sub_stk,
                n_sub_dip,
                delta_x,
                delta_y,
                hyp_stk,
                hyp_dip,
            ) = pl_mng.__unpack_plane_data(segment)
            #
            # Reshape the rupture process
            #
            slip_fault.shape = n_sub_dip, n_sub_stk
            for ny in range(n_sub_dip):
                for nx in range(n_sub_stk):
                    if nx == len_stk and ny == len_dip and i_segment == 0:
                        slip_fault[ny, nx] = max_slip
                    else:
                        slip_fault[ny, nx] = 0
            rake_fault.shape = n_sub_dip, n_sub_stk
            trup_fault.shape = n_sub_dip, n_sub_stk
            trise_fault.shape = n_sub_dip, n_sub_stk
            tfall_fault.shape = n_sub_dip, n_sub_stk
            slip = slip + [slip_fault]
            rake = rake + [rake_fault]
            trup = trup + [trup_fault]
            trise = trise + [trise_fault]
            tfall = tfall + [tfall_fault]

    if option == "Checkerboard":
        min_rise = rise_time["min_rise"]
        i = 0
        for segment, point_source_seg, start, end in zip(
            segments, point_sources, faults_data, headers
        ):
            rake_value = segment["rake"]
            a, b, c, d, e = point_sources[0].shape
            ny = int(c / 2)
            nx = int(d / 2)
            #
            # load FFM input model
            #
            slip_fault = np.array([float(line[0]) for line in jk[start:end]])
            rake_fault = np.array([rake_value for line in jk[start:end]])
            trup_fault = point_source_seg[:, :, ny, nx, 4]
            trise_fault = np.array([min_rise for line in jk[start:end]])
            tfall_fault = np.array([min_rise for line in jk[start:end]])
            (
                n_sub_stk,
                n_sub_dip,
                delta_x,
                delta_y,
                hyp_stk,
                hyp_dip,
            ) = pl_mng.__unpack_plane_data(segment)
            #
            # Reshape the rupture process
            #
            slip_fault.shape = n_sub_dip, n_sub_stk
            len_stk = len_stk if len_stk < 0.5 * n_sub_stk else int(0.45 * n_sub_stk)
            len_dip = len_dip if len_dip < 0.5 * n_sub_dip else int(0.45 * n_sub_dip)
            for ny in range(n_sub_dip):
                for nx in range(n_sub_stk):
                    if (int(nx // len_stk) + int(ny // len_dip) + i) % 2 == 0:
                        slip_fault[ny, nx] = 0
                    else:
                        slip_fault[ny, nx] = max_slip
            rake_fault.shape = n_sub_dip, n_sub_stk
            trup_fault.shape = n_sub_dip, n_sub_stk
            trise_fault.shape = n_sub_dip, n_sub_stk
            tfall_fault.shape = n_sub_dip, n_sub_stk
            slip = slip + [slip_fault]
            rake = rake + [rake_fault]
            trup = trup + [trup_fault]
            trise = trise + [trise_fault]
            tfall = tfall + [tfall_fault]
            i = i + 1

    if option == "Patches":
        min_rise = rise_time["min_rise"]
        for i_segment, (segment, point_source_seg, start, end) in enumerate(
            zip(segments, point_sources, faults_data, headers)
        ):
            rake_value = segment["rake"]
            a, b, c, d, e = point_sources[0].shape
            ny = int(c / 2)
            nx = int(d / 2)
            #
            # load FFM input model
            #
            (
                n_sub_stk,
                n_sub_dip,
                delta_x,
                delta_y,
                hyp_stk,
                hyp_dip,
            ) = pl_mng.__unpack_plane_data(segment)
            slip_fault = np.array([0 for line in jk[start:end]])
            rake_fault = rake_value + 2 * np.random.randn(end - start)  # 2
            trup_fault = point_source_seg[:, :, ny, nx, 4] * (
                np.ones((n_sub_dip, n_sub_stk))
                + 0.2 * np.random.randn(n_sub_dip, n_sub_stk)
            )
            trup_fault = np.maximum(trup_fault, 0)
            trise_fault = np.array([2 * min_rise for line in jk[start:end]])
            tfall_fault = np.array([2 * min_rise for line in jk[start:end]])

            #
            # Reshape the rupture process
            #
            slip_fault.shape = n_sub_dip, n_sub_stk
            patches = 1  # randint(3, 5)

            for patch in range(patches):
                n_sub_xc = randint(0, n_sub_stk - 1)
                n_sub_yc = 0  # randint(0, 2)#n_sub_dip - 1)
                #                if patch == 0:
                #                    n_sub_xc = hyp_stk - 1#randint(0, 2)
                #                    n_sub_yc = hyp_dip - 1#randint(4, 6)
                #                if patch == 1:
                #                    n_sub_xc = hyp_stk + 1#randint(n_sub_stk - 2, n_sub_stk)
                #                    n_sub_yc = hyp_dip + 1#randint(4, 6)
                length = randint(3, 4)  # length = 2
                for ny in range(n_sub_dip):
                    for nx in range(n_sub_stk):
                        dist = (n_sub_yc - ny) ** 2 + (n_sub_xc - nx) ** 2
                        if dist > length**2:
                            continue
                        slip_fault[ny, nx] = (
                            slip_fault[ny, nx]
                            + max_slip * (length - np.sqrt(dist)) / length
                        )
            rake_fault.shape = n_sub_dip, n_sub_stk
            trup_fault.shape = n_sub_dip, n_sub_stk
            trise_fault.shape = n_sub_dip, n_sub_stk
            tfall_fault.shape = n_sub_dip, n_sub_stk
            slip = slip + [slip_fault]
            rake = rake + [rake_fault]
            trup = trup + [trup_fault]
            trise = trise + [trise_fault]
            tfall = tfall + [tfall_fault]

    model = {"slip": slip, "rake": rake, "trup": trup, "trise": trise, "tfall": tfall}
    return model
