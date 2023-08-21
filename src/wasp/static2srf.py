# -*- coding: utf-8 -*-
"""Translate a solution file in static format, to a solution file in FSP format
"""


import datetime
import json
import os
import pathlib
from typing import List, Tuple, Union

import numpy as np
import pandas as pd  # type:ignore
from obspy.geodetics import flinnengdahl  # type:ignore
from pyproj import Geod  # type:ignore

import wasp.fault_plane as pf
import wasp.plane_management as pl_mng
import wasp.seismic_tensor as tensor
from wasp import get_outputs


def static_to_srf(
    tensor_info: dict,
    segments_data: dict,
    used_data: List[str],
    vel_model: dict,
    solution: dict,
    directory: Union[pathlib.Path, str] = pathlib.Path(),
):
    """Write SRF file with the solution of FFM modelling from file Solucion.txt

    :param tensor_info: The moment tensor information
    :type tensor_info: dict
    :param segments_data: The fault segments data
    :type segments_data: dict
    :param used_data: The data types used in modelling
    :type used_data: List[str]
    :param vel_model: The velocity model parameters
    :type vel_model: dict
    :param solution: The solution
    :type solution: dict
    :param directory: Where to read/write files, defaults to pathlib.Path()
    :type directory: Union[pathlib.Path, str], optional
    """
    directory = pathlib.Path(directory)
    print("Writing SRF file output")
    locator = flinnengdahl.FlinnEngdahl()
    segments = segments_data["segments"]
    rise_time = segments_data["rise_time"]
    connections = None
    if "connections" in segments_data:
        connections = segments_data["connections"]
    point_sources = pf.point_sources_param(
        segments, tensor_info, rise_time, connections=connections
    )
    delta_strike = segments[0]["delta_strike"]
    delta_dip = segments[0]["delta_dip"]
    rupture_vel = segments[0]["rupture_vel"]
    subfaults = {"delta_strike": delta_strike, "delta_dip": delta_dip}
    subfaults2 = pf._point_sources_def(rise_time, rupture_vel, subfaults)
    strike_ps = int(subfaults2["strike_ps"] / 2)
    dip_ps = int(subfaults2["dip_ps"] / 2)
    slips = solution["slip"]
    rakes = solution["rake"]
    trup = solution["rupture_time"]
    trise = solution["trise"]
    tfall = solution["tfall"]
    # depths = solution['depth']
    moment = solution["moment"]
    total_moment = 0
    for moment_segment in moment:
        total_moment = total_moment + np.sum(np.array(moment_segment).flatten())

    string = " ---------------------------------- "
    event_lat = tensor_info["lat"]
    event_lon = tensor_info["lon"]
    depth = tensor_info["depth"]
    moment_mag = 2 * np.log10(total_moment) / 3 - 10.7
    location = locator.get_region(event_lon, event_lat)
    date = tensor_info["datetime"]
    tag = date
    now = datetime.datetime.now()

    latitudes = [ps_segment[:, :, dip_ps, strike_ps, 0] for ps_segment in point_sources]
    longitudes = [
        ps_segment[:, :, dip_ps, strike_ps, 1] for ps_segment in point_sources
    ]
    depths = [ps_segment[:, :, dip_ps, strike_ps, 2] for ps_segment in point_sources]
    windows = rise_time["windows"]
    total_segments = len(segments)

    ##########################################################
    ### COUNT UP HOW MANY OF EACH OBSERVATION TYPE WE USED ###
    ##########################################################
    quantity_strong = 0
    strong_phimx = 0
    strong_r = 0
    if "strong_motion" in used_data:
        strong_data = pd.read_json(directory / "strong_motion_waves.json")
        quantity_strong = int(
            len(strong_data) / 3
        )  # divide by 3 for number of stations rather than number of channels
        strong_r = min(strong_data["distance"])
        strong_az = np.sort(strong_data["azimuth"].values.tolist())
        strong_az = np.append(strong_az, strong_az[0] + 360)
        strong_phimx = max(np.diff(strong_az))
    quantity_cgps = 0
    cgps_phimx = 0
    cgps_r = 0
    if "cgps" in used_data:
        cgps_data = pd.read_json(directory / "cgps_waves.json")
        quantity_cgps = int(len(cgps_data) / 3)
        cgps_r = min(cgps_data["distance"])
        cgps_az = np.sort(cgps_data["azimuth"].values.tolist())
        cgps_az = np.append(cgps_az, cgps_az[0] + 360)
        cgps_phimx = max(np.diff(cgps_az))
    quantity_gps = 0
    gps_phimx = 0
    gps_r = 0
    if "gps" in used_data:
        gps_data = pd.read_json(directory / "static_data.json")
        quantity_gps = len(gps_data)
        gps_r = min(gps_data["distance"])
        gps_az = np.sort(gps_data["azimuth"].values.tolist())
        gps_az = np.append(gps_az, gps_az[0] + 360)
        gps_phimx = max(np.diff(gps_az))
    quantity_tele = 0
    tele_phimx = 0
    tele_r = 0
    if "tele_body" in used_data:
        tele_data = pd.read_json(directory / "tele_waves.json")
        quantity_tele = len(tele_data)
        tele_r = min(tele_data["distance"])
        tele_az = np.sort(tele_data["azimuth"].values.tolist())
        tele_az = np.append(tele_az, tele_az[0] + 360)
        tele_phimx = max(np.diff(tele_az))
    quantity_surf = 0
    surf_phimx = 0
    surf_r = 0
    if "surf_tele" in used_data:
        surf_data = pd.read_json(directory / "surf_waves.json")
        quantity_surf = len(surf_data)
        surf_r = min(surf_data["distance"])
        surf_az = np.sort(surf_data["azimuth"].values.tolist())
        surf_az = np.append(surf_az, surf_az[0] + 360)
        surf_phimx = max(np.diff(surf_az))
    quantity_insar_scenes = 0
    quantity_insar_points = 0
    insar_r = 0
    if "insar" in used_data:
        insar_data = pd.read_json(directory / "insar_data.json")
        if "ascending" in insar_data:
            asc_insar = len(insar_data["ascending"])
            quantity_insar_scenes += asc_insar
        if "descending" in insar_data:
            desc_insar = len(insar_data["descending"])
            quantity_insar_scenes += desc_insar
        with open(directory / "insar_data.txt") as f:
            first_line = f.readline().strip("\n")
            quantity_insar_points = int(first_line)
    quantity_dart = 0
    dart_phimx = 0
    dart_r = 0
    if "dart" in used_data:
        dart_data = pd.read_json(directory / "dart_data.json")
        quantity_dart = len(dart_data)
        dart_r = min(dart_data["distance"])
        dart_az = np.sort(dart_data["azimuth"].values.tolist())
        dart_az = np.append(dart_az, dart_az[0] + 360)
        dart_phimx = max(np.diff(dart_az))

    n_layers = len(vel_model["dens"])
    p_vel = [float(v) for v in vel_model["p_vel"]]
    s_vel = [float(v) for v in vel_model["s_vel"]]
    dens = [float(v) for v in vel_model["dens"]]
    thick = [float(v) for v in vel_model["thick"]]
    qp = [float(v) for v in vel_model["qa"]]
    qs = [float(v) for v in vel_model["qb"]]
    depth2 = np.cumsum(np.array(thick))
    depth2 = np.concatenate([[0], depth2])  # - depth2[0]
    zipped = zip(depth2, p_vel, s_vel, dens, qp, qs)
    zipped2 = zip(
        segments,
        point_sources,
        latitudes,
        longitudes,
        depths,
        slips,
        rakes,
        trup,
        trise,
        tfall,
        moment,
    )

    with open(directory / "srf_sol_file.txt", "w") as outfile:
        outfile.write("2.0\n")

        outfile.write(
            "#\n# Data :\tBODY\tSURF\tSTRONG\tcGPS\tGPS\tInSAR\tDART\tTRIL\tLEVEL\tOTHER\n"
        )
        outfile.write(
            "# NoS  :\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t0\t0\t0\n".format(
                quantity_tele,
                quantity_surf,
                quantity_strong,
                quantity_cgps,
                quantity_gps,
                quantity_insar_points,
                quantity_dart,
            )
        )
        outfile.write(
            "# PHImx :\t{:.2f}\t{:.2f}\t{:.2f}\t{:.2f}\t{:.2f}\t{:d}\t{:.2f}\t0.0\t0.0\t0.0\n".format(
                tele_phimx,
                surf_phimx,
                strong_phimx,
                cgps_phimx,
                gps_phimx,
                quantity_insar_scenes,
                dart_phimx,
            )
        )
        outfile.write(
            "# Rmin :\t{:.2f}\t{:.2f}\t{:.2f}\t{:.2f}\t{:.2f}\t--\t{:.2f}\t0.0\t0.0\t0.0\n".format(
                tele_r, surf_r, strong_r, cgps_r, gps_r, insar_r, dart_r  # type:ignore
            )
        )

        outfile.write("#{}FINITE-SOURCE RUPTURE " "MODEL{}\n#\n".format(string, string))
        outfile.write("# Event : {} {} NEIC\n".format(location, date))
        outfile.write("# EventTAG: {}\n#\n".format(tag))
        outfile.write(
            "# Loc  : LAT = {}  LON = {}  DEP = {}\n".format(
                event_lat, event_lon, depth
            )
        )
        outfile.write("#\n#{}{}\n".format(string, string))
        outfile.write("#\n# VELOCITY-DENSITY STRUCTURE\n")
        outfile.write("# No. of layers = {}\n".format(n_layers))
        outfile.write("#\n# DEPTH P_VEL S_VEL DENS QP QS\n")
        outfile.write("# [km] [km/s] [km/s] [g/cm^3]\n")
        for dep, pv, sv, den, qpp, qss in zipped:
            outfile.write(
                "# {:.2f} {:.2f} {:.2f} {:.2f}  {:.1f}  "
                "{:.1f}\n".format(dep, pv, sv, den, qpp, qss)
            )
        outfile.write("#\n#{}{}\n".format(string, string))
        outfile.write(
            "# {}/{}/{} created by degoldberg@usgs.gov"
            "\n".format(now.day, now.month, now.year)
        )
        outfile.write("#\n# SOURCE MODEL PARAMETERS\n")
        outfile.write("PLANE {}\n".format(total_segments))

        ############################################
        ### WRITE PLANE HEADERS FOR EACH SEGMENT ###
        ############################################
        for i_segment, fault_segment_data in enumerate(zipped2):
            segment = fault_segment_data[0]
            ps_seg = fault_segment_data[1]
            lat_fault = fault_segment_data[2].flatten()
            lon_fault = fault_segment_data[3].flatten()
            depth_fault = fault_segment_data[4].flatten()
            strike = segment["strike"]
            dip = segment["dip"]
            plane_info = segments[i_segment]
            (
                stk_subfaults,
                dip_subfaults,
                delta_strike,
                delta_dip,
                hyp_stk,
                hyp_dip,
            ) = pl_mng.__unpack_plane_data(plane_info)
            length = stk_subfaults * delta_strike
            width = dip_subfaults * delta_dip
            min_dep = np.min(ps_seg[:, :, :, :, 2].flatten())
            elon = np.mean(fault_segment_data[3][0])
            elat = np.mean(fault_segment_data[2][0])
            depth_top = min(depth_fault) - np.sin(np.radians(dip)) * (delta_dip / 2)
            dip_hyp_pos = (depth - depth_top) / np.sin(np.radians(dip))

            # Get hypocenter row of subfaults and find coordinates of middle
            i_row = np.argmin(abs(depth_fault - depth))
            test_depth = depth_fault[i_row]
            hypo_row = np.where(depth_fault == test_depth)[0]
            hypo_center_lon = lon_fault[hypo_row].mean()
            hypo_center_lat = lat_fault[hypo_row].mean()
            g = Geod(ellps="WGS84")
            az, baz, dist_from_center = g.inv(
                event_lon, event_lat, hypo_center_lon, hypo_center_lat
            )
            # is dist_from_center postitive or negative?
            if strike > 180:
                strike_rectified = strike - 360
            else:
                strike_rectified = strike
            if np.sign(az) == np.sign(strike_rectified):
                stk_hyp_pos = -dist_from_center / 1000
            else:
                stk_hyp_pos = dist_from_center / 1000
            outfile.write(
                "{:.4f} \t {:.4f} \t {:d} \t {:d} \t {:.2f} \t {:.2f} \n".format(
                    elon,
                    elat,
                    stk_subfaults,
                    dip_subfaults,
                    stk_subfaults * delta_strike,
                    dip_subfaults * delta_dip,
                )
            )
            outfile.write(
                "{:.2f} \t {:.2f} \t {:.2f} \t {:.2f} \t {:.2f} \n".format(
                    strike, dip, depth_top, stk_hyp_pos, dip_hyp_pos
                )
            )

        #########################################
        ### WRITE OUT POINTS FOR EACH SEGMENT ###
        #########################################
        zipped2 = zip(
            segments,
            point_sources,
            latitudes,
            longitudes,
            depths,
            slips,
            rakes,
            trup,
            trise,
            tfall,
            moment,
        )
        stf_dt = 1  # seconds
        vs = -1  # if not known
        density = -1  # if not known
        for i_segment, fault_segment_data in enumerate(zipped2):
            plane_info = segments[i_segment]
            (
                stk_subfaults,
                dip_subfaults,
                delta_strike,
                delta_dip,
                hyp_stk,
                hyp_dip,
            ) = pl_mng.__unpack_plane_data(plane_info)
            subfault_area = delta_strike * delta_dip
            segment = fault_segment_data[0]
            ps_seg = fault_segment_data[1]
            lat_fault = fault_segment_data[2].flatten()
            lon_fault = fault_segment_data[3].flatten()
            depth_fault = fault_segment_data[4].flatten()
            slip_fault = fault_segment_data[5].flatten()
            rake_fault = fault_segment_data[6].flatten()
            trup_fault = fault_segment_data[7].flatten()
            trise_fault = fault_segment_data[8].flatten()
            tfall_fault = fault_segment_data[9].flatten()
            moment_fault = fault_segment_data[10].flatten()

            # calculate total time each subfault can slip:
            max_trise = max(trise_fault)
            max_tfall = max(tfall_fault)
            total_time = max_trise + max_tfall

            outfile.write("POINTS {}\n".format(len(lat_fault)))
            zipped3 = zip(
                lat_fault,
                lon_fault,
                depth_fault,
                slip_fault,
                rake_fault,
                trup_fault,
                trise_fault,
                tfall_fault,
                moment_fault,
            )
            for line in zipped3:
                lat, lon, dep, slip, rake, t_rup, t_ris, t_fal, moment = line
                outfile.write(
                    "{:.4f} {:.4f} {:.2f} {:.2f} {:.2f} {:.2f} {:.2f} {:.1f} {:d} {:d}\n".format(
                        lon,
                        lat,
                        dep,
                        strike,
                        dip,
                        subfault_area,
                        t_rup,
                        stf_dt,
                        vs,
                        density,
                    )
                )
                if slip < 10e-4:
                    print(
                        "Zero slip at: {:.2f}, {:.2f}, {:.2f}km".format(lon, lat, dep)
                    )
                    Nstf = 0  # If zero slip, no slip function
                else:  # Calculate STF
                    tstf, stf = build_source_time_function(
                        t_ris, t_fal, stf_dt, total_time, scale=True
                    )
                    stf_adjust_factor = slip / stf_dt
                    stf = stf * stf_adjust_factor  # now tf is in cm/sec
                    Nstf = len(stf)
                outfile.write("{:.2f} {:.4f} {:d} 0 0 0 0 \n".format(rake, slip, Nstf))
                if slip > 10e-4:
                    # Write stf 6 values per line
                    for kstf in range(Nstf):
                        if kstf == 0:
                            white_space = "  "
                        elif (kstf + 1) % 6 == 0:
                            white_space = "\n"
                        elif (kstf + 1) == Nstf:
                            white_space = "\n"
                        else:
                            white_space = "  "

                        if kstf == 0:
                            pre_white_space = "  "
                        elif (kstf) % 6 == 0:
                            pre_white_space = "  "
                        else:
                            pre_white_space = ""
                        outfile.write(
                            "%s%.6e%s" % (pre_white_space, stf[kstf], white_space)
                        )


def build_source_time_function(
    t_ris: np.ndarray,
    t_fal: np.ndarray,
    dt: float,
    total_time: int,
    scale: bool = True,
    scale_value: float = 1.0,
) -> Tuple[np.ndarray, np.ndarray]:
    """Compute source time function for a given rise time

    :param t_ris: The rise time
    :type t_ris: np.ndarray
    :param t_fal: The fall time
    :type t_fal: np.ndarray
    :param dt: _description_
    :type dt: float
    :param total_time: The total time
    :type total_time: int
    :param scale: Whether to scale, defaults to True
    :type scale: bool, optional
    :param scale_value: The scale value, defaults to 1.0
    :type scale_value: float, optional
    :return: Time and Mdot arrays
    :rtype: Tuple[np.ndarray, np.ndarray]
    """
    from numpy import arange, cos, exp, isnan, pi, roll, sin, where, zeros
    from scipy.integrate import trapz  # type:ignore

    # Initialize outputs
    t = arange(0, total_time + dt, dt)
    Mdot = zeros(t.shape)

    # Up going cosine
    s1 = (1.0 / (t_ris + t_fal)) * (1 - cos((pi * t) / t_ris))
    i = where(t > t_ris)[0]
    s1[i] = 0

    # Down going cosine
    s2 = (1.0 / (t_ris + t_fal)) * (1 + cos((pi * (t - t_ris)) / t_fal))
    i = where(t <= t_ris)[0]
    s2[i] = 0
    i = where(t > t_ris + t_fal)[0]
    s2[i] = 0

    # add the two
    Mdot = s1 + s2

    # Area of STF must be equal to dt
    if scale == True:
        if scale_value is None:  # Re-scale to dt (for traditional convolution)
            target = dt  # this is the target scale value
        else:
            target = scale_value

        area = trapz(Mdot, t)
        Mdot = Mdot * (scale_value / area)
    # Check for errors
    if isnan(Mdot[0]) == True:
        print("ERROR: whoops, STF has nan values!")
        return  # type:ignore

    return t, Mdot


if __name__ == "__main__":
    import argparse
    import errno

    import wasp.manage_parser as mp

    parser = argparse.ArgumentParser()
    parser.add_argument(
        "-f", "--folder", default=os.getcwd(), help="folder where there are input files"
    )
    parser = mp.parser_add_tensor(parser)
    parser = mp.parser_ffm_data(parser)
    args = parser.parse_args()
    os.chdir(args.folder)
    used_data = mp.get_used_data(args)
    if args.gcmt_tensor:
        cmt_file = args.gcmt_tensor
        tensor_info = tensor.get_tensor(cmt_file=cmt_file)
    else:
        tensor_info = tensor.get_tensor()
    segments_data = json.load(open("segments_data.json"))
    segments = segments_data["segments"]
    if not os.path.isfile("velmodel_data.json"):
        raise FileNotFoundError(
            errno.ENOENT, os.strerror(errno.ENOENT), "velmodel_data.json"
        )
    vel_model = json.load(open("velmodel_data.json"))
    solution = get_outputs.read_solution_static_format(segments)
    static_to_srf(tensor_info, segments_data, used_data, vel_model, solution)
