# -*- coding: utf-8 -*-
"""
The routines here create the remaining files in the Downloads section of the NEIC eventpages
"""

import json
import os
import pathlib
import shutil
import warnings
from collections import Counter
from glob import glob
from typing import List, Optional, Union

import cutde.halfspace as half_space  # type:ignore
import matplotlib.pyplot as plt  # type: ignore
import pandas as pd  # type:ignore
import pygmt  # type:ignore
import pyproj  # type: ignore
from numpy import (
    array,
    c_,
    ceil,
    cos,
    deg2rad,
    floor,
    int64,
    linspace,
    meshgrid,
    ones,
    shape,
    sin,
    size,
    zeros,
)
from obspy.imaging.scripts import mopad  # type: ignore

import wasp.management as mng


def temporary_file_reorganization_for_publishing(
    evID: str, directory: Union[pathlib.Path, str] = pathlib.Path()
):
    """Reorganize files for publishing to event pages

    :param evID: The event id
    :type evID: str
    :param directory: The directory where files exist, defaults to pathlib.Path
    :type directory: Union[pathlib.Path, str], optional
    """
    print("Reassigning filenames for use with sendproduct")
    directory = pathlib.Path(directory)
    pub_directory = directory / "PublicationFiles"
    if not os.path.exists(pub_directory):
        os.makedirs(pub_directory)
    ##################
    ### PARAM FILE ###
    ##################
    orig_param = directory / "Solucion.txt"
    shutil.copy(orig_param, pub_directory)
    mv_param = pub_directory / "Solucion.txt"
    pub_param = pub_directory / f"{evID}.param"
    os.rename(mv_param, pub_param)

    ################
    ### FSP FILE ###
    ################
    orig_fsp = directory / "fsp_sol_file.txt"
    shutil.copy(orig_fsp, pub_directory)
    mv_fsp = pub_directory / "fsp_sol_file.txt"
    pub_fsp = pub_directory / f"{evID}.fsp"
    os.rename(mv_fsp, pub_fsp)
    ########################
    ### SHAKEMAP POLYGON ###
    ########################
    shutil.copy(directory / "shakemap_polygon.txt", pub_directory)
    ###########################
    ### SURFACE DEFORMATION ###
    ###########################
    shutil.copy(directory / "surface_deformation.disp", pub_directory)
    #####################
    ### COULOMB INPUT ###
    #####################
    orig_coulomb = directory / "Coulomb.inp"
    shutil.copy(orig_coulomb, pub_directory)
    mv_coulomb = pub_directory / "Coulomb.inp"
    pub_coulomb = pub_directory / f"{evID}_coulomb.inp"
    os.rename(mv_coulomb, pub_coulomb)
    ####################
    ### CMT SOLUTION ###
    ####################
    shutil.copy(directory / "CMTSOLUTION", pub_directory)
    #######################
    ### WAVE PROPERTIES ###
    #######################
    shutil.copy(directory / "wave_properties.json", pub_directory)
    ###################
    ### MOMENT RATE ###
    ###################
    orig_mr = directory / "STF.txt"
    shutil.copy(orig_mr, pub_directory)
    mv_mr = pub_directory / "STF.txt"
    pub_mr = pub_directory / f"{evID}.mr"
    os.rename(mv_mr, pub_mr)
    orig_mr_png = directory / "plots" / "MomentRate.png"
    shutil.copy(orig_mr_png, pub_directory)
    mv_mr_png = pub_directory / "MomentRate.png"
    pub_mr_png = pub_directory / "mr.png"
    os.rename(mv_mr_png, pub_mr_png)
    ###########
    ### MAP ###
    ###########
    orig_map = directory / "plots" / "Map.png"
    shutil.copy(orig_map, pub_directory)
    mv_map = pub_directory / "Map.png"
    pub_map = pub_directory / f"{evID}_basemap.png"
    os.rename(mv_map, pub_map)
    ################################
    ### SLIP DISTRIBUTION FIGURE ###
    ################################
    orig_slip = directory / "plots" / "SlipDist_plane0.png"
    shutil.copy(orig_slip, pub_directory)
    mv_slip = pub_directory / "SlipDist_plane0.png"
    pub_slip = pub_directory / f"{evID}_slip2.png"
    os.rename(mv_slip, pub_slip)
    ############################
    ### DATA/SYNTHETIC PLOTS ###
    ############################
    wave_plots = glob(str(directory / "plots") + "/*_waves.png")
    wave_plots_pub_directory = pub_directory / "fits"
    if not os.path.exists(wave_plots_pub_directory):
        os.makedirs(wave_plots_pub_directory)
    for wp in range(len(wave_plots)):
        orig_wp = wave_plots[wp]
        shutil.copy(orig_wp, wave_plots_pub_directory)
    insar_plots = glob(str(directory / "plots") + "/InSAR_*_fit.png")
    for ip in range(len(insar_plots)):
        orig_ip = insar_plots[ip]
        shutil.copy(orig_ip, wave_plots_pub_directory)

    shutil.make_archive(
        str(pub_directory / "fits"),
        "zip",
        wave_plots_pub_directory,
    )
    #######################
    ### RESAMPLED INSAR ###
    #######################
    resampled_insar_data_asc = glob(str(directory) + "/*ascending.txt")
    resampled_insar_data_desc = glob(str(directory) + "/*descending.txt")
    interferograms_pub_directory = pub_directory / "resampled_interferograms"

    if len(resampled_insar_data_asc) > 0 or len(resampled_insar_data_desc) > 0:
        if not os.path.exists(interferograms_pub_directory):
            os.makedirs(interferograms_pub_directory)
    for asc in range(len(resampled_insar_data_asc)):
        orig_asc = resampled_insar_data_asc[asc]
        shutil.copy(orig_asc, interferograms_pub_directory)
    for desc in range(len(resampled_insar_data_desc)):
        orig_desc = resampled_insar_data_desc[desc]
        shutil.copy(orig_desc, interferograms_pub_directory)
    if os.path.exists(interferograms_pub_directory):
        shutil.make_archive(
            str(pub_directory / "resampled_interferograms"),
            "zip",
            interferograms_pub_directory,
        )


def make_waveproperties_json(directory: Union[pathlib.Path, str] = pathlib.Path()):
    """Write wave_properties.json

    :param directory: The directory where to write the file(s), defaults to pathlib.Path()
    :type directory: Union[pathlib.Path, str], optional
    """
    print("Counting Observations for Waveform Properties JSON")
    directory = pathlib.Path(directory)
    #######################################################
    ### READ IN NUMBER OF EACH OBS TYPE FROM WAVE JSONS ###
    #######################################################
    with open(directory / "tele_waves.json") as f:
        bodywave_channels = json.load(f)
    bodywave_components = [
        bodywave_channel["component"] for bodywave_channel in bodywave_channels
    ]
    bodywave_count = Counter(bodywave_components)
    num_pwaves = int(bodywave_count["BHZ"])
    num_shwaves = int(bodywave_count["BHT"])
    if os.path.exists(directory / "surf_waves.json"):
        with open(directory / "surf_waves.json") as f:
            surfwave_channels = json.load(f)
        num_longwaves = int(len(surfwave_channels))
    else:
        num_longwaves = 0
    #############################################
    ### WRITE RESULTS TO WAVE_PROPERTIES.JSON ###
    #############################################
    with open(directory / "wave_properties.json", "w") as f:
        json.dump(
            {
                "num_longwaves": str(num_longwaves),
                "num_pwaves": str(num_pwaves),
                "num_shwaves": str(num_shwaves),
            },
            f,
        )


def write_CMTSOLUTION_file(
    pdefile: Optional[Union[pathlib.Path, str]] = None,
    directory: Union[pathlib.Path, str] = pathlib.Path(),
):
    """Write the CMTSOLUTION file

    :param pdefile: The path to the cmt file, defaults to None
    :type pdefile: Optional[Union[pathlib.Path, str]], optional
    :param directory: The directory where to write the file(s), defaults to pathlib.Path()
    :type directory: Union[pathlib.Path, str], optional
    """
    directory = pathlib.Path(directory)
    print("Writing CMTSOLUTION to file...")

    if pdefile is None:
        pdefile = pathlib.Path(glob("../../../info/*_cmt_CMT")[0])
    else:
        pdefile = pathlib.Path(pdefile)
    with open(pdefile, "r") as pdef:
        lines = pdef.readlines()
        pde = lines[0]
        eve = lines[1]

    LAT: List[float] = []
    LON: List[float] = []
    DEP: List[float] = []
    SLIP: List[float] = []
    RAKE: List[float] = []
    STRIKE: List[float] = []
    DIP: List[float] = []
    T_RUP: List[float] = []
    T_RIS: List[float] = []
    T_FAL: List[float] = []
    MO: List[float] = []

    with open(directory / "Solucion.txt") as sol:
        for line in sol:
            if "#" in line:  # HEADER LINES
                continue
            if len(array(line.split())) < 4:  # FAULT BOUNDARY LINES
                continue
            else:  # ACTUAL SUBFAULT DETAILS
                (
                    lat,
                    lon,
                    dep,
                    slip,
                    rake,
                    strike,
                    dip,
                    t_rup,
                    t_ris,
                    t_fal,
                    mo,
                ) = [float(v) for v in line.split()]
                # Make rake be -180 to 180 (not 0-360)
                if float(rake) > 180:
                    rake = float(rake) - 360
                LAT.append(float(lat))
                LON.append(float(lon))
                DEP.append(float(dep))
                SLIP.append(float(slip))
                RAKE.append(float(rake))
                STRIKE.append(float(strike))
                DIP.append(float(dip))
                T_RUP.append(float(t_rup))
                T_RIS.append(float(t_ris))
                T_FAL.append(float(t_fal))
                MO.append(float(mo))

    SUBFAULTS = c_[LAT, LON, DEP, SLIP, RAKE, STRIKE, DIP, T_RUP, T_RIS, T_FAL, MO]
    SUBFAULTS = SUBFAULTS[SUBFAULTS[:, 7].argsort()]
    LAT = SUBFAULTS[:, 0]
    LON = SUBFAULTS[:, 1]
    DEP = SUBFAULTS[:, 2]
    SLIP = SUBFAULTS[:, 3]
    RAKE = SUBFAULTS[:, 4]
    STRIKE = SUBFAULTS[:, 5]
    DIP = SUBFAULTS[:, 6]
    T_RUP = SUBFAULTS[:, 7]
    T_RIS = SUBFAULTS[:, 8]
    T_FAL = SUBFAULTS[:, 9]
    MO = SUBFAULTS[:, 10]

    with open(directory / "CMTSOLUTION", "w") as CMTout:
        for ksubf in range(len(SUBFAULTS)):
            MT = mopad.MomentTensor([STRIKE[ksubf], DIP[ksubf], RAKE[ksubf]])
            MT = array(MT.get_M())
            MT = MT * MO[ksubf]

            Mrr = MT[0, 0]
            Mtt = MT[1, 1]
            Mpp = MT[2, 2]
            Mrt = MT[0, 1]
            Mrp = MT[0, 2]
            Mtp = MT[1, 2]

            CMTout.write(pde)
            CMTout.write(eve)
            CMTout.write("time shift:{:>13} \n".format("{:.4f}".format(T_RUP[ksubf])))
            CMTout.write(
                "half duration:{:>10} \n".format("{:.4f}".format(T_RIS[ksubf]))
            )
            CMTout.write("latitude:{:>15} \n".format("{:.4f}".format(LAT[ksubf])))
            CMTout.write("longitude:{:>14} \n".format("{:.4f}".format(LON[ksubf])))
            CMTout.write("depth:{:>18} \n".format("{:.4f}".format(DEP[ksubf])))
            CMTout.write("Mrr: {:>19} \n".format("{:.6e}".format(Mrr)))
            CMTout.write("Mtt: {:>19} \n".format("{:.6e}".format(Mtt)))
            CMTout.write("Mpp: {:>19} \n".format("{:.6e}".format(Mpp)))
            CMTout.write("Mrt: {:>19} \n".format("{:.6e}".format(Mrt)))
            CMTout.write("Mrp: {:>19} \n".format("{:.6e}".format(Mrp)))
            CMTout.write("Mtp: {:>19} \n".format("{:.6e}".format(Mtp)))


def write_Coulomb_file(
    directory: Union[pathlib.Path, str] = pathlib.Path(), eventID: Optional[str] = None
):
    """Write the Coulomb.inp file

    :param directory: The directory where to write the file(s), defaults to pathlib.Path()
    :type directory: Union[pathlib.Path, str], optional
    :param eventID: The event id, defaults to None
    :type eventID: Optional[str], optional
    """
    directory = pathlib.Path(directory)
    print("Writing Coulomb.inp to file...")

    if eventID is None:
        eventID = "EventX"

    ##########################################################################
    ### GRAB RELEVANT SUBFAULT PARAMETER INFO FROM EVENT_MULT.IN, AND .FSP ###
    ##########################################################################

    with open(directory / "fsp_sol_file.txt", "r") as fsp:
        nft = 0
        for line in fsp:
            if line.startswith("% Loc"):
                hypo_lat = float(line.split()[5])
                hypo_lon = float(line.split()[8])
                hypo_dep = float(line.split()[11])
            if line.startswith("% Size"):
                Mw = float(line.split()[13])
                M0 = float(line.split()[16])
            if line.startswith("% Invs : Nx "):
                n_sf_strike = int(line.split()[5])  # Number of subfaults along strike
                n_sf_dip = int(line.split()[8])  # Number of subfaults along dip
            if line.startswith("% Invs : Dx "):
                sf_length_as = float(line.split()[5])  # Subfault length along strike
                sf_length_ad = float(line.split()[9])  # Subfault length along dip
            if line.startswith("% Nsbfs"):
                nft += int(line.split()[3])  # number of subfaults
    ##########################################################################
    ########################### WRITE HEADER INFO ############################
    ##########################################################################
    with open(directory / "Coulomb.inp", "w") as coulomb:
        coulomb.write(
            str(eventID)
            + "_coulomb.inp "
            + "inverted by Dara Goldberg (USGS/NEIC) with Mo = "
            + str("{:.2e}".format(M0))
            + "N-m, and Mw = "
            + str("{:.2f}".format(Mw))
            + "\n"
        )
        coulomb.write(
            "See Hayes(2017), The finite, kinematic rupture properties of great-sized earthquakes since 1990, EPSL 468, 94-100\n"
        )
        coulomb.write("#reg1=  0  #reg2=  0  #fixed= " + str(nft) + "  sym=  1\n")
        coulomb.write(" PR1=       0.250     PR2=       0.250   DEPTH=       10.000\n")
        coulomb.write("  E1=      8.000e+05   E2=      8.000e+05\n")
        coulomb.write("XSYM=       .000     YSYM=       .000\n")
        coulomb.write("FRIC=          0.400\n")
        coulomb.write(
            "S1DR=         19.000 S1DP=         -0.010 S1IN=        100.000 S1GD=          0.000\n"
        )
        coulomb.write(
            "S2DR=         89.990 S2DP=         89.990 S2IN=         30.000 S2GD=          0.000\n"
        )
        coulomb.write(
            "S3DR=        109.000 S3DP=         -0.010 S3IN=          0.000 S3GD=          0.000\n"
        )
        coulomb.write("\n")
        coulomb.write(
            "  #   X-start    Y-start     X-fin      Y-fin   Kode  rake     netslip   dip angle     top      bot\n"
        )
        coulomb.write(
            "xxx xxxxxxxxxx xxxxxxxxxx xxxxxxxxxx xxxxxxxxxx xxx xxxxxxxxxx xxxxxxxxxx xxxxxxxxxx xxxxxxxxxx xxxxxxxxxx\n"
        )

        ##########################################################################
        ########################## SUBFAULT BY SUBFAULT ##########################
        ##########################################################################

        LAT = []
        LON = []
        DEP = []
        SLIP = []
        RAKE = []
        STRIKE = []
        DIP = []

        with open(directory / "Solucion.txt") as param:
            for line in param:
                if "#" in line:  # HEADER LINES
                    continue
                elif len(array(line.split())) < 4:  # FAULT BOUNDARY LINES
                    continue
                else:
                    (
                        lat,
                        lon,
                        dep,
                        slip,
                        rake,
                        strike,
                        dip,
                        t_rup,
                        t_ris,
                        t_fal,
                        mo,
                    ) = line.split()
                    STRIKE.append(float(strike))
                    DIP.append(float(dip))

        with open(directory / "fsp_sol_file.txt", "r") as fsp:
            for line in fsp:
                if line.startswith("%"):
                    continue
                else:
                    (
                        lat,
                        lon,
                        ew,
                        ns,
                        dep,
                        slip,
                        rake,
                        trup,
                        trise,
                        sf_moment,
                    ) = line.split()
                    LAT.append(float(lat))
                    LON.append(float(lon))
                    DEP.append(float(dep))
                    SLIP.append(float(slip))
                    if float(rake) > 360:
                        RAKE.append(float(rake) - 360)
                    else:
                        RAKE.append(float(rake))
        ##########################################################################
        ######### CONVERT COORDINATE SUBFAULT CENTERS TO LOCAL CARTESIAN #########
        ##########################################################################

        Nfaults = len(LAT)
        x = zeros(Nfaults)
        y = zeros(Nfaults)
        g = pyproj.Geod(ellps="WGS84")  # Use WGS84 Ellipsoid
        for nsubfault in range(Nfaults):
            backaz, az, dist_m = pyproj.Geod.inv(
                g, LON[nsubfault], LAT[nsubfault], hypo_lon, hypo_lat
            )
            x[nsubfault] = (dist_m / 1000.0) * sin(deg2rad(az))
            y[nsubfault] = (dist_m / 1000.0) * cos(deg2rad(az))

        top_mid_x = zeros(Nfaults)
        top_mid_y = zeros(Nfaults)
        xstart = zeros(Nfaults)
        ystart = zeros(Nfaults)
        xfin = zeros(Nfaults)
        yfin = zeros(Nfaults)
        ztop = zeros(Nfaults)
        zbot = zeros(Nfaults)

        for ksubfault in range(Nfaults):
            ### ASSUMES .FSP FILE COORDINATES ARE THE **CENTER** OF THE SUBFAULT ###
            top_mid_x[ksubfault] = x[ksubfault] + (
                (sf_length_as / 2) * cos(deg2rad(DIP[ksubfault]))
            ) * sin(deg2rad(STRIKE[ksubfault] - 90))
            top_mid_y[ksubfault] = y[ksubfault] + (
                (sf_length_as / 2) * cos(deg2rad(DIP[ksubfault]))
            ) * cos(deg2rad(STRIKE[ksubfault] - 90))
            xstart[ksubfault] = top_mid_x[ksubfault] + (sf_length_as / 2) * sin(
                deg2rad(STRIKE[ksubfault] - 180)
            )
            ystart[ksubfault] = top_mid_y[ksubfault] + (sf_length_as / 2) * cos(
                deg2rad(STRIKE[ksubfault] - 180)
            )
            xfin[ksubfault] = top_mid_x[ksubfault] + (sf_length_as / 2) * sin(
                deg2rad(STRIKE[ksubfault])
            )
            yfin[ksubfault] = top_mid_y[ksubfault] + (sf_length_as / 2) * cos(
                deg2rad(STRIKE[ksubfault])
            )
            ztop[ksubfault] = DEP[ksubfault] - (sf_length_ad / 2) * sin(
                deg2rad(DIP[ksubfault])
            )
            zbot[ksubfault] = DEP[ksubfault] + (sf_length_ad / 2) * sin(
                deg2rad(DIP[ksubfault])
            )
            out = (
                "  1 %10.4f %10.4f %10.4f %10.4f 100 %10.4f %10.4f %10.4f %10.4f %10.4f FFM%d\n"
                % (
                    xstart[ksubfault],
                    ystart[ksubfault],
                    xfin[ksubfault],
                    yfin[ksubfault],
                    RAKE[ksubfault],
                    SLIP[ksubfault],
                    DIP[ksubfault],
                    ztop[ksubfault],
                    zbot[ksubfault],
                    ksubfault,
                )
            )
            coulomb.write(out)

        ##########################################################################
        ########################### WRITE FOOTER INFO ############################
        ##########################################################################

        deg2km = 111.19
        gsize = 6

        sx = -(gsize / 2.0) * deg2km * cos(deg2rad(hypo_lat))
        ex = (gsize / 2.0) * deg2km * cos(deg2rad(hypo_lat))
        dx = (ex - sx) / 100.0
        sy = -(gsize / 2.0) * deg2km
        ey = gsize / 2.0 * deg2km
        dy = (ey - sy) / 100.0
        min_lat = hypo_lat - (gsize / 2.0)
        max_lat = hypo_lat + (gsize / 2.0)
        min_lon = hypo_lon - (gsize / 2.0)
        max_lon = hypo_lon + (gsize / 2.0)
        dz = 1.0
        eqdep = -hypo_dep

        coulomb.write("\n")
        coulomb.write("     Grid Parameters\n")
        coulomb.write(
            "1 ----------------------------  Start-x =" + "{:16.7f}".format(sx) + "\n"
        )
        coulomb.write(
            "2 ----------------------------  Start-y =" + "{:16.7f}".format(sy) + "\n"
        )
        coulomb.write(
            "3 --------------------------   Finish-x =" + "{:16.7f}".format(ex) + "\n"
        )
        coulomb.write(
            "4 --------------------------   Finish-y =" + "{:16.7f}".format(ey) + "\n"
        )
        coulomb.write(
            "5 -----------------------   x-increment =" + "{:16.7f}".format(dx) + "\n"
        )
        coulomb.write(
            "6 -----------------------   y-increment =" + "{:16.7f}".format(dy) + "\n"
        )
        coulomb.write("\n")

        coulomb.write("     Size Parameters\n")
        coulomb.write(
            "1 --------------------------  Plot size ="
            + "{:16.7f}".format(gsize)
            + "\n"
        )
        coulomb.write(
            "2 --------------  Shade/Color increment ="
            + "{:16.7f}".format(float(1.0))
            + "\n"
        )
        coulomb.write(
            "3 ------  Exaggeration for disp.& dist. ="
            + "{:16.7f}".format(float(10000))
            + "\n"
        )
        coulomb.write("\n")

        coulomb.write("     Cross section default\n")
        coulomb.write(
            "1 ----------------------------  Start-x ="
            + "{:16.7f}".format(min_lon)
            + "\n"
        )
        coulomb.write(
            "2 ----------------------------  Start-y ="
            + "{:16.7f}".format(min_lat)
            + "\n"
        )
        coulomb.write(
            "3 ---------------------------  Finish-x ="
            + "{:16.7f}".format(max_lon)
            + "\n"
        )
        coulomb.write(
            "4 ---------------------------  Finish-y ="
            + "{:16.7f}".format(max_lat)
            + "\n"
        )
        coulomb.write(
            "5 ------------------  Distant-increment =" + "{:16.7f}".format(dx) + "\n"
        )
        coulomb.write(
            "6 ----------------------------- Z-depth ="
            + "{:16.7f}".format(eqdep)
            + "\n"
        )
        coulomb.write(
            "7 ------------------------  Z-increment =" + "{:16.7f}".format(dz) + "\n"
        )
        coulomb.write("\n")

        coulomb.write("     Map info\n")
        coulomb.write(
            "1 ---------------------------  min. lon ="
            + "{:16.7f}".format(min_lon)
            + "\n"
        )
        coulomb.write(
            "2 ---------------------------  max. lon ="
            + "{:16.7f}".format(max_lon)
            + "\n"
        )
        coulomb.write(
            "3 ---------------------------  zero lon ="
            + "{:16.7f}".format(hypo_lon)
            + "\n"
        )
        coulomb.write(
            "4 ---------------------------  min. lat ="
            + "{:16.7f}".format(min_lat)
            + "\n"
        )
        coulomb.write(
            "5 ---------------------------  max. lat ="
            + "{:16.7f}".format(max_lat)
            + "\n"
        )
        coulomb.write(
            "6 ---------------------------  zero lat ="
            + "{:16.7f}".format(hypo_lat)
            + "\n"
        )
        coulomb.write(
            "7 ------------------------  Z-increment =" + "{:16.7f}".format(dz) + "\n"
        )


### OKADA FUNCTIONS ###
def write_Okada_displacements(
    pdefile: Union[pathlib.Path, str],
    directory: Union[pathlib.Path, str] = pathlib.Path(),
):
    """Write out the Okada displacements
    :param pdefile: The path to the CMT file
    :type pdefile: Union[pathlib.Path, str], optional
    :param directory: The directory where to write the file(s), defaults to pathlib.Path()
    :type directory: Union[pathlib.Path, str], optional
    """
    print("Writing Okada Displacement File...")

    directory = pathlib.Path(directory)
    ##########################
    #### FAULT INFORMATION ###
    ##########################
    pdefile = pathlib.Path(pdefile)
    with open(pdefile, "r") as pdef:
        lines = pdef.readlines()
        pde = lines[0]
        hypo_lat = float(pde.split()[7])
        hypo_lon = float(pde.split()[8])
    fault_lat: List[float] = []
    fault_lon: List[float] = []
    fault_depth: List[float] = []
    fault_strike_slip: List[float] = []
    fault_dip_slip: List[float] = []
    fault_tensile_slip: List[float] = []
    fault_strike: List[float] = []
    fault_dip: List[float] = []
    with open(directory / "Solucion.txt", "r") as sol:
        for line in sol:
            if line.startswith("#Fault_segment"):
                sf_length = float(line.split()[7].split("km")[0])
                sf_width = float(line.split()[12].split("km")[0])
                fault_half_length = sf_length / 2
                fault_half_width = sf_width / 2
            elif "#" in line:  # HEADER LINES
                continue
            elif len(array(line.split())) < 4:  # FAULT BOUNDARY LINES
                continue
            else:  # ACTUAL SUBFAULT DETAILS
                (
                    lat,
                    lon,
                    dep,
                    slip,
                    rake,
                    strike,
                    dip,
                    t_rup,
                    t_ris,
                    t_fal,
                    mo,
                ) = [float(v) for v in line.split()]
                # Make rake be -180 to 180 (not 0-360)
                if float(rake) > 180:
                    rake = float(rake) - 360
                # Make Lon 0 to 360 (not -180 to 180)
                if float(lon) < 0:
                    lon = float(lon) + 360
                fault_lat.append(float(lat))
                fault_lon.append(float(lon))
                fault_depth.append(float(dep))
                fault_strike_slip.append(float(slip) / 100.0 * cos(deg2rad(rake)))
                fault_dip_slip.append(float(slip) / 100.0 * sin(deg2rad(rake)))
                fault_tensile_slip.append(0.0)
                fault_strike.append(float(strike))
                fault_dip.append(float(dip))
    Nsubfaults = len(fault_lon)
    obs_n = 40
    obs_x_vec = linspace(-400, 400, obs_n)
    obs_y_vec = linspace(-400, 400, obs_n)
    obs_z_vec = zeros(obs_n)
    obs_x_mat, obs_y_mat = meshgrid(obs_x_vec, obs_y_vec)
    obs_z_mat = zeros((obs_n, obs_n))  # performing calculation at the surface, z = 0

    # Make lon lat grid

    g = pyproj.Geod(ellps="WGS84")  # Use WGS84 Ellipsoid
    _, xLats, _ = pyproj.Geod.fwd(
        g,
        hypo_lon * ones(len(obs_x_vec)),
        hypo_lat * ones(len(obs_x_vec)),
        0 * ones(len(obs_x_vec)),
        obs_x_vec * 1000,
    )
    xLons, _, _ = pyproj.Geod.fwd(
        g,
        hypo_lon * ones(len(obs_y_vec)),
        hypo_lat * ones(len(obs_y_vec)),
        90 * ones(len(obs_y_vec)),
        obs_y_vec * 1000,
    )

    # Make Lon 0 to 360 (not -180 to 180) to avoid plotting issues
    for kLon in range(len(xLons)):
        if xLons[kLon] < 0:
            xLons[kLon] = 360 + xLons[kLon]

    gridLon, gridLat = meshgrid(xLons, xLats)
    fault_mu = 1.0
    fault_poisson_ratio = 0.25
    fault_lmbda = 2 * fault_mu * fault_poisson_ratio / (1 - 2 * fault_poisson_ratio)
    fault_alpha = (fault_lmbda + fault_mu) / (fault_lmbda + 2 * fault_mu)

    # Initialize total
    ux_cutde_total = zeros((obs_n, obs_n))
    uy_cutde_total = zeros((obs_n, obs_n))
    uz_cutde_total = zeros((obs_n, obs_n))

    for ksub in range(Nsubfaults):
        if ksub % 10 == 0:
            print(f"...Subfault: {ksub} of {Nsubfaults}")
        fault_top_depth, fault_bottom_depth = get_top_bottom_from_center(
            fault_depth[ksub], 2 * fault_half_width, fault_dip[ksub]
        )
        # Check for surface breaking
        if fault_top_depth < 0:
            warnings.warn(
                "Warning! Fault might be breaking the surface with top depth of %f"
                % fault_top_depth
            )
            fault_top_depth = 0

        # Check for negative Lon and make Lon 0 to 360 (not -180 to 180)
        if float(fault_lon[ksub]) < 0:
            fault_lon[ksub] = float(fault_lon[ksub]) + 360

        ### Make new observation x and y vectors, which are a grid of distances (in km)  ###
        ###from current subfault (rather than hypocenter) to lon/lat grid defined earlier ###
        fwd_az, b_az, distance_m = pyproj.Geod.inv(
            g,
            fault_lon[ksub] * ones(shape(xLons)),
            fault_lat[ksub] * ones(shape(xLats)),
            xLons,
            xLats,
        )
        distance_km = distance_m / 1000.0
        obs_x_vec_ll = distance_km * sin(deg2rad(fwd_az))
        obs_y_vec_ll = distance_km * cos(deg2rad(fwd_az))

        obs_x_mat, obs_y_mat = meshgrid(obs_x_vec_ll, obs_y_vec_ll)

        # cutde displacement calculation
        ux_okada_cutde = zeros((obs_n, obs_n))
        uy_okada_cutde = zeros((obs_n, obs_n))
        uz_okada_cutde = zeros((obs_n, obs_n))

        rot_x_mat, rot_y_mat = strike_rotation(
            fault_strike[ksub], obs_x_mat, obs_y_mat, direction="fwd"
        )

        # Call the displacement calculation
        (
            ux_okada_cutde,
            uy_okada_cutde,
            uz_okada_cutde,
        ) = get_gridded_okada_displacements_cutde(
            rot_x_mat,
            rot_y_mat,
            obs_z_mat,
            fault_depth[ksub],
            fault_dip[ksub],
            fault_half_length,
            fault_half_width,
            fault_strike_slip[ksub],
            fault_dip_slip[ksub],
            fault_tensile_slip[ksub],
            fault_alpha,
        )

        # un rotate
        unrot_ux_okada_cutde, unrot_uy_okada_cutde = strike_rotation(
            fault_strike[ksub], ux_okada_cutde, uy_okada_cutde, direction="inv"
        )

        ux_cutde_total += unrot_ux_okada_cutde
        uy_cutde_total += unrot_uy_okada_cutde
        uz_cutde_total += uz_okada_cutde

    with open(directory / "surface_deformation.disp", "w") as DISPout:
        DISPout.write(
            "#Longitude, Latitude, Elevation, Easting Displacement (m), Northing Displacement (m), Vertical Displacement (m)\n"
        )
        for kpt in range(len(gridLon.flatten())):
            DISPout.write(
                "%10.4f \t %10.4f \t %10.4f \t %10.4f \t %10.4f \t %10.4f \n"
                % (
                    gridLon.flatten()[kpt],
                    gridLat.flatten()[kpt],
                    0,
                    ux_cutde_total.flatten()[kpt],
                    uy_cutde_total.flatten()[kpt],
                    uz_cutde_total.flatten()[kpt],
                )
            )

    plot_okada_map(
        directory, gridLon, gridLat, ux_cutde_total, uy_cutde_total, uz_cutde_total
    )


def get_top_bottom_from_center(center_depth, width, dip):
    """
    Get the top and bottom depth of a rectangular fault from its width, center, and dip.

    :param center_depth: depth of center of fault plane, in km (positive down)
    :type center_depth: float
    :param width: total downdip width of rectangular fault plane, in km
    :type width: float
    :param dip: dip of fault plane, in degrees (range 0 to 90)
    :type dip: float
    :return: top and bottom of fault plane, in km
    :rtype: float, float
    """
    top = center_depth - (width / 2.0 * sin(deg2rad(dip)))
    bottom = center_depth + (width / 2.0 * sin(deg2rad(dip)))
    return top, bottom


def get_four_corners_eastwest_fault(fault_depth, dip, half_length, half_width):
    """
    Assuming a rectangular fault strikes east-west in a 3D cartesian space,
    and the coordinate (0, 0, depth) is located at its geometric center,
    what are the coordinates of the four corners?

    :param fault_depth: center depth of fault, in km (positive down)
    :type fault_depth: float
    :param dip: dip of fault plane, in degrees (range 0 to 90)
    :type dip: float
    :param half_length: half the along-strike subfault length (in km)
    :type half_length: float
    :param half_width: half the along-dip subfault width (in km)
    :type half_width: float
    :return: corners of subfault
    :rtype: array
    """
    top, bottom = get_top_bottom_from_center(fault_depth, half_width * 2, dip)
    west_point, east_point = -half_length, half_length
    north_point, south_point = half_width * cos(deg2rad(dip)), -half_width * cos(
        deg2rad(dip)
    )
    c1 = (west_point, north_point, top)
    c2 = (east_point, north_point, top)
    c3 = (east_point, south_point, bottom)
    c4 = (west_point, south_point, bottom)
    cutde_fault_pts = array(
        [
            [west_point, north_point, -top],
            [east_point, north_point, -top],
            [east_point, south_point, -bottom],
            [
                west_point,
                south_point,
                -bottom,
            ],  # cutde takes depths as negative downward
        ]
    )
    return cutde_fault_pts


def okada_wrapper_to_cutde(half_length, half_width, depth, dip, ss, ds, ts):
    """
    Organize info needed to treat cutde like former okada-wrapper

    :param half_length: half the along-strike subfault length (in km)
    :type half_length: float
    :param half_width: half the along-dip subfault width (in km)
    :type half_width: float
    :param depth: center depth of fault, in km (positive down)
    :type depth: float
    :param dip: dip of fault plane, in degrees (range 0 to 90)
    :type dip: float
    :param ss: strike-slip component
    :type ss: float
    :param ds: dip-slip component
    :type ds: float
    :param ts: tensile-slip component
    :type ts: float
    :return: fault corners, indices of fault triangles, slip vector
    :rtype: array, array, array
    """
    # Dip is now implemented
    cutde_fault_pts = get_four_corners_eastwest_fault(
        depth, dip, half_length, half_width
    )
    cutde_fault_tris = array([[0, 1, 2], [0, 2, 3]], dtype=int64)
    cutde_slip = array([[ss, ds, ts], [ss, ds, ts]])
    return cutde_fault_pts, cutde_fault_tris, cutde_slip


def dc3dwrapper_cutde(
    fault_alpha,
    coords,
    fault_depth,
    fault_dip,
    fault_half_lengths,
    fault_half_widths,
    slip,
):
    """
    Estimate displacement from subfault

    :param fault_alpha: fault property
    :type fault_alpha: float
    :param coords: subfault corner coordinates
    :type coords: array
    :param fault_depth: center depth of subfault, in km (positive down)
    :type fault_depth: float
    :param fault_dip: dip of fault plane, in degrees (range 0 to 90)
    :type fault_dip: float
    :param fault_half_lengths: half the along-strike subfault length (in km)
    :type fault_half_lengths: float
    :param fault_half_widths: half the along-dip subfault width (in km)
    :type fault_half_widths: flat
    :param slip: strike-slip, dip-slip, tensile-slip components of slip (m)
    :type slip: array
    :return: success flag, displacement array, strain array
    :rtype: int, array, array
    """
    # Material property conversion
    fault_nu = (1 - 2 * fault_alpha) / (-2 * fault_alpha)

    # Geometry and slip conversion
    cutde_fault_pts, cutde_fault_tris, cutde_slip = okada_wrapper_to_cutde(
        fault_half_lengths[1],
        fault_half_widths[1],
        fault_depth,
        fault_dip,
        slip[0],
        -slip[1],
        slip[2],
    )

    obs_pts = (
        array([coords[0], coords[1], -coords[2]]).reshape((3, -1)).T.copy()
    )  # negative sign for depth in cutde

    # Calculate cutde displacements partial derivatives
    u_mat_cutde = half_space.disp_matrix(
        obs_pts=obs_pts,
        tris=cutde_fault_pts[cutde_fault_tris],
        nu=fault_nu,
    )

    # Multiply by displacements
    u_cutde = u_mat_cutde.reshape((-1, 6)).dot(cutde_slip.flatten())

    # Get Displacement Gradients:
    strain_mat_cutde = half_space.strain_matrix(
        obs_pts=obs_pts,
        tris=cutde_fault_pts[cutde_fault_tris],
        nu=fault_nu,
    )
    strain = strain_mat_cutde.reshape((-1, size(cutde_slip))).dot(
        cutde_slip.flatten()
    )  # reshape by len total slip vector
    # strain[:,0] is the xx component of strain, 1 is yy, 2 is zz, 3 is xy, 4 is xz, and 5 is yz.
    strain_tensor = array(
        [
            [strain[0], strain[3], strain[4]],
            [strain[3], strain[1], strain[5]],
            [strain[4], strain[5], strain[2]],
        ]
    )

    success = 1
    grad_u_cutde = 0  # cutde produces the strain tensor directly, not the displacement gradient tensor which must be converted to strain
    return success, u_cutde, strain_tensor


def get_gridded_okada_displacements_cutde(
    obs_x_mat,
    obs_y_mat,
    obs_z_mat,
    fault_depth,
    fault_dip,
    fault_half_length,
    fault_half_width,
    fault_strike_slip,
    fault_dip_slip,
    fault_tensile_slip,
    fault_alpha,
):
    """
    Calculate displacement across area of interest from single subfault
    :param obs_x_mat: matrix of observation points, x-direction
    :type obs_x_mat: array
    :param obs_y_mat: matrix of observation points, y-direction
    :type obs_y_mat: array
    :param obs_z_mat: matrix of observation points, z-directory
    :type obs_z_mat: array
    :param fault_depth: center depth of subfault, in km (positive down)
    :type fault_depth: float
    :param fault_dip: dip of fault plane, in degrees (range 0 to 90)
    :type fault_dip: float
    :param fault_half_length: half the along-strike subfault length (in km)
    :type fault_half_length: float
    :param fault_half_width: half the along-dip subfault width (in km)
    :type fault_half_width: float
    :param fault_strike_slip: slip in the strike-direction
    :type fault_strike_slip: float
    :param fault_dip_slip: slip in the dip-direction
    :type fault_dip_slip: float
    :param fault_tensile_slip: slip in the tensile direction
    :type fault_tensile_slip: float
    :param fault_alpha: fault property
    :type fault_alpha: float
    :return: displacements across observation grid
    :rtype: array, array, array
    """
    ux_okada_cutde = zeros(len(obs_x_mat.flatten()))
    uy_okada_cutde = zeros(len(obs_x_mat.flatten()))
    uz_okada_cutde = zeros(len(obs_x_mat.flatten()))
    for kpt in range(len(obs_x_mat.flatten())):
        _, u, strain_tensor = dc3dwrapper_cutde(
            fault_alpha,
            [
                obs_x_mat.flatten()[kpt],
                obs_y_mat.flatten()[kpt],
                obs_z_mat.flatten()[kpt],
            ],  # variable depths
            fault_depth,
            fault_dip,
            [-fault_half_length, fault_half_length],
            [-fault_half_width, fault_half_width],
            [fault_strike_slip, fault_dip_slip, fault_tensile_slip],
        )
        ux_okada_cutde[kpt] = u[0]
        uy_okada_cutde[kpt] = u[1]
        uz_okada_cutde[kpt] = u[2]
    ux_okada_cutde = ux_okada_cutde.reshape(shape(obs_x_mat))
    uy_okada_cutde = uy_okada_cutde.reshape(shape(obs_x_mat))
    uz_okada_cutde = uz_okada_cutde.reshape(shape(obs_x_mat))
    return ux_okada_cutde, uy_okada_cutde, uz_okada_cutde


def strike_rotation(fault_strike, x_matrix_in, y_matrix_in, direction="fwd"):
    """
    Code assumes E-W fault, this function rotates observation positions by strike angle

    :param fault_strike: strike of fault (degrees, 0-360)
    :type fault_strike: float
    :param x_matrix_in: original observation matrix, x-coordinates
    :type x_matrix_in: array
    :param y_matrix_in: original observation matrix, y-coordinates
    :type y_matrix_in: array
    :param direction: rotate into strike direction (fwd, default) or out of strike direction (inv)
    :type direction: str
    :return: rotated observation matrix, x-coordinates and rotated observation matrix, y-coordinates
    :rtype: array, array
    """
    # Rotate for strike:
    ### Rotation matrices (code assumes E-W fault, must rotate positions by strike angle and then rotate back) ###
    if direction == "fwd":
        theta = fault_strike - 90
    elif direction == "inv":
        theta = -(fault_strike - 90)
    else:
        print("Invalid direction, choose either 'fwd' or 'inv'")
        return
    theta = deg2rad(theta)
    rotation_matrix = array([[cos(theta), -sin(theta)], [sin(theta), cos(theta)]])

    rot_x_mat = zeros(len(x_matrix_in.flatten()))
    rot_y_mat = zeros(len(y_matrix_in.flatten()))

    for kxy in range(len(rot_x_mat)):
        xy = rotation_matrix.dot(
            array([[x_matrix_in.flatten()[kxy]], [y_matrix_in.flatten()[kxy]]])
        )
        rot_x_mat[kxy] = xy[0][0]
        rot_y_mat[kxy] = xy[1][0]
    rot_x_mat = rot_x_mat.reshape(shape(x_matrix_in))
    rot_y_mat = rot_y_mat.reshape(shape(y_matrix_in))
    return rot_x_mat, rot_y_mat


def plot_okada_map(
    directory, gridx, gridy, ux_okada_cutde, uy_okada_cutde, uz_okada_cutde
):
    """
    Plot map of Okada displacements

    :param directory: The directory where to write the file(s), defaults to pathlib.Path() within write_Okada_displacements
    :type directory: Union[pathlib.Path, str], optional
    :param gridx: observation grid, x-coordinate
    :type gridx: array
    :param gridy: observation grid, y-coordinate
    :type gridy: array
    :param ux_okada_cutde: x-coordinate displacement
    :type ux_okada_cutde: array
    :param uy_okada_cutde: y-coordinate displacement
    :type uy_okada_cutde: array
    :param uz_okada_cutde: z-coordinate displacement
    :type uz_okada_cutde: array
    """

    default_dirs = mng.default_dirs()
    min_lon = min(gridx.flatten())
    max_lon = max(gridx.flatten())
    min_lat = min(gridy.flatten())
    max_lat = max(gridy.flatten())
    region = [min_lon, max_lon, min_lat, max_lat]
    ### Fix region to nearest tenth of a degree ###
    region[0] = floor(region[0] * 10) / 10.0
    region[1] = ceil(region[1] * 10) / 10.0
    region[2] = floor(region[2] * 10) / 10.0
    region[3] = ceil(region[3] * 10) / 10.0

    ################################
    ### PLOT BASEMAP/ COASTLINES ###
    ################################

    fig = pygmt.Figure()
    resolution = "03s"
    map_scale_len = "50k"
    if region[1] - region[0] > 5 or region[3] - region[2] > 5:
        resolution = "30s"
        map_scale_len = "100k"
    if region[1] - region[0] > 10 or region[3] - region[2] > 10:
        resolution = "01m"
        map_scale_len = "200k"
    grid = pygmt.datasets.load_earth_relief(resolution=resolution, region=region)
    projection = "M10c"
    map_scale = (
        "g"
        + str(region[0])
        + "/"
        + str(region[2])
        + "+c17.40+w"
        + map_scale_len
        + "+ar+l+jBL+o0.5/0.5+f"
    )
    frame = ["WSen", "xa1f0.5", "ya1f0.5"]

    pygmt.config(
        PS_MEDIA="A0",
        MAP_FRAME_TYPE="plain",
        MAP_FRAME_AXES="WSen",
        FORMAT_GEO_OUT="F",
        FORMAT_GEO_MAP="ddd:mm:ss",
        FONT_ANNOT_PRIMARY="14p,Helvetica,black",
        FONT_ANNOT_SECONDARY="12p,Helvetica,black",
        FONT_HEADING="30p,Helvetica,black",
        FONT_LABEL="14p,Helvetica,black",
        FONT_LOGO="6p,Helvetica,black",
        FONT_SUBTITLE="16p,Helvetica,black",
        FONT_TAG="18p,Helvetica,black",
        FONT_TITLE="18p,Helvetica,black",
        MAP_ANNOT_OFFSET_PRIMARY="3p",
    )

    ### Horizontal Displacement ###
    fig.basemap(
        region=region,
        projection=projection,
        frame=["WSen+tHorizontal Surface Displacement", "xa1f0.5", "ya1f0.5"],
        map_scale=map_scale,
    )
    horizontal_cutde = (ux_okada_cutde**2 + uy_okada_cutde**2) ** 0.5
    max_horiz = max(horizontal_cutde.flatten())
    xx, yy, zz = gridx.flatten(), gridy.flatten(), horizontal_cutde.flatten()
    grid = pygmt.xyz2grd(x=xx, y=yy, z=zz, spacing=(0.3), region=region)
    minmax_horiz = (ceil(max_horiz * 10)) / 10.0  # ceil to the nearest 0.1
    annotation = (floor(minmax_horiz * 10 / 2)) / 10  # contour annotation
    if annotation == 0:
        annotation = 0.1
    pygmt.makecpt(
        cmap=str(default_dirs["root_dir"]) + "/src/wasp/lajolla_white.cpt",
        reverse=True,
        series=[0, minmax_horiz],
    )
    fig.grdimage(grid=grid, projection=projection, frame=frame)
    fig.grdcontour(
        # Pass in the grid made above
        grid=grid,
        # Set the interval for annotated contour lines at 0.5m
        annotation=annotation,
        projection=projection,
        pen="1p,gray40",
    )
    fig.colorbar(frame=["x+lDisplacement (m)"])
    nth = 31
    okada_positions = pd.DataFrame(
        data={
            "x": gridx.flatten()[::nth],
            "y": gridy.flatten()[::nth],
            "east_velocity": ux_okada_cutde.flatten()[::nth] / (max_horiz / 3),
            "north_velocity": uy_okada_cutde.flatten()[::nth] / (max_horiz / 3),
        }
    )

    fig.coast(resolution="h", shorelines="1p,black")
    fig.plot(
        str(default_dirs["root_dir"]) + "/pb2002_boundaries.gmt",
        style="f10/3p",
        region=region,
        pen="2p,white",
    )
    fig.plot(
        str(default_dirs["root_dir"]) + "/pb2002_boundaries.gmt",
        style="f10/3p",
        region=region,
        pen="1p,black",
    )

    fig.velo(
        data=okada_positions,
        region=region,
        projection=projection,
        spec="e1/0",
        pen="1,gray30",
        line=True,
        vector="0.2c+p1p,gray30+eA+a55+n",
    )

    ### Vertical Displacement ###
    fig.shift_origin(xshift="13c")
    fig.basemap(
        region=region,
        projection=projection,
        frame=["WSen+tVertical Surface Displacement", "xa1f0.5", "ya1f0.5"],
        map_scale=map_scale,
    )
    xx, yy, zz = gridx.flatten(), gridy.flatten(), uz_okada_cutde.flatten()
    grid = pygmt.xyz2grd(x=xx, y=yy, z=zz, spacing=(0.3), region=region)
    minmax_uz = (
        ceil(max(abs(uz_okada_cutde.flatten())) * 10)
    ) / 10.0  # ceil to the nearest 0.1
    annotation = (floor(minmax_uz * 10 / 2)) / 10  # contour annotation
    if annotation == 0:
        annotation = 0.1
    pygmt.makecpt(cmap="polar", series=[-minmax_uz, minmax_uz])
    fig.grdimage(grid=grid, projection=projection, frame=frame)
    fig.grdcontour(
        # Pass in the grid made above
        grid=grid,
        # Set the interval for annotated contour lines at 0.5m
        annotation=annotation,
        projection=projection,
        pen="1p,gray40",
    )
    fig.colorbar(frame=["x+lDisplacement (m)"])
    fig.coast(resolution="h", shorelines="1p,black")
    fig.plot(
        str(default_dirs["root_dir"]) + "/pb2002_boundaries.gmt",
        style="f10/3p",
        region=region,
        pen="2p,white",
    )
    fig.plot(
        str(default_dirs["root_dir"]) + "/pb2002_boundaries.gmt",
        style="f10/3p",
        region=region,
        pen="1p,black",
    )

    fig.savefig(directory / "Okada_Displacement.png")
