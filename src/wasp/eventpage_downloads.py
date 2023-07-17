# -*- coding: utf-8 -*-
"""
The routines here create the remaining files in the Downloads section of the NEIC eventpages
"""

import json
import os
import pathlib
import shutil
from collections import Counter
from glob import glob
from typing import List, Optional, Union

import matplotlib.pyplot as plt  # type: ignore
import pyproj  # type: ignore
from numpy import array, c_, cos, deg2rad, linspace, meshgrid, ones, shape, sin, zeros
from obspy.imaging.scripts import mopad  # type: ignore
from okada_wrapper import dc3dwrapper  # type: ignore


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
    num_shwaves = int(bodywave_count["SH"])
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


def write_Okada_displacements(directory: Union[pathlib.Path, str] = pathlib.Path()):
    """Write out the Okada displacements

    :param directory: The directory where to write the file(s), defaults to pathlib.Path()
    :type directory: Union[pathlib.Path, str], optional
    """
    print("Writing Okada Displacement File...")

    directory = pathlib.Path(directory)
    ##########################
    #### FAULT INFORMATION ###
    ##########################
    with open(directory / "fsp_sol_file.txt", "r") as fsp:
        for line in fsp:
            if line.startswith("% Loc"):
                hypo_lat = float(line.split()[5])
                hypo_lon = float(line.split()[8])
    with open(directory / "Solucion.txt", "r") as sol:
        for line in sol:
            if line.startswith("#Fault_segment"):
                sf_length_as = float(line.split()[7].split("km")[0])
                sf_length_ad = float(line.split()[12].split("km")[0])
    alpha = 2.0 / 3
    #########################################
    ### COORDINATES OF OBSERVATION POINTS ###
    #########################################
    x = linspace(-400, 400, 40)  # in km
    y = linspace(-400, 400, 40)  # in km
    g = pyproj.Geod(ellps="WGS84")  # Use WGS84 Ellipsoid
    # Convert to grid of lon/lat around epicenter
    _, xLats, _ = pyproj.Geod.fwd(
        g, hypo_lon * ones(len(x)), hypo_lat * ones(len(x)), 0 * ones(len(x)), x * 1000
    )
    xLons, _, _ = pyproj.Geod.fwd(
        g, hypo_lon * ones(len(y)), hypo_lat * ones(len(y)), 90 * ones(len(y)), y * 1000
    )
    # Make Lon 0 to 360 (not -180 to 180) to avoid plotting issues
    for kLon in range(len(xLons)):
        if xLons[kLon] < 0:
            xLons[kLon] = 360 + xLons[kLon]

    gridLon, gridLat = meshgrid(xLons, xLats)

    LAT: List[float] = []
    LON: List[float] = []
    DEP: List[float] = []
    SLIP: List[float] = []
    RAKE: List[float] = []
    STRIKE: List[float] = []
    DIP: List[float] = []

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
                # Make Lon 0 to 360 (not -180 to 180)
                if float(lon) < 0:
                    lon = float(lon) + 360
                LAT.append(float(lat))
                LON.append(float(lon))
                DEP.append(float(dep))
                SLIP.append(float(slip))
                RAKE.append(float(rake))
                STRIKE.append(float(strike))
                DIP.append(float(dip))

    ###############################
    ### CALCULATE DISPLACEMENTS ###
    ###############################
    Nfaults = len(LAT)
    ### Initialize northing, easting, and vertical displacements ###
    ux = zeros(len(xLats) * len(xLons))
    uy = zeros(len(xLats) * len(xLons))
    uz = zeros(len(xLats) * len(xLons))

    # for ksub in range(1):
    for ksub in range(Nfaults):
        if ksub % 100 == 0:
            print("...Subfault: " + str(ksub))
        sf_lon = LON[ksub]
        sf_lat = LAT[ksub]
        strike = STRIKE[ksub]
        dip = DIP[ksub]
        slip = SLIP[ksub]
        rake = RAKE[ksub]
        depth = DEP[ksub]
        ss_slip = (slip / 100.0) * cos(deg2rad(rake))
        ds_slip = (slip / 100.0) * sin(deg2rad(rake))

        ### Make new x,y vectors, grid of distances in km from current subfault to lon/lat grid defined earlier ###
        fwd_az, b_az, distance_m = pyproj.Geod.inv(
            g,
            sf_lon * ones(shape(gridLon)),
            sf_lat * ones(shape(gridLat)),
            gridLon,
            gridLat,
        )
        distance_km = distance_m / 1000.0
        x = distance_km * sin(deg2rad(fwd_az))
        y = distance_km * cos(deg2rad(fwd_az))
        ### Rotation matrices (code assumes N-S fault, must rotate positions by strike angle and then rotate back) ###
        theta = strike - 90
        theta = deg2rad(theta)
        R_fwd = array([[cos(theta), -sin(theta)], [sin(theta), cos(theta)]])
        R_inv = array([[cos(-theta), -sin(-theta)], [sin(-theta), cos(-theta)]])

        k = 0
        for kxy in range(len(x) * len(y)):
            # Calculate on rotated position
            xy = R_fwd.dot(array([[x.flatten()[kxy]], [y.flatten()[kxy]]]))
            success, u, grad_u = dc3dwrapper(
                alpha,
                [xy[0], xy[1], 0.0],
                depth,
                dip,
                [-sf_length_as / 2, sf_length_as / 2],
                [-sf_length_ad / 2, sf_length_ad / 2],
                [ss_slip, ds_slip, 0.0],
            )
            # "un-rotate"
            urot = R_inv.dot(array([[u[0]], [u[1]]]))
            u[0] = urot[0]
            u[1] = urot[1]

            # populate output vector
            ux[k] = ux[k] + u[0]
            uy[k] = uy[k] + u[1]
            uz[k] = uz[k] + u[2]
            k += 1
    with open(directory / "surface_deformation.disp", "w") as DISPout:
        DISPout.write(
            "#Longitude, Latitude, Elevation, Easting Displacement (m), Northing Displacement (m), Vertical Displacement (m)\n"
        )
        ko = 0
        # for klon in range(len(xLons)):
        #    for klat in range(len(xLats)):
        for klat in range(len(xLats)):
            for klon in range(len(xLons)):
                DISPout.write(
                    "%10.4f \t %10.4f \t %10.4f \t %10.4f \t %10.4f \t %10.4f \n"
                    % (xLons[klon], xLats[klat], 0, ux[ko], uy[ko], uz[ko])
                )
                ko += 1
    ############
    ### PLOT ###
    ###########
    plt.figure(figsize=(8, 6))
    horizontal = (ux**2 + uy**2) ** 0.5
    horizontal_grid = horizontal.reshape((shape(gridLon)))
    horiz = abs(horizontal)
    horiz.sort()
    minmax = horiz[-2]
    plt.scatter(
        gridLon,
        gridLat,
        marker="s",
        c=horizontal * 100,
        lw=0,
        s=100,
        vmin=0.0,
        vmax=minmax * 100,
        cmap="rainbow",
    )
    cb = plt.colorbar()
    cb.set_label("Horizontal Surface Displacement (cm)")
    i = 0
    for klat in range(len(xLats)):
        for klon in range(len(xLons)):
            if i % 15 == 0:
                plt.quiver(
                    xLons[klon],
                    xLats[klat],
                    ux[i] / horizontal[i],
                    uy[i] / horizontal[i],
                    pivot="mid",
                    linewidths=0.01,
                    width=0.005,
                    color="k",
                    scale=30,
                )
            i += 1
    plt.xlim([xLons.min(), xLons.max()])
    plt.ylim([xLats.min(), xLats.max()])
    plt.ylabel("Latitude")
    plt.xlabel("Longitude")
    plt.savefig(directory / "Horizontal_Surface_Displacement.png", dpi=300)
    plt.close()

    plt.figure(figsize=(8, 6))
    vertical_grid = uz.reshape((shape(gridLon)))
    abs_uz = abs(uz)
    abs_uz.sort()
    minmax = abs_uz[-2]
    plt.scatter(
        gridLon,
        gridLat,
        marker="s",
        c=uz * 100,
        lw=0,
        s=100,
        vmin=-minmax * 100,
        vmax=minmax * 100,
        cmap="coolwarm",
    )
    cb = plt.colorbar()
    cb.set_label("Vertical Surface Displacement (cm)")
    plt.xlim([xLons.min(), xLons.max()])
    plt.ylim([xLats.min(), xLats.max()])
    plt.ylabel("Latitude")
    plt.xlabel("Longitude")
    plt.savefig(directory / "Vertical_Surface_Displacement.png", dpi=300)
    plt.close()


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
    print(nft)
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
        # T_RUP=[]
        # T_RIS=[]
        # T_FAL=[]
        # MO=[]
        # SS_SLIP = []
        # DS_SLIP = []

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
