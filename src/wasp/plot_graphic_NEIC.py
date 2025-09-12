# -*- coding: utf-8 -*-
"""The routines here allow to plot the solution model FFM modelling, as well
as moment rate function, and waveform fits.
"""

import argparse
import collections
import errno
import glob
import json
import os
import pathlib
from datetime import datetime
from shutil import move
from typing import Any, List, Optional, Tuple, Union

import cartopy.crs as ccrs  # type: ignore
import cartopy.feature as cf  # type: ignore
import cartopy.io.shapereader as shpreader  # type: ignore
import matplotlib  # type: ignore
import numpy as np  # type: ignore
import pandas as pd  # type: ignore
import pygmt  # type: ignore
from matplotlib import ticker  # type: ignore
from matplotlib import colormaps, colors, gridspec  # type: ignore
from matplotlib import pyplot as plt
from matplotlib.axes import Axes
from matplotlib.colors import ListedColormap  # type: ignore
from matplotlib.image import AxesImage  # type: ignore
from obspy.imaging.beachball import beach, beachball  # type: ignore
from pyproj import Geod  # type: ignore
from pyrocko import moment_tensor as pmt  # type: ignore
from pyrocko import plot  # type: ignore
from pyrocko.plot import beachball  # type: ignore
from scipy.interpolate import griddata  # type: ignore

#
# local modules
#
import wasp.fault_plane as pf
import wasp.plane_management as pl_mng
import wasp.seismic_tensor as tensor
import wasp.shakemap_tools as shakemap
import wasp.velocity_models as mv
from wasp import get_outputs, load_ffm_model
from wasp.plot_maps_NEIC import plot_map, set_map_cartopy
from wasp.static2fsp import static_to_fsp
from wasp.static2srf import static_to_srf
from wasp.waveform_plots_NEIC import plot_waveform_fits

"""
Set colorbar for slip
"""
rm = 100  # amount of lines to remove on black end of magma_r
ad = 50  # how much at the zero end should be *just* white before transitioning to meet colors
magma_cpt = colormaps.get_cmap("magma_r")  # start with magma_r
white_bit = np.array([255 / 256, 250 / 256, 250 / 256, 1])  # create array of white
slip_cpt = magma_cpt(np.linspace(0, 1, 512))  # initialize slip_cpt
slip_cpt[rm:, :] = slip_cpt[0:-rm, :]  # move beginning up to remove black end
r_s = np.linspace(
    white_bit[0], slip_cpt[rm][0], rm - ad
)  # gradient from white to beginning of new magma
g_s = np.linspace(white_bit[1], slip_cpt[rm][1], rm - ad)
b_s = np.linspace(white_bit[2], slip_cpt[rm][2], rm - ad)
slip_cpt[ad:rm, :][:, 0] = r_s
slip_cpt[ad:rm, :][:, 1] = g_s
slip_cpt[ad:rm, :][:, 2] = b_s
slip_cpt[:ad, :] = white_bit
slipcpt = ListedColormap(slip_cpt)


def plot_ffm_sol(
    tensor_info: dict,
    segments_data: dict,
    point_sources: list,
    shear: list,
    solution: dict,
    default_dirs: dict,
    autosize: bool = False,
    max_val: Optional[float] = None,
    mr_time: int = False,
    files_str: Optional[dict] = None,
    stations_gps: Optional[zip] = None,
    stations_cgps: Optional[str] = None,
    legend_len: Optional[float] = None,
    scale: Optional[float] = None,
    limits: List[Optional[float]] = [None, None, None, None],
    separate_planes: bool = False,
    label_stations: bool = False,
    directory: Union[pathlib.Path, str] = pathlib.Path(),
):
    """Execute plotting routines

    :param tensor_info: The tensor information
    :type tensor_info: dict
    :param segments: The segment properties
    :type segments: dict
    :param point_sources: The point source locations
    :type point_sources: list
    :param shear: The shear moduli
    :type shear: list
    :param solution: The kinematic model read from Solution.txt
    :type solution: dict
    :param default_dirs: The location of default directories
    :type default_dirs: dict
    :param autosize: Autosize the plot, defaults to False
    :type autosize: bool, optional
    :param files_str: The stations file properties, defaults to None
    :type files_str:Optional[dict], optional
    :param mr_time: Add moment rate time, defaults to False
    :type mr_time: int, optional
    :param files_str: The stations file properties, defaults to None
    :type files_str:Optional[dict], optional
    :param stations_gps: The gps stations description, defaults to None
    :type stations_gps:  Optional[zip], optional
    :param stations_cgps: The cgps stations description, defaults to None
    :type stations_cgps: Optional[str], optional
    :param legend_len: The length of the legend, defaults to None
    :type legend_len: Optional[float], optional
    :param scale: The scale, defaults to None
    :type scale: Optional[int], optional
    :param limits: The extent of the map, defaults to [None, None, None, None]
    :type limits: List[Optional[float]], optional
    :param separate_planes: Whether separate files, defaults to False
    :type separate_planes: bool, optional
    :param label_stations: Whether to label the stations, defaults to False
    :type label_stations: bool, optional
    :param directory: The location of data files, defaults to pathlib.Path()
    :type directory: Union[pathlib.Path, str], optional
    """
    directory = pathlib.Path(directory)
    segments = segments_data["segments"]

    path_velmodel = os.path.join(directory, "velmodel_data.json")
    if not os.path.exists(path_velmodel):
        vel_model = mv.select_velmodel(
            tensor_info=tensor_info,
            default_dirs=default_dirs,
            directory=directory,
        )
    else:
        with open(path_velmodel) as v:
            vel_model = json.load(v)

    _plot_vel_model(vel_model, directory=directory)
    plot_moment_rate_function(
        segments_data,
        shear,
        point_sources,
        mr_time=mr_time,
        separate_planes=separate_planes,
        directory=directory,
    )
    PlotSlipDistribution(
        segments,
        point_sources,
        solution,
        autosize=autosize,
        max_val=max_val,
        directory=directory,
    )
    PlotMap(
        tensor_info,
        segments,
        point_sources,
        solution,
        default_dirs,
        files_str=files_str,
        stations_gps=stations_gps,
        stations_cgps=stations_cgps,
        max_slip=max_val,
        legend_len=legend_len,
        scale=scale,
        limits=limits,
        label_stations=label_stations,
        directory=directory,
    )
    PlotSlipTimes(segments, point_sources, solution, directory=directory)


def plot_misfit(
    used_data_type: List[str], directory: Union[pathlib.Path, str] = pathlib.Path()
):
    """Plot misfit of observed and synthetic data

    :param used_data_type: The list of data types used
    :type used_data_type: List[str]
    :param directory: The location of data files, defaults to pathlib.Path()
    :type directory: Union[pathlib.Path, str], optional
    :raises FileNotFoundError: When a data type's json file is not found
    """
    directory = pathlib.Path(directory)
    if "body" in used_data_type:
        if not os.path.isfile(directory / "tele_waves.json"):
            raise FileNotFoundError(
                errno.ENOENT, os.strerror(errno.ENOENT), "tele_waves.json"
            )
        with open(directory / "tele_waves.json") as t:
            traces_info = json.load(t)
        traces_info = get_outputs.get_data_dict(
            traces_info, syn_file="synthetics_body.txt", directory=directory
        )
        values = [["BHZ"], ["BHT"]]
        for components in values:
            plot_waveform_fits(
                traces_info, components, "body", plot_directory=directory
            )
    if "surf" in used_data_type:
        if not os.path.isfile(directory / "surf_waves.json"):
            raise FileNotFoundError(
                errno.ENOENT, os.strerror(errno.ENOENT), "surf_waves.json"
            )
        with open(directory / "surf_waves.json") as t:
            traces_info = json.load(t)
        traces_info = get_outputs.get_data_dict(
            traces_info, syn_file="synthetics_surf.txt", margin=0, directory=directory
        )
        values = [["BHZ"], ["BHT"]]
        for components in values:
            plot_waveform_fits(
                traces_info, components, "surf", plot_directory=directory
            )
    if "strong" in used_data_type:
        if not os.path.isfile(directory / "strong_motion_waves.json"):
            raise FileNotFoundError(
                errno.ENOENT, os.strerror(errno.ENOENT), "strong_motion_waves.json"
            )
        with open(directory / "strong_motion_waves.json") as t:
            traces_info = json.load(t)
        traces_info = get_outputs.get_data_dict(
            traces_info, syn_file="synthetics_strong.txt", directory=directory
        )
        values = [["HLZ", "HNZ"], ["HLE", "HNE"], ["HLN", "HNN"]]
        plot_waveform_fits(
            traces_info,
            values,
            "strong",
            start_margin=10,
            plot_directory=directory,
        )
    if "cgps" in used_data_type:
        if not os.path.isfile(directory / "cgps_waves.json"):
            raise FileNotFoundError(
                errno.ENOENT, os.strerror(errno.ENOENT), "cgps_waves.json"
            )
        with open(directory / "cgps_waves.json") as t:
            traces_info = json.load(t)
        traces_info = get_outputs.get_data_dict(
            traces_info, syn_file="synthetics_cgps.txt", directory=directory
        )
        values = [["LXZ", "LHZ", "LYZ"], ["LXE", "LHE", "LYE"], ["LXN", "LHN", "LYN"]]
        plot_waveform_fits(
            traces_info, values, "cgps", start_margin=10, plot_directory=directory
        )
    return


def _plot_vel_model(
    velmodel: dict,
    directory: Union[pathlib.Path, str] = pathlib.Path(),
):
    """Plot the seismic velocity model as a function of depth

    :param velmodel: The velocity model
    :type velmodel: dict
    :param directory: The location where to write plots, defaults to pathlib.Path()
    :type directory: Union[pathlib.Path, str], optional
    """
    directory = pathlib.Path(directory)
    print("Creating Velocity Model Plot...")
    p_vel = np.array(velmodel["p_vel"]).astype(float)
    sh_vel = np.array(velmodel["s_vel"]).astype(float)
    thick = np.array(velmodel["thick"]).astype(float)

    depths = np.zeros(len(thick) + 1)

    depths[1:] = np.cumsum(thick)
    depths = np.array([depth for depth in depths if depth < 70])
    depths = np.append([depths], [70])
    plt.plot((p_vel[0], p_vel[0]), (depths[0], depths[1]), "b-", label="P")
    plt.plot((sh_vel[0], sh_vel[0]), (depths[0], depths[1]), "r-", label="SH")
    j = len(depths) - 2
    for i in range(j):
        plt.plot((p_vel[i], p_vel[i]), (depths[i], depths[i + 1]), "b-")
        plt.plot((p_vel[i], p_vel[i + 1]), (depths[i + 1], depths[i + 1]), "b-")
        plt.plot((sh_vel[i], sh_vel[i]), (depths[i], depths[i + 1]), "r-")
        plt.plot((sh_vel[i], sh_vel[i + 1]), (depths[i + 1], depths[i + 1]), "r-")

    plt.plot((p_vel[j], p_vel[j]), (depths[j], depths[j + 1]), "b-")
    plt.plot((sh_vel[j], sh_vel[j]), (depths[j], depths[j + 1]), "r-")

    plt.title("Crust model")
    plt.xlabel("Body wave velocity $(km/s)$")
    plt.ylabel("Depth $(km)$")
    plt.legend(loc="upper right")

    ax = plt.gca()
    ax.invert_yaxis()
    plt.savefig(directory / "crust_body_wave_vel_model.png", bbox_inches="tight")
    plt.close()


def _PlotRuptTime(
    segments: dict,
    point_sources: list,
    solution: dict,
    directory: Union[pathlib.Path, str] = pathlib.Path(),
):
    """Plot time distribution based on the FFM solution model

    :param segments: The segment properties
    :type segments: dict
    :param point_sources: The point source locations
    :type point_sources: list
    :param solution: The kinematic solution read from Solution.txt
    :type solution: dict
    :param directory: The location where to write plots, defaults to pathlib.Path()
    :type directory: Union[pathlib.Path, str], optional
    """
    directory = pathlib.Path(directory)
    print("Creating Rupture Time Plot...")
    rupt_time = solution["rupture_time"]
    max_rupt_time = [np.max(rupt_time_seg.flatten()) for rupt_time_seg in rupt_time]
    max_rupt_time = np.max(max_rupt_time)
    slip = solution["slip"]
    max_slip = [np.max(slip_seg.flatten()) for slip_seg in slip]
    max_slip = np.max(max_slip)
    x_label = "Distance along strike $(km)$"
    y_label = "Distance along dip $(km)$"
    for i_segment, (segment, rupt_time_seg, slip_seg, ps_seg) in enumerate(
        zip(segments, rupt_time, slip, point_sources)
    ):
        #
        # Plot the slip distribution
        #
        fig = plt.figure(figsize=(15, 8))
        ax = fig.add_subplot(111)
        ax.set_ylabel(y_label)
        ax.set_xlabel(x_label)
        fig.subplots_adjust(right=0.75)
        if i_segment == 0:
            ax.plot(0, 0, "w*", ms=20)
        ax, im = __several_axes(
            rupt_time_seg, segment, ps_seg, ax, max_val=max_rupt_time
        )
        cbar_ax = fig.add_axes((0.85, 0.15, 0.05, 0.7))
        cb = fig.colorbar(im, cax=cbar_ax)
        cb.set_label("Rupt_time (s)")
        plt.savefig(
            directory / "RuptTime_plane{}.png".format(i_segment), bbox_inches="tight"
        )
        plt.close()
    return


def _PlotRiseTime(
    segments: dict,
    point_sources: list,
    solution: dict,
    directory: Union[pathlib.Path, str] = pathlib.Path(),
):
    """Plot rise time distribution based on the FFM solution model

    :param segments: The segment properties
    :type segments: dict
    :param point_sources: The point source locations
    :type point_sources: list
    :param solution: The kinematic solution read from Solution.txt
    :type solution: dict
    :param directory: The location where to write plots, defaults to pathlib.Path()
    :type directory: Union[pathlib.Path, str], optional
    """
    directory = pathlib.Path(directory)
    print("Creating Rise Time Plot...")
    rise_time = solution["trise"]
    max_trise = [np.max(trise_seg.flatten()) for trise_seg in rise_time]
    max_trise = np.max(max_trise)
    fall_time = solution["tfall"]
    max_tfall = [np.max(tfall_seg.flatten()) for tfall_seg in fall_time]
    max_tfall = np.max(max_tfall)
    slip = solution["slip"]
    max_slip: Union[float, list] = [np.max(slip_seg.flatten()) for slip_seg in slip]
    max_slip = np.max(max_slip)
    x_label = "Distance along strike $(km)$"
    y_label = "Distance along dip $(km)$"
    for i_segment, (segment, trise_seg, tfall_seg, slip_seg, ps_seg) in enumerate(
        zip(segments, rise_time, fall_time, slip, point_sources)
    ):
        #
        # Plot the slip distribution
        #
        indexes = np.where(slip_seg < 0.1 * max_slip)  # type:ignore
        trise_seg[indexes] = 0
        tfall_seg[indexes] = 0
        fig, axes = plt.subplots(1, 2, figsize=(20, 10), sharex=True, sharey=True)
        ax1: Axes = axes[0]  # type: ignore
        ax2: Axes = axes[1]  # type: ignore
        fig.subplots_adjust(bottom=0.15)
        ax1.set_ylabel(y_label)
        ax1.set_xlabel(x_label)
        if i_segment == 0:
            ax1.plot(0, 0, "w*", ms=20)
        ax1, im = __several_axes(
            trise_seg, segment, ps_seg, ax1, max_val=max_trise, autosize=False
        )
        ax2.set_ylabel(y_label)
        ax2.set_xlabel(x_label)
        if i_segment == 0:
            ax2.plot(0, 0, "w*", ms=20)
        ax2, im = __several_axes(
            tfall_seg, segment, ps_seg, ax2, max_val=max_tfall, autosize=False
        )
        cbar_ax = fig.add_axes((0.1, 0.05, 0.8, 0.05))
        cb = fig.colorbar(im, cax=cbar_ax, orientation="horizontal")
        cb.set_label("Rise_time (s)")
        plt.savefig(
            directory / "RiseTime_plane{}.png".format(i_segment), bbox_inches="tight"
        )
        plt.close()
    return


def _PlotMultiSlipDist(
    segments: dict,
    point_sources: list,
    solution: dict,
    autosize: bool = False,
    max_val: Optional[float] = None,
    directory: Union[pathlib.Path, str] = pathlib.Path(),
):
    """Plot slip distribution based on the FFM solution model

    :param segments: The segment properties
    :type segments: dict
    :param point_sources: The point source locations
    :type point_sources: list
    :param solution: The kinematic solution read from Solution.txt
    :type solution: dict
    :param autosize: Automatically size the figure, defaults to False
    :type autosize: bool, optional
    :param max_val: Add a maximum value, defaults to None
    :type max_val: Optional[float], optional
    :param directory: The location where to write plots, defaults to pathlib.Path()
    :type directory: Union[pathlib.Path, str], optional
    """
    directory = pathlib.Path(directory)
    print("Creating Slip Distribution Plot...")
    slip = solution["slip"]
    rake = solution["rake"]
    rupt_time = solution["rupture_time"]
    max_slip: Union[float, list] = [np.max(slip_seg.flatten()) for slip_seg in slip]
    max_slip = np.max(max_slip)
    print("Max Slip of Solution: {:.2f} cm".format(max_slip))  # type:ignore
    if max_val != None:
        max_slip = max_val  # type:ignore
    print("Max Slip for Plotting: {:.2f} cm".format(max_slip))  # type:ignore
    x_label = "Distance Along Strike (km)"
    y_label = "Distance Along Dip (km)"
    for i_segment, (segment, slip_seg, rake_seg, rupttime_seg, ps_seg) in enumerate(
        zip(segments, slip, rake, rupt_time, point_sources)
    ):
        max_slip_seg = np.max(slip_seg.flatten())
        max_slip_seg = max_slip
        u = slip_seg * np.cos(rake_seg * np.pi / 180.0) / max_slip_seg
        v = slip_seg * np.sin(rake_seg * np.pi / 180.0) / max_slip_seg
        #
        # Plot the slip distribution
        #
        plt.rc("axes", titlesize=20)
        plt.rc("axes", labelsize=20)
        plt.rc("xtick", labelsize=16)
        plt.rc("ytick", labelsize=16)
        plt.rc("font", size=20)
        fig = plt.figure(figsize=(9, 7))
        ax = fig.add_subplot(111)
        ax.set_ylabel(y_label, fontsize=16)
        ax.set_xlabel(x_label, fontsize=16)
        ax.xaxis.tick_top()
        ax.xaxis.set_label_position("top")
        ax.spines["bottom"].set_linewidth(3)
        ax.spines["top"].set_linewidth(3)
        ax.spines["left"].set_linewidth(3)
        ax.spines["right"].set_linewidth(3)
        fig.subplots_adjust(right=0.85, top=0.85, bottom=0.3)
        (
            stk_subfaults,
            dip_subfaults,
            delta_strike,
            delta_dip,
            hyp_stk,
            hyp_dip,
        ) = pl_mng.__unpack_plane_data(segment)
        x = np.arange(stk_subfaults) * delta_strike - hyp_stk * delta_strike
        y = np.arange(dip_subfaults) * delta_dip - hyp_dip * delta_dip
        X = np.linspace(min(x), max(x), 20 * len(x))
        Y = np.linspace(min(y), max(y), 20 * len(y))
        z = (u**2 + v**2) ** 0.5
        xcols, yrows = np.meshgrid(x, y)
        XCOLS, YROWS = np.meshgrid(X, Y)
        orig_grid = np.transpose(np.array([xcols.flatten(), yrows.flatten()]))
        new_grid = np.transpose(np.array([XCOLS.flatten(), YROWS.flatten()]))
        grid_rupttime = griddata(
            orig_grid, rupttime_seg.flatten(), new_grid, method="linear"
        )
        grid_rupttime_reshape = grid_rupttime.reshape((np.shape(XCOLS)))
        contplot = ax.contour(
            XCOLS,
            YROWS,
            grid_rupttime_reshape,
            colors="0.75",
            linestyles="dashed",
            levels=range(10, 800, 10),
            linewidths=1.0,
        )
        plt.clabel(contplot, fmt="%.0f", inline=True, fontsize=14, colors="k")
        grid_z = griddata(orig_grid, z.flatten(), new_grid, method="linear")
        grid_z_reshape = grid_z.reshape((np.shape(XCOLS)))
        ax.quiver(x, y, u, v, scale=20.0, width=0.002, color="0.5", clip_on=False)
        ax.plot(0, 0, "w*", ms=15, markeredgewidth=1.5, markeredgecolor="k")
        ax, im = __several_axes(z, segment, ps_seg, ax, max_val=1.0, autosize=autosize)
        cbar_ax = fig.add_axes((0.125, 0.15, 0.5, 0.07))
        sm = plt.cm.ScalarMappable(
            cmap=slipcpt,
            norm=plt.Normalize(vmin=0.0, vmax=max_slip / 100.0),  # type:ignore
        )
        cb = fig.colorbar(sm, cax=cbar_ax, orientation="horizontal")
        cb.outline.set_linewidth(3)  # type:ignore
        cb.set_label("Slip (m)", fontsize=18)
        strk = segment["strike"]
        dip = segment["dip"]
        ax.text(
            -0.1,
            1.22,
            "Strike = " + str(int(strk)),
            fontsize=15,
            fontweight="bold",
            transform=ax.transAxes,
            va="top",
            ha="left",
        )
        ax.text(
            -0.1,
            1.13,
            "Dip = " + str(int(dip)),
            fontsize=15,
            fontweight="bold",
            transform=ax.transAxes,
            va="top",
            ha="left",
        )
        ax.text(
            0,
            -0.04,
            "Rupture Front Contours Plotted Every 10 s",
            fontsize=15,
            fontweight="bold",
            transform=ax.transAxes,
            va="top",
            ha="left",
        )
        plt.savefig(directory / "SlipDist_plane{}.png".format(i_segment), dpi=300)
        plt.savefig(directory / "SlipDist_plane{}.ps".format(i_segment))
        plt.close()
    return


def PlotSlipDistribution(
    segments: dict,
    point_sources: list,
    solution: dict,
    autosize: bool = False,
    max_val: Optional[float] = None,
    directory: Union[pathlib.Path, str] = pathlib.Path(),
):
    """Plot slip distribution based on the FFM solution model

    :param segments: The segment properties
    :type segments: dict
    :param point_sources: The point source locations
    :type point_sources: list
    :param solution: The kinematic solution read from Solution.txt
    :type solution: dict
    :param autosize: Automatically size the figure, defaults to False
    :type autosize: bool, optional
    :param max_val: Add a maximum value, defaults to None
    :type max_val: Optional[float], optional
    :param directory: The location where to write plots, defaults to pathlib.Path()
    :type directory: Union[pathlib.Path, str], optional
    """
    directory = pathlib.Path(directory)
    print("Creating Slip Distribution Plot...")
    slip = solution["slip"]
    rake = solution["rake"]
    rupt_time = solution["rupture_time"]
    max_slip = [np.max(slip_seg.flatten()) for slip_seg in slip]
    max_slip = np.max(max_slip)
    print("Max Slip of Solution: {:.2f} cm".format(max_slip))  # type:ignore
    if max_val != None:
        max_slip = max_val  # type:ignore
    print("Max Slip for Plotting: {:.2f} cm".format(max_slip))  # type:ignore
    x_label = "Distance Along Strike (km)"
    y_label = "Distance Along Dip (km)"
    for i_segment, (segment, slip_seg, rake_seg, rupttime_seg, ps_seg) in enumerate(
        zip(segments, slip, rake, rupt_time, point_sources)
    ):
        max_slip_seg = np.max(slip_seg.flatten())
        max_slip_seg = max_slip
        u = slip_seg * np.cos(rake_seg * np.pi / 180.0) / max_slip_seg
        v = slip_seg * np.sin(rake_seg * np.pi / 180.0) / max_slip_seg
        #
        # Plot the slip distribution
        #
        plt.rc("axes", titlesize=20)
        plt.rc("axes", labelsize=20)
        plt.rc("xtick", labelsize=16)
        plt.rc("ytick", labelsize=16)
        plt.rc("font", size=20)
        fig = plt.figure(figsize=(9, 7))
        ax = fig.add_subplot(111)
        ax.set_ylabel(y_label, fontsize=16)
        ax.set_xlabel(x_label, fontsize=16)
        ax.xaxis.tick_top()
        ax.xaxis.set_label_position("top")
        ax.spines["bottom"].set_linewidth(3)
        ax.spines["top"].set_linewidth(3)
        ax.spines["left"].set_linewidth(3)
        ax.spines["right"].set_linewidth(3)
        fig.subplots_adjust(right=0.85, top=0.85, bottom=0.3)
        (
            stk_subfaults,
            dip_subfaults,
            delta_strike,
            delta_dip,
            hyp_stk,
            hyp_dip,
        ) = pl_mng.__unpack_plane_data(segment)
        x = np.arange(stk_subfaults) * delta_strike - hyp_stk * delta_strike
        y = np.arange(dip_subfaults) * delta_dip - hyp_dip * delta_dip
        X = np.linspace(min(x), max(x), 20 * len(x))
        Y = np.linspace(min(y), max(y), 20 * len(y))
        z = (u**2 + v**2) ** 0.5
        xcols, yrows = np.meshgrid(x, y)
        XCOLS, YROWS = np.meshgrid(X, Y)
        orig_grid = np.transpose(np.array([xcols.flatten(), yrows.flatten()]))
        new_grid = np.transpose(np.array([XCOLS.flatten(), YROWS.flatten()]))
        grid_rupttime = griddata(
            orig_grid, rupttime_seg.flatten(), new_grid, method="linear"
        )
        grid_rupttime_reshape = grid_rupttime.reshape((np.shape(XCOLS)))
        contplot = ax.contour(
            XCOLS,
            YROWS,
            grid_rupttime_reshape,
            colors="0.75",
            linestyles="dashed",
            levels=range(10, 800, 10),
            linewidths=1.0,
        )
        plt.clabel(contplot, fmt="%.0f", inline=True, fontsize=14, colors="k")
        grid_z = griddata(orig_grid, z.flatten(), new_grid, method="linear")
        grid_z_reshape = grid_z.reshape((np.shape(XCOLS)))
        if (
            autosize == True
        ):  # make sure that rake vectors are same scale for each segment
            scale = (1.0 / 10) * (
                stk_subfaults * delta_strike
            )  # One-tenth of segment width
            width = 1 / (scale * 25)  # quiver width
            ax.quiver(x, y, u, v, scale=scale, width=width, color="0.5", clip_on=False)
        else:  # will be plotted to default aspect ratio
            ax.quiver(x, y, u, v, scale=20.0, width=0.002, color="0.5", clip_on=False)
        ax.plot(0, 0, "w*", ms=15, markeredgewidth=1.5, markeredgecolor="k")
        ax, im = __several_axes(
            grid_z_reshape, segment, ps_seg, ax, max_val=1.0, autosize=autosize
        )
        cbar_ax = fig.add_axes((0.125, 0.15, 0.5, 0.07))
        sm = plt.cm.ScalarMappable(
            cmap=slipcpt,
            norm=plt.Normalize(vmin=0.0, vmax=max_slip / 100.0),  # type:ignore
        )
        cb = fig.colorbar(sm, cax=cbar_ax, orientation="horizontal")
        cb.outline.set_linewidth(3)  # type:ignore
        cb.set_label("Slip (m)", fontsize=18)
        strk = segment["strike"]
        dip = segment["dip"]
        ax.text(
            -0.1,
            1.22,
            "Strike = " + str(int(strk)),
            fontsize=15,
            fontweight="bold",
            transform=ax.transAxes,
            va="top",
            ha="left",
        )
        ax.text(
            -0.1,
            1.13,
            "Dip = " + str(int(dip)),
            fontsize=15,
            fontweight="bold",
            transform=ax.transAxes,
            va="top",
            ha="left",
        )
        ax.text(
            0,
            -0.04,
            "Rupture Front Contours Plotted Every 10 s",
            fontsize=15,
            fontweight="bold",
            transform=ax.transAxes,
            va="top",
            ha="left",
        )
        plt.savefig(directory / "SlipDist_plane{}.png".format(i_segment), dpi=300)
        plt.savefig(directory / "SlipDist_plane{}.ps".format(i_segment))
        plt.close()
    return


def PlotSlipTimes(
    segments: dict,
    point_sources: list,
    solution: dict,
    directory: Union[pathlib.Path, str] = pathlib.Path(),
):
    """Plot slip distribution based on the FFM solution model

    :param segments: The segment properties
    :type segments: dict
    :param point_sources: The point source locations
    :type point_sources: list
    :param solution: The kinematic solution read from Solution.txt
    :type solution: dict
    :param directory: The location where to write plots, defaults to pathlib.Path()
    :type directory: Union[pathlib.Path, str], optional
    """
    directory = pathlib.Path(directory)
    print("Creating Slip Times Plot...")
    slip = solution["slip"]
    rake = solution["rake"]
    rupt_time = solution["rupture_time"]
    rise_time = solution["trise"]
    fall_time = solution["tfall"]
    max_slip = [np.max(slip_seg.flatten()) for slip_seg in slip]
    max_slip = np.max(max_slip)
    x_label = "Distance Along Strike (km)"
    y_label = "Dist Along Dip (km)"
    for i_segment, (
        segment,
        slip_seg,
        rake_seg,
        rupttime_seg,
        risetime_seg,
        falltime_seg,
        ps_seg,
    ) in enumerate(
        zip(segments, slip, rake, rupt_time, rise_time, fall_time, point_sources)
    ):
        max_slip_seg = np.max(slip_seg.flatten())
        u = slip_seg * np.cos(rake_seg * np.pi / 180.0) / max_slip_seg
        v = slip_seg * np.sin(rake_seg * np.pi / 180.0) / max_slip_seg
        #
        # Plot the slip distribution
        #
        plt.rc("axes", titlesize=16)
        plt.rc("axes", labelsize=16)
        plt.rc("xtick", labelsize=14)
        plt.rc("ytick", labelsize=14)
        plt.rc("font", size=16)
        fig = plt.figure(figsize=(8, 13))

        ax1 = fig.add_subplot(311)
        ax1.set_ylabel(y_label, fontsize=16)
        ax1.set_xlabel(x_label, fontsize=16)
        ax1.xaxis.tick_top()
        ax1.xaxis.set_label_position("top")
        ax1.spines["bottom"].set_linewidth(3)
        ax1.spines["top"].set_linewidth(3)
        ax1.spines["left"].set_linewidth(3)
        ax1.spines["right"].set_linewidth(3)
        fig.subplots_adjust(right=0.75, hspace=0.38)
        (
            stk_subfaults,
            dip_subfaults,
            delta_strike,
            delta_dip,
            hyp_stk,
            hyp_dip,
        ) = pl_mng.__unpack_plane_data(segment)
        x = np.arange(stk_subfaults) * delta_strike - hyp_stk * delta_strike
        y = np.arange(dip_subfaults) * delta_dip - hyp_dip * delta_dip
        X = np.linspace(min(x), max(x), 20 * len(x))
        Y = np.linspace(min(y), max(y), 20 * len(y))
        z = (u**2 + v**2) ** 0.5
        xcols, yrows = np.meshgrid(x, y)
        XCOLS, YROWS = np.meshgrid(X, Y)
        orig_grid = np.transpose(np.array([xcols.flatten(), yrows.flatten()]))
        new_grid = np.transpose(np.array([XCOLS.flatten(), YROWS.flatten()]))
        grid_rupttime = griddata(
            orig_grid, rupttime_seg.flatten(), new_grid, method="linear"
        )
        grid_rupttime_reshape = grid_rupttime.reshape((np.shape(XCOLS)))
        contplot = ax1.contour(
            XCOLS,
            YROWS,
            grid_rupttime_reshape,
            colors="0.75",
            linestyles="dashed",
            levels=range(10, 500, 10),
            linewidths=1.0,
        )
        ax1.quiver(x, y, u, v, scale=15.0, width=0.002, color="0.5", clip_on=False)
        if i_segment == 0:
            ax1.plot(0, 0, "w*", ms=15, markeredgewidth=1.5, markeredgecolor="k")
        ax1, im = __several_axes(
            rupttime_seg, segment, ps_seg, ax1, autosize=False, cmap="plasma"
        )
        plt.clabel(contplot, fmt="%.0f", inline=True, fontsize=11, colors="k")
        cbar_ax1 = fig.add_axes((0.85, 0.678, 0.03, 0.2))
        sm = plt.cm.ScalarMappable(
            cmap="plasma", norm=plt.Normalize(vmin=0.0, vmax=max(grid_rupttime))
        )
        cb = fig.colorbar(sm, cax=cbar_ax1, orientation="vertical")
        cb.outline.set_linewidth(3)  # type:ignore
        cb.set_label("Rupture Time (s)", fontsize=16)
        strk = segment["strike"]
        dip = segment["dip"]
        ax1.text(
            -0.1,
            1.25,
            "Strike = " + str(int(strk)),
            fontsize=13,
            fontweight="bold",
            transform=ax1.transAxes,
            va="top",
            ha="left",
        )
        ax1.text(
            -0.1,
            1.19,
            "Dip = " + str(int(dip)),
            fontsize=13,
            fontweight="bold",
            transform=ax1.transAxes,
            va="top",
            ha="left",
        )
        ax1.text(
            0,
            -0.04,
            "Rupture Front Contours Plotted Every 10 s",
            fontsize=13,
            fontweight="bold",
            transform=ax1.transAxes,
            va="top",
            ha="left",
        )

        ax2 = fig.add_subplot(312)
        ax2.set_ylabel(y_label, fontsize=16)
        ax2.set_xlabel(x_label, fontsize=16)
        ax2.xaxis.tick_top()
        ax2.xaxis.set_label_position("top")
        ax2.spines["bottom"].set_linewidth(3)
        ax2.spines["top"].set_linewidth(3)
        ax2.spines["left"].set_linewidth(3)
        ax2.spines["right"].set_linewidth(3)
        vrup_ref = segment["rupture_vel"]

        rupt_vel = rupttime_seg
        idx0 = np.where(rupt_vel < 0.1)
        rupt_vel[idx0] = 1.0
        diag_dist = np.sqrt(xcols**2 + yrows**2)
        rupt_vel = diag_dist / rupt_vel
        rupt_vel[idx0] = vrup_ref
        idx10 = np.where(slip_seg > 0.05 * max_slip)  # type:ignore
        mean_vrup = np.mean(rupt_vel[idx10].flatten())
        contplot = ax2.contour(
            XCOLS,
            YROWS,
            grid_rupttime_reshape,
            colors="0.9",
            linestyles="dashed",
            levels=range(10, 500, 10),
            linewidths=1.0,
        )
        ax2.quiver(x, y, u, v, scale=15.0, width=0.002, color="0.5", clip_on=False)
        if i_segment == 0:
            ax2.plot(0, 0, "w*", ms=15, markeredgewidth=1.5, markeredgecolor="k")
        ax2, im = __several_axes(
            rupt_vel,
            segment,
            ps_seg,
            ax2,
            autosize=False,
            cmap="viridis_r",
            min_val=min(rupt_vel.flatten()),
            max_val=max(rupt_vel.flatten()),
        )
        plt.clabel(contplot, fmt="%.0f", inline=True, fontsize=11, colors="k")
        sm = plt.cm.ScalarMappable(
            cmap="viridis_r", norm=plt.Normalize(vmin=0.0, vmax=max(rupt_vel.flatten()))
        )
        cbar_ax2 = fig.add_axes((0.85, 0.396, 0.03, 0.2))
        cb = fig.colorbar(sm, cax=cbar_ax2, orientation="vertical")
        cb.outline.set_linewidth(3)  # type:ignore
        cb.set_label("Rupture Velocity (km/s)", fontsize=16)
        ax2.text(
            0,
            -0.04,
            "Average* Rupture Velocity to Subfault: "
            + "{:.2f}".format(mean_vrup)
            + " km/s",
            fontsize=13,
            fontweight="bold",
            transform=ax2.transAxes,
            va="top",
            ha="left",
        )

        ax3 = fig.add_subplot(313)
        ax3.set_ylabel(y_label, fontsize=16)
        ax3.set_xlabel(x_label, fontsize=16)
        ax3.xaxis.tick_top()
        ax3.xaxis.set_label_position("top")
        ax3.spines["bottom"].set_linewidth(3)
        ax3.spines["top"].set_linewidth(3)
        ax3.spines["left"].set_linewidth(3)
        ax3.spines["right"].set_linewidth(3)
        contplot = ax3.contour(
            XCOLS,
            YROWS,
            grid_rupttime_reshape,
            colors="0.5",
            linestyles="dashed",
            levels=range(10, 500, 10),
            linewidths=1.0,
        )
        ax3.quiver(x, y, u, v, scale=15.0, width=0.002, color="0.5", clip_on=False)
        if i_segment == 0:
            ax3.plot(0, 0, "w*", ms=15, markeredgewidth=1.5, markeredgecolor="k")
        slip_duration = risetime_seg + falltime_seg
        ax3, im = __several_axes(
            slip_duration, segment, ps_seg, ax3, autosize=False, cmap="cividis_r"
        )
        plt.clabel(contplot, fmt="%.0f", inline=True, fontsize=11, colors="k")
        cbar_ax3 = fig.add_axes((0.85, 0.115, 0.03, 0.2))
        sm = plt.cm.ScalarMappable(
            cmap="cividis_r",
            norm=plt.Normalize(vmin=0.0, vmax=max(slip_duration.flatten())),
        )
        cb = fig.colorbar(sm, cax=cbar_ax3, orientation="vertical")
        cb.outline.set_linewidth(3)  # type:ignore
        cb.set_label("Slip Duration (s)", fontsize=16)
        mean_rise = np.mean(slip_duration[idx10])
        ax3.text(
            0,
            -0.04,
            "Average* Slip Duration (Rise Time): " + "{:.2f}".format(mean_rise) + " s",
            fontsize=13,
            fontweight="bold",
            transform=ax3.transAxes,
            va="top",
            ha="left",
        )
        ax3.text(
            0,
            -0.12,
            "*Includes subfaults with >5% maximum slip",
            fontsize=13,
            fontweight="bold",
            transform=ax3.transAxes,
            va="top",
            ha="left",
        )
        plt.savefig(directory / "SlipTime_plane{}.png".format(i_segment), dpi=300)
        plt.close()
    return


def _PlotCumulativeSlip(
    segments: dict,
    point_sources: list,
    solution: dict,
    tensor_info: dict,
    evID: Optional[str] = None,
    directory: Union[pathlib.Path, str] = pathlib.Path(),
) -> Tuple[str, str, str, str]:
    """Plot the cumulative slip

    :param segments: The segment properties
    :type segments: dict
    :param point_sources: The point source locations
    :type point_sources: list
    :param solution: The kinematic solution read from Solution.txt
    :type solution: dict
    :param tensor_info: The tensor information
    :type tensor_info: dict
    :param evID: The event name/ID, defaults to None
    :type evID: Optional[str], optional
    :param directory: The location where plots should be written, defaults to pathlib.Path()
    :type directory: Union[pathlib.Path, str], optional
    :return: The location of corners
    :rtype: Tuple[str, str, str, str]
    """
    directory = pathlib.Path(directory)
    print("Creating Along Strike/Dip Plot...")
    slip = solution["slip"]
    max_slip = [np.max(slip_seg.flatten()) for slip_seg in slip]
    max_slip = np.max(max_slip)
    x_label = "Distance Along Strike (km)"
    y_label = "Dist Along Dip (km)"
    for i_segment, (segment, slip_seg, ps_seg) in enumerate(
        zip(segments, slip, point_sources)
    ):
        max_slip_seg = np.max(slip_seg.flatten())
        #
        # Plot the slip distribution
        #
        plt.rc("axes", titlesize=11)
        plt.rc("axes", labelsize=11)
        plt.rc("xtick", labelsize=11)
        plt.rc("ytick", labelsize=11)
        plt.rc("font", size=11)
        fig = plt.figure(figsize=(8, 8))

        (
            stk_subfaults,
            dip_subfaults,
            delta_strike,
            delta_dip,
            hyp_stk,
            hyp_dip,
        ) = pl_mng.__unpack_plane_data(segment)
        x = np.arange(stk_subfaults) * delta_strike - hyp_stk * delta_strike
        y = np.arange(dip_subfaults) * delta_dip - hyp_dip * delta_dip

        ### PLOT IT ###
        grid = plt.GridSpec(4, 4, hspace=0.2, wspace=0.2)
        main_ax = fig.add_subplot(grid[1:3, 1:4])
        main_ax.spines["bottom"].set_linewidth(3)
        main_ax.spines["top"].set_linewidth(3)
        main_ax.spines["left"].set_linewidth(3)
        main_ax.spines["right"].set_linewidth(3)
        if i_segment == 0:
            main_ax.plot(0, 0, "w*", ms=15, markeredgewidth=1.5, markeredgecolor="k")
        main_ax, im = __several_axes(
            slip_seg,
            segment,
            ps_seg,
            main_ax,
            max_val=max_slip,
            autosize=False,
            depth_ax=False,
        )

        ### Slip_threshold to ignore low-slip region ###
        slip_threshold = 0.1 * (max_slip / 100.0)  # type:ignore
        if slip_threshold < 1.0:
            slip_threshold = 1.0
        if max_slip / 100.0 < 3.0:  # type:ignore
            slip_threshold = 0.5
        if max_slip / 100.0 < 1.0:  # type:ignore
            slip_threshold = 0.2
        slip_seg0 = slip_seg
        cumslip_seg = slip_seg
        idxslip = np.where(cumslip_seg < 100.0 * slip_threshold)
        cumslip_seg[idxslip] = 0.0
        along_strike = np.sum(cumslip_seg, axis=0)
        along_dip = np.sum(cumslip_seg, axis=1)
        ### Calculate equivalent slip length in along-strike and along-dip directions ###
        eq_len_AS = shakemap.equivalent_slip_length(along_strike / 100, delta_strike)
        left_edge_AS = shakemap.locate_equivalent_slip(
            along_strike / 100, delta_strike, eq_len_AS  # type:ignore
        )
        left_edge_AS = left_edge_AS + x[0]
        eq_len_AD = shakemap.equivalent_slip_length(
            along_dip / 100, delta_dip  # type:ignore
        )
        left_edge_AD = shakemap.locate_equivalent_slip(
            along_dip / 100, delta_dip, eq_len_AD
        )
        left_edge_AD = left_edge_AD + y[0]

        main_ax.hlines(
            left_edge_AD, left_edge_AS, left_edge_AS + eq_len_AS, "k", ls="--"
        )
        main_ax.hlines(
            left_edge_AD + eq_len_AD,
            left_edge_AS,
            left_edge_AS + eq_len_AS,
            "k",
            ls="--",
        )
        main_ax.vlines(
            left_edge_AS, left_edge_AD, left_edge_AD + eq_len_AD, "k", ls="--"
        )
        main_ax.vlines(
            left_edge_AS + eq_len_AS,
            left_edge_AD,
            left_edge_AD + eq_len_AD,
            "k",
            ls="--",
        )
        main_ax.set_xlabel(x_label, fontsize=11)
        main_ax.set_ylabel(y_label, fontsize=11)
        main_ax.xaxis.tick_bottom()
        main_ax.xaxis.set_label_position("bottom")
        main_ax.yaxis.tick_right()
        main_ax.yaxis.set_label_position("right")
        main_ax.text(
            1,
            1.05,
            "*Fault Dimensions Not To Scale",
            fontsize=9,
            fontweight="bold",
            transform=main_ax.transAxes,
            va="top",
            ha="right",
        )

        ax1 = fig.add_subplot(grid[0, 1:4], sharex=main_ax)
        ax1.xaxis.tick_top()
        ax1.xaxis.set_label_position("top")
        ax1.spines["bottom"].set_linewidth(3)
        ax1.spines["top"].set_linewidth(3)
        ax1.spines["left"].set_linewidth(3)
        ax1.spines["right"].set_linewidth(3)
        ax1.plot(x, along_strike / 100, "r", lw=2)
        ax1.scatter(x, along_strike / 100, c="0.5", s=25, edgecolor="k", lw=1, zorder=5)
        ax1.set_xlabel("Distance Along Strike (km)", fontsize=11)
        ax1.set_ylabel("Cumulative Slip (m)", fontsize=11)
        ax1.invert_yaxis()
        ax1.yaxis.tick_right()
        ax1.yaxis.set_label_position("right")
        ax1.set_ylim(0, max(along_strike / 100) + 0.1 * max(along_strike / 100))
        ax1.vlines(
            left_edge_AS,
            0,
            max(along_strike / 100) + 0.1 * max(along_strike / 100),
            "k",
            ls="--",
        )
        ax1.vlines(
            left_edge_AS + eq_len_AS,
            0,
            max(along_strike / 100) + 0.1 * max(along_strike / 100),
            "k",
            ls="--",
        )

        ax2 = fig.add_subplot(grid[1:3, 0], sharey=main_ax)
        ax2.yaxis.tick_left()
        ax2.xaxis.set_label_position("bottom")
        ax2.yaxis.set_label_position("left")
        ax2.spines["bottom"].set_linewidth(3)
        ax2.spines["top"].set_linewidth(3)
        ax2.spines["left"].set_linewidth(3)
        ax2.spines["right"].set_linewidth(3)
        ax2.plot(along_dip / 100, y, "r", lw=2)
        ax2.scatter(along_dip / 100, y, c="0.5", s=25, edgecolor="k", lw=1, zorder=5)
        ax2.set_ylabel("Distance Along Dip (km)", fontsize=11)  # , rotation=270)
        ax2.set_xlabel("Cumulative Slip (m)", fontsize=11)
        ax2.set_xlim(0, max(along_dip / 100) + 0.1 * max(along_dip / 100))
        ax2.hlines(
            left_edge_AD,
            0,
            max(along_dip / 100) + 0.1 * max(along_dip / 100),
            "k",
            ls="--",
        )
        ax2.hlines(
            left_edge_AD + eq_len_AD,
            0,
            max(along_dip / 100) + 0.1 * max(along_dip / 100),
            "k",
            ls="--",
        )
        ax2.invert_xaxis()
        ###LEGEND###
        cbar_ax = fig.add_axes((0.124, 0.712, 0.169, 0.025))
        sm = plt.cm.ScalarMappable(
            cmap=slipcpt,
            norm=plt.Normalize(vmin=0.0, vmax=max_slip / 100.0),  # type:ignore
        )
        cb = fig.colorbar(sm, cax=cbar_ax, orientation="horizontal")
        cbar_ax.xaxis.set_ticks_position("top")
        cbar_ax.xaxis.set_label_position("top")
        cb.outline.set_linewidth(3)  # type:ignore
        cb.set_label("Slip (m)", fontsize=10, fontweight="bold")
        strk = segment["strike"]
        dip = segment["dip"]
        ax1.text(
            -0.2,
            0.85,
            "Strike = " + str(int(strk)),
            fontsize=11,
            transform=ax1.transAxes,
            va="top",
            ha="center",
        )
        ax1.text(
            -0.2,
            0.7,
            "Dip = " + str(int(dip)),
            fontsize=11,
            transform=ax1.transAxes,
            va="top",
            ha="center",
        )
        if evID != None:
            ax1.text(
                -0.2,
                1.1,
                evID,  # type:ignore
                fontsize=14,
                fontweight="bold",
                transform=ax1.transAxes,
                va="top",
                ha="center",
            )

        plt.savefig(directory / "CumulativeSlip_plane{}.png".format(i_segment), dpi=300)
        plt.savefig(directory / "CumulativeSlip_plane{}.ps".format(i_segment))
        plt.close()

        ### recover edge lat/lon/dep from _AS and _AD info: ###
        hyp_lon = tensor_info["lon"]
        hyp_lat = tensor_info["lat"]
        hyp_dep = tensor_info["depth"]
        corner_1, corner_2, corner_3, corner_4 = shakemap.translate_xy_to_latlondep(
            segment,
            hyp_lon,
            hyp_lat,
            hyp_dep,
            eq_len_AS,
            eq_len_AD,
            left_edge_AS,
            left_edge_AD,
        )
    return corner_1, corner_2, corner_3, corner_4


def PlotSlipDist_Compare(
    segments: dict,
    point_sources: list,
    input_model: dict,
    solution: dict,
    max_val: Optional[float] = None,
    directory: Union[pathlib.Path, str] = pathlib.Path(),
):
    """Plot slip distribution based on the FFM solution model

    :param segments: The segment properties
    :type segments: dict
    :param point_sources: The point source locations
    :type point_sources: list
    :param input_model: The model
    :type input_model: dict
    :param solution: The kinematic solution read from Solution.txt
    :type solution: dict
    :param max_val: Specify a max value, defaults to None
    :type max_val: Optional[float], optional
    :param directory: The location where to write plots, defaults to pathlib.Path()
    :type directory: Union[pathlib.Path,str], optional
    """
    directory = pathlib.Path(directory)
    print("Creating Checkerboard Comparison Plot")
    slip = solution["slip"]
    rake = solution["rake"]
    slip2 = input_model["slip"]
    rake2 = input_model["rake"]
    max_slip = [np.max(slip_seg.flatten()) for slip_seg in slip]
    max_slip = np.max(max_slip)
    max_slip2 = [np.max(slip_seg2.flatten()) for slip_seg2 in slip2]
    max_slip2 = np.max(max_slip2)
    max_slip = np.maximum(max_slip, max_slip2)  # type:ignore
    if max_val == None:
        max_val = max_slip  # type:ignore
    x_label = "Distance along strike $(km)$"
    y_label = "Distance along dip $(km)$"
    zipped = zip(segments, slip, rake, slip2, rake2, point_sources)
    for i_segment, (
        segment,
        slip_seg,
        rake_seg,
        slip_seg2,
        rake_seg2,
        ps_seg,
    ) in enumerate(zipped):
        max_slip_seg = np.max(slip_seg.flatten())
        max_slip_seg2 = np.max(slip_seg2.flatten())
        max_slip_seg = np.maximum(max_slip_seg, max_slip_seg2)
        u = slip_seg * np.cos(rake_seg * np.pi / 180.0) / max_slip_seg
        v = slip_seg * np.sin(rake_seg * np.pi / 180.0) / max_slip_seg
        u2 = slip_seg2 * np.cos(rake_seg2 * np.pi / 180.0) / max_slip_seg
        v2 = slip_seg2 * np.sin(rake_seg2 * np.pi / 180.0) / max_slip_seg
        #
        # Plot the slip distribution
        #
        fig, axes = plt.subplots(1, 2, figsize=(30, 8))
        ax0: Axes = axes[0]  # type: ignore
        ax1: Axes = axes[1]  # type: ignore
        ax0.set_ylabel(y_label)
        ax0.set_xlabel(x_label)
        ax1.set_ylabel(y_label)
        ax1.set_xlabel(x_label)
        fig.subplots_adjust(right=0.75)
        (
            stk_subfaults,
            dip_subfaults,
            delta_strike,
            delta_dip,
            hyp_stk,
            hyp_dip,
        ) = pl_mng.__unpack_plane_data(segment)
        x = np.arange(stk_subfaults) * delta_strike - hyp_stk * delta_strike
        y = np.arange(dip_subfaults) * delta_dip - hyp_dip * delta_dip
        ax0.quiver(x, y, u, v, scale=15.0, width=0.003)
        ax1.quiver(x, y, u2, v2, scale=15.0, width=0.003)
        if i_segment == 0:
            ax0.plot(0, 0, "w*", ms=20)
            ax1.plot(0, 0, "w*", ms=20)
        ax0, im = __several_axes(
            slip_seg, segment, ps_seg, ax0, max_val=max_val, autosize=False
        )
        ax1, im = __several_axes(
            slip_seg2, segment, ps_seg, ax1, max_val=max_val, autosize=False
        )
        ax0.set_title("Inverted model", fontsize=20)
        ax1.set_title("Original model", fontsize=20)
        cbar_ax = fig.add_axes((0.85, 0.15, 0.05, 0.7))
        cb = fig.colorbar(im, cax=cbar_ax)
        cb.set_label("Slip (cm)")
        plt.savefig(
            directory / "Checkerboard_SlipDist_plane{}.png".format(i_segment),
            bbox_inches="tight",
        )
        plt.close()
    return


def PlotMap(
    tensor_info: dict,
    segments: List[dict],
    point_sources: list,
    solution: dict,
    default_dirs: dict,
    files_str: Optional[dict] = None,
    stations_gps: Optional[zip] = None,
    stations_cgps: Optional[str] = None,
    max_slip: Optional[float] = None,
    legend_len: Optional[float] = None,
    scale: Optional[float] = None,
    limits: List[Optional[float]] = [None, None, None, None],
    label_stations: bool = False,
    directory: Union[pathlib.Path, str] = pathlib.Path(),
):
    """Plot map

    :param tensor_info: The tensor information
    :type tensor_info: dict
    :param segments: The segment properties
    :type segments: List[dict]
    :param point_sources: The point source locations
    :type point_sources: list
    :param solution: The kinematic solution read from Solution.txt
    :type solution: dict
    :param default_dirs: The location of default directories
    :type default_dirs: dict
    :param files_str: The stations file properties, defaults to None
    :type files_str:Optional[dict], optional
    :param stations_gps: The gps stations description, defaults to None
    :type stations_gps:  Optional[zip], optional
    :param stations_cgps: The cgps stations description, defaults to None
    :type stations_cgps: Optional[str], optional
    :param max_slip: Specify maximum slip, defaults to None
    :type max_slip: Optional[float], optional
    :param legend_len: The length of the legend, defaults to None
    :type legend_len: Optional[int], optional
    :param scale: The scale, defaults to None
    :type scale: Optional[float], optional
    :param limits: The extent of the map, defaults to [None, None, None, None]
    :type limits: List[Optional[float]], optional
    :param label_stations: Label the stations, defaults to False
    :type label_stations: bool, optional
    :param directory: the location where to write plots, defaults to pathlib.Path()
    :type directory: Union[pathlib.Path, str], optional
    """
    directory = pathlib.Path(directory)
    g = Geod(ellps="WGS84")
    ################################
    ### GET DESIRED PLOT REGION ####
    ################################

    lat0 = tensor_info["lat"]
    lon0 = tensor_info["lon"]

    segments_lats, segments_lons, segments_deps = __redefine_lat_lon(
        segments, point_sources
    )
    min_lats = [min(segment_lat.flatten()) for segment_lat in segments_lats]
    max_lats = [max(segment_lat.flatten()) for segment_lat in segments_lats]
    min_lons = [min(segment_lon.flatten()) for segment_lon in segments_lons]
    max_lons = [max(segment_lon.flatten()) for segment_lon in segments_lons]
    min_lat = np.min(min_lats)
    max_lat = np.max(max_lats)
    min_lon = np.min(min_lons)
    max_lon = np.max(max_lons)

    if files_str is not None:
        for file in files_str:
            name = file["name"]
            latp, lonp = file["location"]
            min_lat = min(min_lat, latp)
            max_lat = max(max_lat, latp)
            min_lon = min(min_lon, lonp)
            max_lon = max(max_lon, lonp)
    if stations_cgps is not None:
        for file in stations_cgps:
            name = file["name"]
            latp, lonp = file["location"]
            min_lat = min(min_lat, latp)
            max_lat = max(max_lat, latp)
            min_lon = min(min_lon, lonp)
            max_lon = max(max_lon, lonp)
    if stations_gps is not None:
        max_obs = np.zeros(3)
        stations_gps2: list = []
        for name, sta_lat, sta_lon, obs, syn, error in stations_gps:
            min_lat = min(min_lat, sta_lat)
            max_lat = max(max_lat, sta_lat)
            min_lon = min(min_lon, sta_lon)
            max_lon = max(max_lon, sta_lon)
            stations_gps2 = stations_gps2 + [[name, sta_lat, sta_lon, obs, syn, error]]
            max_obs = np.maximum([abs(float(v)) for v in obs], max_obs)
        max_obs = np.max(max_obs)  # type:ignore
        if legend_len == None:
            if max_obs < 5:
                legend_len = 1
            elif max_obs < 10:
                legend_len = 5
            elif max_obs < 20:
                legend_len = 10
            elif max_obs < 50:
                legend_len = 20
            elif max_obs < 100:
                legend_len = 50
            else:
                legend_len = 100
        if scale == None:
            scale = 2
        max_obs = max_obs / scale  # type:ignore

    if limits == [None, None, None, None]:
        region = [min_lon - 0.5, max_lon + 0.5, min_lat - 0.5, max_lat + 0.5]
    else:
        region = [limits[0], limits[1], limits[2], limits[3]]
    ### Fix region to nearest tenth of a degree ###
    region[0] = np.floor(region[0] * 10) / 10.0
    region[1] = np.ceil(region[1] * 10) / 10.0
    region[2] = np.floor(region[2] * 10) / 10.0
    region[3] = np.ceil(region[3] * 10) / 10.0

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

    pygmt.config(
        PS_MEDIA="A0",
        MAP_FRAME_TYPE="plain",
        MAP_FRAME_AXES="WseN",
        FORMAT_GEO_OUT="F",
        FORMAT_GEO_MAP="ddd:mm:ss",
        FONT_ANNOT_PRIMARY="10p,Helvetica,black",
        FONT_ANNOT_SECONDARY="6p,Helvetica,black",
        FONT_HEADING="30p,Helvetica,black",
        FONT_LABEL="10p,Helvetica,black",
        FONT_LOGO="6p,Helvetica,black",
        FONT_SUBTITLE="16p,Helvetica,black",
        FONT_TAG="18p,Helvetica,black",
        FONT_TITLE="22p,Helvetica,black",
        MAP_ANNOT_OFFSET_PRIMARY="3p",
    )
    fig.basemap(region=region, projection=projection, frame=["WseN", "afg"])
    fig.grdimage(grid=grid, cmap="oleron", shading=True, transparency=20)
    fig.shift_origin(xshift="-0.45c", yshift="0.75c")
    fig.colorbar(
        position="n0.05/-0.1+jBL+w3c/8%+h",
        frame="x+lElevation (km)",
        scale=0.001,
    )
    fig.shift_origin(xshift="0.45c", yshift="-0.75c")
    fig.plot(
        str(default_dirs["root_dir"]) + "/pb2002_boundaries.gmt",
        style="f10/3p",
        region=region,
        pen="2p,black",
    )
    fig.plot(
        str(default_dirs["root_dir"]) + "/pb2002_boundaries.gmt",
        style="f10/3p",
        region=region,
        pen="1p,white",
    )
    fig.basemap(region=region, projection=projection, frame="ag1")
    fig.coast(resolution="h", shorelines=True)

    ###############################
    ### PLOT FINITE FAULT MODEL ###
    ###############################
    with open(directory / "segments_data.json") as s:
        segments_data = json.load(s)
    segments = segments_data["segments"]
    solution = get_outputs.read_solution_static_format(segments, data_dir=directory)
    plane_info = segments[0]
    (
        stk_subfaults,
        dip_subfaults,
        delta_strike,
        delta_dip,
        hyp_stk,
        hyp_dip,
    ) = pl_mng.__unpack_plane_data(plane_info)

    segments_lats, segments_lons, segments_deps = __redefine_lat_lon(
        segments, point_sources
    )
    slip = solution["slip"]

    ### MAKE CPT OF SLIP ###
    if max_slip == None:
        maxslip = 0
        for segment in range(len(segments_lats)):
            slips = slip[segment].flatten()
            maxslip = np.max([maxslip, np.array(slips).max() / 100])
    else:
        maxslip = max_slip  # type:ignore
    pygmt.makecpt(
        cmap=str(default_dirs["root_dir"]) + "/src/wasp/fault2.cpt", series=[0, maxslip]
    )

    plane_info = segments[0]
    strike = plane_info["strike"]
    dip = plane_info["dip"]

    rupture_vel = segments[0]["rupture_vel"]
    subfaults = {"delta_strike": delta_strike, "delta_dip": delta_dip}
    rise_time = segments_data["rise_time"]
    subfaults2 = pf._point_sources_def(rise_time, rupture_vel, subfaults)
    strike_ps = int(subfaults2["strike_ps"] / 2)
    dip_ps = int(subfaults2["dip_ps"] / 2)
    latitudes = [ps_segment[:, :, dip_ps, strike_ps, 0] for ps_segment in point_sources]
    longitudes = [
        ps_segment[:, :, dip_ps, strike_ps, 1] for ps_segment in point_sources
    ]

    for segment in range(len(segments_lats)):
        plane_info = segments[segment]
        strike = plane_info["strike"]
        dip = plane_info["dip"]
        lons = longitudes[segment].flatten()
        lats = latitudes[segment].flatten()
        slips = slip[segment].flatten()
        fig.plot(
            x=lons,
            y=lats,
            style="J"
            + str(strike)
            + "/"
            + str(delta_strike)
            + "/"
            + str(delta_dip * np.cos(np.radians(dip))),
            cmap=True,
            fill=slips / 100.0,
        )
    fig.shift_origin(xshift="3c", yshift="0.75c")
    fig.colorbar(
        position="n0.05/-0.1+jBL+w3c/8%+h",
        frame="x+lSlip (m)",
    )
    fig.shift_origin(xshift="-3c", yshift="-0.75c")
    # Just in case fault plane is going under coastline, plot that coast again #
    fig.coast(resolution="h", shorelines=True)

    #################################
    ### PLOT AFTERSHOCKS OVER TOP ###
    #################################
    aftershocks = glob.glob(str(directory) + "/*aftershock*")
    if len(aftershocks) > 0:
        for kafter in range(len(aftershocks)):
            print("...Adding aftershocks from: " + str(aftershocks[kafter]))
            aftershock = np.genfromtxt(
                aftershocks[kafter], delimiter="\t", skip_header=1
            )
            aftershock_lat = aftershock[:, 3]
            aftershock_lon = aftershock[:, 4]
            aftershock_mag = aftershock[:, 6]
            fig.plot(
                x=aftershock_lon,
                y=aftershock_lat,
                style="cp",
                size=aftershock_mag,
                fill="grey",
                pen="black",
                transparency=60,
            )

    #################################
    # outline the segments in black #
    #################################
    for segment in range(len(segments_lats)):
        plane_info = segments[segment]
        (
            stk_subfaults,
            dip_subfaults,
            delta_strike,
            delta_dip,
            hyp_stk,
            hyp_dip,
        ) = pl_mng.__unpack_plane_data(plane_info)
        strike = plane_info["strike"]
        dip = plane_info["dip"]
        depths = segments_deps[segment].flatten()
        min_lats_idx = np.where(segments_lats[segment].flatten() == min_lats[segment])[
            0
        ][0]
        cornerA = [
            segments_lons[segment].flatten()[min_lats_idx],
            min_lats[segment],
            depths[min_lats_idx],
            0,
        ]
        min_lons_idx = np.where(segments_lons[segment].flatten() == min_lons[segment])[
            0
        ][0]
        cornerB = [
            min_lons[segment],
            segments_lats[segment].flatten()[min_lons_idx],
            depths[min_lons_idx],
            1,
        ]
        max_lats_idx = np.where(segments_lats[segment].flatten() == max_lats[segment])[
            0
        ][-1]
        cornerC = [
            segments_lons[segment].flatten()[max_lats_idx],
            max_lats[segment],
            depths[max_lats_idx],
            2,
        ]
        max_lons_idx = np.where(segments_lons[segment].flatten() == max_lons[segment])[
            0
        ][-1]
        cornerD = [
            max_lons[segment],
            segments_lats[segment].flatten()[max_lons_idx],
            depths[max_lons_idx],
            3,
        ]
        # find outside corner of the subfault (rather than middle)
        angle = np.degrees(
            np.arctan((delta_dip * np.cos(np.radians(dip))) / delta_strike)
        )
        corner_offset = 1000 * (
            (
                (0.5 * delta_dip * np.cos(np.radians(dip))) ** 2
                + (0.5 * delta_strike) ** 2
            )
            ** 0.5
        )
        letter = ["A", "B", "C", "D"]
        corners = np.c_[cornerA, cornerB, cornerC, cornerD]
        max_dep = max(corners[2, :])
        idx_max = np.where(np.array(corners[2, :]) == max_dep)[0][0]
        min1 = min2 = corners[:, idx_max]
        for i in range(len(corners)):
            if corners[2, i] < min1[2]:
                min2 = min1
                min1 = corners[:, i]
            elif corners[2, i] < min2[2] and collections.Counter(
                corners[:, i]
            ) != collections.Counter(min1):
                min2 = corners[:, i]
        updip = np.c_[min1, min2]
        corners = np.c_[corners, cornerA]
        fig.plot(x=corners[0, :], y=corners[1, :], pen="1p,black")
        fig.plot(x=updip[0, :], y=updip[1, :], pen="1p,red")

        ######################
        # plot hypocenter(s) #
        ######################
        fig.plot(x=lon0, y=lat0, style="a7p", fill="white", pen="black")

    ###########################
    ### PLOT LOCAL STATIONS ###
    ###########################

    ### STRONG-MOTION ACCELEROMETER ###
    if files_str is not None:
        print("...Plotting Strong Motion Stations")
        for file in files_str:
            name = file["name"]
            latp, lonp = file["location"]
            weight = file["trace_weight"]
            comp = file["component"]
            if comp[-1] == "E":
                if weight == 0:
                    fig.plot(x=lonp, y=latp, style="i7p", fill="lightgrey", pen="black")
                else:
                    fig.plot(x=lonp, y=latp, style="i10p", fill="white", pen="black")
                if label_stations == True:
                    fig.text(x=lonp, y=latp, text=name, justify="ML", font="7p")
        ### ADD TO LEGEND ###
        fig.shift_origin(xshift="-2.7c", yshift="-2.3c")
        fig.plot(
            x=region[1],
            y=region[2],
            no_clip=True,
            style="i10p",
            fill="white",
            pen="black",
        )
        fig.shift_origin(xshift="2.7c", yshift="2.3c")
        fig.shift_origin(xshift="-2.5c", yshift="-2.3c")
        fig.text(
            x=region[1],
            y=region[2],
            text="Accelerometer",
            no_clip=True,
            justify="ML",
        )
        fig.shift_origin(xshift="2.5c", yshift="2.3c")

    ### HIGH-RATE GNSS ###
    if stations_cgps is not None:
        print("...Plotting cGPS Stations")
        for file in stations_cgps:
            name = file["name"]
            latp, lonp = file["location"]
            weight = file["trace_weight"]
            comp = file["component"]
            if comp[-1] == "E":
                if weight == 0:
                    fig.plot(
                        x=lonp, y=latp, style="t7p", fill="lightsteelblue", pen="black"
                    )
                else:
                    fig.plot(x=lonp, y=latp, style="t10p", fill="navy", pen="black")
                if label_stations == True:
                    fig.text(x=lonp, y=latp, text=name, justify="TR", font="7p")
        ### ADD TO LEGEND ###
        fig.shift_origin(xshift="0.2c", yshift="-2.3c")
        fig.plot(
            x=region[1],
            y=region[2],
            no_clip=True,
            style="t10p",
            fill="navy",
            pen="black",
        )
        fig.shift_origin(xshift="-0.2c", yshift="2.3c")
        fig.shift_origin(xshift="0.4c", yshift="-2.3c")
        fig.text(
            x=region[1],
            y=region[2],
            text="HR GNSS",
            no_clip=True,
            justify="ML",
            offset=0 / 10,
        )
        fig.shift_origin(xshift="-0.4c", yshift="2.3c")

    ### STATIC GNSS ####
    if stations_gps is not None:
        print("...Plotting Static GPS Stations")
        for name, sta_lat, sta_lon, obs, syn, error in stations_gps2:
            gps_z_obs, gps_n_obs, gps_e_obs = obs
            gps_z_syn, gps_n_syn, gps_e_syn = syn
            err_z, err_n, err_e = error
            if label_stations == True:
                fig.text(x=sta_lon, y=sta_lat, text=name, justify="TR", font="7p")
            staticv_obs = pd.DataFrame(
                data={
                    "x": [sta_lon],
                    "y": [sta_lat],
                    "east_velocity": [float(gps_e_obs) / max_obs],
                    "north_velocity": [float(gps_n_obs) / max_obs],
                    "east_sigma": [float(err_e) / max_obs],
                    "north_sigma": [float(err_n) / max_obs],
                    "correlation_EN": [0.0],
                }
            )
            v_obs = "0.3c+a45+p0.5p,grey+e+h0.5+gblack"
            # Plot thick white arrow behind, to get white outline on black arrow
            fig.velo(
                data=staticv_obs,
                pen="0.07c,grey",
                line="grey",
                fill="grey",
                spec="e1/0",
                vector=v_obs,
            )
            fig.velo(
                data=staticv_obs,
                pen="0.05c,black",
                line="BLACK",
                fill="BLACK",
                spec="e1/0.34",
                vector=v_obs,
            )

            staticv_syn = pd.DataFrame(
                data={
                    "x": [sta_lon],
                    "y": [sta_lat],
                    "east_velocity": [float(gps_e_syn) / max_obs],
                    "north_velocity": [float(gps_n_syn) / max_obs],
                    "east_sigma": [0.0],
                    "north_sigma": [0.0],
                    "correlation_EN": [0.0],
                }
            )
            v_syn = "0.3c+p0.5p,black+e+h0.5+gred"
            # Plot thick black arrow behind, to get black outline on red arrow
            fig.velo(
                data=staticv_syn,
                pen=".05c,black",
                line="black",
                fill="black",
                spec="e1/0",
                vector=v_syn,
            )
            # overlay thin red arrow
            fig.velo(
                data=staticv_syn,
                pen=".03c,red",
                line="red",
                fill="red",
                spec="e1/0",
                vector=v_syn,
            )
        ### ADD TO LEGEND ###
        fig.shift_origin(xshift="-2.9c", yshift="-1.45c")
        fig.text(
            x=region[1],
            y=region[2],
            text="Observed GNSS",
            font="10p,Helvetica,black",
            no_clip=True,
            justify="ML",
        )
        fig.shift_origin(xshift="2.9c", yshift="1.45c")
        fig.shift_origin(xshift="-2.9c", yshift="-1.8c")
        fig.text(
            x=region[1],
            y=region[2],
            text="Synthetic  GNSS",
            font="10p,Helvetica,black",
            no_clip=True,
            justify="ML",
        )
        fig.shift_origin(xshift="2.9c", yshift="1.8c")
        static_legend = pd.DataFrame(
            data={
                "x": [region[1]],
                "y": [region[2]],
                "east_velocity": [legend_len / max_obs],  # type:ignore
                "north_velocity": [0],
                "east_sigma": [legend_len / max_obs / 10],  # type:ignore
                "north_sigma": [legend_len / max_obs / 10],  # type:ignore
                "correlation_EN": [0],
            }
        )
        # Plot thick white arrow behind, to get white outline on black arrow
        fig.shift_origin(xshift="0c", yshift="-1.45c")
        fig.velo(
            data=static_legend,
            pen="0.07c,grey",
            line="grey",
            fill="grey",
            spec="e1/0",
            vector=v_obs,
            no_clip=True,
        )
        fig.velo(
            data=static_legend,
            pen="0.05c,black",
            line="BLACK",
            fill="BLACK",
            spec="e1/0.34",
            vector=v_obs,
            no_clip=True,
        )
        fig.shift_origin(xshift="0c", yshift="1.45c")
        # Plot thick black arrow behind, to get black outline on red arrow
        fig.shift_origin(xshift="0c", yshift="-1.8c")
        fig.velo(
            data=static_legend,
            pen=".05c,black",
            line="black",
            fill="black",
            spec="e1/0",
            vector=v_syn,
            no_clip=True,
        )
        # overlay thin red arrow
        fig.velo(
            data=static_legend,
            pen=".03c,red",
            line="red",
            fill="red",
            spec="e1/0",
            vector=v_syn,
            no_clip=True,
        )
        fig.shift_origin(xshift="0c", yshift="1.8c")
        fig.shift_origin(xshift="-0.2c", yshift="-1.1c")
        if legend_len <= 10:  # type:ignore
            fig.text(
                x=region[1],
                y=region[2],
                text=str(legend_len * 10)  # type:ignore
                + "+/-"
                + str(legend_len)
                + " mm",  # type:ignore
                no_clip=True,
                justify="ML",
            )
        else:
            fig.text(
                x=region[1],
                y=region[2],
                text=str(legend_len)
                + "+/-"
                + str(legend_len / 10)  # type:ignore
                + " cm",  # type:ignore
                no_clip=True,
                justify="ML",
            )
        fig.shift_origin(xshift="0.2c", yshift="1.1c")
    ##################################
    ### PLOT FAULT TRACES OVER TOP ###
    ##################################
    faults = glob.glob(str(directory) + "/*.fault")
    if len(faults) > 0:
        for kfault in range(len(faults)):
            print("...Adding fault trace from: " + str(faults[kfault]))
            fault_trace = np.genfromtxt(faults[kfault])
            fig.plot(x=fault_trace[:, 0], y=fault_trace[:, 1], pen="1p,black")
            fig.plot(x=fault_trace[:, 0], y=fault_trace[:, 1], pen="0.5p,grey")

    ###################
    ### INSET GLOBE ###
    ###################
    with fig.inset(position="jTR+w100p+o-50p", margin=0):
        fig.coast(
            region="g",
            projection="G" + str(lon0) + "/" + str(lat0) + "/?",
            land="gray",
            water="white",
            frame="g",
        )
        fig.plot(x=lon0, y=lat0, style="a7p", fill="gold", pen="black")

    fig.savefig(directory / "Map.eps")
    fig.savefig(directory / "Map.png")
    fig.savefig(directory / "Map.pdf")


def PlotInsar(
    tensor_info: dict,
    segments: List[dict],
    point_sources: list,
    solution: dict,
    insar_points: List[dict],
    scene: str,
    los: str = "ascending",
    directory: Union[pathlib.Path, str] = pathlib.Path(),
):
    """Plot InSAR observation and model fits

    :param tensor_info: The tensor information
    :type tensor_info: dict
    :param segments: The segment properties
    :type segments: List[dict]
    :param point_sources: The point source locations
    :type point_sources: list
    :param solution: The kinematic solution read from Solution.txt
    :type solution: dict
    :param insar_points: List of insar data
    :type insar_points: List[dict]
    :param scene: The scene name or number
    :type scene: str
    :param los: The direction of the path, defaults to "ascending"
    :type los: str, optional
    :param directory: The location where plots should be written, defaults to pathlib.Path()
    :type directory: Union[pathlib.Path, str], optional
    """
    directory = pathlib.Path(directory)
    print(f"Creating InSAR plots: {los}")
    plane_info = segments[0]
    (
        stk_subfaults,
        dip_subfaults,
        delta_strike,
        delta_dip,
        hyp_stk,
        hyp_dip,
    ) = pl_mng.__unpack_plane_data(plane_info)
    segments_lats, segments_lons, segments_deps = __redefine_lat_lon(
        segments, point_sources
    )
    min_lats = [np.min(segment_lat.flatten()) for segment_lat in segments_lats]
    max_lats = [np.max(segment_lat.flatten()) for segment_lat in segments_lats]
    min_lons = [np.min(segment_lon.flatten()) for segment_lon in segments_lons]
    max_lons = [np.max(segment_lon.flatten()) for segment_lon in segments_lons]
    min_lat = np.min(min_lats)  # - 0.5
    max_lat = np.max(max_lats)  # + 0.5
    min_lon = np.min(min_lons)  # - 0.5
    max_lon = np.max(max_lons)  # + 0.5

    max_diff = -1
    min_diff = 1
    lats = [point["lat"] for point in insar_points]
    lons = [point["lon"] for point in insar_points]
    min_lat = min(min_lat, np.min(lats))
    max_lat = max(max_lat, np.max(lats))
    min_lon = min(min_lon, np.min(lons))
    max_lon = max(max_lon, np.max(lons))
    observed = [point["observed"] for point in insar_points]
    synthetic = [point["synthetic"] for point in insar_points]
    ramp = [point["ramp"] for point in insar_points]
    diffs = [obs - syn for obs, syn in zip(observed, synthetic)]
    obs_no_ramp = [obs - ramp for obs, ramp in zip(observed, ramp)]
    syn_no_ramp = [syn - ramp for syn, ramp in zip(synthetic, ramp)]

    max_obs_abs = np.max(np.abs(observed))
    min_obs_abs = -max_obs_abs

    fig = pygmt.Figure()
    pygmt.config(
        PS_MEDIA="50ix50i",
        MAP_FRAME_TYPE="plain",
        MAP_FRAME_AXES="WSen",
        MAP_TITLE_OFFSET="-0.3",
        FORMAT_GEO_OUT="F",
        FORMAT_GEO_MAP="ddd:mm:ss",
        FONT_ANNOT_PRIMARY="10p,Helvetica,black",
        FONT_ANNOT_SECONDARY="6p,Helvetica,black",
        FONT_HEADING="30p,Helvetica,black",
        FONT_LABEL="12p,Helvetica,black",
        FONT_LOGO="6p,Helvetica,black",
        FONT_SUBTITLE="16p,Helvetica,black",
        FONT_TAG="18p,Helvetica,black",
        FONT_TITLE="16p,Helvetica,black",
        MAP_ANNOT_OFFSET_PRIMARY="3p",
    )

    region = [min_lon - 0.2, max_lon + 0.2, min_lat - 0.2, max_lat + 0.2]
    region[0] = np.floor(region[0] * 10) / 10.0
    region[1] = np.ceil(region[1] * 10) / 10.0
    region[2] = np.floor(region[2] * 10) / 10.0
    region[3] = np.ceil(region[3] * 10) / 10.0
    resolution = "03s"
    map_scale_len = "50k"
    if region[1] - region[0] > 5 or region[3] - region[2] > 5:
        resolution = "30s"
        map_scale_len = "100k"
    if region[1] - region[0] > 10 or region[3] - region[2] > 10:
        resolution = "01m"
        map_scale_len = "200k"
    grid = pygmt.datasets.load_earth_relief(resolution=resolution, region=region)

    map_width = 6
    map_height = int(
        np.ceil(map_width / ((region[1] - region[0]) / (region[3] - region[2])))
    )

    projection = "M" + str(map_width) + "c"
    map_scale = (
        "g"
        + str(region[1])
        + "/"
        + str(region[3])
        + "+c17.40+w"
        + map_scale_len
        + "+ar+l+jTR+o0.8/0.5+f"
    )

    titles = [
        "Observed",
        "Synthetic",
        "Misfit",
        "Observed - Ramp",
        "Synthetic - Ramp",
        "Ramp",
    ]
    legends = [
        "Observed LOS Displacement (cm)",
        "Modeled LOS Displacement (cm)",
        "Residual (cm)",
        "Observed LOS Displacement (cm)",
        "Modeled LOS Displacement (cm)",
        "Modeled Ramp (cm)",
    ]
    values = [observed, synthetic, diffs, obs_no_ramp, syn_no_ramp, ramp]
    xshift = [
        0,
        map_width + 2,
        map_width + 2,
        -2 * (map_width + 2),
        map_width + 2,
        map_width + 2,
    ]
    yshift = [23, 0, 0, -(5 + map_height), 0, 0]

    sub = 0
    for subplot in range(6):
        title = "+t" + titles[subplot]
        print(f"...Subplot {titles[subplot]}")
        fig.shift_origin(xshift=xshift[sub], yshift=yshift[sub])
        fig.basemap(
            region=region,
            projection=projection,
            frame=["a", title],
            # map_scale=map_scale,
        )
        fig.coast(resolution="h", shorelines=True)
        fig.grdimage(grid=grid, cmap="oleron", shading=True, transparency=80)

        ##################################
        ### PLOT FAULT TRACES OVER TOP ###
        ##################################
        faults = glob.glob(str(directory) + "/*.fault")
        if len(faults) > 0:
            for kfault in range(len(faults)):
                print("...Adding fault trace from: " + str(faults[kfault]))
                fault_trace = np.genfromtxt(faults[kfault])
                fig.plot(x=fault_trace[:, 0], y=fault_trace[:, 1], pen="1p,black")

        pygmt.makecpt(cmap="roma", series=[min_obs_abs, max_obs_abs])
        fig.plot(
            x=lons, y=lats, style="c0.1c", cmap=True, fill=values[sub], pen="0.5p,black"
        )

        fig.colorbar(
            position="jBC",
            frame="x+l" + legends[subplot],
            # box="+p2p,black+ggray80"
        )

        #################################
        # outline the segments in black #
        #################################
        for segment in range(len(segments_lats)):
            plane_info = segments[segment]
            (
                stk_subfaults,
                dip_subfaults,
                delta_strike,
                delta_dip,
                hyp_stk,
                hyp_dip,
            ) = pl_mng.__unpack_plane_data(plane_info)
            depths = segments_deps[segment].flatten()
            min_lats_idx = np.where(
                segments_lats[segment].flatten() == min_lats[segment]
            )[0][0]
            cornerA = [
                segments_lons[segment].flatten()[min_lats_idx],
                min_lats[segment],
                depths[min_lats_idx],
                0,
            ]
            min_lons_idx = np.where(
                segments_lons[segment].flatten() == min_lons[segment]
            )[0][0]
            cornerB = [
                min_lons[segment],
                segments_lats[segment].flatten()[min_lons_idx],
                depths[min_lons_idx],
                1,
            ]
            max_lats_idx = np.where(
                segments_lats[segment].flatten() == max_lats[segment]
            )[0][-1]
            cornerC = [
                segments_lons[segment].flatten()[max_lats_idx],
                max_lats[segment],
                depths[max_lats_idx],
                2,
            ]
            max_lons_idx = np.where(
                segments_lons[segment].flatten() == max_lons[segment]
            )[0][-1]
            cornerD = [
                max_lons[segment],
                segments_lats[segment].flatten()[max_lons_idx],
                depths[max_lons_idx],
                3,
            ]
            corners = np.c_[cornerA, cornerB, cornerC, cornerD]
            max_dep = max(corners[2, :])
            idx_max = np.where(np.array(corners[2, :]) == max_dep)[0][0]
            min1 = min2 = corners[:, idx_max]
            for i in range(len(corners)):
                if corners[2, i] < min1[2]:
                    min2 = min1
                    min1 = corners[:, i]
                elif corners[2, i] < min2[2] and collections.Counter(
                    corners[:, i]
                ) != collections.Counter(min1):
                    min2 = corners[:, i]
            updip = np.c_[min1, min2]
            corners = np.c_[corners, cornerA]
            fig.plot(x=corners[0, :], y=corners[1, :], pen="1p,black")
            fig.plot(x=updip[0, :], y=updip[1, :], pen="1p,red")
        sub += 1
    fig.savefig(directory / "InSAR_{}_fit_{}.png".format(los, scene))
    fig.savefig(directory / "InSAR_{}_fit_{}.pdf".format(los, scene))
    plt.close()
    return


def PlotComparisonMap(
    tensor_info: dict,
    segments: List[dict],
    point_sources: list,
    input_model: dict,
    solution: dict,
    max_val: Optional[float] = None,
    directory: Union[pathlib.Path, str] = pathlib.Path(),
):
    """Plot the slip map

    :param tensor_info: The tensor information
    :type tensor_info: dict
    :param segments: The segment properties
    :type segments: List[dict]
    :param point_sources: The point source locations
    :type point_sources: list
    :param input_model: the model
    :type input_model: dict
    :param solution: The kinematic solution read from Solution.txt
    :type solution: dict
    :param max_val: Specify a maximum value, defaults to None
    :type max_val: Optional[float], optional
    :param directory: Where to write the plots, defaults to pathlib.Path()
    :type directory: Union[pathlib.Path,str], optional
    """
    directory = pathlib.Path(directory)
    print("Plotting Comparison Map")
    input_slip = input_model["slip"]
    plane_info = segments[0]
    (
        stk_subfaults,
        dip_subfaults,
        delta_strike,
        delta_dip,
        hyp_stk,
        hyp_dip,
    ) = pl_mng.__unpack_plane_data(plane_info)
    if stk_subfaults * dip_subfaults == 1:
        return
    slip = solution["slip"]
    #
    # accurate plot coordinates
    #
    segments_lats, segments_lons, segments_deps = __redefine_lat_lon(
        segments, point_sources
    )

    margin = min(1.5 * (stk_subfaults * delta_strike) / 111.19, 10)  # 3
    lat0 = tensor_info["lat"]
    lon0 = tensor_info["lon"]
    margins = [lon0 - margin, lon0 + margin, lat0 - margin, lat0 + margin]
    dictn = {"projection": ccrs.PlateCarree(), "facecolor": "#eafff5"}
    shpfilename = shpreader.natural_earth(
        resolution="10m", category="cultural", name="admin_0_countries"
    )
    countries = cf.ShapelyFeature(
        shpreader.Reader(shpfilename).geometries(),
        ccrs.PlateCarree(),
        edgecolor="black",
        facecolor="lightgray",
    )

    fig, axes = plt.subplots(1, 2, figsize=(30, 15), subplot_kw=dictn)
    ax1: Axes = axes[0]  # type: ignore
    ax2: Axes = axes[1]  # type: ignore
    ax1.set_title("Inverted model", fontsize=22)
    ax2.set_title("Original model", fontsize=22)
    fig.subplots_adjust(hspace=0, wspace=0.1, top=0.9, bottom=0.3)
    for ax in [ax1, ax2]:
        ax = set_map_cartopy(ax, margins, countries=countries)
        ax.plot(lon0, lat0, "w*", markersize=15, transform=ccrs.PlateCarree(), zorder=4)
    max_slip = max([np.amax(slip_fault) for slip_fault in slip])
    max_slip2 = max([np.amax(input_slip2) for input_slip2 in input_slip])
    max_slip = max(max_slip, max_slip2)
    #
    # plot slip map
    #
    if max_val == None:
        max_val = max_slip
    ax1, cs1 = plot_map(
        ax1,
        segments_lats,  # type:ignore
        segments_lons,  # type:ignore
        slip,  # type:ignore
        max_val=max_val,
        transform=dictn["projection"],
        cmap=slipcpt,
    )
    ax2, cs2 = plot_map(
        ax2,
        segments_lats,  # type:ignore
        segments_lons,  # type:ignore
        input_slip,  # type:ignore
        max_val=max_val,
        transform=dictn["projection"],
        cmap=slipcpt,
    )
    fig.subplots_adjust(bottom=0.15)
    cbar_ax = fig.add_axes((0.1, 0.05, 0.8, 0.05))
    cb = fig.colorbar(cs2, cax=cbar_ax, orientation="horizontal")
    cb.set_label("Slip (cm)")
    plt.savefig(directory / "Checkerboard_Map_Comparison.png", bbox_inches="tight")
    plt.close()
    return


def __redefine_lat_lon(
    segments: List[dict], point_sources: list
) -> Tuple[list, list, list]:
    """Redefine the lat/lon

    :param segments: The segments
    :type segments: List[dict]
    :param point_sources: The point sources
    :type point_sources: list
    :return: segments latitudes, segments longitudes, and segments depths
    :rtype: Tuple[list, list, list]
    """
    segments_lats: list = [[]] * len(segments)
    segments_lons: list = [[]] * len(segments)
    segments_deps: list = [[]] * len(segments)
    for i, point_sources_seg in enumerate(point_sources):
        lat = point_sources_seg[:, :, :, :, 0]
        lon = point_sources_seg[:, :, :, :, 1]
        dep = point_sources_seg[:, :, :, :, 2]
        ny, nx, a, b = lat.shape
        new_lat = np.zeros((ny + 1, nx + 1))
        new_lon = np.zeros((ny + 1, nx + 1))
        new_dep = np.zeros((ny + 1, nx + 1))
        for j in range(ny):
            for k in range(nx):
                new_lat[j, k] = lat[j, k, 0, 0]
                new_lon[j, k] = lon[j, k, 0, 0]
                new_dep[j, k] = dep[j, k, 0, 0]
        for k in range(nx):
            new_lat[-1, k] = lat[-1, k, -1, 0]
            new_lon[-1, k] = lon[-1, k, -1, 0]
            new_dep[-1, k] = dep[-1, k, -1, 0]
        for j in range(ny):
            new_lat[j, -1] = lat[j, -1, 0, -1]
            new_lon[j, -1] = lon[j, -1, 0, -1]
            new_dep[j, -1] = dep[j, -1, 0, -1]
        new_lat[-1, -1] = lat[-1, -1, -1, -1]
        new_lon[-1, -1] = lon[-1, -1, -1, -1]
        new_dep[-1, -1] = dep[-1, -1, -1, -1]
        segments_lats[i] = new_lat
        segments_lons[i] = new_lon
        segments_deps[i] = new_dep
    return segments_lats, segments_lons, segments_deps


def plot_moment_rate_function(
    segments_data: dict,
    shear: list,
    point_sources: list,
    mr_time: Optional[int] = None,
    separate_planes: bool = False,
    directory: Union[pathlib.Path, str] = pathlib.Path(),
):
    """Plot moment rate function

    :param segments_data: The segment properties
    :type segments_data: dict
    :param shear: The shear moduli
    :type shear: list
    :param point_sources: The point source locations
    :type point_sources: list
    :param mr_time: The moment rate time, defaults to None
    :type mr_time: Optional[int], optional
    :param separate_planes: Whether there are separate planes, defaults to False
    :type separate_planes: bool, optional
    :param directory: Where to write plots, defaults to pathlib.Path()
    :type directory: Union[pathlib.Path,str], optional
    """
    directory = pathlib.Path(directory)
    print("Creating Moment Rate Plot...")
    segments = segments_data["segments"]
    plane_info = segments[0]
    dt = 0.01
    model = load_ffm_model.load_ffm_model(
        segments_data, point_sources, directory=directory
    )
    slip = model["slip"]
    trup = model["trup"]
    tl = model["trise"]
    tr = model["tfall"]
    properties = pl_mng.__unpack_plane_data(plane_info)
    delta_strike, delta_dip = [properties[2], properties[3]]
    stk_subfaults, dip_subfaults = [properties[0], properties[1]]
    tmax = 1.5 * np.max(
        [
            np.amax(trup_seg + tl_seg + tr_seg)
            for trup_seg, tl_seg, tr_seg in zip(trup, tl, tr)
        ]
    )
    tmax = (
        tmax if stk_subfaults * dip_subfaults > 1 else (tl[0][0, 0] + tr[0][0, 0]) * 8.0
    )
    nmax = int((tmax / dt + 1))
    mr = np.zeros(nmax)
    mr_seg = np.zeros(nmax)
    mr_seg_all = np.zeros((nmax, len(segments)))
    seismic_moment = 0
    seg_num = 0
    list = [item for sublist in tl for item in sublist]
    flat_list = [item for sublist in list for item in sublist]
    if all(flat_list) == 0:
        print("This looks like a static solution. Skipping moment rate plot.")
    else:
        ## Set up plot ##
        fig = plt.figure(figsize=(10, 8))
        plt.rc("axes", titlesize=16)
        plt.rc("axes", labelsize=16)
        plt.rc("xtick", labelsize=14)
        plt.rc("ytick", labelsize=14)
        plt.rc("font", size=16)
        ax = fig.add_subplot(111)
        ax.spines["bottom"].set_linewidth(3)
        ax.spines["top"].set_linewidth(3)
        ax.spines["left"].set_linewidth(3)
        ax.spines["right"].set_linewidth(3)
        ax.set_xlabel("Time (s)")
        ax.set_ylabel("Relative Moment Rate (Nm/s)")

        for (
            segment,
            slip_seg,
            trup_seg,
            trise_seg,
            tfall_seg,
            shear_seg,
            point_sources_seg,
        ) in zip(segments, slip, trup, tl, tr, shear, point_sources):
            dip_subfaults, stk_subfaults = np.shape(slip_seg)
            moment_rate = np.zeros(nmax)
            for iy in range(dip_subfaults):
                for ix in range(stk_subfaults):
                    rupt_vel = segment["rupture_vel"]
                    rise_time = np.zeros(nmax)
                    tfall = tfall_seg[iy, ix]
                    trise = trise_seg[iy, ix]
                    array1 = np.arange(0, trise, dt)
                    tmid = len(array1)
                    rise_time[:tmid] = (1 - np.cos(np.pi * array1 / trise)) / (
                        trise + tfall
                    )
                    array2 = np.arange(0, tfall, dt)
                    tend = tmid + len(array2)
                    rise_time[tmid:tend] = (1 + np.cos(np.pi * array2 / tfall)) / (
                        tfall + trise
                    )
                    duration = int(max(delta_strike, delta_dip) / dt / rupt_vel)
                    source_dur = np.ones(duration) / duration
                    start_index = max(0, int(trup_seg[iy, ix] / dt))
                    product = slip_seg[iy, ix] * shear_seg[iy, ix] / 100 / 10
                    sub_rise_time = rise_time[:tend]
                    convolve = np.convolve(source_dur, sub_rise_time)
                    moment_rate[start_index : start_index + len(convolve)] = (
                        moment_rate[start_index : start_index + len(convolve)]
                        + convolve * product
                    )

            seismic_moment = seismic_moment + np.sum(
                (slip_seg / 100)
                * (shear_seg / 10)
                * (delta_strike * 1000)
                * (delta_dip * 1000)
            )
            #
            # find moment rate function
            #
            for i in range(nmax):
                time = i * dt
                mr[i] = mr[i] + moment_rate[i] * (delta_strike * 1000) * (
                    delta_dip * 1000
                )
                mr_seg[i] = moment_rate[i] * (delta_strike * 1000) * (delta_dip * 1000)
                mr_seg_all[i, seg_num] = mr_seg[i]
            t_seg = np.arange(nmax) * dt
            with open(directory / f"STF_{seg_num}.txt", "w") as outsegf:
                outsegf.write("dt: {}\n".format(dt))
                outsegf.write("Time[s]    Moment_Rate [Nm/s]\n")
                for t, val in zip(t_seg, mr_seg):
                    outsegf.write("{:8.2f}\t\t{:8.4e}\n".format(t, val))
            seg_num += 1
        time = np.arange(nmax) * dt  # type:ignore
        with open(directory / "STF.txt", "w") as outf:
            outf.write("dt: {}\n".format(dt))
            outf.write("Time[s]     Moment_Rate [Nm/s]\n")
            for t, val in zip(time, mr):  # type:ignore
                outf.write("{:8.2f}\t\t{:8.4e}\n".format(t, val))

        seismic_moment = np.trapz(mr, dx=0.01)
        magnitude = 2.0 * (np.log10(seismic_moment * 10**7) - 16.1) / 3.0
        rel_mr = mr / (max(mr))
        plt.text(
            0.99 * max(time),  # type:ignore
            0.95 * max(rel_mr),
            "Max Mr: {:.2E} Nm/sec".format(max(mr)),  # type:ignore
            ha="right",
            fontweight="bold",
        )
        plt.text(
            0.99 * max(time),  # type:ignore
            0.90 * max(rel_mr),
            "M$_0$: {:.2E} Nm".format(seismic_moment),
            ha="right",
            fontweight="bold",
        )
        plt.text(
            0.99 * max(time),  # type:ignore
            0.85 * max(rel_mr),
            "M$_w$: {:.2f}".format(magnitude),
            ha="right",
            fontweight="bold",
        )
        plt.grid(visible=True, ls="dotted")
        ax.fill_between(time, rel_mr, color="0.9")
        ax.plot(time, rel_mr, "k", lw=2, label="Total")
        ### Plot multiple planes separately? ###
        if separate_planes == True:
            for seg in range(np.shape(mr_seg_all)[1]):
                ax.plot(time, mr_seg_all[:, seg] / max(mr), label="Segment " + str(seg))
            plt.legend(loc="lower right")
        if mr_time == None:
            tenth = np.ones(len(time)) * 0.1  # type:ignore
            idx = np.argwhere(np.diff(np.sign(rel_mr - tenth))).flatten()
            ax.vlines(
                time[idx][-1],  # type:ignore
                0,
                1,
                "r",
                linestyle="dashed",
                lw=2,
                dashes=(9, (5, 3)),  # type:ignore
            )
        else:
            ax.vlines(
                mr_time,  # type:ignore
                0,
                1,
                "r",
                linestyle="dashed",
                lw=2,
                dashes=(9, (5, 3)),  # type:ignore
            )
        ax.set_ylim((0, 1))
        ax.set_xlim((0, max(time)))  # type:ignore
        plt.grid(which="minor", linestyle="dotted", color="0.5")
        plt.grid(which="major", linestyle="dotted", color="0.5")
        ax.yaxis.set_major_locator(matplotlib.ticker.MultipleLocator(0.2))
        ax.yaxis.set_minor_locator(matplotlib.ticker.MultipleLocator(0.1))
        ax.xaxis.set_major_locator(matplotlib.ticker.MultipleLocator(20))
        ax.xaxis.set_minor_locator(matplotlib.ticker.MultipleLocator(10))
        plt.savefig(directory / "MomentRate.png")
        plt.close()
        return


def _PlotSnapshotSlip(
    segments: dict, solution: dict, directory: Union[pathlib.Path, str] = pathlib.Path()
):
    """Plot snapshots of the rupture process

    :param segments: The segment properties
    :type segments: dict
    :param solution: The kinematic solution read from Solution.txt
    :type solution: dict
    :param directory: Where to write the plots, defaults to pathlib.Path()
    :type directory: Union[pathlib.Path,str], optional
    """
    directory = pathlib.Path(directory)
    plane_info = segments[0]
    dt = 0.01
    slip = solution["slip"]
    trup = solution["rupture_time"]
    tl = solution["t_rise"]
    tr = solution["t_fall"]
    (
        stk_subfaults,
        dip_subfaults,
        delta_strike,
        delta_dip,
        hyp_stk,
        hyp_dip,
    ) = pl_mng.__unpack_plane_data(plane_info)
    if [stk_subfaults, dip_subfaults] == [1, 1]:
        return
    tmid = trup + tl
    tstop = trup + tl + tr
    srmax = slip / (tr + tl) * 2.0
    #
    # Define the vector field of the slip and plotting parameters.
    #
    x = np.arange(stk_subfaults) * delta_strike - hyp_stk * delta_strike
    y = np.arange(dip_subfaults) * delta_dip - hyp_dip * delta_dip
    ratio = max(
        1, (9 / 16) * ((stk_subfaults * delta_strike) / (dip_subfaults * delta_dip))
    )
    vmax = np.amax(slip)
    tmax = np.amax(trup)
    step = int((tmax / dt + 1) / 9.0)
    #
    # snapshots of rupture process
    #
    fig, axes = plt.subplots(nrows=3, ncols=3, figsize=(9, 9))
    i = 1

    for ax in axes.flat:  # type: ignore
        time = i * step * dt  # type:ignore
        srate, cslip, broken = __rupture_process(
            time, slip, srmax, trup, tl, tr, tmid, tstop  # type:ignore
        )
        ax.set_yticklabels([])
        ax.set_xticklabels([])
        if np.max(broken) > np.min(broken):
            ax.contour(-broken, 1, colors="k", extent=__extent_plot(plane_info))
        ax.contourf(
            cslip, cmap="jet", vmin=0, vmax=vmax, extent=__extent_plot(plane_info)
        )
        ax.plot(0, 0, "r*", ms=20)
        ax.invert_yaxis()
        ax.set_aspect(ratio)
        ax.set_title("Time: {0:.2f} s".format(time))
        if i == 9:
            im = ax.contourf(x, y, cslip, cmap="jet")
        i = i + 1

    fig.subplots_adjust(right=0.8)
    cbar_ax = fig.add_axes((0.85, 0.15, 0.05, 0.7))
    cb = fig.colorbar(im, cax=cbar_ax)
    cb.set_label("Slip (cm)")
    fig.text(
        0.1,
        0.1,
        "CSN Automatic \nSolution",
        fontsize=50,
        color="gray",
        ha="left",
        va="bottom",
        alpha=0.5,
        wrap=True,
    )
    plt.savefig(directory / "SlipSnapshot.png", bbox_inches="tight")
    plt.close()
    return


def shakemap_polygon(
    segments: dict,
    point_sources: list,
    solution: dict,
    tensor_info: dict,
    evID: str,
    directory: Union[pathlib.Path, str] = pathlib.Path(),
):
    """Write the shakemap polygon

    :param segments: The segment properties
    :type segments: dict
    :param point_sources: The point source locations
    :type point_sources: list
    :param solution: The kinematic solution read from Solution.txt
    :type solution: dict
    :param tensor_info: The tensor information
    :type tensor_info: dict

    :param evID: The event name/ID
    :type evID: str
    :param directory: Where to write the file, defaults to pathlib.Path()
    :type directory: Union[pathlib.Path,str], optional
    """
    directory = pathlib.Path(directory)
    corner_1, corner_2, corner_3, corner_4 = _PlotCumulativeSlip(
        segments, point_sources, solution, tensor_info, evID=evID, directory=directory
    )
    # translate above output into lon/lat/dep for txt output #
    print("Writing shakemap_polygon.txt to file...")
    ### Create shakemap_polygon.txt file ###
    if evID == None:
        evID = "None Provided"
    now = datetime.utcnow().strftime("%Y-%m-%d %H:%M:%S")
    corners: list = []
    with open(directory / "shakemap_polygon.txt", "w") as txtout:
        txtout.write("#Source: USGS NEIC Rapid Finite Fault \n")
        txtout.write("#Event ID: " + evID + "\n")
        txtout.write("#Model created: " + now + "\n")
        txtout.write(corner_1 + "\n")
        txtout.write(corner_2 + "\n")
        txtout.write(corner_3 + "\n")
        txtout.write(corner_4 + "\n")
        txtout.write(corner_1 + "\n")


def calculate_cumulative_moment_tensor(
    solution: dict, directory: Union[pathlib.Path, str] = pathlib.Path()
):
    """Calculate the cumulative moment tensor

    :param solution: The kinematic solution read from Solution.txt
    :type solution: dict
    :param directory: Whwere the segments file is located and where to write to,
                    defaults to pathlib.Path()
    :type directory: Union[pathlib.Path, str], optional
    """
    directory = pathlib.Path(directory)
    print("Calculating Cumulative Moment Tensor...")

    with open(directory / "segments_data.json") as s:
        segments_data = json.load(s)
    segments = segments_data["segments"]

    slip = solution["slip"]
    rake = solution["rake"]
    rupture_time = solution["rupture_time"]
    trise = solution["trise"]
    tfall = solution["tfall"]
    lat = solution["lat"]
    lon = solution["lon"]
    depth = solution["depth"]
    moment = solution["moment"]

    num_seg = len(solution["slip"])

    Mm6 = [0, 0, 0, 0, 0, 0]
    total_moment = 0
    fig = plt.figure(figsize=(num_seg + 2.0, 2.0), dpi=300)
    fig.subplots_adjust(left=0.0, right=1.0, bottom=0.0, top=1.0)

    for kseg in range(num_seg):
        print("Segment: ", kseg)
        axes = fig.add_subplot(1, num_seg + 1, kseg + 1)
        axes.set_axis_off()
        strk = segments[kseg]["strike"]
        dip = segments[kseg]["dip"]
        Mm6_seg = [0, 0, 0, 0, 0, 0]
        seg_moment = sum(moment[kseg].flatten())
        seg_Mw = (2.0 / 3) * (np.log10(seg_moment * 1e-7) - 9.1)
        print("   Mw:", seg_Mw)
        total_moment += seg_moment
        for subfault in range(len(slip[kseg].flatten())):
            Mmt_seg = pmt.MomentTensor(
                strike=strk,
                dip=dip,
                rake=rake[kseg].flatten()[subfault],
                scalar_moment=moment[kseg].flatten()[subfault],
            )
            Mm6_seg = [
                Mm6_seg[0] + Mmt_seg.mnn,
                Mm6_seg[1] + Mmt_seg.mee,
                Mm6_seg[2] + Mmt_seg.mdd,
                Mm6_seg[3] + Mmt_seg.mne,
                Mm6_seg[4] + Mmt_seg.mnd,
                Mm6_seg[5] + Mmt_seg.med,
            ]
        print("   MT:", Mm6_seg)
        beachball.plot_beachball_mpl(
            pmt.as_mt(Mm6_seg),
            axes,
            size=60,
            position=(0.5, 0.5),
            beachball_type="full",
            color_t=plot.mpl_color("blue"),
            linewidth=1.0,
        )
        Mm6 = [
            Mm6[0] + Mm6_seg[0],
            Mm6[1] + Mm6_seg[1],
            Mm6[2] + Mm6_seg[2],
            Mm6[3] + Mm6_seg[3],
            Mm6[4] + Mm6_seg[4],
            Mm6[5] + Mm6_seg[5],
        ]
        (s1, d1, r1), (s2, d2, r2) = pmt.as_mt(Mm6_seg).both_strike_dip_rake()
        axes.text(0.1, 0.8, "Segment " + str(kseg) + "\n Mw" + "{:0.2f}".format(seg_Mw))
        axes.text(0.1, 0.1, "s1/d1/r1: %d, %d, %d" % (s1, d1, r1), fontsize=5)
        axes.text(0.1, 0.05, "s2/d2/r2: %d, %d, %d" % (s2, d2, r2), fontsize=5)
        if kseg < num_seg - 1:
            axes.text(1.05, 0.5, "+")
        else:
            axes.text(1.05, 0.5, "=")

    axes = fig.add_subplot(1, num_seg + 1, num_seg + 1)
    axes.set_axis_off()
    beachball.plot_beachball_mpl(
        pmt.as_mt(Mm6),
        axes,
        size=60,
        position=(0.5, 0.5),
        beachball_type="full",
        color_t=plot.mpl_color("blue"),
        linewidth=1.0,
    )
    total_Mw = (2.0 / 3) * (np.log10(total_moment * 1e-7) - 9.1)
    axes.text(0.1, 0.8, "Total MT \n Mw" + "{:0.2f}".format(total_Mw))
    (s1, d1, r1), (s2, d2, r2) = pmt.as_mt(Mm6).both_strike_dip_rake()
    axes.text(0.1, 0.1, "s1/d1/r1: %d, %d, %d" % (s1, d1, r1), fontsize=5)
    axes.text(0.1, 0.05, "s2/d2/r2: %d, %d, %d" % (s2, d2, r2), fontsize=5)
    print("Total Mw:", total_Mw)
    print("Total MT:", Mm6)
    plt.savefig(directory / "Cumulative_Moment_Tensor.png")
    plt.close()
    return


def plot_beachball(
    segments: dict,
    files: Optional[List[dict]] = None,
    phase: Optional[str] = None,
    directory: Union[pathlib.Path, str] = pathlib.Path(),
):
    """Plot the beachball for the event. Optionally, Add the
    location of teleseismic data used in the FFM modelling

    :param segments: The segment properties
    :type segments: dict
    :param files: The file(s) information, defaults to None
    :type files: Optional[List[dict]], optional
    :param phase: The phase, defaults to None
    :type phase: Optional[str], optional
    :param directory: Where to write the plot, defaults to pathlib.Path()
    :type directory: Union[pathlib.Path,str], optional
    """
    directory = pathlib.Path(directory)
    #
    # Get the focal mechanism
    #
    plane_info = segments[0]
    strike = plane_info["strike"]
    dip = plane_info["dip"]
    rake = plane_info["rake"]
    fm = [strike, dip, rake]
    #
    # Plot the beach ball
    #
    fig = plt.figure(figsize=(7, 7))
    bb = beach(fm, width=3.0, zorder=1)
    ax = plt.gca()
    ax.add_collection(bb)
    plt.plot(0, 0, "b*", markersize=20)
    ax.set_aspect("equal")
    #
    # Load the station informations
    #
    if files:
        plt.plot(0, 0, "r*", markersize=20)
        for file in files:
            comp = file["component"]
            if comp == "BHZ":
                comp = "P"
            dist = file["distance"]
            az = file["azimuth"]
            nam = file["name"]
            if comp == phase:
                r = dist * (3.0 / 2.0) / 90.0
                x1 = np.sin(az * np.pi / 180.0) * r
                y1 = np.cos(az * np.pi / 180.0) * r
                plt.plot(x1, y1, "ro")
                plt.text(x1 + 0.01, y1 + 0.01, nam)

    fig.patch.set_visible(False)
    plt.gca().axis("off")
    name_plot = "{}_azimuthcover.png".format(phase) if phase else "Tensor.png"
    plt.savefig(directory / name_plot)
    plt.close()
    return


def _plot_waveforms(
    files: List[dict],
    components: List[str],
    type_str: str,
    start_margin: int = 10,
    forward: bool = False,
    directory: Union[pathlib.Path, str] = pathlib.Path(),
):
    """Plot the observed and synthetic data for a set of stations and
    channels

    :param files: The file(s) information, defaults to None
    :type files: Optional[List[str]], optional
    :param components: List of components
    :type components: List[str]
    :param type_str: The data type
    :type type_str: str
    :param start_margin: Where to start, defaults to 10
    :type start_margin: int, optional
    :param forward: Whether model is result of kinematic modelling or not,
                    defaults to False
    :type forward: bool, optional
    :param directory: Where to write the files, defaults to pathlib.Path()
    :type directory: Union[pathlib.Path,str], optional
    """
    directory = pathlib.Path(directory)
    files = [file for file in files if file["component"] in components]
    files = sorted(files, key=lambda k: k["azimuth"])
    azimuth = [file["azimuth"] for file in files]
    fig = plt.figure(figsize=(13, 9))
    numrows_phase = len(files) // 4 + 1
    gs = gridspec.GridSpec(max(4, numrows_phase), 4)
    for file in files:
        dt = file["dt"]
        nstart = file["start_signal"]
        margin = int(start_margin / dt) if nstart > int(start_margin / dt) else 0
        obs = np.array(file["observed"])
        if nstart >= 0:
            obs = np.concatenate((np.zeros(nstart), obs))
        if type_str in ["body", "strong", "cgps"]:
            obs = obs[nstart - margin :]
        if type_str == "cgps":
            obs = obs - obs[10]
        if type_str == "surf":
            if nstart >= 0:
                obs = obs[nstart:]
            else:
                obs = np.concatenate((np.zeros(-nstart), obs))

        obs = obs if not forward else 0 * obs
        syn = np.array(file["synthetic"])
        syn = np.concatenate((np.zeros(margin), syn))
        length = min(len(obs), len(syn), file["duration"])
        length = min(length, int(950 / dt)) if not type_str in ["surf"] else length
        obs = np.array(obs[:length])
        syn = np.array(syn[:length])
        dt = file["dt"]
        az = file["azimuth"]
        dist = file["distance"]
        name = file["name"]
        time = (
            np.arange(-margin, length - margin) * dt
            if not type_str == "dart"
            else np.arange(-margin, length - margin)
        )
        jj = azimuth.index(az)
        weight = file["trace_weight"]
        alpha = 1 if weight > 0 else 0.1
        ax = fig.add_subplot(gs[jj % numrows_phase, jj // numrows_phase])
        ax.plot(time, obs, "k", linewidth=0.8, alpha=alpha)
        ax.plot(time, syn, "r", linewidth=0.8, alpha=alpha)
        min_val = min(np.min(obs), np.min(syn))
        max_val = max(np.max(obs), np.max(syn))
        ax.vlines(0, min_val, max_val)
        ax.text(
            0.1,
            0.9,
            "{:0.1f}".format(az),
            ha="center",
            va="center",
            transform=ax.transAxes,
        )
        ax.text(
            0.1,
            0.1,
            "{:0.1f}".format(dist),
            ha="center",
            va="center",
            transform=ax.transAxes,
        )
        ax.text(0.9, 0.9, name, ha="center", va="center", transform=ax.transAxes)
        ax.set_xlim((np.min(time), np.max(time)))
        ax.set_ylim((min_val, max_val))
        ax.xaxis.set_major_locator(ticker.MaxNLocator(nbins=3, min_n_ticks=3))
        ax.yaxis.set_major_locator(ticker.MaxNLocator(nbins=3, min_n_ticks=3))

    if type_str == "cgps":
        if "LXZ" in components:
            plot_name = "LXZ_cgps_waves.png"
        if "LXN" in components:
            plot_name = "LXN_cgps_waves.png"
        if "LXE" in components:
            plot_name = "LXE_cgps_waves.png"

    if type_str == "strong":
        if "HNZ" in components:
            plot_name = "HNZ_strong_motion_waves.png"
        if "HNN" in components:
            plot_name = "HNN_strong_motion_waves.png"
        if "HNE" in components:
            plot_name = "HNE_strong_motion_waves.png"

    if type_str == "body":
        if "BHZ" in components:
            plot_name = "P_body_waves.png"
        if "BHT" in components:
            plot_name = "SH_body_waves.png"

    if type_str == "surf":
        if "BHZ" in components:
            plot_name = "Rayleigh_surf_waves.png"
        if "BHT" in components:
            plot_name = "Love_surf_waves.png"

    if type_str == "dart":
        plot_name = "Dart_waves.png"

    plt.savefig(directory / plot_name, bbox_inches="tight")
    plt.close()
    return


def __extent_plot(plane_info: dict) -> List[float]:
    """Get the extent of the plot

    :param plane_info: The plane information
    :type plane_info: dict
    :return: The extent
    :rtype: List[float]
    """
    (
        stk_subfaults,
        dip_subfaults,
        delta_strike,
        delta_dip,
        hyp_stk,
        hyp_dip,
    ) = pl_mng.__unpack_plane_data(plane_info)
    return [
        -hyp_stk * delta_strike,
        (stk_subfaults - hyp_stk) * delta_strike,
        -hyp_dip * delta_dip,
        (dip_subfaults - hyp_dip) * delta_dip,
    ]


def __several_axes(
    data: np.ndarray,
    segment: dict,
    point_source_seg: np.ndarray,
    ax: plt.Axes,
    min_val: Optional[Union[list, float]] = None,
    max_val: Optional[Union[list, float]] = None,
    autosize: bool = True,
    cmap: Any = slipcpt,
    depth_ax: bool = True,
) -> Tuple[plt.Axes, AxesImage]:
    """Setup several axes

    :param data: The data to plot
    :type data: np.ndarray
    :param segment: The segment
    :type segment: dict
    :param point_source_seg: The point sources in the segment
    :type point_source_seg: np.ndarray
    :param ax: The axes
    :type ax: plt.Axes
    :param min_val: Specify a minimum value, defaults to None
    :type min_val: Optional[float], optional
    :param max_val: Specify a maximum value, defaults to None
    :type max_val: Optional[float], optional
    :param autosize: Autosize the plot, defaults to True
    :type autosize: bool, optional
    :param cmap: Specify a colormap, defaults to slipcpt
    :type cmap: Any, optional
    :param depth_ax: Create a depth axes, defaults to True
    :type depth_ax: bool, optional
    :return: Return axes
    :rtype: Tuple[plt.Axes, AxesImage]
    """
    (
        stk_subfaults,  # type:ignore
        dip_subfaults,  # type:ignore
        delta_strike,
        delta_dip,
        hyp_stk,
        hyp_dip,
    ) = pl_mng.__unpack_plane_data(segment)
    min_dist = -(hyp_dip + 0.5) * delta_dip
    max_dist = (dip_subfaults - hyp_dip - 0.5) * delta_dip
    min_strike = -(hyp_stk + 0.5) * delta_strike
    max_strike = (stk_subfaults - hyp_stk - 0.5) * delta_strike
    dep = point_source_seg[:, :, :, :, 2]
    min_depth = dep[0, 0, 0, 0]
    max_depth = dep[-1, 0, -1, 0]
    if stk_subfaults * dip_subfaults == 1:
        dip = segment["dip"]
        delta_z = delta_dip * np.sin(dip * np.pi / 180.0) / 2
        min_depth = min_depth - delta_z
        max_depth = max_depth + delta_z
    max_val = np.max(data.flatten()) if not max_val else max_val
    min_val = np.min(data.flatten()) if not min_val else min_val
    im = ax.imshow(
        data,
        cmap=cmap,
        origin="lower",
        vmin=min_val,  # type:ignore
        vmax=max_val,  # type:ignore
        aspect="auto",
        extent=(min_strike, max_strike, min_dist, max_dist),
    )
    ax.set(adjustable="datalim")
    if autosize:
        ax.figure.set_size_inches(  # type:ignore
            4 * stk_subfaults * delta_strike / dip_subfaults / delta_dip, 4
        )
    if depth_ax:
        ax2 = ax.twinx()
        ax2.set_xlim((min_strike, max_strike))
        ax2.set_ylim((min_depth, max_depth))
        ax2.set(adjustable="datalim")
        ax2.set_ylabel("Depth (km)", fontsize=16)
        ax2.invert_yaxis()
        if autosize:
            ax2.figure.set_size_inches(  # type:ignore
                4 * stk_subfaults * delta_strike / dip_subfaults / delta_dip, 4
            )
    ax.invert_yaxis()
    return ax, im


def __rupture_process(
    time: np.ndarray,
    slip: np.ndarray,
    trup: np.ndarray,
    tmid: np.ndarray,
    tstop: np.ndarray,
    rise_time: np.ndarray,
    dt: float,
) -> Tuple[np.ndarray, np.ndarray, np.ndarray]:
    """Give slip rate, rupture front, and accumulated slip at a certain
    time ``time``

    :param time: The time array
    :type time: np.ndarray
    :param slip: The slip array
    :type slip: np.ndarray
    :param trup: The trup array
    :type trup: np.ndarray
    :param tmid: The tmid array
    :type tmid: np.ndarray
    :param tstop: The tstop array
    :type tstop: np.ndarray
    :param rise_time: The rise time array
    :type rise_time: np.ndarray
    :param dt: The time delta
    :type dt: float
    :return: srate, cslip, and broken
    :rtype: Tuple[np.ndarray,np.ndarray, np.ndarray]
    """
    dip_subfaults, stk_subfaults = np.shape(slip)
    srate = np.zeros((dip_subfaults, stk_subfaults))
    cslip = np.zeros((dip_subfaults, stk_subfaults))
    broken = np.ones((dip_subfaults, stk_subfaults))

    for i in range(dip_subfaults):
        for j in range(stk_subfaults):
            convolve = rise_time[i, j, :]
            index = int(time / dt - trup[i, j] / dt)
            if time < trup[i, j]:
                broken[i, j] = 0.0
            elif index < len(convolve):
                srate[i, j] = convolve[index] * slip[i, j]
            if trup[i, j] < time <= tmid[i, j]:
                cslip[i, j] = (time - trup[i, j]) * srate[i, j] / 2.0
            if tmid[i, j] < time <= tstop[i, j]:
                cslip[i, j] = slip[i, j] - (tstop[i, j] - time) * srate[i, j] / 2.0
            if time > tstop[i, j]:
                cslip[i, j] = slip[i, j]
    return srate, cslip, broken


def __add_watermark(fig: plt.Figure) -> plt.Figure:
    """Add a watermark to the plot

    :param fig: The figure
    :type fig: plt.Figure
    :return: The updated figure
    :rtype: plt.Figure
    """
    fig.text(
        0.1,
        0.1,
        "CSN Automatic \nSolution",
        fontsize=50,
        color="gray",
        ha="left",
        va="bottom",
        alpha=0.5,
        wrap=True,
    )
    return fig
