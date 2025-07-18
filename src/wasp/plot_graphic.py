# -*- coding: utf-8 -*-
"""The routines here allow to plot the solution model FFM modelling, as well
as moment rate function, and waveform fits.
"""


import argparse
import errno
import glob
import json
import os
import pathlib
from shutil import move
from typing import List, Optional, Tuple, Union

import cartopy.crs as ccrs  # type: ignore
import cartopy.feature as cf  # type: ignore
import cartopy.io.shapereader as shpreader  # type: ignore
import numpy as np  # type: ignore
from matplotlib import colors, patches  # type: ignore
from matplotlib import pyplot as plt  # type: ignore
from matplotlib.axes import Axes
from matplotlib.image import AxesImage  # type: ignore
from obspy.imaging.beachball import beach  # type: ignore

#
# local modules
#
import wasp.fault_plane as pf
import wasp.plane_management as pl_mng
import wasp.seismic_tensor as tensor
import wasp.velocity_models as mv
from wasp import get_outputs, load_ffm_model
from wasp.many_events import select_waveforms_event
from wasp.plot_maps import plot_borders, plot_map, set_map_cartopy
from wasp.waveform_plots import plot_waveform_fits


def plot_ffm_sol(
    tensor_info: dict,
    segments_data: dict,
    point_sources: list,
    shear: list,
    solution: dict,
    vel_model: dict,
    default_dirs: dict,
    use_waveforms: bool = True,
    event: Optional[str] = None,
    directory: Union[pathlib.Path, str] = pathlib.Path(),
):
    """Main routine. Coordinates execution of different plotting routines

    :param tensor_info: The tensor information
    :type tensor_info: dict
    :param segments: The segment properties
    :type segments: dict
    :param point_sources: The point source locations
    :type point_sources: list
    :param shear: The shear moduli
    :type shear: list
    :param solution: The kinematic model read from Solucion.txt
    :type solution: dict
    :param vel_model: The velocity model
    :type vel_model: dict
    :param default_dirs: The location of default directories
    :type default_dirs: dict
    :param use_waveforms: Use the waveforms, defaults to True
    :type use_waveforms: bool, optional
    :param event: The event, defaults to None
    :type event: Optional[str], optional
    :param directory: Where to read/write file(s), defaults to pathlib.Path()
    :type directory: Union[pathlib.Path,str], optional

    .. rubric:: Example:

    First, we load necessary modules.

    >>> import json
    >>> import get_outputs # Allows us to get properties of inverted model
    >>> import management as mng # Allows us to load location of plotting files
    >>> import fault_plane as pf
    >>> import plane_management as pl_mng
    Next, we load necessary data for plots.
    >>> vel_model = json.load(open('velmodel_data.json')) # Assume vel_model stored in file 'velmodel_data.json'
    >>> segments, rise_time, point_sources = pl_mng.__read_planes_info() # Loads point sources and segments information
    >>> solution = get_outputs.read_solution_static_format(segments, point_sources)
    >>> shear = pf.shear_modulous(point_sources, velmodel=vel_model)
    >>> tensor_info = {
            'moment_mag': 7 * 10 ** 27,
            'date_origin': UTCDateTime(2014, 04, 01, 23, 46, 47)
            'lat': -19.5,
            'lon': -70.5,
            'depth': 25,
            'time_shift': 44,
            'half_duration': 40,
            'centroid_lat': -21,
            'centroid_lon': -70,
            'centroid_depth': 35
        }
    Next, we plot solution
    >>> default_dirs = mng.default_dirs()
    >>> plot_ffm_sol(tensor_info, segments, point_sources, shear, solution,
    >>>              vel_model, default_dirs)

    .. note::

        To plot the results of the FFM modelling, we need to run this code
        in a folder whih contains files Solucion.txt, Fault.time, Fault.pos,
        Event_mult.in, and some among the files synm.tele, synm.str_low,
        synm.str and synm.cgps.

    .. note::

        When running this code manually, it is good idea to check if
        the information among files Solucion.txt, Fault.pos, Fault.time,
        and Event_mult.in is consistent.

    """
    directory = pathlib.Path(directory)
    segments = segments_data["segments"]
    _plot_vel_model(vel_model, point_sources, directory=directory)
    if use_waveforms:
        plot_moment_rate_function(
            segments_data, shear, point_sources, event=event, directory=directory
        )
        _PlotRiseTime(segments, point_sources, solution, directory=directory)
        _PlotRuptTime(segments, point_sources, solution, directory=directory)
    PlotSlipDistribution(segments, point_sources, solution, directory=directory)
    PlotMap(
        tensor_info,
        segments,
        point_sources,
        solution,
        default_dirs,
        event=event,
        directory=directory,
    )


def plot_beachballs(
    segments: dict,
    data_type: List[str],
    directory: Union[pathlib.Path, str] = pathlib.Path(),
):
    """Plot the beachball for the event

    :param segments: The segment properties
    :type segments: dict
    :param data_type: The available data types
    :type data_type: List[str]
    :param directory: Where to write the plot, defaults to pathlib.Path()
    :type directory: Union[pathlib.Path, str], optional
    :raises FileNotFoundError: If tele_waves.json cannot be found
    """
    directory = pathlib.Path(directory)
    if "body" in data_type:
        if not os.path.isfile(directory / "tele_waves.json"):
            raise FileNotFoundError(
                errno.ENOENT, os.strerror(errno.ENOENT), "tele_waves.json"
            )
        with open(directory / "tele_waves.json") as tw:
            traces_info = json.load(tw)
        plot_beachball(segments, files=traces_info, phase="P", directory=directory)
        plot_beachball(segments, files=traces_info, phase="SH", directory=directory)


def plot_misfit(
    used_data_type: List[str],
    event: Optional[int] = None,
    directory: Union[pathlib.Path, str] = pathlib.Path(),
):
    """Plot misfit of observed and synthetic data

    :param used_data_type: The data types used
    :type used_data_type: List[str]
    :param event: The event, defaults to None
    :type event: Optional[int], optional
    :param directory: Where to write the plot, defaults to pathlib.Path()
    :type directory: Union[pathlib.Path, str], optional
    :raises FileNotFoundError: If one of the input json files cannot be found
    """
    directory = pathlib.Path(directory)
    if "dart" in used_data_type:
        if not os.path.isfile(directory / "dart_waves.json"):
            raise FileNotFoundError(
                errno.ENOENT, os.strerror(errno.ENOENT), "dart_waves.json"
            )
        with open(directory / "dart_waves.json") as dw:
            traces_info = json.load(dw)
        traces_info = get_outputs.get_data_dict(
            traces_info,
            obs_file="waveforms_dart.txt",
            syn_file="synthetics_dart.txt",
            directory=directory,
        )
        values = [["dart"]]
        for components in values:
            plot_waveform_fits(
                traces_info,
                components,
                "dart",
                start_margin=0,
                plot_directory=directory,
            )
    if "body" in used_data_type:
        if not os.path.isfile(directory / "tele_waves.json"):
            raise FileNotFoundError(
                errno.ENOENT, os.strerror(errno.ENOENT), "tele_waves.json"
            )
        with open(directory / "tele_waves.json") as tw:
            traces_info = json.load(tw)
        traces_info = get_outputs.get_data_dict(
            traces_info,
            syn_file="synthetics_body.txt",
            directory=directory,
        )
        if event:
            traces_info = select_waveforms_event(
                traces_info,
                event,
            )
        values = [["BHZ"], ["BHT"]]
        for components in values:
            plot_waveform_fits(
                traces_info,
                components,
                "body",
                event=event,
                plot_directory=directory,
            )
    if "surf" in used_data_type:
        if not os.path.isfile(directory / "surf_waves.json"):
            raise FileNotFoundError(
                errno.ENOENT, os.strerror(errno.ENOENT), "surf_waves.json"
            )
        with open(directory / "surf_waves.json") as sw:
            traces_info = json.load(sw)
        traces_info = get_outputs.get_data_dict(
            traces_info,
            syn_file="synthetics_surf.txt",
            margin=0,
            directory=directory,
        )
        if event:
            traces_info = select_waveforms_event(
                traces_info,
                event,
            )
        values = [["BHZ"], ["BHT"]]
        for components in values:
            plot_waveform_fits(
                traces_info,
                components,
                "surf",
                event=event,
                plot_directory=directory,
            )
    if "strong" in used_data_type:
        if not os.path.isfile(directory / "strong_motion_waves.json"):
            raise FileNotFoundError(
                errno.ENOENT, os.strerror(errno.ENOENT), "strong_motion_waves.json"
            )
        with open(directory / "strong_motion_waves.json") as smw:
            traces_info = json.load(smw)
        traces_info = get_outputs.get_data_dict(
            traces_info,
            syn_file="synthetics_strong.txt",
            directory=directory,
        )
        if event:
            traces_info = select_waveforms_event(
                traces_info,
                event,
            )
        values = [["HLZ", "HNZ"], ["HLE", "HNE"], ["HLN", "HNN"]]
        for components in values:
            plot_waveform_fits(
                traces_info,
                components,
                "strong",
                event=event,
                plot_directory=directory,
            )
    print(1)
    if "cgps" in used_data_type:
        if not os.path.isfile(directory / "cgps_waves.json"):
            raise FileNotFoundError(
                errno.ENOENT, os.strerror(errno.ENOENT), "cgps_waves.json"
            )
        print(2)
        with open(directory / "cgps_waves.json") as cw:
            traces_info = json.load(cw)
        print(3)
        traces_info = get_outputs.get_data_dict(
            traces_info,
            syn_file="synthetics_cgps.txt",
            directory=directory,
        )
        print(4)
        if event:
            traces_info = select_waveforms_event(
                traces_info,
                event,
            )
        print(5)
        values = [["LXZ", "LHZ", "LYZ"], ["LXE", "LHE", "LYE"], ["LXN", "LHN", "LYN"]]
        print(6)
        for components in values:
            plot_waveform_fits(
                traces_info,
                components,
                "cgps",
                event=event,
                plot_directory=directory,
            )
    return


def _plot_vel_model(
    velmodel: dict,
    point_sources: list,
    directory: Union[pathlib.Path, str] = pathlib.Path(),
):
    """Plot the seismic velocity model as a function of depth

    :param velmodel: The velocity model
    :type velmodel: dict
    :param point_sources: The point source locations
    :type point_sources: list
    :param directory: Where to read/write file(s), defaults to pathlib.Path()
    :type directory: Union[pathlib.Path,str], optional
    """
    directory = pathlib.Path(directory)
    max_depth = [
        max(ps_segment[:, :, :, :, 2].flatten()) for ps_segment in point_sources
    ]
    max_depth = max(max_depth)
    p_vel = np.array(velmodel["p_vel"]).astype(float)
    sh_vel = np.array(velmodel["s_vel"]).astype(float)
    thick = np.array(velmodel["thick"]).astype(float)

    depths = np.zeros(len(thick) + 1)

    depths[1:] = np.cumsum(thick)
    depths = np.array([depth for depth in depths if depth < 70])
    depths = np.append([depths], [70])  # [max_depth])
    plt.plot((p_vel[0], p_vel[0]), (depths[0], depths[1]), "b-", label="P")
    plt.plot((sh_vel[0], sh_vel[0]), (depths[0], depths[1]), "r-", label="SH")
    j = len(depths) - 3  # 2
    for i in range(j):
        plt.plot((p_vel[i], p_vel[i]), (depths[i], depths[i + 1]), "b-")
        plt.plot((p_vel[i], p_vel[i + 1]), (depths[i + 1], depths[i + 1]), "b-")
        plt.plot((sh_vel[i], sh_vel[i]), (depths[i], depths[i + 1]), "r-")
        plt.plot((sh_vel[i], sh_vel[i + 1]), (depths[i + 1], depths[i + 1]), "r-")

    plt.plot((p_vel[j], p_vel[j]), (depths[j], depths[j + 1]), "b-")
    plt.plot((sh_vel[j], sh_vel[j]), (depths[j], depths[j + 1]), "r-")

    plt.title("Crust model for north of Chile")  #'Body wave velocity model')
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
    :param solution: The kinematic solution read from Solucion.txt
    :type solution: dict
    :param directory: The location where to write plots, defaults to pathlib.Path()
    :type directory: Union[pathlib.Path, str], optional
    """
    directory = pathlib.Path(directory)
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
        indexes = np.where(slip_seg < 0.1 * max_slip)  # type:ignore
        rupt_time_seg[indexes] = 0
        fig = plt.figure(figsize=(15, 8))
        ax = fig.add_subplot(111)
        ax.set_ylabel(y_label)
        ax.set_xlabel(x_label)
        fig.subplots_adjust(right=0.75)
        if i_segment == 0:
            ax.plot(0, 0, "w*", ms=20)
        ax, im = __several_axes(
            rupt_time_seg, segment, ps_seg, ax, max_val=max_rupt_time  # type: ignore
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
    :param solution: The kinematic solution read from Solucion.txt
    :type solution: dict
    :param directory: The location where to write plots, defaults to pathlib.Path()
    :type directory: Union[pathlib.Path, str], optional
    """
    directory = pathlib.Path(directory)
    rise_time = solution["trise"]
    max_trise = [np.max(trise_seg.flatten()) for trise_seg in rise_time]
    max_trise = np.max(max_trise)
    fall_time = solution["tfall"]
    max_tfall = [np.max(tfall_seg.flatten()) for tfall_seg in fall_time]
    max_tfall = np.max(max_tfall)
    slip = solution["slip"]
    max_slip = [np.max(slip_seg.flatten()) for slip_seg in slip]
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
        fig.subplots_adjust(bottom=0.15)
        a1: Axes = axes[0]  # type: ignore
        a2: Axes = axes[1]  # type: ignore
        a1.set_ylabel(y_label)
        a1.set_xlabel(x_label)
        if i_segment == 0:
            a1.plot(0, 0, "w*", ms=20)
        a1, im = __several_axes(
            trise_seg,
            segment,
            ps_seg,
            a1,
            max_val=max_trise,  # type:ignore
            autosize=False,
        )
        a2.set_ylabel(y_label)
        a2.set_xlabel(x_label)
        if i_segment == 0:
            a2.plot(0, 0, "w*", ms=20)
        a2, im = __several_axes(
            tfall_seg,
            segment,
            ps_seg,
            a2,
            max_val=max_tfall,  # type:ignore
            autosize=False,
        )
        cbar_ax = fig.add_axes((0.1, 0.05, 0.8, 0.05))
        cb = fig.colorbar(im, cax=cbar_ax, orientation="horizontal")
        cb.set_label("Rise_time (s)")
        plt.savefig(
            directory / "RiseTime_plane{}.png".format(i_segment), bbox_inches="tight"
        )
        plt.close()
    return


def PlotSlipDistribution(
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
    :param solution: The kinematic solution read from Solucion.txt
    :type solution: dict
    :param directory: The location where to write plots, defaults to pathlib.Path()
    :type directory: Union[pathlib.Path, str], optional
    """
    directory = pathlib.Path(directory)
    slip = solution["slip"]
    rake = solution["rake"]
    max_slip = [np.max(slip_seg.flatten()) for slip_seg in slip]
    max_slip = np.max(max_slip)
    x_label = "Distance along strike $(km)$"
    y_label = "Distance along dip $(km)$"
    for i_segment, (segment, slip_seg, rake_seg, ps_seg) in enumerate(
        zip(segments, slip, rake, point_sources)
    ):
        max_slip_seg = np.max(slip_seg.flatten())
        u = slip_seg * np.cos(rake_seg * np.pi / 180.0) / max_slip_seg
        v = slip_seg * np.sin(rake_seg * np.pi / 180.0) / max_slip_seg
        #
        # Plot the slip distribution
        #
        fig = plt.figure(figsize=(15, 8))
        ax = fig.add_subplot(111)
        ax.set_ylabel(y_label)
        ax.set_xlabel(x_label)
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
        ax.quiver(x, y, u, v, scale=15.0, width=0.003)
        if i_segment == 0:
            ax.plot(0, 0, "w*", ms=20)
        ax, im = __several_axes(
            slip_seg, segment, ps_seg, ax, max_val=max_slip  # type:ignore
        )
        cbar_ax = fig.add_axes((0.85, 0.15, 0.05, 0.7))
        cb = fig.colorbar(im, cax=cbar_ax)
        cb.set_label("Slip (cm)")
        plt.savefig(
            directory / "SlipDist_plane{}.png".format(i_segment), bbox_inches="tight"
        )
        plt.close()
    return


def PlotSlipDist_Compare(
    segments: dict,
    point_sources: list,
    input_model: dict,
    solution: dict,
    directory: Union[pathlib.Path, str] = pathlib.Path(),
):
    """Plot slip distribution based on the FFM solution model and compare
    inverted model to an input model

    :param segments: The segment properties
    :type segments: dict
    :param point_sources: The point source locations
    :type point_sources: list
    :param input_model: The model
    :type input_model: dict
    :param solution: The kinematic solution read from Solucion.txt
    :type solution: dict
    :param directory: The location where to write plots, defaults to pathlib.Path()
    :type directory: Union[pathlib.Path,str], optional
    """
    directory = pathlib.Path(directory)
    slip = solution["slip"]
    rake = solution["rake"]
    slip2 = input_model["slip"]
    rake2 = input_model["rake"]
    max_slip = [np.max(slip_seg.flatten()) for slip_seg in slip]
    max_slip = np.max(max_slip)
    max_slip2 = [np.max(slip_seg2.flatten()) for slip_seg2 in slip2]
    max_slip2 = np.max(max_slip2)
    max_slip = np.maximum(max_slip, max_slip2)  # type:ignore
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
            slip_seg,
            segment,
            ps_seg,
            ax0,
            max_val=max_slip,  # type:ignore
            autosize=False,
        )
        ax1, im = __several_axes(
            slip_seg2,
            segment,
            ps_seg,
            ax1,
            max_val=max_slip,  # type:ignore
            autosize=False,
        )
        ax0.set_title("Inverted model", fontsize=20)
        ax1.set_title("Original model", fontsize=20)
        cbar_ax = fig.add_axes((0.85, 0.15, 0.05, 0.7))
        cb = fig.colorbar(im, cax=cbar_ax)
        cb.set_label("Slip (cm)")
        plt.savefig(
            directory / "SlipDist_Compare_plane{}.png".format(i_segment),
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
    max_slip: Optional[float] = None,
    event: Optional[str] = None,
    directory: Union[pathlib.Path, str] = pathlib.Path(),
):
    """Plot slip map

    :param tensor_info: The tensor information
    :type tensor_info: dict
    :param segments: The segment properties
    :type segments: List[dict]
    :param point_sources: The point source locations
    :type point_sources: list
    :param solution: The kinematic solution read from Solucion.txt
    :type solution: dict
    :param default_dirs: The location of default directories
    :type default_dirs: dict
    :param files_str: The stations file properties, defaults to None
    :type files_str:Optional[dict], optional
    :param stations_gps: The gps stations description, defaults to None
    :type stations_gps:  Optional[zip], optional
    :param max_slip: Specify maximum slip, defaults to None
    :type max_slip: Optional[float], optional
    :param event: The event, defaults to None
    :type event: Optional[str], optional
    :param directory: the location where to write plots, defaults to pathlib.Path()
    :type directory: Union[pathlib.Path, str], optional
    """
    directory = pathlib.Path(directory)
    plane_info = segments[0]
    (
        stk_subfaults,
        dip_subfaults,
        delta_strike,
        delta_dip,
        hyp_stk,
        hyp_dip,
    ) = pl_mng.__unpack_plane_data(plane_info)
    slip = solution["slip"]

    point_sources2 = point_sources.copy()
    segments2 = segments.copy()
    if event is not None:
        zipped = zip(segments, slip)
        slip = [slip_seg for segment, slip_seg in zipped if segment["event"] == event]
        zipped = zip(segments, point_sources)
        point_sources2 = [
            ps_seg
            for segment, ps_seg in zipped
            if segment["event"] == event  # type:ignore
        ]
        segments2 = [segment for segment in segments if segment["event"] == event]
    #
    # accurate plot coordinates
    #
    segments_lats, segments_lons = __redefine_lat_lon(segments2, point_sources2)
    min_lats = [np.min(segment_lat.flatten()) for segment_lat in segments_lats]
    max_lats = [np.max(segment_lat.flatten()) for segment_lat in segments_lats]
    min_lons = [np.min(segment_lon.flatten()) for segment_lon in segments_lons]
    max_lons = [np.max(segment_lon.flatten()) for segment_lon in segments_lons]
    min_lat = np.min(min_lats)  # - 0.5
    max_lat = np.max(max_lats)  # + 0.5
    min_lon = np.min(min_lons)  # - 0.5
    max_lon = np.max(max_lons)  # + 0.5

    margin = 1.3 * (stk_subfaults * delta_strike) / 111.19
    lat0 = tensor_info["lat"]
    lon0 = tensor_info["lon"]
    tectonic = "{}.shp".format(default_dirs["trench_graphics"])
    dictn = {"projection": ccrs.PlateCarree(), "facecolor": "#eafff5"}

    fig, ax = plt.subplots(1, 1, figsize=(15, 15), subplot_kw=dictn)
    fig.subplots_adjust(hspace=0, wspace=0, top=0.9, bottom=0.1, right=0.8)
    tectonic = cf.ShapelyFeature(
        shpreader.Reader(tectonic).geometries(),
        ccrs.PlateCarree(),
        edgecolor="red",
        facecolor=(198 / 255.0, 236 / 255.0, 253 / 255.0),
    )
    shpfilename = shpreader.natural_earth(
        resolution="10m", category="cultural", name="admin_0_countries"
    )
    countries = cf.ShapelyFeature(
        shpreader.Reader(shpfilename).geometries(),
        ccrs.PlateCarree(),
        edgecolor="black",
        facecolor="lightgray",
    )

    if files_str is not None:
        for file in files_str:
            name = file["name"]
            latp, lonp = file["location"]
            min_lat = min(min_lat, latp)
            max_lat = max(max_lat, latp)
            min_lon = min(min_lon, lonp)
            max_lon = max(max_lon, lonp)
            distance = max(np.abs(latp - lat0), np.abs(lonp - lon0))
            margin = max(margin, 1.2 * distance)
            ax.plot(
                lonp, latp, "wo", markersize=10, transform=ccrs.PlateCarree(), zorder=4
            )
            ax.text(
                lonp + 0.1,
                latp + 0.1,
                "{}".format(name),
                transform=ccrs.PlateCarree(),
                zorder=4,
            )
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
            distance = max(np.abs(sta_lat - lat0), np.abs(sta_lon - lon0))
            margin = max(margin, 1.2 * distance)
        max_obs = np.max(max_obs)  # type:ignore
        plt.text(
            lon0 + margin - 2,
            lat0 + margin - 0.25,
            "{:.2f} cm".format(max_obs),  # type:ignore
            transform=ccrs.PlateCarree(),
        )
        plt.text(
            lon0 + margin - 2,
            lat0 + margin - 0.45,
            "{:.2f} cm".format(max_obs),  # type:ignore
            transform=ccrs.PlateCarree(),
        )
        for name, sta_lat, sta_lon, obs, syn, error in stations_gps2:
            plt.plot(
                sta_lon, sta_lat, "ks", transform=ccrs.PlateCarree(), markersize=14
            )
            gps_z, gps_n, gps_e = syn
            east_west = float(gps_e) / max_obs
            north_south = float(gps_n) / max_obs
            plt.arrow(
                sta_lon,
                sta_lat,
                east_west,  # type:ignore
                north_south,  # type:ignore
                color="r",
                zorder=3,
                linewidth=2,
                head_width=0.05,
                head_length=0.05,
                transform=ccrs.PlateCarree(),
            )
            up_down = float(gps_z) / max_obs
            plt.arrow(
                sta_lon,
                sta_lat,
                0.0,
                up_down,  # type:ignore
                color="r",
                zorder=3,
                linewidth=2,
                head_width=0.05,
                head_length=0.05,
                transform=ccrs.PlateCarree(),
            )
            gps_z, gps_n, gps_e = obs
            east_west = float(gps_e) / max_obs
            north_south = float(gps_n) / max_obs
            plt.arrow(
                sta_lon,
                sta_lat,
                east_west,  # type:ignore
                north_south,  # type:ignore
                zorder=3,
                linewidth=2,
                head_width=0.05,
                head_length=0.05,
                transform=ccrs.PlateCarree(),
            )
            up_down = float(gps_z) / max_obs
            plt.arrow(
                sta_lon,
                sta_lat,
                0.0,
                up_down,  # type:ignore
                zorder=3,
                linewidth=2,
                head_width=0.05,
                head_length=0.05,
                transform=ccrs.PlateCarree(),
            )
            plt.text(
                sta_lon + 0.1,
                sta_lat + 0.1,
                "{}".format(name),
                transform=ccrs.PlateCarree(),
            )
            err_z, err_n, err_e = error
            width = float(err_e) / max_obs  # / 100
            height = float(err_n) / max_obs  # / 100
            ellipse = patches.Ellipse(
                (sta_lon + east_west, sta_lat + north_south),
                width,  # type:ignore
                height,  # type:ignore
                zorder=4,
                color="k",
                linewidth=10,
                transform=ccrs.PlateCarree(),
            )
            plt.gca().add_patch(ellipse)
        plt.arrow(
            lon0 + margin - 0.2,
            lat0 + margin - 0.2,
            -1,
            0,
            color="r",
            zorder=3,
            linewidth=2,
            head_width=0.05,
            head_length=0.05,
            transform=ccrs.PlateCarree(),
        )
        plt.arrow(
            lon0 + margin - 0.2,
            lat0 + margin - 0.4,
            -1,
            0,
            color="k",
            zorder=3,
            linewidth=2,
            head_width=0.05,
            head_length=0.05,
            transform=ccrs.PlateCarree(),
        )
    max_slip = (
        max([np.amax(slip_fault) for slip_fault in slip]) if not max_slip else max_slip
    )
    margins = [min_lon - 0.5, max_lon + 0.5, min_lat - 0.5, max_lat + 0.5]
    ax = set_map_cartopy(ax, margins, tectonic=tectonic, countries=countries)
    ax.plot(lon0, lat0, "w*", markersize=15, transform=ccrs.PlateCarree(), zorder=4)
    #
    # plot slip map
    #
    ax, cs = plot_map(
        ax,
        segments_lats,  # type:ignore
        segments_lons,  # type:ignore
        slip,
        max_val=max_slip,
        transform=dictn["projection"],
    )
    cbar_ax = fig.add_axes((0.85, 0.1, 0.05, 0.8))
    cbar = plt.colorbar(cs, cax=cbar_ax)
    cbar.set_label("Slip (cm)", size=15)
    cbar.ax.yaxis.set_ticks_position("left")
    ax.set_aspect("auto", adjustable=None)
    plot_name = "Map"
    if event is not None:
        plot_name = "{}_event{}".format(plot_name, event)
    plt.savefig(directory / plot_name, bbox_inches="tight")
    plt.close()
    return


def PlotInsar(
    tensor_info: dict,
    segments: List[dict],
    point_sources: list,
    default_dirs: dict,
    insar_points: List[dict],
    los: str = "ascending",
    directory: Union[pathlib.Path, str] = pathlib.Path(),
):
    """Compare insar data with insar produced by the inverted earthquake model

    :param tensor_info: The tensor information
    :type tensor_info: dict
    :param segments: The segment properties
    :type segments: List[dict]
    :param point_sources: The point source locations
    :type point_sources: list
    :param default_dirs: The location of default directories
    :type default_dirs: dict
    :param insar_points: List of insar data
    :type insar_points: List[dict]
    :param los: The direction of the path, defaults to "ascending"
    :type los: str, optional
    :param directory: The location where plots should be written, defaults to pathlib.Path()
    :type directory: Union[pathlib.Path, str], optional
    """
    directory = pathlib.Path(directory)
    plane_info = segments[0]
    (
        stk_subfaults,
        dip_subfaults,
        delta_strike,
        delta_dip,
        hyp_stk,
        hyp_dip,
    ) = pl_mng.__unpack_plane_data(plane_info)
    #
    # accurate plot coordinates
    #
    segments_lats, segments_lons = __redefine_lat_lon(segments, point_sources)
    min_lats = [np.min(segment_lat.flatten()) for segment_lat in segments_lats]
    max_lats = [np.max(segment_lat.flatten()) for segment_lat in segments_lats]
    min_lons = [np.min(segment_lon.flatten()) for segment_lon in segments_lons]
    max_lons = [np.max(segment_lon.flatten()) for segment_lon in segments_lons]
    min_lat = np.min(min_lats)
    max_lat = np.max(max_lats)
    min_lon = np.min(min_lons)
    max_lon = np.max(max_lons)

    margin = 1.3 * (stk_subfaults * delta_strike) / 111.19
    lat0 = tensor_info["lat"]
    lon0 = tensor_info["lon"]
    tectonic = "{}.shp".format(default_dirs["trench_graphics"])
    dictn = {"projection": ccrs.PlateCarree(), "facecolor": "#eafff5"}

    fig, axes = plt.subplots(2, 3, figsize=(30, 15), subplot_kw=dictn)
    tectonic = cf.ShapelyFeature(
        shpreader.Reader(tectonic).geometries(),
        ccrs.PlateCarree(),
        edgecolor="red",
        facecolor=(198 / 255.0, 236 / 255.0, 253 / 255.0),
    )
    shpfilename = shpreader.natural_earth(
        resolution="10m", category="cultural", name="admin_0_countries"
    )
    countries = cf.ShapelyFeature(
        shpreader.Reader(shpfilename).geometries(),
        ccrs.PlateCarree(),
        edgecolor="black",
        facecolor="lightgray",
    )

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

    margins = [min_lon - 0.5, max_lon + 0.5, min_lat - 0.5, max_lat + 0.5]

    values = [observed, synthetic, diffs, obs_no_ramp, syn_no_ramp, ramp]
    titles = [
        "Observed",
        "Synthetic",
        "Misfit",
        "Observed - Ramp",
        "Synthetic - Ramp",
        "Ramp",
    ]
    labels = [
        "Observed LOS (m)",
        "Modeled LOS (m)",
        "Residual (m)",
        "Observed LOS (m)",
        "Modeled LOS (m)",
        "Modeled Ramp (m)",
    ]
    rows = [0, 0, 0, 1, 1, 1]
    cols = [0, 1, 2, 0, 1, 2]
    zipped = zip(values, titles, labels, rows, cols)
    for value, title, label, row, col in zipped:
        axes[row][col].set_title(title, fontdict={"fontsize": 20})  # type: ignore
        max_abs = np.max(np.abs(value))
        vmin = -max_abs - 0.001  # if not title == 'Misfit' else -40
        vmax = max_abs + 0.001  # if not title == 'Misfit' else 40
        norm = colors.TwoSlopeNorm(vmin=vmin, vcenter=0, vmax=vmax)
        cs = axes[row][col].scatter(  # type: ignore
            lons,
            lats,
            zorder=4,
            c=value,
            cmap="bwr",
            norm=norm,
            transform=ccrs.PlateCarree(),
        )
        ax = set_map_cartopy(
            axes[row][col], margins, tectonic=tectonic, countries=countries  # type: ignore
        )
        ax.plot(lon0, lat0, "y*", markersize=15, transform=ccrs.PlateCarree(), zorder=4)
        ax = plot_borders(
            ax,
            segments_lats,  # type:ignore
            segments_lons,  # type:ignore
            transform=dictn["projection"],
        )
        fig.colorbar(cs, ax=ax, orientation="horizontal")
        ax.set_aspect("auto", adjustable=None)

    fig.tight_layout()
    plt.savefig(directory / "Insar_{}_fit.png".format(los), bbox_inches="tight")
    plt.close()
    return


def PlotComparisonMap(
    tensor_info: dict,
    segments: List[dict],
    point_sources: list,
    input_model: dict,
    solution: dict,
    directory: Union[pathlib.Path, str] = pathlib.Path(),
):
    """Plot slip map and compare with map from input model

    :param tensor_info: The tensor information
    :type tensor_info: dict
    :param segments: The segment properties
    :type segments: List[dict]
    :param point_sources: The point source locations
    :type point_sources: list
    :param input_model: the model
    :type input_model: dict
    :param solution: The kinematic solution read from Solucion.txt
    :type solution: dict
    :param directory: Where to write the plots, defaults to pathlib.Path()
    :type directory: Union[pathlib.Path,str], optional
    """
    directory = pathlib.Path(directory)
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
    segments_lats, segments_lons = __redefine_lat_lon(segments, point_sources)

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
    ax1, cs1 = plot_map(
        ax1,
        segments_lats,  # type:ignore
        segments_lons,  # type:ignore
        slip,
        max_val=max_slip,
        transform=dictn["projection"],
    )
    ax2, cs2 = plot_map(
        ax2,
        segments_lats,  # type:ignore
        segments_lons,  # type:ignore
        input_slip,
        max_val=max_slip,
        transform=dictn["projection"],
    )
    fig.subplots_adjust(bottom=0.15)
    cbar_ax = fig.add_axes((0.1, 0.05, 0.8, 0.05))
    cb = fig.colorbar(cs2, cax=cbar_ax, orientation="horizontal")
    cb.set_label("Slip (cm)")
    plt.savefig(directory / "Comparison.png", bbox_inches="tight")
    plt.close()
    return


def __redefine_lat_lon(segments: List[dict], point_sources: list) -> Tuple[list, list]:
    """Redefine the lat/lon

    :param segments: The segments
    :type segments: List[dict]
    :param point_sources: The point sources
    :type point_sources: list
    :return: segments latitudes, segments longitudes
    :rtype: Tuple[list, list]
    """
    segments_lats: list = [[]] * len(segments)
    segments_lons: list = [[]] * len(segments)
    for i, point_sources_seg in enumerate(point_sources):
        lat = point_sources_seg[:, :, :, :, 0]
        lon = point_sources_seg[:, :, :, :, 1]
        ny, nx, a, b = lat.shape
        new_lat = np.zeros((ny + 1, nx + 1))
        new_lon = np.zeros((ny + 1, nx + 1))
        for j in range(ny):
            for k in range(nx):
                new_lat[j, k] = lat[j, k, 0, 0]
                new_lon[j, k] = lon[j, k, 0, 0]
        for k in range(nx):
            new_lat[-1, k] = lat[-1, k, -1, 0]
            new_lon[-1, k] = lon[-1, k, -1, 0]
        for j in range(ny):
            new_lat[j, -1] = lat[j, -1, 0, -1]
            new_lon[j, -1] = lon[j, -1, 0, -1]
        new_lat[-1, -1] = lat[-1, -1, -1, -1]
        new_lon[-1, -1] = lon[-1, -1, -1, -1]
        segments_lats[i] = new_lat
        segments_lons[i] = new_lon
    return segments_lats, segments_lons


def plot_moment_rate_function(
    segments_data: dict,
    shear: list,
    point_sources: list,
    event: Optional[str] = None,
    directory: Union[pathlib.Path, str] = pathlib.Path(),
):
    """Plot moment rate function

    :param segments_data: The segment properties
    :type segments_data: dict
    :param shear: The shear moduli
    :type shear: list
    :param point_sources: The point source locations
    :type point_sources: list
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
    seismic_moment = 0

    shear2 = shear.copy()
    point_sources2 = point_sources.copy()
    segments2 = segments.copy()
    if event is not None:
        zipped = zip(segments, slip)
        slip = [slip_seg for segment, slip_seg in zipped if segment["event"] == event]
        zipped = zip(segments, trup)
        trup = [trup_seg for segment, trup_seg in zipped if segment["event"] == event]
        zipped = zip(segments, tl)
        tl = [trise_seg for segment, trise_seg in zipped if segment["event"] == event]
        zipped = zip(segments, tr)
        tr = [tfall_seg for segment, tfall_seg in zipped if segment["event"] == event]
        zipped = zip(segments, shear)
        shear2 = [
            shear_seg for segment, shear_seg in zipped if segment["event"] == event
        ]
        zipped = zip(segments, point_sources)
        point_sources2 = [
            ps_seg
            for segment, ps_seg in zipped
            if segment["event"] == event  # type:ignore
        ]
        segments2 = [segment for segment in segments if segment["event"] == event]

    for (
        segment,
        slip_seg,
        trup_seg,
        trise_seg,
        tfall_seg,
        shear_seg,
        point_sources_seg,
    ) in zip(segments2, slip, trup, tl, tr, shear2, point_sources2):
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
            mr[i] = mr[i] + moment_rate[i] * (delta_strike * 1000) * (delta_dip * 1000)

    time = np.arange(nmax) * dt  # type:ignore
    with open(directory / "STF.txt", "w") as outf:
        outf.write("dt: {}\n".format(dt))
        outf.write("Time[s]     Moment_Rate [Nm]\n")
        for t, val in zip(time, mr):  # type:ignore
            outf.write("{:8.2f}:   {:8.4e}\n".format(t, val))

    seismic_moment = np.trapz(mr, dx=0.01)
    magnitude = 2.0 * (np.log10(seismic_moment * 10**7) - 16.1) / 3.0
    plt.xlabel("Time $(s)$")
    plt.ylabel("Moment rate $(Nm/s)$")
    plt.text(
        0.5 * max(time),  # type:ignore
        0.95 * max(mr),
        "$M_0$: {:.2E} $Nm$".format(seismic_moment),
    )
    plt.text(
        0.5 * max(time),  # type:ignore
        0.85 * max(mr),
        "$M_w$: {:.2f}".format(magnitude),
    )
    plt.grid(visible=True)
    plt.fill_between(time, mr)
    plot_name = "MomentRate"
    if event is not None:
        plot_name = "{}_event{}".format(plot_name, event)
    plt.savefig(directory / plot_name, bbox_inches="tight")
    plt.close()
    return


def _PlotSnapshotSlip(
    segments: dict, solution: dict, directory: Union[pathlib.Path, str] = pathlib.Path()
):
    """Plot snapshots of the rupture process

    :param segments: The segment properties
    :type segments: dict
    :param solution: The kinematic solution read from Solucion.txt
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
        time = i * step * dt
        srate, cslip, broken = __rupture_process(  # type:ignore
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
    segment = segments[0]
    #
    # Get the focal mechanism
    #
    strike = segment["strike"]
    dip = segment["dip"]
    rake = segment["rake"]
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
    max_val: Optional[float] = None,
    autosize: bool = True,
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
    :param max_val: Specify a maximum value, defaults to None
    :type max_val: Optional[float], optional
    :param autosize: Autosize the plot, defaults to True
    :type autosize: bool, optional
    :return: Return axes
    :rtype: Tuple[plt.Axes, AxesImage]
    """
    (
        stk_subfaults,
        dip_subfaults,
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
    im = ax.imshow(
        data,
        cmap="jet",
        origin="lower",
        vmax=max_val,
        aspect="auto",
        extent=(min_strike, max_strike, min_dist, max_dist),
    )
    ax2 = ax.twinx()
    ax2.set_xlim((min_strike, max_strike))
    ax2.set_ylim((min_depth, max_depth))
    ax.set(adjustable="datalim")
    ax2.set(adjustable="datalim")
    if autosize:
        ax.figure.set_size_inches(  # type:ignore
            4 * stk_subfaults * delta_strike / dip_subfaults / delta_dip, 4
        )
        ax2.figure.set_size_inches(  # type:ignore
            4 * stk_subfaults * delta_strike / dip_subfaults / delta_dip, 4
        )
    ax2.set_ylabel("Depth $(km)$")
    ax.invert_yaxis()
    ax2.invert_yaxis()
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
                srate[i, j] = convolve[index] * slip[i, j]  # srmax[iy, ix]
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
