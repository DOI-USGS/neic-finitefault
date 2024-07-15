#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Jul 11 11:21:45 2022

@author: degoldberg
"""
import collections
import pathlib
from glob import glob
from typing import List, Optional, Tuple, Union

import cartopy  # type: ignore
import cartopy.crs as ccrs  # type: ignore
import numpy as np
from cartopy.crs import Projection  # type: ignore
from cartopy.feature import Feature  # type: ignore
from cartopy.mpl.geoaxes import GeoAxes  # type: ignore
from matplotlib import colormaps  # type: ignore
from matplotlib import pyplot as plt  # type: ignore
from matplotlib.colors import ListedColormap  # type: ignore
from matplotlib.patches import Rectangle  # type: ignore

import wasp.plane_management as pl_mng
from wasp.plot_graphic_NEIC import __redefine_lat_lon
from wasp.plot_maps_NEIC import plot_map

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
#####


def _write_KML(
    segments: List[dict],
    point_sources: list,
    evID: str,
    margins: list,
    directory: Union[pathlib.Path, str] = pathlib.Path(),
):
    """Write the kml file

    :param segments: The segment properties
    :type segments: List[dict]
    :param point_sources: The point source locations
    :type point_sources: list
    :param evID: The event id/name
    :type evID: str
    :param margins: The extent of the map
    :type margins: list
    :param directory: The directory to read/write to, defaults to pathlib.Path()
    :type directory: Union[pathlib.Path, str], optional
    """
    print("Writing KML file...")
    directory = pathlib.Path(directory)

    kmlfile = str(evID) + ".kml"
    with open(directory / kmlfile, "w") as kml:
        kml.write('<?xml version="1.0" encoding="UTF-8"?>\n')
        kml.write('<kml xmlns="http://earth.google.com/kml/2.2">\n')
        kml.write("<Document>\n")
        kml.write("	<name>" + str(evID) + "</name>\n")
        kml.write('	<StyleMap id="msn_ylw-pushpin">\n')
        kml.write("	<Pair>\n")
        kml.write("	<key>normal</key>\n")
        kml.write("	<styleUrl>#sn_ylw-pushpin</styleUrl>\n")
        kml.write("		</Pair>\n")
        kml.write("		<Pair>\n")
        kml.write("		<key>highlight</key>\n")
        kml.write("	<styleUrl>#sh_ylw-pushpin</styleUrl>\n")
        kml.write("		</Pair>\n")
        kml.write("	</StyleMap>\n")
        kml.write('	<Style id="sn_ylw-pushpin">\n')
        kml.write("		<IconStyle>\n")
        kml.write("			<scale>1.1</scale>\n")
        kml.write("			<Icon>\n")
        kml.write(
            "				<href>http://maps.google.com/mapfiles/kml/pushpin/ylw-pushpin.png</href>\n"
        )
        kml.write("			</Icon>\n")
        kml.write('			<hotSpot x="20" y="2" xunits="pixels" yunits="pixels"/>\n')
        kml.write("		</IconStyle>\n")
        kml.write("		<LineStyle>\n")
        kml.write("		<color>ff00ffff</color>\n")
        kml.write("	<width>5</width>\n")
        kml.write("		</LineStyle>\n")
        kml.write("	</Style>\n")
        kml.write('	<Style id="sh_ylw-pushpin">\n')
        kml.write("		<IconStyle>\n")
        kml.write("		<scale>1.3</scale>\n")
        kml.write("	<Icon>\n")
        kml.write(
            "				<href>http://maps.google.com/mapfiles/kml/pushpin/ylw-pushpin.png</href>\n"
        )
        kml.write("		</Icon>\n")
        kml.write('			<hotSpot x="20" y="2" xunits="pixels" yunits="pixels"/>\n')
        kml.write("		</IconStyle>\n")
        kml.write("	<LineStyle>\n")
        kml.write("	<color>ff00ffff</color>\n")
        kml.write("			<width>5</width>\n")
        kml.write("	</LineStyle>\n")
        kml.write("	</Style>\n")

        segments_lats, segments_lons, segments_deps = __redefine_lat_lon(
            segments, point_sources
        )
        min_lats = [min(segment_lat.flatten()) for segment_lat in segments_lats]
        max_lats = [max(segment_lat.flatten()) for segment_lat in segments_lats]
        min_lons = [min(segment_lon.flatten()) for segment_lon in segments_lons]
        max_lons = [max(segment_lon.flatten()) for segment_lon in segments_lons]

        # find outline of segments #
        segment_number = 0
        for segment in range(len(segments_lats)):
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

            corners = np.c_[corners, cornerA]
            # write outline of segments to kml file: #
            kml.write(" <Placemark>\n")
            kml.write("		<name>Fault Segment " + str(segment_number) + "</name>\n")
            kml.write("		<styleUrl>#msn_ylw-pushpin</styleUrl>\n")
            kml.write("		<LineString>\n")
            kml.write("			<tessellate>1</tessellate>\n")
            kml.write("			<coordinates>\n")
            kml.write(
                " "
                + str(corners[0, 0])
                + ","
                + str(corners[1, 0])
                + ","
                + str(corners[2, 0])
                + "   "
                + str(corners[0, 1])
                + ","
                + str(corners[1, 1])
                + ","
                + str(corners[2, 1])
                + "   "
                + str(corners[0, 2])
                + ","
                + str(corners[1, 2])
                + ","
                + str(corners[2, 2])
                + "   "
                + str(corners[0, 3])
                + ","
                + str(corners[1, 3])
                + ","
                + str(corners[2, 3])
                + "   "
                + str(corners[0, 4])
                + ","
                + str(corners[1, 4])
                + ","
                + str(corners[2, 4])
                + "   </coordinates>\n"
            )
            kml.write("		</LineString>\n")
            kml.write("	</Placemark>\n")
            segment_number += 1

        kml.write("     <GroundOverlay>\n")
        kml.write("             <name>" + str(evID) + " FFM Image Overlay</name>\n")
        kml.write("             <Icon>\n")
        kml.write("                     <href>Map_kml.png</href>\n")
        kml.write("                     <viewBoundScale>0.75</viewBoundScale>\n")
        kml.write("             </Icon>\n")
        kml.write("		<color>bfffffff</color>\n")
        kml.write("             <LatLonBox>\n")
        kml.write("                     <north>" + str(margins[3]) + "</north>\n")
        kml.write("                     <south>" + str(margins[2]) + "</south>\n")
        kml.write("                     <east>" + str(margins[1]) + "</east>\n")
        kml.write("                     <west>" + str(margins[0]) + "</west>\n")
        kml.write("             </LatLonBox>\n")
        kml.write("     </GroundOverlay>\n")

        kml.write(" </Document>\n")
        kml.write(" </kml>\n")
    return


def PlotMap_KML(
    tensor_info: dict,
    segments: List[dict],
    point_sources: list,
    solution: dict,
    default_dirs: dict,
    stations_str: Optional[dict] = None,
    stations_gps: Optional[zip] = None,
    stations_cgps: Optional[str] = None,
    max_slip: Optional[float] = None,
    legend_len: Optional[float] = None,
    scale: Optional[float] = None,
    limits: List[Optional[float]] = [None, None, None, None],
    evID: Optional[str] = None,
    directory: Union[pathlib.Path, str] = pathlib.Path(),
):
    """Write the slip map to KML

    :param tensor_info: The tensor information
    :type tensor_info: dict
    :param segments: The segment properties
    :type segments: List[dict]
    :param point_sources: The point source locations
    :type point_sources: list
    :param solution: The solution read from Solucion.txt
    :type solution: dict
    :param default_dirs: The location of default directories
    :type default_dirs: dict
    :param stations_str: The stations file properties, defaults to None
    :type stations_str:Optional[dict], optional
    :param stations_gps: The gps stations description, defaults to None
    :type stations_gps:  Optional[zip], optional
    :param stations_cgps: The cgps stations description, defaults to None
    :type stations_cgps: Optional[str], optional
    :param max_slip: A specified maximum slip, defaults to None
    :type max_slip: Optional[float], optional
    :param legend_len: The length of the legend, defaults to None
    :type legend_len: Optional[float], optional
    :param scale: The scale, defaults to None
    :type scale: Optional[float], optional
    :param limits: The extent of the map, defaults to [None, None, None, None]
    :type limits: List[Optional[float]], optional
    :param evID: The event name/ID, defaults to None
    :type evID: Optional[str], optional
    :param directory: The directory to read/write to, defaults to pathlib.Path()
    :type directory: Union[pathlib.Path, str], optional
    """
    directory = pathlib.Path(directory)

    print("Creating Slip Map...")
    plane_info = segments[0]
    (
        stk_subfaults,
        dip_subfaults,
        delta_strike,
        delta_dip,
        hyp_stk,
        hyp_dip,
    ) = pl_mng.__unpack_plane_data(plane_info)
    x = np.arange(stk_subfaults) * delta_strike - hyp_stk * delta_strike
    y = np.arange(dip_subfaults) * delta_dip - hyp_dip * delta_dip
    slip = solution["slip"]
    #
    # accurate plot coordinates
    #
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

    margin = 1.3 * (stk_subfaults * delta_strike) / 111.19
    lat0 = tensor_info["lat"]
    lon0 = tensor_info["lon"]
    tectonic: Optional[str] = "{}.shp".format(default_dirs["trench_graphics"])
    dictn = {
        "projection": ccrs.AlbersEqualArea(
            lon0, lat0, standard_parallels=(lat0 - 3.0, lat0 + 3)
        ),
        "transform": ccrs.PlateCarree(),
        "facecolor": "None",
    }

    fig = plt.figure(figsize=(15, 15))
    ax = fig.add_subplot(
        111, projection=dictn["projection"], facecolor=dictn["facecolor"], frameon=False
    )
    ax.spines["geo"].set_linewidth(2)
    gl = ax.gridlines()  # type:ignore
    gl.top_labels = False
    gl.bottom_labels = False
    gl.left_labels = False
    gl.right_labels = False
    fig.subplots_adjust(hspace=0, wspace=0, top=0.9, bottom=0.1, right=0.8)
    tectonic = None
    shpfilename = None
    countries = None

    if stations_str is not None:
        for file in stations_str:
            name = file["name"]
            latp, lonp = file["location"]
            min_lat = min(min_lat, latp)
            max_lat = max(max_lat, latp)
            min_lon = min(min_lon, lonp)
            max_lon = max(max_lon, lonp)
            distance = max(np.abs(latp - lat0), np.abs(lonp - lon0))
            margin = max(margin, 1.2 * distance)
            ax.plot(
                lonp,
                latp,
                "w",
                marker="v",
                markersize=15,
                lw=0.3,
                markeredgecolor="k",
                transform=dictn["transform"],
                zorder=4,
            )
    if stations_cgps is not None:
        for file in stations_cgps:
            name = file["name"]
            latp, lonp = file["location"]
            min_lat = min(min_lat, latp)
            max_lat = max(max_lat, latp)
            min_lon = min(min_lon, lonp)
            max_lon = max(max_lon, lonp)
            distance = max(np.abs(latp - lat0), np.abs(lonp - lon0))
            margin = max(margin, 1.2 * distance)
            ax.plot(
                lonp,
                latp,
                "b",
                marker="^",
                markersize=10,
                lw=0.3,
                markeredgecolor="k",
                transform=dictn["transform"],
                zorder=4,
            )
    if stations_gps is not None:
        max_obs = np.zeros(3)
        stations_gps2: List[List[Union[float, np.ndarray, str]]] = []
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
        for name, sta_lat, sta_lon, obs, syn, error in stations_gps2:  # type:ignore
            if scale == None:
                scale = 2  # bigger number here makes arrows look longer
            scale2 = 2
            gps_z, gps_n, gps_e = syn
            east_west = float(gps_e) / max_obs / (1.0 / scale)  # type:ignore
            north_south = float(gps_n) / max_obs / (1.0 / scale)  # type:ignore
            plt.arrow(
                sta_lon,
                sta_lat,
                east_west,
                north_south,
                facecolor="r",
                zorder=7,
                linewidth=0.5,
                width=0.02 * scale2,
                head_width=0.05 * scale2,
                head_length=0.05 * scale2,
                transform=dictn["transform"],
                edgecolor="k",
            )
            gps_z, gps_n, gps_e = obs
            east_west = float(gps_e) / max_obs / (1.0 / scale)  # type:ignore
            north_south = float(gps_n) / max_obs / (1.0 / scale)  # type:ignore
            plt.arrow(
                sta_lon,
                sta_lat,
                east_west,
                north_south,
                facecolor="k",
                zorder=5,
                linewidth=0.5,
                width=0.02 * scale2,
                head_width=0.05 * scale2,
                head_length=0.05 * scale2,
                transform=dictn["transform"],
            )
        if legend_len == None:
            legend_len = 20
        plt.arrow(
            0.83 * (max_lon - min_lon) + (min_lon + 0.5),
            0.08 * (max_lat - min_lat) + (min_lat - 0.5),
            legend_len / max_obs / (1.0 / scale),  # type:ignore
            0,
            facecolor="k",
            edgecolor="k",
            zorder=10,
            linewidth=0.5,
            width=0.02 * scale2,
            head_width=0.05 * scale2,
            head_length=0.05 * scale2,
            transform=dictn["transform"],
        )
        plt.arrow(
            0.83 * (max_lon - min_lon) + (min_lon + 0.5),
            0.04 * (max_lat - min_lat) + (min_lat - 0.5),
            legend_len / max_obs / (1.0 / scale),  # type:ignore
            0,
            facecolor="r",
            edgecolor="k",
            zorder=10,
            linewidth=0.5,
            width=0.02 * scale2,
            head_width=0.05 * scale2,
            head_length=0.05 * scale2,
            transform=dictn["transform"],
        )
        plt.text(
            0.9,
            0.09,
            f"{int(legend_len)} cm",  # type:ignore
            horizontalalignment="center",
            verticalalignment="center",
            transform=ax.transAxes,
        )
        plt.text(
            0.82,
            0.06,
            "Observed GNSS",
            horizontalalignment="right",
            verticalalignment="center",
            transform=ax.transAxes,
        )
        plt.text(
            0.82,
            0.03,
            "Synthetic GNSS",
            horizontalalignment="right",
            verticalalignment="center",
            transform=ax.transAxes,
        )
    max_slip = (
        max([np.amax(slip_fault) for slip_fault in slip]) if not max_slip else max_slip
    )
    if limits == [None, None, None, None]:
        margins = [min_lon - 0.5, max_lon + 0.5, min_lat - 0.5, max_lat + 0.5]
    else:
        margins = [
            limits[0],
            limits[1],
            limits[2],
            limits[3],
        ]

    ax.add_patch(
        Rectangle(
            (margins[0] + 0.05, margins[2] + 0.05),
            0.975 * (margins[1] - margins[0]),
            0.1 * (margins[3] - margins[2]),
            edgecolor="k",
            facecolor="0.5",
            alpha=0.5,
            transform=dictn["transform"],
            zorder=3,
        )
    )
    ax.set_extent(margins)  # type:ignore
    ax = set_KML_map_cartopy(
        ax,
        margins,
        tectonic=tectonic,
        countries=countries,
        bathymetry=None,
        faults=True,
        aftershocks=True,
        transform=dictn["transform"],
        directory=directory,
    )
    ax.plot(
        lon0,
        lat0,
        "*",
        markersize=15,
        transform=dictn["transform"],
        zorder=5,
        markerfacecolor="None",
        markeredgecolor="k",
        markeredgewidth=1.5,
    )

    # outline the segments in black #
    for segment in range(len(segments_lats)):
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

        ax.plot(
            corners[0, :],
            corners[1, :],
            "k",
            lw=2,
            zorder=4,
            transform=dictn["transform"],
        )
        ax.plot(
            updip[0, :],
            updip[1, :],
            "r",
            lw=1.5,
            zorder=4,
            transform=dictn["transform"],
        )

        # plot multiple segment hypocenters if any
        if "hypocenter" in segments[segment]:
            hyp = segments[segment]["hypocenter"]
            ax.plot(
                hyp["lon"],
                hyp["lat"],
                "*",
                markersize=15,
                transform=dictn["transform"],
                zorder=5,
                markerfacecolor="None",
                markeredgecolor="k",
                markeredgewidth=1.5,
            )
    #
    # plot slip map
    #
    ax, cs = plot_map(
        ax,
        segments_lats,
        segments_lons,
        slip,
        max_val=max_slip,
        transform=dictn["transform"],
    )
    sm = plt.cm.ScalarMappable(
        cmap=slipcpt, norm=plt.Normalize(vmin=0.0, vmax=max_slip / 100.0)
    )

    cb_ax = ax.inset_axes((0.06, 0.03, 0.3, 0.02))
    cbar = plt.colorbar(sm, cax=cb_ax, orientation="horizontal")
    cbar.outline.set_linewidth(3)  # type:ignore
    cbar.set_label("Slip (m)")
    cbar.ax.xaxis.set_ticks_position("top")
    cbar.ax.xaxis.set_label_position("top")

    plt.savefig(directory / "Map_kml.png", dpi=300, bbox_inches="tight")
    plt.savefig(directory / "Map_kml.ps")
    plt.close()

    _write_KML(segments, point_sources, str(evID), margins, directory=directory)
    return


def set_KML_map_cartopy(
    ax: GeoAxes,
    margins: list,
    tectonic: Optional[Feature] = None,
    countries: Optional[Feature] = None,
    bathymetry: Optional[Feature] = None,
    faults: bool = True,
    aftershocks: bool = True,
    transform: Optional[Projection] = None,
    directory: Union[pathlib.Path, str] = pathlib.Path(),
) -> GeoAxes:
    """Setup the cartopy map

    :param ax: The axes of the plot object
    :type ax: GeoAxes
    :param margins: The spatial margins of the plot
                    ([lon_min, lon_max, lat_min, lat_max])
    :type margins: list
    :param tectonic: An optional tectonic feature to add, defaults to None
    :type tectonic: Optional[Feature], optional
    :param countries: An optional country feature to add, defaults to None
    :type countries: Optional[Feature], optional
    :param bathymetry: An optional bathymetry feature, defaults to None
    :type bathymetry: Optional[Feature], optional
    :param faults: Whether to plot faults, defaults to True
    :type faults: bool, optional
    :param aftershocks: Whether to plot aftershocks, defaults to True
    :type aftershocks: bool, optional
    :param transform: Coordinate transform to use, defaults to None
    :type transform: Optional[Projection], optional
    :param directory: The directory to read/write from, defaults to "."
    :type directory: Union[pathlib.Path, str], optional
    :return: The updated axes
    :rtype: GeoAxes
    """
    directory = pathlib.Path(directory)
    ax.set_extent(margins)
    ax.spines["bottom"].set_linewidth(10)
    ax.spines["top"].set_linewidth(10)
    ax.spines["left"].set_linewidth(10)
    ax.spines["right"].set_linewidth(10)
    gl = ax.gridlines(linewidth=1, color="black", alpha=0.3, draw_labels=True)
    gl.top_labels = False
    gl.right_labels = False
    gl.left_labels = False
    gl.bottom_labels = False
    ocean = cartopy.feature.NaturalEarthFeature(
        "physical",
        "ocean",
        scale="110m",
        edgecolor="none",
        facecolor=cartopy.feature.COLORS["water"],
    )
    if tectonic:
        ax.add_feature(tectonic)
    if countries:
        ax.add_feature(countries)
    if bathymetry:
        ax.add_feature(bathymetry)
    if faults == True:
        fault_files = glob(f"{str(directory)}/*.fault")
        if len(fault_files) > 0:
            for kfault in range(len(fault_files)):
                print("...Adding fault trace from: " + str(fault_files[kfault]))
                fault_trace = np.genfromtxt(fault_files[kfault])
                ax.plot(
                    fault_trace[:, 0],
                    fault_trace[:, 1],
                    "k",
                    zorder=100,
                    transform=transform,
                )
    if aftershocks == True:
        aftershock_files = glob(f"{str(directory)}/*aftershock*")
        if len(aftershock_files) > 0:
            for kafter in range(len(aftershock_files)):
                print("...Adding aftershocks from: " + str(aftershock_files[kafter]))
                aftershock = np.genfromtxt(aftershock_files[kafter], delimiter="\t")
                aftershock_lat = aftershock[:, 3]
                aftershock_lon = aftershock[:, 4]
                aftershock_mag = aftershock[:, 6]
                ax.scatter(
                    aftershock_lon,
                    aftershock_lat,
                    s=aftershock_mag * 10,
                    c="0.5",
                    zorder=200,
                    transform=transform,
                )

    min_lon, max_lon, min_lat, max_lat = margins
    return ax
