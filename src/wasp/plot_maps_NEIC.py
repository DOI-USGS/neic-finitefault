import pathlib
from glob import glob
from typing import List, Optional, Union

import cartopy  # type: ignore
import numpy as np
from cartopy.crs import Projection  # type: ignore
from cartopy.feature import Feature  # type: ignore
from cartopy.io.img_tiles import Stamen  # type: ignore
from cartopy.mpl.geoaxes import GeoAxes  # type: ignore
from matplotlib import colormaps, patches  # type: ignore
from matplotlib.colors import Colormap, ListedColormap  # type: ignore

"""
Set colorbar for slip
"""

rm = 100  # amount of lines to remove on black end of magma_r
magma_cpt = colormaps.get_cmap("magma_r")  # start with magma_r
white_bit = np.array([255 / 256, 250 / 256, 250 / 256, 1])  # create array of white
slip_cpt = magma_cpt(np.linspace(0, 1, 512))  # initialize slip_cpt
slip_cpt[rm:, :] = slip_cpt[0:-rm, :]  # move beginning up to remove black end
r_s = np.linspace(
    white_bit[0], slip_cpt[rm][0], rm
)  # gradient from white to beginning of new magma
g_s = np.linspace(white_bit[1], slip_cpt[rm][1], rm)
b_s = np.linspace(white_bit[2], slip_cpt[rm][2], rm)
slip_cpt[:rm, :][:, 0] = r_s
slip_cpt[:rm, :][:, 1] = g_s
slip_cpt[:rm, :][:, 2] = b_s
slipcpt = ListedColormap(slip_cpt)


def set_map_cartopy(
    ax: GeoAxes,
    margins: List[float],
    tectonic: Optional[Feature] = None,
    countries: Optional[Feature] = None,
    bathymetry: Optional[Feature] = None,
    faults: bool = True,
    aftershocks: bool = True,
    transform: Optional[Projection] = None,
    directory: Union[pathlib.Path, str] = pathlib.Path(),
) -> GeoAxes:
    """Initialize map with some features

    :param ax: The axes of the plot object
    :type ax: GeoAxes
    :param margins: The spatial margins of the plot
                    ([lon_min, lon_max, lat_min, lat_max])
    :type margins: List[float]
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
    :param directory: Where to write plots to, defaults to "."
    :type directory: Union[pathlib.Path, str], optional
    :return: The updated axes
    :rtype: GeoAxes
    """
    directory = pathlib.Path(directory)
    ax.set_extent(margins)
    ax.coastlines(resolution="10m", zorder=3)
    ax.spines["bottom"].set_linewidth(10)
    ax.spines["top"].set_linewidth(10)
    ax.spines["left"].set_linewidth(10)
    ax.spines["right"].set_linewidth(10)
    tiler = Stamen("terrain-background")
    ax.add_image(tiler, 10)
    gl = ax.gridlines(linewidth=1, color="black", alpha=0.3, draw_labels=True)
    gl.top_labels = False
    gl.right_labels = False
    # ax.add_feature(cartopy.feature.OCEAN)
    ocean = cartopy.feature.NaturalEarthFeature(
        "physical",
        "ocean",
        scale="110m",
        edgecolor="none",
        facecolor=cartopy.feature.COLORS["water"],
    )
    #    ax.add_feature(ocean)
    if tectonic:
        ax.add_feature(tectonic)
    if countries:
        ax.add_feature(countries)
    if bathymetry:
        ax.add_feature(bathymetry)
    if faults == True:
        faults_files = glob(str(directory) + "/*.fault")
        if len(faults_files) > 0:
            for kfault in range(len(faults_files)):
                print("...Adding fault trace from: " + str(faults_files[kfault]))
                fault_trace = np.genfromtxt(faults_files[kfault])
                ax.plot(
                    fault_trace[:, 0],
                    fault_trace[:, 1],
                    "k",
                    zorder=100,
                    transform=transform,
                )
    if aftershocks == True:
        aftershock_files = glob(str(directory) + "/*aftershock*")
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


def plot_map(
    ax: GeoAxes,
    latitudes: list,
    longitudes: list,
    values: list,
    min_val: Optional[float] = None,
    max_val: Optional[float] = None,
    transform: Optional[Projection] = None,
    cmap: Colormap = slipcpt,
) -> GeoAxes:
    """Plot data values in a mesh on a map

    :param ax: The axes of the plot object
    :type ax: GeoAxes
    :param latitudes: Array of latitudes
    :type latitudes: list
    :param longitudes: Array of longitudes
    :type longitudes: list
    :param values: Array of values
    :type values: list
    :param min_val: Minimum value of the colorbar, defaults to None
    :type min_val: Optional[float], optional
    :param max_val: Maximum value of the colorbar, defaults to None
    :type max_val: Optional[float], optional
    :param transform: Coordinate transform to use, defaults to None
    :type transform: Optional[Projection], optional
    :param cmap: The colormap to use, defaults to slipcpt
    :type cmap: Colormap, optional
    :return: The updated axes and mesh
    :rtype: GeoAxes
    """
    min_val = min([np.amin(value) for value in values]) if not min_val else min_val
    max_val = max([np.amax(value) for value in values]) if not max_val else max_val
    zipped = zip(longitudes, latitudes, values)
    for longitude, latitude, value in zipped:
        if np.prod(longitude.shape) > 1:
            cs = ax.pcolormesh(
                longitude,
                latitude,
                value,
                zorder=3,
                vmin=min_val,
                cmap=cmap,
                vmax=max_val,
                edgecolor=None,
                transform=transform,
            )
            # vmax=max_val, edgecolor='0.5', lw=0.5, transform=transform)
    return ax, cs


def plot_borders(
    ax: GeoAxes,
    latitudes: np.ndarray,
    longitudes: np.ndarray,
    transform: Optional[Projection] = None,
) -> GeoAxes:
    """Plot borders of a grid on a map

    :param ax: The axes of the plot object
    :type ax: GeoAxes
    :param latitudes: Array of latitudes
    :type latitudes: np.ndarray
    :param longitudes: Array of longitudes
    :type longitudes: np.ndarray
    :param transform: Coordinate transform to use, defaults to None
    :type transform: Optional[Projection], optional
    :return: The updated axes
    :rtype: GeoAxes
    """
    zipped = zip(longitudes, latitudes)
    for longitude, latitude in zipped:
        edge1 = [longitude[0, 0], latitude[0, 0]]
        edge2 = [longitude[-1, 0], latitude[-1, 0]]
        edge3 = [longitude[0, -1], latitude[0, -1]]
        edge4 = [longitude[-1, -1], latitude[-1, -1]]
        poly = patches.Polygon(
            [edge1, edge2, edge4, edge3, edge1], facecolor="0.9", edgecolor="1"
        )
        if np.prod(longitude.shape) > 1:
            ax.add_patch(
                patches.Polygon(
                    poly.get_xy(),
                    closed=True,
                    ec="k",
                    lw=3,
                    fill=False,
                    transform=transform,
                    zorder=3,
                )
            )
    return ax
