from typing import List, Optional, Tuple

import numpy as np
from cartopy.crs import Projection  # type: ignore
from cartopy.feature import Feature  # type: ignore
from cartopy.mpl.geoaxes import GeoAxes  # type: ignore
from matplotlib import patches  # type: ignore
from matplotlib.collections import QuadMesh  # type: ignore


def set_map_cartopy(
    ax: GeoAxes,
    margins: List[float],
    tectonic: Optional[Feature] = None,
    countries: Optional[Feature] = None,
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
    :return: The updated axes
    :rtype: GeoAxes
    """
    ax.coastlines(resolution="10m", zorder=3)
    gl = ax.gridlines(linewidth=1, color="black", alpha=0.3, draw_labels=True)
    gl.bottom_labels = False
    gl.right_labels = False
    if tectonic:
        ax.add_feature(tectonic)
    if countries:
        ax.add_feature(countries)
    min_lon, max_lon, min_lat, max_lat = margins
    ax.set_xlim(min_lon, max_lon)
    ax.set_ylim(min_lat, max_lat)
    return ax


def plot_map(
    ax: GeoAxes,
    latitudes: np.ndarray,
    longitudes: np.ndarray,
    values: np.ndarray,
    min_val: Optional[float] = None,
    max_val: Optional[float] = None,
    transform: Optional[Projection] = None,
) -> Tuple[GeoAxes, QuadMesh]:
    """Plot data values in a mesh on a map

    :param ax: The axes of the plot object
    :type ax: GeoAxes
    :param latitudes: Array of latitudes
    :type latitudes: np.ndarray
    :param longitudes: Array of longitudes
    :type longitudes: np.ndarray
    :param values: Array of values
    :type values: np.ndarray
    :param min_val: Minimum value of the colorbar, defaults to None
    :type min_val: Optional[float], optional
    :param max_val: Maximum value of the colorbar, defaults to None
    :type max_val: Optional[float], optional
    :param transform: Coordinate transform to use, defaults to None
    :type transform: Optional[Projection], optional
    :return: The updated axes and mesh
    :rtype: Tuple[GeoAxes, QuadMesh]
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
                cmap="jet",
                vmax=max_val,
                edgecolor="none",
                transform=transform,
            )
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
        print(longitude)
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
