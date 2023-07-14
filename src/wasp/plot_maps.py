import os
import numpy as np
#import matplotlib
#matplotlib.use('Qt5Agg')
from matplotlib import pyplot as plt
from matplotlib import gridspec, ticker, patches
import cartopy.crs as ccrs
import cartopy.io.shapereader as shpreader
import cartopy.feature as cf


def set_map_cartopy(ax, margins, tectonic=None, countries=None):
    """Initialize cartopy map with some properties

    :param ax: object where to plot cartopy map
    :param margins: (lat, lon) margins of plot
    :param tectonic: tectonic features to plot
    :param countries: plot countries boundaries
    :type ax: Axes
    :type margins: list
    :type tectonic: shapefile, optional
    :type countries: shapefile, optional
    """
    ax.coastlines(resolution='10m', zorder=3)
    gl = ax.gridlines(linewidth=1, color='black', alpha=0.3, draw_labels=True)
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


def plot_map(ax, latitudes, longitudes, values, min_val=None, max_val=None,
             transform=None):
    """Plot data values in a (lat, lon) grid in a cartopy map

    :param ax: object where to plot cartopy map
    :param latitudes: latitudes of the grid
    :param longitudes: longitudes of the grid
    :param values: values over the (lat, lon) grid
    :param min_val: minimum value in the colorbar of the plot
    :param max_val: maximum value in the colorbar of the plot
    :param transform: coordinate transform to use
    :type ax: Axes
    :type latitudes: array
    :type longitudes: array
    :type values: array
    :type min_val: float, optional
    :type max_val: float, optional
    :type transform: CRS
    """
    min_val = min([np.amin(value) for value in values]) if not min_val else min_val
    max_val = max([np.amax(value) for value in values]) if not max_val else max_val
    zipped = zip(longitudes, latitudes, values)
    for longitude, latitude, value in zipped:
        if np.prod(longitude.shape) > 1:
            cs = ax.pcolormesh(
                longitude, latitude, value, zorder=3, vmin=min_val, cmap='jet',
                vmax=max_val, edgecolor='none', transform=transform)
    return ax, cs


def plot_borders(ax, latitudes, longitudes, transform=None):
    """Plot borders of a (lat, lon) grid in a cartopy map.

    :param ax: object where to plot cartopy map
    :param latitudes: latitudes of the grid
    :param longitudes: longitudes of the grid
    :param transform: coordinate transform to use
    :type ax: Axes
    :type latitudes: array
    :type longitudes: array
    :type transform: CRS
    """
    zipped = zip(longitudes, latitudes)
    for longitude, latitude in zipped:
        edge1 = [longitude[0, 0], latitude[0, 0]]
        edge2 = [longitude[-1, 0], latitude[-1, 0]]
        edge3 = [longitude[0, -1], latitude[0, -1]]
        edge4 = [longitude[-1, -1], latitude[-1, -1]]
        poly = patches.Polygon(
            [edge1, edge2, edge4, edge3, edge1], facecolor='0.9', edgecolor='1')
        if np.prod(longitude.shape) > 1:
            ax.add_patch(
                patches.Polygon(
                    poly.get_xy(), closed=True, ec='k', lw=3,
                    fill=False, transform=transform, zorder=3))
    return ax


if __name__ == '__main__':
    print(1)
