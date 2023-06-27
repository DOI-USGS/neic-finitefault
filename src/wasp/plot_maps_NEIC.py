import os
import numpy as np
import matplotlib
#matplotlib.use('Qt5Agg')
from matplotlib import pyplot as plt
from matplotlib import gridspec, ticker, patches
import cartopy.crs as ccrs
import cartopy
import cartopy.io.shapereader as shpreader
import cartopy.feature as cf
from matplotlib import cm
from matplotlib.colors import ListedColormap
from glob import glob
from cartopy.io.img_tiles import Stamen

"""
Set colorbar for slip
"""

rm = 100 # amount of lines to remove on black end of magma_r
magma_cpt = cm.get_cmap('magma_r', 512) #start with magma_r
white_bit = np.array([255/256, 250/256, 250/256, 1]) #create array of white
slip_cpt = magma_cpt(np.linspace(0,1,512)) #initialize slip_cpt
slip_cpt[rm:,:] = slip_cpt[0:-rm,:] #move beginning up to remove black end
r_s = np.linspace(white_bit[0],slip_cpt[rm][0],rm) # gradient from white to beginning of new magma
g_s = np.linspace(white_bit[1],slip_cpt[rm][1],rm)
b_s = np.linspace(white_bit[2],slip_cpt[rm][2],rm)
slip_cpt[:rm,:][:,0] = r_s
slip_cpt[:rm,:][:,1] = g_s
slip_cpt[:rm,:][:,2] = b_s
slipcpt = ListedColormap(slip_cpt)

def set_map_cartopy(ax, margins, tectonic=None, countries=None, bathymetry=None, faults=True, aftershocks=True, transform=None):
    """
    """
    ax.set_extent(margins)
    ax.coastlines(resolution='10m', zorder=3)
    ax.spines["bottom"].set_linewidth(10)
    ax.spines["top"].set_linewidth(10)
    ax.spines["left"].set_linewidth(10)
    ax.spines["right"].set_linewidth(10)
    tiler = Stamen('terrain-background')
    ax.add_image(tiler, 10)
    gl = ax.gridlines(linewidth=1, color='black', alpha=0.3, draw_labels=True)
    gl.top_labels = False
    gl.right_labels = False
    #ax.add_feature(cartopy.feature.OCEAN)
    ocean = cartopy.feature.NaturalEarthFeature('physical', 'ocean', \
        scale='110m', edgecolor='none', facecolor=cartopy.feature.COLORS['water'])
#    ax.add_feature(ocean)
    if tectonic:
        ax.add_feature(tectonic)
    if countries:
        ax.add_feature(countries)
    if bathymetry:
        ax.add_feature(bathymetry)
    if faults==True:
        faults = glob('*.fault')
        if len(faults) > 0:
            for kfault in range(len(faults)):
                print('...Adding fault trace from: '+str(faults[kfault]))
                fault_trace = np.genfromtxt(faults[kfault])
                ax.plot(fault_trace[:,0],fault_trace[:,1],'k',zorder=100,transform=transform)
    if aftershocks==True:
        aftershocks = glob('*aftershock*')
        if len(aftershocks) > 0:
            for kafter in range(len(aftershocks)):
                print('...Adding aftershocks from: '+str(aftershocks[kafter]))
                aftershock = np.genfromtxt(aftershocks[kafter], delimiter="\t")
                aftershock_lat = aftershock[:,3]
                aftershock_lon = aftershock[:,4]
                aftershock_mag = aftershock[:,6]
                ax.scatter(aftershock_lon, aftershock_lat, s=aftershock_mag*10, c='0.5', zorder = 200,transform=transform)

    min_lon, max_lon, min_lat, max_lat = margins
    return ax


def plot_map(ax, latitudes, longitudes, values, min_val=None, max_val=None,
             transform=None, cmap=slipcpt):
    """
    """
    min_val = min([np.amin(value) for value in values]) if not min_val else min_val
    max_val = max([np.amax(value) for value in values]) if not max_val else max_val
    zipped = zip(longitudes, latitudes, values)
    for longitude, latitude, value in zipped:
        if np.prod(longitude.shape) > 1:
            cs = ax.pcolormesh(
                longitude, latitude, value, zorder=3, vmin=min_val, cmap=cmap,
                vmax=max_val, edgecolor=None, transform=transform)
                #vmax=max_val, edgecolor='0.5', lw=0.5, transform=transform)
    return ax, cs

def plot_borders(ax, latitudes, longitudes, transform=None):
    """
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
