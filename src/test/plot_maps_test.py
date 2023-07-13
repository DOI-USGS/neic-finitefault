import pathlib
import shutil
from tempfile import mkdtemp

import cartopy.crs as ccrs
import cartopy.feature as cfeature
import matplotlib.pyplot as plt
import numpy as np

from wasp.plot_maps import plot_borders, plot_map, set_map_cartopy


def test_plot_borders():
    tempdir = mkdtemp()
    try:
        tempfile = pathlib.Path(tempdir) / "test.png"
        ax = plt.axes(projection=ccrs.PlateCarree())
        ax = set_map_cartopy(ax, [46, 50, 10, 13], cfeature.OCEAN, cfeature.BORDERS)
        updated_ax = plot_borders(
            ax,
            [
                np.array([[46.909, 48.359], [46.909, 48.359]]),
                np.array([[47.217, 49.062], [47.217, 49.062]]),
            ],
            [
                np.array([[11.727, 12.650], [12.650, 11.727]]),
                np.array([[10.022, 11.317], [11.317, 10.022]]),
            ],
            transform=None,
        )
        plt.savefig(tempfile)
        plt.close()
    finally:
        print("Cleaning up test directory")
        shutil.rmtree(tempdir)


def test_set_map_cartopy():
    tempdir = mkdtemp()
    try:
        tempfile = pathlib.Path(tempdir) / "test.png"
        ax = plt.axes(projection=ccrs.PlateCarree())
        margins = [11.5, 21.7, 44.1, 49.7]
        updated_ax = set_map_cartopy(ax, margins, cfeature.RIVERS, cfeature.BORDERS)
        plt.savefig(tempfile)
        plt.close()
    finally:
        print("Cleaning up test directory")
        shutil.rmtree(tempdir)


def test_plot_map():
    tempdir = mkdtemp()
    try:
        tempfile = pathlib.Path(tempdir) / "test.png"
        ax = plt.axes(projection=ccrs.PlateCarree())
        ax = set_map_cartopy(ax, [-119, -116, 33, 34], cfeature.OCEAN, cfeature.BORDERS)
        x = np.linspace(-118.5, -117.5, 3)
        y = np.linspace(33.6, 33.8, 3)
        X, Y = np.meshgrid(x, y)
        Z = [[[1, 2, 3], [-1, 0, 1], [20, 11, -3]]]
        updated_ax = plot_map(
            ax,
            np.array([Y]),
            np.array([X]),
            Z,
            transform=None,
        )
        plt.savefig(tempfile)
        plt.close()
    finally:
        print("Cleaning up test directory")
        shutil.rmtree(tempdir)
