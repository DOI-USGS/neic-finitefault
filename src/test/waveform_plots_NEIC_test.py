import json
import os
import pathlib
import shutil
import tempfile
from copy import deepcopy

import matplotlib.pyplot as plt
import numpy as np
import pytest

from wasp.waveform_plots_NEIC import (
    add_metadata,
    filt_waveform,
    plot_spectra,
    plot_waveform_fits,
    plot_waveforms,
)

from .testutils import update_manager_file_locations

RESULTS = pathlib.Path(__file__).parent / "data" / "end_to_end" / "results"


with open(RESULTS / "NP1" / "tele_waves.json", "r") as f:
    TELE_WAVES = update_manager_file_locations(json.load(f)[0:3], RESULTS / "data")


def test_plot_waveforms():
    ax = plt.axes()
    time = np.arange(0, 4 * np.pi, np.pi / 8)
    waveform1 = np.cos(time)
    waveform2 = np.sin(time)
    ax = plot_waveforms(
        [ax, ax, ax, ax],
        [time, time, time, time],
        [waveform1, waveform2, waveform2, waveform2],
        [1, 0, 1, 1],
        ["body", "surf", "strong", "body"],
        ["BHZ", "SH", "BHZ", "SH"],
        custom="fill",
    )
    plt.close()


def test_add_metadata():
    ax = plt.axes()
    ax = add_metadata(
        [ax, ax, ax, ax],
        azimuths=[30, 10, 20, 5],
        comps="BHZ",
        distances=[5, 100, 20, 10],
        names=["name1", "name2", "name3", "name4"],
        type_str="body",
        weights=[0, 1, 2, 3],
    )
    plt.close()


def test_plot_waveform_fits():
    files = [
        {
            "file": "test_file1",
            "name": "test_name1",
            "component": "HN1",
            "azimuth": 1,
            "weight": 2,
            "trace_weight": 0,
            "dt": 1,
            "distance": 3,
            "start_signal": 1,
            "observed": np.arange(0, 100),
            "synthetic": np.arange(0, 100),
        },
        {
            "file": "test_file2",
            "name": "test_name2",
            "component": "BHZ",
            "dt": 1,
            "azimuth": 1,
            "weight": 2,
            "distance": 3,
            "trace_weight": 1,
            "start_signal": 1,
            "observed": np.arange(0, 100),
            "synthetic": np.arange(0, 100),
        },
    ]
    tempdir = tempfile.mkdtemp()
    try:
        plot_waveform_fits(
            files,
            ["LXZ", "LXN", "LXE"],
            "cgnss",
            start_margin=10,
            plot_directory=tempdir,
        )
        plot_waveform_fits(
            files,
            ["HNZ", "HNN", "HNE"],
            "strong",
            start_margin=10,
            plot_directory=tempdir,
        )
        plot_waveform_fits(
            files,
            ["BHZ", "SH"],
            "body",
            start_margin=10,
            plot_directory=tempdir,
        )
        plot_waveform_fits(
            files,
            ["BHZ", "SH"],
            "surf",
            start_margin=10,
            plot_directory=tempdir,
        )
    finally:
        shutil.rmtree(tempdir)


def test_filt_waveform():
    file = deepcopy(TELE_WAVES[0])
    tempdir = tempfile.mkdtemp()
    try:
        fdict = filt_waveform(file, 0.9)
        np.testing.assert_almost_equal(
            fdict["synthetic"][0:100],
            [
                52.77696617,
                58.55397582,
                63.99139388,
                69.01639714,
                73.55929998,
                77.55494286,
                80.94394265,
                83.67386233,
                85.70037786,
                86.98848255,
                87.51369777,
                87.26318012,
                86.23654801,
                84.44623629,
                81.91726843,
                78.68650192,
                74.80155314,
                70.31962511,
                65.3063265,
                59.83440427,
                53.98226517,
                47.83226693,
                41.46890433,
                34.97705834,
                28.44039299,
                21.93987278,
                15.55234535,
                9.34920141,
                3.39520542,
                -2.25239844,
                -7.54446844,
                -12.44028384,
                -16.90768674,
                -20.92310967,
                -24.47149749,
                -27.54609393,
                -30.14805705,
                -32.28588464,
                -33.97465142,
                -35.23508439,
                -36.09253831,
                -36.5759704,
                -36.71701878,
                -36.54923915,
                -36.10747003,
                -35.42723903,
                -34.54413301,
                -33.49311577,
                -32.30782634,
                -31.01989251,
                -29.65826808,
                -28.24860016,
                -26.81267306,
                -25.36802796,
                -23.92786976,
                -22.50132284,
                -21.09400936,
                -19.70884657,
                -18.34693005,
                -17.00839218,
                -15.6931708,
                -14.40165701,
                -13.13519608,
                -11.89640533,
                -10.68927604,
                -9.51905936,
                -8.39198999,
                -7.31494547,
                -6.2951406,
                -5.33990705,
                -4.4565315,
                -3.65207389,
                -2.93310065,
                -2.30533557,
                -1.77329391,
                -1.3399715,
                -1.00661994,
                -0.7726147,
                -0.63545343,
                -0.59096672,
                -0.63380331,
                -0.75814372,
                -0.95847571,
                -1.23023114,
                -1.57016,
                -1.97642982,
                -2.44850135,
                -2.98682755,
                -3.59239949,
                -4.26615625,
                -5.00828427,
                -5.81743732,
                -6.68992335,
                -7.61894467,
                -8.59402772,
                -9.60078781,
                -10.62111717,
                -11.63379009,
                -12.61540103,
                -13.54150947,
            ],
            decimal=8,
        )
    finally:
        shutil.rmtree(tempdir)


@pytest.mark.skipif(
    os.getenv("CI_REGISTRY") is not None,
    reason="Build runner does not have the resources to run",
)
def test_plot_spectra():
    file = deepcopy(TELE_WAVES[0])

    tempdir = tempfile.mkdtemp()
    try:
        fdict = filt_waveform(file, 0.9)
        plot_spectra([fdict], 0.2, tempdir)
    finally:
        shutil.rmtree(tempdir)
