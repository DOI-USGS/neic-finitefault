import pathlib
import shutil
import tempfile

import matplotlib.pyplot as plt
import numpy as np

from wasp.waveform_plots_NEIC import (
    add_metadata,
    filt_waveform,
    plot_spectra,
    plot_waveform_fits,
    plot_waveforms,
)


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
        ["tele_body", "surf_tele", "strong_motion", "tele_body"],
        ["BHZ", "SH", "BHZ", "SH"],
        custom="fill",
    )
    plt.close()


def test_add_metadata():
    ax = plt.axes()
    ax = add_metadata(
        [ax, ax, ax, ax],
        azimuths=[30, 10, 20, 5],
        comps=["BHZ", "SH", "BHZ", "SH"],
        distances=[5, 100, 20, 10],
        names=["name1", "name2", "name3", "name4"],
        type_str=["tele_body", "surf_tele", "strong_motion", "tele_body"],
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
            "cgps",
            start_margin=10,
            plot_directory=tempdir,
        )
        plot_waveform_fits(
            files,
            ["HNZ", "HNN", "HNE"],
            "strong_motion",
            start_margin=10,
            plot_directory=tempdir,
        )
        plot_waveform_fits(
            files,
            ["BHZ", "SH"],
            "tele_body",
            start_margin=10,
            plot_directory=tempdir,
        )
        plot_waveform_fits(
            files,
            ["BHZ", "SH"],
            "surf_tele",
            start_margin=10,
            plot_directory=tempdir,
        )
    finally:
        shutil.rmtree(tempdir)


def test_filt_waveform():
    file = {
        "file": pathlib.Path(__file__).parent
        / "data"
        / "end_to_end"
        / "teleseismic"
        / "G_CCD__BHE00.sac",
        "name": "test_name1",
        "component": "BHE",
        "azimuth": 1,
        "weight": 2,
        "trace_weight": 0,
        "dt": 1,
        "distance": 3,
        "start_signal": 1,
        "observed": np.arange(0, 100),
        "synthetic": np.arange(0, 100),
    }
    tempdir = tempfile.mkdtemp()
    try:
        fdict = filt_waveform(file, 0.9)
        np.testing.assert_almost_equal(
            fdict["synthetic"][0:7],
            [
                -201.69425217,
                -198.23057874,
                -194.73799219,
                -191.23556024,
                -187.74207612,
                -184.27608824,
                -180.85579492,
            ],
            decimal=6,
        )
    finally:
        shutil.rmtree(tempdir)


def test_plot_spectra():
    file = {
        "file": pathlib.Path(__file__).parent
        / "data"
        / "end_to_end"
        / "teleseismic"
        / "G_CCD__BHE00.sac",
        "name": "test_name1",
        "component": "BHE",
        "azimuth": 1,
        "weight": 2,
        "trace_weight": 0,
        "dt": 1,
        "distance": 3,
        "start_signal": 1,
        "observed": np.arange(0, 100),
        "synthetic": np.arange(0, 100),
    }

    tempdir = tempfile.mkdtemp()
    try:
        fdict = filt_waveform(file, 0.9)
        plot_spectra([fdict], 0.2, tempdir)
    finally:
        shutil.rmtree(tempdir)
