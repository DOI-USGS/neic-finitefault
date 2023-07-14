import pathlib
import shutil
import tempfile

import matplotlib.pyplot as plt
import numpy as np

from wasp.waveform_plots import add_metadata, plot_waveform_fits, plot_waveforms


def test_plot_waveforms():
    ax = plt.axes()
    time = np.arange(0, 4 * np.pi, np.pi / 8)
    waveform1 = np.cos(time)
    waveform2 = np.sin(time)
    ax = plot_waveforms(
        [ax, ax, ax, ax],
        [time, time, None, None],
        [waveform1, waveform2, [], None],
        custom="fill",
    )
    plt.close()


def test_add_metadata():
    ax = plt.axes()
    ax = add_metadata(
        [ax, ax, ax],
        azimuths=[30, 10, None],
        distances=[5, 100, None],
        names=["name1", "name2", None],
        weights=[0, 1, None],
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
            event="test_event",
            plot_directory=tempdir,
        )
        plot_waveform_fits(
            files,
            ["HNZ", "HNN", "HNE"],
            "strong_motion",
            start_margin=10,
            event="test_event",
            plot_directory=tempdir,
        )
        plot_waveform_fits(
            files,
            ["BHZ", "SH"],
            "tele_body",
            start_margin=10,
            event="test_event",
            plot_directory=tempdir,
        )
        plot_waveform_fits(
            files,
            ["BHZ", "SH"],
            "surf_tele",
            start_margin=10,
            event="test_event",
            plot_directory=tempdir,
        )
    finally:
        shutil.rmtree(tempdir)
