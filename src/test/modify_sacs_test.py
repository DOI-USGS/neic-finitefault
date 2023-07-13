import json
import pathlib
import shutil
import tempfile
from copy import deepcopy
from glob import glob
from unittest import mock

import numpy as np
from obspy import read

from wasp.modify_sacs import __is_number, correct_waveforms, plot_channels

DATA_DIR = pathlib.Path(__file__).parent / "data" / "end_to_end" / "teleseismic"

TELE_WAVES = [
    {
        "name": "G_CCD__BHZ00",
        "file": str(DATA_DIR / "G_CCD__BHE00.sac"),
        "component": "BHZ",
        "station": "G_CCD",
        # made up a start signal
        "start_signal": 120,
        "dt": 0.2,
        "duration": 60000,
    },
    {
        "name": "IUTUC__BH200",
        "file": str(DATA_DIR / "IUTUC__BH200.sac"),
        "component": "BH2",
        "station": "IUTUC",
        # made up a start signal
        "start_signal": 100,
        # made up a dt
        "dt": 0.2,
        "duration": 60000,
    },
]


def test_is_number():
    assert __is_number("INVALID") == False
    assert __is_number("999") == True


def test_plot_channels():
    tempdir = tempfile.mkdtemp()
    try:
        with open(pathlib.Path(tempdir) / "tele_waves.json", "w") as f:
            json.dump(TELE_WAVES, f)
        with open(pathlib.Path(tempdir) / "strong_motion_waves.json", "w") as f:
            json.dump(TELE_WAVES, f)
        with open(pathlib.Path(tempdir) / "surf_waves.json", "w") as f:
            json.dump(TELE_WAVES, f)

        plot_channels("tele_waves.json", tempdir, tempdir)
        plot_channels("strong_motion_waves.json", tempdir, tempdir)
        plot_channels("surf_waves.json", tempdir, tempdir)
        written_plots = [f.replace(f"{tempdir}/", "") for f in glob(f"{tempdir}/*/*")]
        assert len(written_plots) == 6
        assert "review_strong/IUTUC__BH200_BH2.png" in written_plots
        assert "review_strong/G_CCD__BHZ00_BHZ.png" in written_plots
        assert "review_surf/IUTUC__BH200_BH2.png" in written_plots
        assert "review_surf/G_CCD__BHZ00_BHZ.png" in written_plots
        assert "review_tele/IUTUC__BH200_BH2.png" in written_plots
        assert "review_tele/G_CCD__BHZ00_BHZ.png" in written_plots
    finally:
        shutil.rmtree(tempdir)


@mock.patch(
    "builtins.input",
    side_effect=[
        "G_CCD__BHZ00",
        "BHZ",
        -1,
        "",
        "exit",
        "IUTUC__BH200",
        "BH2",
        "",
        -10,
        "exit",
    ],
)
def test_correct_waveforms_input(mock_input):
    tempdir = tempfile.mkdtemp()
    try:
        tele_waves = deepcopy(TELE_WAVES)
        for idx, t in enumerate(TELE_WAVES):
            new_file = tempdir + "/" + t["file"].split("/")[-1]
            shutil.copyfile(t["file"], new_file)
            tele_waves[idx]["file"] = new_file
            stream_before = read(new_file)

        print("channels_before", tele_waves)
        with open(pathlib.Path(tempdir) / "tele_waves.json", "w") as f:
            json.dump(tele_waves, f)
        correct_waveforms(pathlib.Path(tempdir) / "tele_waves.json")
        correct_waveforms(pathlib.Path(tempdir) / "tele_waves.json")

        # check baseline shift
        stream = read(pathlib.Path(tempdir) / "IUTUC__BH200.sac")
        assert np.max(stream[0].data) == 1090667.0

        # check time shift
        with open(pathlib.Path(tempdir) / "tele_waves.json", "r") as f:
            updated_tele_waves = json.load(f)
            print("channels_after", updated_tele_waves)
        assert updated_tele_waves[0]["start_signal"] == 140
    finally:
        shutil.rmtree(tempdir)


def test_correct_waveforms_dict():
    tempdir = tempfile.mkdtemp()
    try:
        tele_waves = deepcopy(TELE_WAVES)
        for idx, t in enumerate(TELE_WAVES):
            new_file = tempdir + "/" + t["file"].split("/")[-1]
            shutil.copyfile(t["file"], new_file)
            tele_waves[idx]["file"] = new_file
            stream_before = read(new_file)

        print("channels_before", tele_waves)
        with open(pathlib.Path(tempdir) / "tele_waves.json", "w") as f:
            json.dump(tele_waves, f)
        correct_waveforms(
            pathlib.Path(tempdir) / "tele_waves.json",
            input=False,
            station_dict={
                "G_CCD__BHZ00": {
                    "BHZ": {
                        "time_correction": -1,
                        "baseline_correction": None,
                    }
                },
                "IUTUC__BH200": {
                    "BH2": {
                        "time_correction": None,
                        "baseline_correction": -10,
                    }
                },
            },
        )

        # check baseline shift
        stream = read(pathlib.Path(tempdir) / "IUTUC__BH200.sac")
        assert np.max(stream[0].data) == 1090667.0

        # check time shift
        with open(pathlib.Path(tempdir) / "tele_waves.json", "r") as f:
            updated_tele_waves = json.load(f)
            print("channels_after", updated_tele_waves)
        assert updated_tele_waves[0]["start_signal"] == 140
    finally:
        shutil.rmtree(tempdir)
