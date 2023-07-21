import json
import pathlib
import shutil
import tempfile
from copy import deepcopy
from glob import glob
from os import mkdir
from unittest import mock

import numpy as np
from obspy import read

from wasp.modify_sacs import __is_number, correct_waveforms, plot_channels

from .testutils import update_manager_file_locations

RESULTS = pathlib.Path(__file__).parent / "data" / "end_to_end" / "results"


with open(RESULTS / "NP1" / "tele_waves.json", "r") as f:
    TELE_WAVES = update_manager_file_locations(json.load(f)[0:3], RESULTS / "data")
with open(RESULTS / "NP1" / "surf_waves.json", "r") as f:
    SURF_WAVES = update_manager_file_locations(json.load(f)[0:3], RESULTS / "data")
with open(RESULTS / "NP1" / "strong_motion_waves.json", "r") as f:
    STRONG_WAVES = update_manager_file_locations(json.load(f)[0:3], RESULTS / "data")


def test_is_number():
    assert __is_number("INVALID") == False
    assert __is_number("999") == True


def test_plot_channels():
    tempdir = tempfile.mkdtemp()
    try:
        with open(pathlib.Path(tempdir) / "tele_waves.json", "w") as f:
            json.dump(TELE_WAVES, f)
        with open(pathlib.Path(tempdir) / "strong_motion_waves.json", "w") as f:
            json.dump(STRONG_WAVES, f)
        with open(pathlib.Path(tempdir) / "surf_waves.json", "w") as f:
            json.dump(SURF_WAVES, f)

        plot_channels("tele_waves.json", tempdir, tempdir)
        plot_channels("strong_motion_waves.json", tempdir, tempdir)
        plot_channels("surf_waves.json", tempdir, tempdir)
        written_plots = [f.replace(f"{tempdir}/", "") for f in glob(f"{tempdir}/*/*")]
        assert len(written_plots) == 9
        target_plots = [
            "review_strong/VA03_HNZ.png",
            "review_strong/GO04_HNN.png",
            "review_strong/CO03_HNZ.png",
            "review_surf/MACI_BHZ.png",
            "review_surf/MPG_BHZ.png",
            "review_surf/RCBR_BHZ.png",
            "review_tele/MACI_BHZ.png",
            "review_tele/MPG_BHZ.png",
            "review_tele/RCBR_BHZ.png",
        ]
        for t in target_plots:
            assert t in written_plots
    finally:
        shutil.rmtree(tempdir)


@mock.patch(
    "builtins.input",
    side_effect=[
        "MACI",
        "BHZ",
        -100,
        "",
        "exit",
        "RCBR",
        "BHZ",
        "",
        -10,
        "exit",
    ],
)
def test_correct_waveforms_input(mock_input):
    tempdir = tempfile.mkdtemp()
    try:
        new_tele_waves = update_manager_file_locations(
            TELE_WAVES, tempdir, replace_dir=str(RESULTS / "data")
        )
        mkdir(pathlib.Path(tempdir) / "P")
        for o, n in zip(TELE_WAVES, new_tele_waves):
            shutil.copyfile(o["file"], n["file"])
        with open(pathlib.Path(tempdir) / "tele_waves.json", "w") as f:
            json.dump(new_tele_waves, f)
        correct_waveforms(pathlib.Path(tempdir) / "tele_waves.json")
        correct_waveforms(pathlib.Path(tempdir) / "tele_waves.json")

        # check baseline shift
        stream = read(pathlib.Path(tempdir) / "P" / "final_IU_RCBR_BHZ.sac")
        assert np.max(stream[0].data) == 356.4468994140625

        # check time shift
        with open(pathlib.Path(tempdir) / "tele_waves.json", "r") as f:
            updated_tele_waves = json.load(f)
            print("channels_after", updated_tele_waves)
        assert updated_tele_waves[0]["start_signal"] == 595
    finally:
        shutil.rmtree(tempdir)


def test_correct_waveforms_dict():
    tempdir = tempfile.mkdtemp()
    try:
        new_tele_waves = update_manager_file_locations(
            TELE_WAVES, tempdir, replace_dir=str(RESULTS / "data")
        )
        mkdir(pathlib.Path(tempdir) / "P")
        for o, n in zip(TELE_WAVES, new_tele_waves):
            shutil.copyfile(o["file"], n["file"])
        with open(pathlib.Path(tempdir) / "tele_waves.json", "w") as f:
            json.dump(new_tele_waves, f)

        correct_waveforms(
            pathlib.Path(tempdir) / "tele_waves.json",
            input=False,
            station_dict={
                "MACI": {
                    "BHZ": {
                        "time_correction": -100,
                        "baseline_correction": None,
                    }
                },
                "RCBR": {
                    "BHZ": {
                        "time_correction": None,
                        "baseline_correction": -10,
                    }
                },
            },
        )

        # check baseline shift
        stream = read(pathlib.Path(tempdir) / "P" / "final_IU_RCBR_BHZ.sac")
        assert np.max(stream[0].data) == 356.4468994140625

        # check time shift
        with open(pathlib.Path(tempdir) / "tele_waves.json", "r") as f:
            updated_tele_waves = json.load(f)
            print("channels_after", updated_tele_waves)
        assert updated_tele_waves[0]["start_signal"] == 595
    finally:
        shutil.rmtree(tempdir)
