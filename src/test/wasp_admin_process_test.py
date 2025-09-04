import json
import os
import pathlib
import shutil
import tempfile
from unittest import mock

import numpy as np
from obspy import read  # type:ignore
from typer.testing import CliRunner

from .testutils import (
    DATA_DIR,
    END_TO_END_DIR,
    HOME,
    RESULTS_DIR,
    get_strong_motion_json,
    get_tele_waves_json,
    update_manager_file_locations,
)

runner = CliRunner()


@mock.patch(target="wasp.green_functions.gf_retrieve")
def test_greens(p1):
    from wasp.wasp_admin.process import app

    tempdir = pathlib.Path(tempfile.mkdtemp())
    try:
        shutil.copyfile(
            RESULTS_DIR / "NP1" / "sampling_filter.json",
            tempdir / "sampling_filter.json",
        )
        shutil.copyfile(
            RESULTS_DIR / "NP1" / "tensor_info.json",
            tempdir / "tensor_info.json",
        )
        with open(DATA_DIR / "config.ini") as f:
            config = f.read().replace("/home/user/neic-finitefault", str(HOME))
        with open(tempdir / "config.ini", "w") as wf:
            wf.write(config)
        with open(tempdir / "GF_cgnss", "w"):
            pass
        with open(tempdir / "GF_strong", "w"):
            pass
        os.mkdir(tempdir / "logs")
        shutil.copy(
            RESULTS_DIR / "NP1" / "waveforms_body.txt", tempdir / "waveforms_body.txt"
        )
        shutil.copy(
            RESULTS_DIR / "NP1" / "tensor_info.json", tempdir / "tensor_info.json"
        )
        # test greens
        result = runner.invoke(
            app,
            [
                "greens",
                str(tempdir),
                "-t",
                "cgnss",
                "-t",
                "gnss",
                "-t",
                "strong",
                "-t",
                "body",
                "-c",
                str(tempdir / "config.ini"),
            ],
        )
        print(result.exception)
        assert result.exit_code == 0
        with open(tempdir / "strong_motion_gf.json") as f:
            strong = json.load(f)
        with open(tempdir / "cgnss_gf.json") as f:
            cgnss = json.load(f)
        assert cgnss == {
            "location": str(tempdir / "GF_cgnss"),
            "min_depth": 1,
            "max_depth": 49.8,
            "min_dist": 0,
            "max_dist": 1000,
            "dt": 0.4,
            "time_corr": 25,
        }
        assert strong == {
            "location": str(tempdir / "GF_strong"),
            "min_depth": 1,
            "max_depth": 49.8,
            "min_dist": 0,
            "max_dist": 1000,
            "dt": 0.4,
            "time_corr": 10,
        }
    finally:
        print("Cleaning up test directory.")
        shutil.rmtree(tempdir)


def test_process_tele():
    from wasp.wasp_admin.process import app

    tempdir = pathlib.Path(tempfile.mkdtemp())
    try:
        tele_waves = get_tele_waves_json()
        new_tele_waves = update_manager_file_locations(
            tele_waves, tempdir, replace_dir=str(RESULTS_DIR / "data")
        )
        shutil.copyfile(
            RESULTS_DIR / "NP1" / "tensor_info.json",
            tempdir / "tensor_info.json",
        )
        shutil.copyfile(
            RESULTS_DIR / "NP1" / "sampling_filter.json",
            tempdir / "sampling_filter.json",
        )
        os.mkdir(pathlib.Path(tempdir) / "P")
        for o, n in zip(tele_waves, new_tele_waves):
            shutil.copyfile(o["file"], n["file"])
        with open(pathlib.Path(tempdir) / "tele_waves.json", "w") as f:
            json.dump(new_tele_waves, f)

        # test base process tele
        result = runner.invoke(
            app,
            [
                "process-all",
                str(tempdir),
                "body",
            ],
        )
        assert result.exit_code == 0

        np.testing.assert_array_almost_equal(
            read(new_tele_waves[0]["file"])[0].data[:10],
            [
                -38.251144,
                -38.221546,
                -38.199234,
                -38.189747,
                -38.194893,
                -38.213146,
                -38.240723,
                -38.27213,
                -38.30074,
                -38.320812,
            ],
            decimal=5,
        )
        np.testing.assert_array_almost_equal(
            read(new_tele_waves[1]["file"])[0].data[:10],
            [
                -5.4753532,
                -5.5141807,
                -5.5615206,
                -5.61239,
                -5.6537414,
                -5.675031,
                -5.6726656,
                -5.6468496,
                -5.599341,
                -5.537252,
            ],
            decimal=5,
        )
    finally:
        print("Cleaning up test directory.")
        shutil.rmtree(tempdir)


def test_process_missing_file():
    from wasp.wasp_admin.process import app

    result = runner.invoke(
        app,
        [
            "process-all",
            ".",
            "body",
        ],
    )
    assert result.exit_code == 1
    assert isinstance(result.exception, FileNotFoundError)
    assert "does not exist!" in str(result.exception)


def test_remove_baseline():
    from wasp.wasp_admin.process import app

    tempdir = pathlib.Path(tempfile.mkdtemp())
    try:
        strong_motion = get_strong_motion_json()
        new_strong_motion = update_manager_file_locations(
            strong_motion, tempdir, replace_dir=str(RESULTS_DIR)
        )
        shutil.copyfile(
            RESULTS_DIR / "NP1" / "sampling_filter.json",
            tempdir / "sampling_filter.json",
        )
        os.mkdir(pathlib.Path(tempdir) / "data")
        os.mkdir(pathlib.Path(tempdir) / "data" / "STR")
        for o, n in zip(strong_motion, new_strong_motion):
            shutil.copyfile(o["file"], n["file"])
        with open(pathlib.Path(tempdir) / "strong_motion_waves.json", "w") as f:
            json.dump(new_strong_motion, f)
        shutil.copyfile(
            RESULTS_DIR / "NP1" / "tensor_info.json",
            tempdir / "tensor_info.json",
        )
        # test remove baseline
        result = runner.invoke(
            app,
            [
                "remove-baseline",
                str(tempdir / "data"),
            ],
        )
        assert result.exit_code == 0
    finally:
        print("Cleaning up test directory.")
        shutil.rmtree(tempdir)


def test_shift_match():
    from wasp.wasp_admin.process import app

    tempdir = pathlib.Path(tempfile.mkdtemp())
    try:
        tele_waves = get_tele_waves_json(all=True)
        new_tele_waves = update_manager_file_locations(
            tele_waves, tempdir, replace_dir=str(RESULTS_DIR)
        )
        shutil.copyfile(
            RESULTS_DIR / "NP1" / "sampling_filter.json",
            tempdir / "sampling_filter.json",
        )
        shutil.copyfile(
            RESULTS_DIR / "NP1" / "tensor_info.json",
            tempdir / "tensor_info.json",
        )
        os.mkdir(pathlib.Path(tempdir) / "data")
        os.mkdir(pathlib.Path(tempdir) / "data" / "P")
        os.mkdir(pathlib.Path(tempdir) / "data" / "SH")
        for o, n in zip(tele_waves, new_tele_waves):
            shutil.copyfile(o["file"], n["file"])
        with open(pathlib.Path(tempdir) / "tele_waves.json", "w") as f:
            json.dump(new_tele_waves, f)
        shutil.copyfile(
            RESULTS_DIR / "NP1" / "synthetics_body.txt", tempdir / "synthetics_body.txt"
        )
        # test base process tele
        result = runner.invoke(
            app,
            [
                "shift-match",
                str(tempdir),
                "body",
            ],
        )
        result = runner.invoke(
            app,
            [
                "shift-match",
                str(tempdir),
                "body",
                "-o",
                "manual",
            ],
        )
        print(result.exception)
        assert result.exit_code == 0
    finally:
        print("Cleaning up test directory.")
        shutil.rmtree(tempdir)
