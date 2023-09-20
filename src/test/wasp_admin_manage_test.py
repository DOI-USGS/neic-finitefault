import json
import os
import pathlib
import shutil
import tempfile
from unittest import mock

import numpy as np
from obspy import read
from typer.testing import CliRunner

from .testutils import (
    END_TO_END_DIR,
    RESULTS_DIR,
    get_tele_waves_json,
    update_manager_file_locations,
)

runner = CliRunner()
CHANNELS = [
    {
        "component": "BHZ",
        "name": "EFGH",
        "trace_weight": 1.0,
    },
    {
        "component": "BH2",
        "name": "EFGH",
        "trace_weight": 1.0,
    },
    {
        "component": "BH1",
        "name": "EFGH",
        "trace_weight": 1.0,
    },
    {
        "component": "SH",
        "name": "ABCD",
        "trace_weight": 1.0,
    },
]


@mock.patch("wasp.data_acquisition.acquisition", return_value=None)
def test_acquire(p1):
    from wasp.wasp_admin.manage import app

    tempdir = pathlib.Path(tempfile.mkdtemp())
    try:
        shutil.copyfile(
            END_TO_END_DIR / "info" / "20003k7a_cmt_CMT", tempdir / "20003k7a_cmt_CMT"
        )
        result = runner.invoke(
            app,
            ["acquire", str(tempdir), str(tempdir / "20003k7a_cmt_CMT")],
        )
    finally:
        print("Cleaning up test directory.")

        shutil.rmtree(tempdir)


def test_fill_dicts():
    from wasp.wasp_admin.manage import app

    tempdir = pathlib.Path(tempfile.mkdtemp())
    try:
        shutil.copyfile(
            END_TO_END_DIR / "info" / "20003k7a_cmt_CMT",
            tempdir / "20003k7a_cmt_CMT",
        )
        shutil.copyfile(
            RESULTS_DIR / "NP1" / "sampling_filter.json",
            tempdir / "sampling_filter.json",
        )
        result = runner.invoke(
            app,
            [
                "fill-dicts",
                str(tempdir),
                str(tempdir / "20003k7a_cmt_CMT"),
                "-t",
                "body",
            ],
        )
        print(result.stdout)
        assert result.exit_code == 0
    finally:
        print("Cleaning up test directory.")
        shutil.rmtree(tempdir)


def test_fill_dicts_missing_file():
    from wasp.wasp_admin.manage import app

    tempdir = pathlib.Path(tempfile.mkdtemp())
    try:
        shutil.copyfile(
            END_TO_END_DIR / "info" / "20003k7a_cmt_CMT",
            tempdir / "20003k7a_cmt_CMT",
        )
        result = runner.invoke(
            app,
            [
                "fill-dicts",
                str(tempdir),
                str(tempdir / "20003k7a_cmt_CMT"),
                "-t",
                "body",
            ],
        )
        assert result.exit_code == 1
        assert isinstance(result.exception, FileNotFoundError)
        assert "does not exist!" in str(result.exception)
    finally:
        print("Cleaning up test directory.")
        shutil.rmtree(tempdir)


def test_modify_dicts():
    from wasp.wasp_admin.manage import app

    tempdir = pathlib.Path(tempfile.mkdtemp())
    try:
        channels = tempdir / "tele_waves.json"
        with open(channels, "w") as f:
            json.dump(CHANNELS, f)

        # test downweight
        result = runner.invoke(
            app,
            [
                "modify-dicts",
                str(tempdir),
                "downweight",
                "body",
                "-sc",
                "ABCD:SH",
                "-sc",
                "EFGH:BH1,BHZ",
            ],
        )
        assert result.exit_code == 0
        with open(channels, "r") as f:
            downweighted = json.load(f)
        target = [
            {
                "component": "BHZ",
                "name": "EFGH",
                "trace_weight": 0,
            },
            {
                "component": "BH2",
                "name": "EFGH",
                "trace_weight": 1.0,
            },
            {
                "component": "BH1",
                "name": "EFGH",
                "trace_weight": 0,
            },
            {
                "component": "SH",
                "name": "ABCD",
                "trace_weight": 0,
            },
        ]
        for idx, t in enumerate(target):
            assert downweighted[idx] == t

        # test delete
        result = runner.invoke(
            app,
            [
                "modify-dicts",
                str(tempdir),
                "delete",
                "body",
                "-sc",
                "ABCD:SH",
                "-sc",
                "EFGH:BH2,BHZ",
            ],
        )
        assert result.exit_code == 0
        with open(channels, "r") as f:
            deleted = json.load(f)
        target = [
            {
                "component": "BH1",
                "name": "EFGH",
                "trace_weight": 0,
            }
        ]
        for idx, t in enumerate(target):
            assert deleted[idx] == t
    finally:
        print("Cleaning up test directory.")
        shutil.rmtree(tempdir)


def test_modify_dicts_missing_file():
    from wasp.wasp_admin.manage import app

    result = runner.invoke(
        app,
        [
            "modify-dicts",
            ".",
            "downweight",
            "body",
            "-sc",
            "ABCD:SH",
            "-sc",
            "EFGH:BH2,BHZ",
        ],
    )
    assert result.exit_code == 1
    assert isinstance(result.exception, FileNotFoundError)
    assert "does not exist!" in str(result.exception)


def test_modify_sacs():
    from wasp.wasp_admin.manage import app

    tempdir = pathlib.Path(tempfile.mkdtemp())
    try:
        tele_waves = get_tele_waves_json()
        new_tele_waves = update_manager_file_locations(
            tele_waves, tempdir, replace_dir=str(RESULTS_DIR / "data")
        )
        os.mkdir(pathlib.Path(tempdir) / "P")
        for o, n in zip(tele_waves, new_tele_waves):
            shutil.copyfile(o["file"], n["file"])
        with open(pathlib.Path(tempdir) / "tele_waves.json", "w") as f:
            json.dump(new_tele_waves, f)

        # test downweight
        result = runner.invoke(
            app,
            [
                "modify-sacs",
                str(tempdir),
                "body",
                "-b",
                "RCBR:BHZ=-10",
                "-t",
                "MACI:BHZ=-100",
            ],
        )
        print(result.exception)
        assert result.exit_code == 0
        # check baseline shift
        stream = read(pathlib.Path(tempdir) / "P" / "final_IU_RCBR_BHZ.sac")
        assert np.max(stream[0].data) == 356.4468994140625

        # check time shift
        with open(pathlib.Path(tempdir) / "tele_waves.json", "r") as f:
            updated_tele_waves = json.load(f)
            print("channels_after", updated_tele_waves)
        assert updated_tele_waves[0]["start_signal"] == 595
    finally:
        print("Cleaning up test directory.")
        shutil.rmtree(tempdir)


def test_modify_sacs_missing_file():
    from wasp.wasp_admin.manage import app

    result = runner.invoke(
        app,
        ["modify-sacs", ".", "body", "-b", "ABCD:SH=10", "-p"],
    )
    print(result.stdout, result.exception)
    assert result.exit_code == 1
    assert isinstance(result.exception, FileNotFoundError)
    assert "does not exist!" in str(result.exception)
