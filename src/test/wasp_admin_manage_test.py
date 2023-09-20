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
    DATA_DIR,
    END_TO_END_DIR,
    HOME,
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


def test_create_ff():
    from wasp.wasp_admin.manage import app

    tempdir = pathlib.Path(tempfile.mkdtemp())
    try:
        shutil.copyfile(
            END_TO_END_DIR / "info" / "20003k7a_cmt_CMT",
            tempdir / "20003k7a_cmt_CMT",
        )
        shutil.copyfile(
            RESULTS_DIR / "NP1" / "Event_mult.in",
            tempdir / "Event_mult.in",
        )
        result = runner.invoke(
            app,
            [
                "create-ff",
                str(tempdir),
                str(tempdir / "20003k7a_cmt_CMT"),
                "6.613912311529926",
                "19.280827965117993",
                "109.27817171619564",
                "-t",
                "body",
            ],
        )
        assert result.exit_code == 0
        with open(tempdir / "segments_data.json") as d:
            data = json.load(d)
        with open(RESULTS_DIR / "NP1" / "segments_data.json") as t:
            target = json.load(t)
        assert data == target
    finally:
        print("Cleaning up test directory.")
        shutil.rmtree(tempdir)


def test_create_ff():
    from wasp.wasp_admin.manage import app

    tempdir = pathlib.Path(tempfile.mkdtemp())
    try:
        shutil.copyfile(
            END_TO_END_DIR / "info" / "20003k7a_cmt_CMT",
            tempdir / "20003k7a_cmt_CMT",
        )
        shutil.copyfile(
            RESULTS_DIR / "NP1" / "Event_mult.in",
            tempdir / "Event_mult.in",
        )
        result = runner.invoke(
            app,
            [
                "create-ff",
                str(tempdir),
                str(tempdir / "20003k7a_cmt_CMT"),
                "6.613912311529926",
                "19.280827965117993",
                "109.27817171619564",
                "-t",
                "body",
            ],
        )
        assert result.exit_code == 0
        with open(tempdir / "segments_data.json") as d:
            data = json.load(d)
        with open(RESULTS_DIR / "NP1" / "segments_data.json") as t:
            target = json.load(t)
        assert data == target
    finally:
        print("Cleaning up test directory.")
        shutil.rmtree(tempdir)


def test_eventmult_to_json():
    from wasp.wasp_admin.manage import app

    tempdir = pathlib.Path(tempfile.mkdtemp())
    try:
        shutil.copyfile(
            END_TO_END_DIR / "info" / "20003k7a_cmt_CMT",
            tempdir / "20003k7a_cmt_CMT",
        )
        shutil.copyfile(
            RESULTS_DIR / "NP1" / "Event_mult.in",
            tempdir / "Event_mult.in",
        )
        result = runner.invoke(
            app,
            [
                "eventmult-to-json",
                str(tempdir),
            ],
        )
        assert result.exit_code == 0
        with open(tempdir / "segments_data.json") as d:
            data = json.load(d)
        with open(RESULTS_DIR / "NP1" / "segments_data.json") as t:
            target = json.load(t)
        assert data == target
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


def test_sampling_filtering():
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
                "sampling-filtering",
                str(tempdir),
                str(tempdir / "20003k7a_cmt_CMT"),
            ],
        )
        assert result.exit_code == 0
        with open(RESULTS_DIR / "NP1" / "sampling_filter.json") as t:
            target = json.load(t)
        with open(tempdir / "sampling_filter.json") as d:
            data = json.load(d)
        assert data == target
    finally:
        print("Cleaning up test directory.")
        shutil.rmtree(tempdir)


def test_static_srf():
    from wasp.wasp_admin.manage import app

    tempdir = pathlib.Path(tempfile.mkdtemp())
    try:
        shutil.copyfile(RESULTS_DIR / "NP1" / "Solucion.txt", tempdir / "Solucion.txt")
        shutil.copyfile(
            RESULTS_DIR / "NP1" / "strong_motion_waves.json",
            tempdir / "strong_motion_waves.json",
        )
        shutil.copyfile(
            RESULTS_DIR / "NP1" / "cgps_waves.json",
            tempdir / "cgps_waves.json",
        )
        shutil.copyfile(
            RESULTS_DIR / "NP1" / "static_data.json",
            tempdir / "static_data.json",
        )
        shutil.copyfile(
            RESULTS_DIR / "NP1" / "tele_waves.json",
            tempdir / "tele_waves.json",
        )
        shutil.copyfile(
            RESULTS_DIR / "NP1" / "surf_waves.json",
            tempdir / "surf_waves.json",
        )
        shutil.copyfile(
            RESULTS_DIR / "NP1" / "insar_data.json",
            tempdir / "insar_data.json",
        )
        shutil.copyfile(
            RESULTS_DIR / "NP1" / "insar_data.txt",
            tempdir / "insar_data.txt",
        )
        shutil.copyfile(
            END_TO_END_DIR / "info" / "20003k7a_cmt_CMT", tempdir / "20003k7a_cmt_CMT"
        )
        shutil.copyfile(
            RESULTS_DIR / "NP1" / "velmodel_data.json", tempdir / "velmodel_data.json"
        )
        shutil.copyfile(
            RESULTS_DIR / "NP1" / "segments_data.json", tempdir / "segments_data.json"
        )
        shutil.copyfile(
            RESULTS_DIR / "NP1" / "strong_motion_waves.json",
            tempdir / "strong_motion_waves.json",
        )
        shutil.copyfile(
            RESULTS_DIR / "NP1" / "tele_waves.json",
            tempdir / "tele_waves.json",
        )
        shutil.copyfile(
            RESULTS_DIR / "NP1" / "surf_waves.json",
            tempdir / "surf_waves.json",
        )
        shutil.copyfile(
            RESULTS_DIR / "NP1" / "cgps_waves.json",
            tempdir / "cgps_waves.json",
        )

        result = runner.invoke(
            app,
            [
                "static-to-srf",
                str(tempdir),
                str(tempdir / "20003k7a_cmt_CMT"),
                "-t",
                "cgps",
                "-t",
                "gps",
                "-t",
                "insar",
                "-t",
                "strong",
                "-t",
                "surf",
                "-t",
                "body",
            ],
        )
        result = runner.invoke(
            app,
            [
                "static-to-srf",
                str(tempdir),
                str(tempdir / "20003k7a_cmt_CMT"),
                "-t",
                "cgps",
                "-t",
                "gps",
                "-t",
                "insar",
                "-t",
                "strong",
                "-t",
                "surf",
                "-t",
                "body",
                "-s",
                str(tempdir / "segments_data.json"),
                "-v",
                str(tempdir / "velmodel_data.json"),
            ],
        )
        print(result.stdout, result.exception)
        assert result.exit_code == 0
        with open(tempdir / "srf_sol_file.txt") as f:
            data = f.readlines()[:40]
        target = [
            "2.0\n",
            "#\n",
            "# Data :\tBODY\tSURF\tSTRONG\tcGPS\tGPS\tInSAR\tDART\tTRIL\tLEVEL\tOTHER\n",
            "# NoS  :\t19\t19\t3\t3\t10\t2166\t0\t0\t0\t0\n",
            "# PHImx :\t186.81\t186.81\t246.51\t312.71\t183.06\t2\t0.00\t0.0\t0.0\t0.0\n",
            "# Rmin :\t40.80\t40.80\t1.11\t0.89\t0.36\t--\t0.00\t0.0\t0.0\t0.0\n",
            "# ---------------------------------- FINITE-SOURCE RUPTURE MODEL ---------------------------------- \n",
            "#\n",
            "# Event : NEAR COAST OF CENTRAL CHILE 2015-09-16T22:54:32 NEIC\n",
            "# EventTAG: 2015-09-16T22:54:32\n",
            "#\n",
            "# Loc  : LAT = -31.57  LON = -71.67  DEP = 22.4\n",
            "#\n",
            "# ----------------------------------  ---------------------------------- \n",
            "#\n",
            "# VELOCITY-DENSITY STRUCTURE\n",
            "# No. of layers = 6\n",
            "#\n",
            "# DEPTH P_VEL S_VEL DENS QP QS\n",
            "# [km] [km/s] [km/s] [g/cm^3]\n",
            "# 0.00 3.35 1.44 2.04  1200.0  600.0\n",
            "# 0.07 6.23 3.61 2.71  1200.0  600.0\n",
            "# 12.15 6.75 3.87 2.83  1200.0  600.0\n",
            "# 25.09 7.65 4.36 2.97  1200.0  600.0\n",
            "# 40.90 8.08 4.47 3.38  1200.0  500.0\n",
            "# 236.90 8.59 4.66 3.45  360.0  140.0\n",
            "#\n",
            "# ----------------------------------  ---------------------------------- \n",
            "# 21/8/2023 created by degoldberg@usgs.gov\n",
            "#\n",
            "# SOURCE MODEL PARAMETERS\n",
            "PLANE 1\n",
            "-72.1958 \t -31.0308 \t 23 \t 9 \t 412.33 \t 134.32 \n",
            "6.61 \t 19.28 \t 0.22 \t -53.67 \t 67.16 \n",
            "POINTS 207\n",
            "-72.4357 -32.7937 2.69 6.61 19.28 267.56 123.94 1.0 -1 -1\n",
            "121.47 67.7587 31 0 0 0 0 \n",
            "  0.000000e+00  4.935634e-02  1.952683e-01  4.313587e-01  7.473094e-01  1.129312e+00\n",
            "  1.560670e+00  2.022533e+00  2.494714e+00  2.956577e+00  3.387935e+00  3.769938e+00\n",
            "  4.085888e+00  4.321979e+00  4.467891e+00  4.517247e+00  4.467891e+00  4.321979e+00\n",
        ]
        for idx, target_line in enumerate(target):
            if "@usgs.gov" in target_line:
                continue
            assert data[idx] == target_line
    finally:
        shutil.rmtree(tempdir)


def test_static_fsp():
    from wasp.wasp_admin.manage import app

    tempdir = pathlib.Path(tempfile.mkdtemp())
    try:
        shutil.copyfile(RESULTS_DIR / "NP1" / "Solucion.txt", tempdir / "Solucion.txt")
        shutil.copyfile(
            RESULTS_DIR / "NP1" / "strong_motion_waves.json",
            tempdir / "strong_motion_waves.json",
        )
        shutil.copyfile(
            RESULTS_DIR / "NP1" / "cgps_waves.json",
            tempdir / "cgps_waves.json",
        )
        shutil.copyfile(
            RESULTS_DIR / "NP1" / "static_data.json",
            tempdir / "static_data.json",
        )
        shutil.copyfile(
            RESULTS_DIR / "NP1" / "tele_waves.json",
            tempdir / "tele_waves.json",
        )
        shutil.copyfile(
            RESULTS_DIR / "NP1" / "surf_waves.json",
            tempdir / "surf_waves.json",
        )
        shutil.copyfile(
            RESULTS_DIR / "NP1" / "insar_data.json",
            tempdir / "insar_data.json",
        )
        shutil.copyfile(
            RESULTS_DIR / "NP1" / "insar_data.txt",
            tempdir / "insar_data.txt",
        )
        shutil.copyfile(
            END_TO_END_DIR / "info" / "20003k7a_cmt_CMT", tempdir / "20003k7a_cmt_CMT"
        )
        shutil.copyfile(
            RESULTS_DIR / "NP1" / "velmodel_data.json", tempdir / "velmodel_data.json"
        )
        shutil.copyfile(
            RESULTS_DIR / "NP1" / "segments_data.json", tempdir / "segments_data.json"
        )
        shutil.copyfile(
            RESULTS_DIR / "NP1" / "strong_motion_waves.json",
            tempdir / "strong_motion_waves.json",
        )
        shutil.copyfile(
            RESULTS_DIR / "NP1" / "tele_waves.json",
            tempdir / "tele_waves.json",
        )
        shutil.copyfile(
            RESULTS_DIR / "NP1" / "surf_waves.json",
            tempdir / "surf_waves.json",
        )
        shutil.copyfile(
            RESULTS_DIR / "NP1" / "cgps_waves.json",
            tempdir / "cgps_waves.json",
        )

        result = runner.invoke(
            app,
            [
                "static-to-fsp",
                str(tempdir),
                str(tempdir / "20003k7a_cmt_CMT"),
                "-t",
                "cgps",
                "-t",
                "gps",
                "-t",
                "insar",
                "-t",
                "strong",
                "-t",
                "surf",
                "-t",
                "body",
            ],
        )
        result = runner.invoke(
            app,
            [
                "static-to-fsp",
                str(tempdir),
                str(tempdir / "20003k7a_cmt_CMT"),
                "-t",
                "cgps",
                "-t",
                "gps",
                "-t",
                "insar",
                "-t",
                "strong",
                "-t",
                "surf",
                "-t",
                "body",
                "-s",
                str(tempdir / "segments_data.json"),
                "-v",
                str(tempdir / "velmodel_data.json"),
            ],
        )
        assert result.exit_code == 0
        with open(tempdir / "fsp_sol_file.txt") as f:
            data = f.readlines()[41:]
        with open(RESULTS_DIR / "NP1" / "fsp_sol_file.txt") as f:
            target = f.readlines()[41:]
        for l, t in zip(data, target):
            assert l == t
    finally:
        shutil.rmtree(tempdir)


def test_velmodel_from_tensor():
    from wasp.wasp_admin.manage import app

    tempdir = pathlib.Path(tempfile.mkdtemp())
    try:
        shutil.copyfile(RESULTS_DIR / "NP1" / "Solucion.txt", tempdir / "Solucion.txt")
        shutil.copyfile(
            RESULTS_DIR / "NP1" / "strong_motion_waves.json",
            tempdir / "strong_motion_waves.json",
        )
        shutil.copyfile(
            RESULTS_DIR / "NP1" / "cgps_waves.json",
            tempdir / "cgps_waves.json",
        )
        shutil.copyfile(
            RESULTS_DIR / "NP1" / "static_data.json",
            tempdir / "static_data.json",
        )
        shutil.copyfile(
            RESULTS_DIR / "NP1" / "tele_waves.json",
            tempdir / "tele_waves.json",
        )
        shutil.copyfile(
            RESULTS_DIR / "NP1" / "surf_waves.json",
            tempdir / "surf_waves.json",
        )
        shutil.copyfile(
            RESULTS_DIR / "NP1" / "insar_data.json",
            tempdir / "insar_data.json",
        )
        shutil.copyfile(
            RESULTS_DIR / "NP1" / "insar_data.txt",
            tempdir / "insar_data.txt",
        )
        shutil.copyfile(
            END_TO_END_DIR / "info" / "20003k7a_cmt_CMT", tempdir / "20003k7a_cmt_CMT"
        )
        shutil.copyfile(
            RESULTS_DIR / "NP1" / "velmodel_data.json", tempdir / "velmodel_data.json"
        )
        shutil.copyfile(
            RESULTS_DIR / "NP1" / "segments_data.json", tempdir / "segments_data.json"
        )
        shutil.copyfile(
            RESULTS_DIR / "NP1" / "strong_motion_waves.json",
            tempdir / "strong_motion_waves.json",
        )
        shutil.copyfile(
            RESULTS_DIR / "NP1" / "tele_waves.json",
            tempdir / "tele_waves.json",
        )
        shutil.copyfile(
            RESULTS_DIR / "NP1" / "surf_waves.json",
            tempdir / "surf_waves.json",
        )
        shutil.copyfile(
            RESULTS_DIR / "NP1" / "cgps_waves.json",
            tempdir / "cgps_waves.json",
        )

        result = runner.invoke(
            app,
            [
                "static-to-fsp",
                str(tempdir),
                str(tempdir / "20003k7a_cmt_CMT"),
                "-t",
                "cgps",
                "-t",
                "gps",
                "-t",
                "insar",
                "-t",
                "strong",
                "-t",
                "surf",
                "-t",
                "body",
            ],
        )
        result = runner.invoke(
            app,
            [
                "static-to-fsp",
                str(tempdir),
                str(tempdir / "20003k7a_cmt_CMT"),
                "-t",
                "cgps",
                "-t",
                "gps",
                "-t",
                "insar",
                "-t",
                "strong",
                "-t",
                "surf",
                "-t",
                "body",
                "-s",
                str(tempdir / "segments_data.json"),
                "-v",
                str(tempdir / "velmodel_data.json"),
            ],
        )
        assert result.exit_code == 0
        with open(tempdir / "fsp_sol_file.txt") as f:
            data = f.readlines()[41:]
        with open(RESULTS_DIR / "NP1" / "fsp_sol_file.txt") as f:
            target = f.readlines()[41:]
        for l, t in zip(data, target):
            assert l == t
    finally:
        shutil.rmtree(tempdir)


def test_tensor_from_gcmt():
    from wasp.wasp_admin.manage import app

    tempdir = pathlib.Path(tempfile.mkdtemp())
    try:
        shutil.copyfile(
            END_TO_END_DIR / "info" / "20003k7a_cmt_CMT", tempdir / "20003k7a_cmt_CMT"
        )

        result = runner.invoke(
            app,
            [
                "tensor-from-gcmt",
                str(tempdir / "20003k7a_cmt_CMT"),
                "-d",
                str(tempdir),
            ],
        )
        assert result.exit_code == 0
        with open(tempdir / "tensor_info.json") as f:
            data = json.load(f)
        with open(RESULTS_DIR / "NP1" / "tensor_info.json") as f:
            target = json.load(f)
        del target["timedelta"]
        del data["timedelta"]
        assert data == target
    finally:
        shutil.rmtree(tempdir)


def test_velmodel_from_tensor():
    from wasp.wasp_admin.manage import app

    tempdir = pathlib.Path(tempfile.mkdtemp())
    try:
        shutil.copyfile(RESULTS_DIR / "NP1" / "Solucion.txt", tempdir / "Solucion.txt")
        shutil.copyfile(
            RESULTS_DIR / "NP1" / "strong_motion_waves.json",
            tempdir / "strong_motion_waves.json",
        )
        shutil.copyfile(
            RESULTS_DIR / "NP1" / "cgps_waves.json",
            tempdir / "cgps_waves.json",
        )
        shutil.copyfile(
            RESULTS_DIR / "NP1" / "static_data.json",
            tempdir / "static_data.json",
        )
        shutil.copyfile(
            RESULTS_DIR / "NP1" / "tele_waves.json",
            tempdir / "tele_waves.json",
        )
        shutil.copyfile(
            RESULTS_DIR / "NP1" / "surf_waves.json",
            tempdir / "surf_waves.json",
        )
        shutil.copyfile(
            RESULTS_DIR / "NP1" / "insar_data.json",
            tempdir / "insar_data.json",
        )
        shutil.copyfile(
            RESULTS_DIR / "NP1" / "insar_data.txt",
            tempdir / "insar_data.txt",
        )
        shutil.copyfile(
            END_TO_END_DIR / "info" / "20003k7a_cmt_CMT", tempdir / "20003k7a_cmt_CMT"
        )
        shutil.copyfile(
            RESULTS_DIR / "NP1" / "velmodel_data.json", tempdir / "velmodel_data.json"
        )
        shutil.copyfile(
            RESULTS_DIR / "NP1" / "segments_data.json", tempdir / "segments_data.json"
        )
        shutil.copyfile(
            RESULTS_DIR / "NP1" / "strong_motion_waves.json",
            tempdir / "strong_motion_waves.json",
        )
        shutil.copyfile(
            RESULTS_DIR / "NP1" / "tele_waves.json",
            tempdir / "tele_waves.json",
        )
        shutil.copyfile(
            RESULTS_DIR / "NP1" / "surf_waves.json",
            tempdir / "surf_waves.json",
        )
        shutil.copyfile(
            RESULTS_DIR / "NP1" / "cgps_waves.json",
            tempdir / "cgps_waves.json",
        )

        result = runner.invoke(
            app,
            [
                "static-to-fsp",
                str(tempdir),
                str(tempdir / "20003k7a_cmt_CMT"),
                "-t",
                "cgps",
                "-t",
                "gps",
                "-t",
                "insar",
                "-t",
                "strong",
                "-t",
                "surf",
                "-t",
                "body",
            ],
        )
        result = runner.invoke(
            app,
            [
                "static-to-fsp",
                str(tempdir),
                str(tempdir / "20003k7a_cmt_CMT"),
                "-t",
                "cgps",
                "-t",
                "gps",
                "-t",
                "insar",
                "-t",
                "strong",
                "-t",
                "surf",
                "-t",
                "body",
                "-s",
                str(tempdir / "segments_data.json"),
                "-v",
                str(tempdir / "velmodel_data.json"),
            ],
        )
        assert result.exit_code == 0
        with open(tempdir / "fsp_sol_file.txt") as f:
            data = f.readlines()[41:]
        with open(RESULTS_DIR / "NP1" / "fsp_sol_file.txt") as f:
            target = f.readlines()[41:]
        for l, t in zip(data, target):
            assert l == t
    finally:
        shutil.rmtree(tempdir)


def test_velmodel_to_json():
    from wasp.wasp_admin.manage import app

    tempdir = pathlib.Path(tempfile.mkdtemp())
    try:
        with open(tempdir / "velmodel", "w") as v:
            v.write(
                "6\n"
                "3.35 1.44 2.0447 0.074 1200.0 600.0\n"
                "6.23 3.61 2.7074199 12.076 1200.0 600.0\n"
                "6.75 3.87 2.83123 12.945 1200.0 600.0\n"
                "7.65 4.36 2.96935 15.8029995 1200.0 600.0\n"
                "8.08 4.473 3.3754 196.0 1200.0 500.0\n"
                "8.594 4.657 3.4465 36.0 360.0 140.0\n"
            )

        result = runner.invoke(
            app,
            ["velmodel-to-json", str(tempdir), str(tempdir / "velmodel")],
        )

        assert result.exit_code == 0
        with open(tempdir / "velmodel_data.json") as v:
            data = json.load(v)

        target = {
            "dens": ["2.0447", "2.7074199", "2.83123", "2.96935", "3.3754", "3.4465"],
            "p_vel": ["3.35", "6.23", "6.75", "7.65", "8.08", "8.594"],
            "qa": ["1200", "1200", "1200", "1200", "1200", "360"],
            "qb": ["600", "600", "600", "600", "500", "140"],
            "s_vel": ["1.44", "3.61", "3.87", "4.36", "4.473", "4.657"],
            "thick": ["0.074", "12.076", "12.945", "15.8029995", "196.0", "36.0"],
        }
        for t in target:
            assert data[t] == target[t]
    finally:
        print("Cleaning up test directory.")
        shutil.rmtree(tempdir)
