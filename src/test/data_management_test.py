import json
import os
import pathlib
import shutil
import tempfile

import numpy as np
from obspy.io.sac import SACTrace

from wasp.data_management import (
    __failsafe,
    __is_number,
    __s2nr,
    __used_stations,
    _dict_trace,
    cgps_traces,
    get_traces_files,
    insar_data,
    select_tele_stations,
    static_data,
    strong_motion_traces,
    tele_body_traces,
    tele_surf_traces,
)

from .testutils import (
    RESULTS_DIR,
    get_cgps_json,
    get_insar_json,
    get_sampling_filter,
    get_static_json,
    get_strong_motion_json,
    get_surf_waves_json,
    get_tele_waves_json,
    get_tensor_info,
    update_manager_file_locations,
)

CGPS_WAVES = get_cgps_json()
CMT = get_tensor_info()
INSAR_DATA = get_insar_json()
SAMPLING_FILTER = get_sampling_filter()
STATIC_DATA = get_static_json()
STRONG_WAVES = get_strong_motion_json()
SURF_WAVES = get_surf_waves_json()
TELE_WAVES = get_tele_waves_json()


def test_cgps_traces():
    tempdir = pathlib.Path(tempfile.mkdtemp())
    try:
        new_cgps_waves = update_manager_file_locations(
            CGPS_WAVES, tempdir, replace_dir=str(RESULTS_DIR / "data")
        )
        os.mkdir(tempdir / "cGPS")
        for o, n in zip(CGPS_WAVES, new_cgps_waves):
            shutil.copyfile(o["file"], n["file"])

        cgps_info = cgps_traces(
            [f["file"] for f in new_cgps_waves],
            CMT,
            SAMPLING_FILTER,
            tempdir,
        )
        with open(tempdir / "cgps_waves.json", "r") as f:
            cgps_json = json.load(f)
        for idx, d in enumerate(new_cgps_waves):
            assert cgps_info[idx] == cgps_json[idx]
            assert cgps_info[idx] == d
    finally:
        shutil.rmtree(tempdir)


def test_get_traces_files():
    assert [
        f.split("/")[-1] for f in get_traces_files("tele_body", RESULTS_DIR / "data")
    ] == ["final_G_MPG_BHZ.sac", "final_IU_RCBR_BHZ.sac", "final_IU_MACI_BHZ.sac"]
    assert [
        f.split("/")[-1] for f in get_traces_files("surf_tele", RESULTS_DIR / "data")
    ] == ["final_G_MPG_BHZ.sac", "final_IU_RCBR_BHZ.sac", "final_IU_MACI_BHZ.sac"]
    assert [
        f.split("/")[-1]
        for f in get_traces_files("strong_motion", RESULTS_DIR / "data")
    ] == [
        "final_vel_GO04_HNN_C.sac",
        "final_vel_VA03_HNZ_C1.sac",
        "final_vel_CO03_HNZ_C1.sac",
    ]
    assert [
        f.split("/")[-1] for f in get_traces_files("cgps", RESULTS_DIR / "data")
    ] == ["final_ovll_LXZ.sac", "final_pedr_LXZ.sac", "final_pfrj_LXE.sac"]


def test_failsafe():
    # test cgps
    __failsafe(
        SAMPLING_FILTER["cgps_filter"], SACTrace.read(CGPS_WAVES[0]["file"]), True
    )

    # test not cgps
    __failsafe(SAMPLING_FILTER["tele_filter"], SACTrace.read(TELE_WAVES[0]["file"]))

    # test bad match
    try:
        __failsafe(
            {"low_freq": 1000000, "high_freq": 1000000},
            SACTrace.read(CGPS_WAVES[0]["file"]),
            True,
        )
        raise Exception("This should have failed!")
    except Exception as e:
        assert "Selected filter doesn't match filter applied to data" in str(e)
    try:
        __failsafe(
            {"low_freq": 1000000, "high_freq": 1000000},
            SACTrace.read(TELE_WAVES[0]["file"]),
        )
        raise Exception("This should have failed!")
    except Exception as e:
        assert "Selected filter doesn't match filter applied to data" in str(e)


def test_insar_data():
    tempdir = pathlib.Path(tempfile.mkdtemp())
    try:
        new_insar = update_manager_file_locations(
            INSAR_DATA, tempdir, replace_dir=str(RESULTS_DIR / "NP1"), file_key="name"
        )

        for key in ["ascending", "descending"]:
            for o, n in zip(INSAR_DATA[key], new_insar[key]):
                shutil.copyfile(o["name"], n["name"])

        insar_info = insar_data(
            [f["name"] for f in new_insar["ascending"]],
            [f["name"] for f in new_insar["descending"]],
            directory=tempdir,
        )
        with open(tempdir / "insar_data.json", "r") as f:
            insar_json = json.load(f)
        for idx, d in enumerate(new_insar["ascending"]):
            assert insar_info["ascending"][idx] == insar_json["ascending"][idx]
            assert insar_info["ascending"][idx] == d
        for idx, d in enumerate(new_insar["descending"]):
            assert insar_info["descending"][idx] == insar_json["descending"][idx]
            assert insar_info["descending"][idx] == d
    finally:
        shutil.rmtree(tempdir)


def test_is_number():
    assert __is_number("INVALID") == False
    assert __is_number("999") == True


def test_dict_trace():
    assert _dict_trace(
        "test_file",
        "station_name",
        "test_channel",
        20,
        100,
        10,
        50,
        2,
        1.0,
        0.0,
        "test_synthetic",
        "test_observed",
        [39, -110],
        True,
    ) == {
        "file": "test_file",
        "dt": 10,
        "duration": 50,
        "wavelet_weight": 0,
        "start_signal": 2,
        "trace_weight": 1.0,
        "synthetic": "test_synthetic",
        "observed": "test_observed",
        "name": "station_name",
        "component": "test_channel",
        "azimuth": 20,
        "distance": 100,
        "location": [39, -110],
        "derivative": True,
    }


def test_select_tele_stations():
    file_dicts = TELE_WAVES + SURF_WAVES
    files = [f["file"] for f in file_dicts]
    love = [f.split("/")[-1] for f in select_tele_stations(files, "Love", CMT)]
    p = [f.split("/")[-1] for f in select_tele_stations(files, "P", CMT)]
    rayleigh = [f.split("/")[-1] for f in select_tele_stations(files, "Rayleigh", CMT)]
    sh = [f.split("/")[-1] for f in select_tele_stations(files, "SH", CMT)]
    assert love == [
        "final_G_MPG_BHZ.sac",
        "final_IU_MACI_BHZ.sac",
        "final_IU_RCBR_BHZ.sac",
    ]
    assert p == [
        "final_G_MPG_BHZ.sac",
        "final_IU_MACI_BHZ.sac",
        "final_IU_RCBR_BHZ.sac",
    ]
    assert rayleigh == [
        "final_G_MPG_BHZ.sac",
        "final_IU_MACI_BHZ.sac",
        "final_IU_RCBR_BHZ.sac",
    ]
    assert sh == [
        "final_G_MPG_BHZ.sac",
        "final_IU_MACI_BHZ.sac",
        "final_IU_RCBR_BHZ.sac",
    ]


def test_s2nr():
    file_dicts = TELE_WAVES + SURF_WAVES
    snrs = []
    for window, phase in zip([1500, 100, 1500, 200], ["Love", "P", "Rayleigh", "SH"]):
        for f in file_dicts:
            fpath = f["file"]
            snrs += [__s2nr(fpath, phase, window)]
    np.testing.assert_array_almost_equal(
        snrs,
        [
            20,
            20,
            20,
            5.405919,
            4.297204,
            6.7423224,
            13.786471,
            20,
            8.615324,
            1.1287758,
            1.3141725,
            1.4124081,
            20,
            20,
            20,
            5.405919,
            4.297204,
            6.7423224,
            20,
            20,
            20,
            1.0934643,
            5.136194,
            1.6105795,
        ],
        decimal=6,
    )


# def test_static_data():
#     tempdir = pathlib.Path(tempfile.mkdtemp())
#     try:
#         new_cgps_waves = update_manager_file_locations(
#             CGPS_WAVES, tempdir, replace_dir=str(RESULTS_DIR / "data")
#         )
#         os.mkdir(tempdir / "cGPS")
#         for o, n in zip(CGPS_WAVES, new_cgps_waves):
#             shutil.copyfile(o["file"], n["file"])

#         with open(tempdir / "cgps_waves.json", "w") as f:
#             json.dump(new_cgps_waves, f)

#         static_info1 = static_data(
#             CMT,
#             directory=tempdir,
#         )
#         with open(tempdir / "static_data.json", "r") as f:
#             static_json1 = json.load(f)
#         target_no_gps_data = [
#             {
#                 "azimuth": None,
#                 "component": None,
#                 "data_error": ["2.3578446", 0, 0],
#                 "derivative": False,
#                 "distance": None,
#                 "dt": None,
#                 "duration": 1,
#                 "file": "None",
#                 "location": [-30.6037540435791, -71.20388793945312],
#                 "name": "ovll",
#                 "observed": ["-11.523461", 0, 0],
#                 "start_signal": 10,
#                 "synthetic": None,
#                 "trace_weight": ["0.5", 0, 0],
#                 "wavelet_weight": None,
#             },
#             {
#                 "azimuth": None,
#                 "component": None,
#                 "data_error": [0, 0, "2.5028799"],
#                 "derivative": False,
#                 "distance": None,
#                 "dt": None,
#                 "duration": 1,
#                 "file": "None",
#                 "location": [-30.674747467041016, -71.63543701171875],
#                 "name": "pfrj",
#                 "observed": [0, 0, "-140.64612"],
#                 "start_signal": 10,
#                 "synthetic": None,
#                 "trace_weight": [0, 0, "1.0"],
#                 "wavelet_weight": None,
#             },
#             {
#                 "azimuth": None,
#                 "component": None,
#                 "data_error": ["2.9627585", 0, 0],
#                 "derivative": False,
#                 "distance": None,
#                 "dt": None,
#                 "duration": 1,
#                 "file": "None",
#                 "location": [-30.838970184326172, -70.68913269042969],
#                 "name": "pedr",
#                 "observed": ["-3.6849258", 0, 0],
#                 "start_signal": 10,
#                 "synthetic": None,
#                 "trace_weight": ["0.5", 0, 0],
#                 "wavelet_weight": None,
#             },
#         ]
#         for idx, d in enumerate(target_no_gps_data):
#             print(static_json1[idx], static_info1[idx], d)
#             assert static_json1[idx] == static_info1[idx] == d

#         shutil.copyfile(RESULTS_DIR / "NP1" / "gps_data", tempdir / "gps_data")

#         static_info2 = static_data(
#             CMT,
#             directory=tempdir,
#         )
#         with open(tempdir / "static_data.json", "r") as f:
#             static_json2 = json.load(f)
#         for idx, d in enumerate(STATIC_DATA):
#             for key in [
#                 "duration",
#                 "start_signal",
#                 "trace_weight",
#                 "observed",
#                 "name",
#                 "component",
#                 "azimuth",
#                 "distance",
#                 "location",
#                 "derivative",
#                 "data_error",
#             ]:
#                 assert static_info2[idx][key] == static_json2[idx][key]
#                 assert static_info2[idx][key] == d[key]
#     finally:
#         shutil.rmtree(tempdir)


def test_strong_motion_traces():
    tempdir = pathlib.Path(tempfile.mkdtemp())
    try:
        new_strong_waves = update_manager_file_locations(
            STRONG_WAVES, tempdir, replace_dir=str(RESULTS_DIR / "data")
        )
        os.mkdir(tempdir / "STR")
        for o, n in zip(STRONG_WAVES, new_strong_waves):
            shutil.copyfile(o["file"], n["file"])

        strong_info = strong_motion_traces(
            [f["file"] for f in new_strong_waves],
            CMT,
            SAMPLING_FILTER,
            tempdir,
        )
        with open(tempdir / "strong_motion_waves.json", "r") as f:
            strong_json = json.load(f)
        for idx, d in enumerate(new_strong_waves):
            assert strong_info[idx] == strong_json[idx]
            assert strong_info[idx] == d
        with open(tempdir / "outlier_strong_motion_waves.json", "r") as f:
            outlier_json = json.load(f)
        assert outlier_json == []
    finally:
        shutil.rmtree(tempdir)


def test_tele_body_traces():
    tempdir = pathlib.Path(tempfile.mkdtemp())
    try:
        new_tele_waves = update_manager_file_locations(
            TELE_WAVES, tempdir, replace_dir=str(RESULTS_DIR / "data")
        )
        os.mkdir(tempdir / "P")
        for o, n in zip(TELE_WAVES, new_tele_waves):
            shutil.copyfile(o["file"], n["file"])

        tele_info = tele_body_traces(
            [f["file"] for f in new_tele_waves],
            CMT,
            SAMPLING_FILTER,
            tempdir,
        )
        with open(tempdir / "tele_waves.json", "r") as f:
            tele_json = json.load(f)
        for idx, d in enumerate(new_tele_waves):
            assert tele_info[idx] == tele_json[idx]
            assert tele_info[idx] == d
    finally:
        shutil.rmtree(tempdir)


def test_tele_surf_traces():
    tempdir = pathlib.Path(tempfile.mkdtemp())
    try:
        new_surf_waves = update_manager_file_locations(
            SURF_WAVES, tempdir, replace_dir=str(RESULTS_DIR / "data")
        )
        os.mkdir(tempdir / "LONG")
        print(new_surf_waves)
        for o, n in zip(SURF_WAVES, new_surf_waves):
            shutil.copyfile(o["file"], n["file"])

        surf_info = tele_surf_traces(
            [f["file"] for f in new_surf_waves],
            CMT,
            SAMPLING_FILTER,
            tempdir,
        )
        with open(tempdir / "surf_waves.json", "r") as f:
            surf_json = json.load(f)
        for idx, d in enumerate(new_surf_waves):
            assert surf_info[idx] == surf_json[idx]
            assert surf_info[idx] == d
    finally:
        shutil.rmtree(tempdir)


def test_used_stations():
    file_dicts = TELE_WAVES + SURF_WAVES + STRONG_WAVES
    files = [f["file"] for f in file_dicts]
    assert __used_stations(2, files, CMT) == 5
    assert __used_stations(4, files, CMT) == 5
    assert __used_stations(8, files, CMT) == 5
    assert __used_stations(10, files, CMT) == 4
    assert __used_stations(12, files, CMT) == 4
