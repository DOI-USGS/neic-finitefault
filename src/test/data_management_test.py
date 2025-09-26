import json
import os
import pathlib
import shutil
import tempfile

import numpy as np
from obspy.io.sac import SACTrace  # type:ignore

from wasp.data_management import (
    __failsafe,
    __is_number,
    __s2nr,
    __used_stations,
    __wavelets_dart,
    _dict_trace,
    cgnss_traces,
    duration_dart,
    duration_strong_motion,
    duration_tele_waves,
    filling_data_dicts,
    get_traces_files,
    imagery_data,
    select_tele_stations,
    static_data,
    strong_motion_traces,
    tele_body_traces,
    tele_surf_traces,
    wavelets_body_waves,
    wavelets_strong_motion,
    wavelets_surf_tele,
)
from wasp.management import _distazbaz

from .testutils import (
    RESULTS_DIR,
    get_cgnss_json,
    get_imagery_json,
    get_sampling_filter,
    get_static_json,
    get_strong_motion_json,
    get_surf_waves_json,
    get_tele_waves_json,
    get_tensor_info,
    update_manager_file_locations,
)

CGNSS_WAVES = get_cgnss_json()
CMT = get_tensor_info()
IMAGERY_DATA = get_imagery_json()
SAMPLING_FILTER = get_sampling_filter()
STATIC_DATA = get_static_json()
STRONG_WAVES = get_strong_motion_json()
SURF_WAVES = get_surf_waves_json()
TELE_WAVES = get_tele_waves_json()


def test_cgnss_traces():
    tempdir = pathlib.Path(tempfile.mkdtemp())
    try:
        new_cgnss_waves = update_manager_file_locations(
            CGNSS_WAVES, tempdir, replace_dir=str(RESULTS_DIR / "data")
        )
        os.mkdir(tempdir / "cGNSS")
        for o, n in zip(CGNSS_WAVES, new_cgnss_waves):
            shutil.copyfile(o["file"], n["file"])

        cgnss_info = cgnss_traces(
            [f["file"] for f in new_cgnss_waves],
            CMT,
            SAMPLING_FILTER,
            tempdir,
        )
        with open(tempdir / "cgnss_waves.json", "r") as f:
            cgnss_json = json.load(f)
        print(json.dumps(cgnss_json, indent=4))
        for idx, d in enumerate(new_cgnss_waves):
            for key, target_values in d.items():
                test_values = cgnss_info[idx][key]
                if isinstance(target_values, float):
                    np.testing.assert_almost_equal(
                        target_values, test_values, decimal=5
                    )
                else:
                    assert test_values == target_values
    finally:
        shutil.rmtree(tempdir)


def test_duration_dart():
    depth = CMT["depth"]
    # dont have dart records to made up some distances
    distances = [20, 30, 40]
    arrivals = [np.sqrt(dist**2 + depth**2) / 5 for dist in distances]
    duration = duration_dart(distances, arrivals, CMT, 0.2)
    assert duration == 950


def test_duration_tele_waves():
    for wavelet in TELE_WAVES:
        duration = duration_tele_waves(CMT, wavelet["dt"])
        assert duration == wavelet["duration"]


def test_duration_strong_motion():
    depth = CMT["depth"]
    distances = [
        _distazbaz(w["location"][0], w["location"][1], CMT["lat"], CMT["lon"])[0]
        for w in STRONG_WAVES
    ]
    arrivals = [np.sqrt(dist**2 + depth**2) / 5 for dist in distances]
    for wavelet in STRONG_WAVES:
        duration = duration_strong_motion(distances, arrivals, CMT, 0.4)
        assert duration == wavelet["duration"]


def test_filling_data_dicts():
    tempdir = pathlib.Path(tempfile.mkdtemp())
    try:
        new_surf_waves = update_manager_file_locations(
            SURF_WAVES, tempdir / "data", replace_dir=str(RESULTS_DIR / "data")
        )
        new_tele_waves = update_manager_file_locations(
            TELE_WAVES, tempdir / "data", replace_dir=str(RESULTS_DIR / "data")
        )
        new_strong_waves = update_manager_file_locations(
            STRONG_WAVES,
            tempdir / "data",
            replace_dir=str(RESULTS_DIR / "data"),
        )
        new_cgnss_waves = update_manager_file_locations(
            CGNSS_WAVES, tempdir / "data", replace_dir=str(RESULTS_DIR / "data")
        )
        new_imagery = update_manager_file_locations(
            IMAGERY_DATA, tempdir, replace_dir=str(RESULTS_DIR / "NP1"), file_key="name"
        )
        os.mkdir(tempdir / "data")
        os.mkdir(tempdir / "data" / "STR")
        os.mkdir(tempdir / "data" / "cGNSS")
        os.mkdir(tempdir / "data" / "SH")
        os.mkdir(tempdir / "data" / "P")
        os.mkdir(tempdir / "data" / "LOVE")
        os.mkdir(tempdir / "data" / "RAYLEIGH")
        for a, b, c, d, e, f, g, h in zip(
            SURF_WAVES,
            new_surf_waves,
            TELE_WAVES,
            new_tele_waves,
            STRONG_WAVES,
            new_strong_waves,
            CGNSS_WAVES,
            new_cgnss_waves,
        ):
            shutil.copyfile(a["file"], b["file"])
            shutil.copyfile(c["file"], d["file"])
            shutil.copyfile(e["file"], f["file"])
            shutil.copyfile(g["file"], h["file"])
        for key in ["insar_ascending", "insar_descending"]:
            for o, n in zip(IMAGERY_DATA[key], new_imagery[key]):
                shutil.copyfile(o["name"], n["name"])
        shutil.copyfile(RESULTS_DIR / "NP1" / "gnss_data", tempdir / "gnss_data")
        imagery_files = [w["name"] for w in new_imagery["insar_ascending"]]
        imagery_files += [w["name"] for w in new_imagery["insar_descending"]]
        filling_data_dicts(
            CMT,
            ["cgnss", "gnss", "imagery", "strong", "surf", "body"],
            SAMPLING_FILTER,
            tempdir / "data",
            imagery_files,
            None,
            working_directory=tempdir,
        )
        for fname, target in zip(
            [
                "surf_waves.json",
                "tele_waves.json",
                "strong_motion_waves.json",
                "cgnss_waves.json",
                "imagery_data.json",
                "static_data.json",
            ],
            [
                new_surf_waves,
                new_tele_waves,
                new_strong_waves,
                new_cgnss_waves,
                new_imagery,
                STATIC_DATA,
            ],
        ):
            with open(tempdir / fname, "r") as f:
                props = json.load(f)
            if "imagery" in fname:
                for idx, d in enumerate(target["insar_ascending"]):
                    assert props["insar_ascending"][idx] == d
                for idx, d in enumerate(target["insar_descending"]):
                    assert props["insar_descending"][idx] == d
            elif "static" in fname:
                for idx, d in enumerate(target):
                    for key in [
                        "duration",
                        "start_signal",
                        "trace_weight",
                        "observed",
                        "name",
                        "component",
                        "location",
                        "derivative",
                        "data_error",
                    ]:
                        assert props[idx][key] == d[key]
                    for key in [
                        "azimuth",
                        "distance",
                    ]:
                        np.testing.assert_almost_equal(props[idx][key], d[key], 5)
            else:
                for idx, t in enumerate(target):
                    match = False
                    for idy, d in enumerate(props):
                        if t["name"] == d["name"] and t["component"] == d["component"]:
                            match = True
                            for key, target_values in t.items():
                                test_values = d[key]
                                if isinstance(target_values, float):
                                    np.testing.assert_almost_equal(
                                        test_values, target_values, decimal=5
                                    )
                                else:
                                    assert test_values == target_values
                    if not match:
                        raise Exception(
                            f"No match found for target {t['name']} {t['component']}"
                        )
    finally:
        shutil.rmtree(tempdir)


def test_get_traces_files():
    tele_files = [
        f.split("/")[-1] for f in get_traces_files("body", RESULTS_DIR / "data")
    ]
    assert len(tele_files) == 20
    for target in [
        "processed_G_MPG_BHZ.sac",
        "processed_IU_RCBR_BHZ.sac",
        "processed_IU_MACI_BHZ.sac",
    ]:
        assert target in tele_files
    surf_files = [
        f.split("/")[-1] for f in get_traces_files("surf", RESULTS_DIR / "data")
    ]
    assert len(surf_files) == 20
    for target in [
        "processed_G_MPG_BHZ.sac",
        "processed_IU_RCBR_BHZ.sac",
        "processed_IU_MACI_BHZ.sac",
    ]:
        assert target in surf_files
    sm_files = [
        f.split("/")[-1] for f in get_traces_files("strong", RESULTS_DIR / "data")
    ]
    assert len(sm_files) == 9
    for target in [
        "processed_vel_GO04_HNN_C.sac",
        "processed_vel_VA03_HNZ_C1.sac",
        "processed_vel_CO03_HNZ_C1.sac",
    ]:
        assert target in sm_files
    cgnss_files = [
        f.split("/")[-1] for f in get_traces_files("cgnss", RESULTS_DIR / "data")
    ]
    assert len(cgnss_files) == 9
    for target in [
        "processed_ovll_LXZ.sac",
        "processed_pedr_LXZ.sac",
        "processed_pfrj_LXE.sac",
    ]:
        assert target in cgnss_files


def test_failsafe():
    # test cgnss
    __failsafe(
        SAMPLING_FILTER["cgnss_filter"], SACTrace.read(CGNSS_WAVES[0]["file"]), True
    )

    # test not cgnss
    __failsafe(SAMPLING_FILTER["tele_filter"], SACTrace.read(TELE_WAVES[0]["file"]))

    # test bad match
    try:
        __failsafe(
            {"low_freq": 1000000, "high_freq": 1000000},
            SACTrace.read(CGNSS_WAVES[0]["file"]),
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


def test_imagery_data():
    tempdir = pathlib.Path(tempfile.mkdtemp())
    try:
        new_imagery = update_manager_file_locations(
            IMAGERY_DATA, tempdir, replace_dir=str(RESULTS_DIR / "NP1"), file_key="name"
        )

        for key in ["insar_ascending", "insar_descending"]:
            for o, n in zip(IMAGERY_DATA[key], new_imagery[key]):
                shutil.copyfile(o["name"], n["name"])

        imagery_files = [f["name"] for f in new_imagery["insar_ascending"]]
        imagery_files += [f["name"] for f in new_imagery["insar_descending"]]

        imagery_info = imagery_data(
            imagery_files,
            directory=tempdir,
        )
        with open(tempdir / "imagery_data.json", "r") as f:
            imagery_json = json.load(f)
        for idx, d in enumerate(new_imagery["insar_ascending"]):
            assert (
                imagery_info["insar_ascending"][idx]
                == imagery_json["insar_ascending"][idx]
            )
            assert imagery_info["insar_ascending"][idx] == d
        for idx, d in enumerate(new_imagery["insar_descending"]):
            assert (
                imagery_info["insar_descending"][idx]
                == imagery_json["insar_descending"][idx]
            )
            assert imagery_info["insar_descending"][idx] == d
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
        "processed_G_MPG_BHZ.sac",
        "processed_IU_MACI_BHZ.sac",
        "processed_IU_RCBR_BHZ.sac",
    ]
    assert p == [
        "processed_G_MPG_BHZ.sac",
        "processed_IU_MACI_BHZ.sac",
        "processed_IU_RCBR_BHZ.sac",
    ]
    assert rayleigh == [
        "processed_G_MPG_BHZ.sac",
        "processed_IU_MACI_BHZ.sac",
        "processed_IU_RCBR_BHZ.sac",
    ]
    assert sh == [
        "processed_G_MPG_BHZ.sac",
        "processed_IU_MACI_BHZ.sac",
        "processed_IU_RCBR_BHZ.sac",
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
        decimal=5,
    )


def test_static_data():
    tempdir = pathlib.Path(tempfile.mkdtemp())
    try:
        new_cgnss_waves = update_manager_file_locations(
            CGNSS_WAVES, tempdir, replace_dir=str(RESULTS_DIR / "data")
        )
        os.mkdir(tempdir / "cGNSS")
        for o, n in zip(CGNSS_WAVES, new_cgnss_waves):
            shutil.copyfile(o["file"], n["file"])

        with open(tempdir / "cgnss_waves.json", "w") as f:
            json.dump(new_cgnss_waves, f)

        static_info1 = static_data(
            CMT,
            directory=tempdir,
        )
        with open(tempdir / "static_data.json", "r") as f:
            static_json1 = json.load(f)
        target_no_gnss_data = [
            {
                "azimuth": None,
                "component": None,
                "data_error": ["2.3578446", 0, 0],
                "derivative": False,
                "distance": None,
                "dt": None,
                "duration": 1,
                "file": "None",
                "location": [-30.6037540435791, -71.20388793945312],
                "name": "ovll",
                "observed": ["-11.523461", 0, 0],
                "start_signal": 10,
                "synthetic": None,
                "trace_weight": ["0.5", 0, 0],
                "wavelet_weight": None,
            },
            {
                "azimuth": None,
                "component": None,
                "data_error": [0, 0, "2.5028799"],
                "derivative": False,
                "distance": None,
                "dt": None,
                "duration": 1,
                "file": "None",
                "location": [-30.674747467041016, -71.63543701171875],
                "name": "pfrj",
                "observed": [0, 0, "-140.64612"],
                "start_signal": 10,
                "synthetic": None,
                "trace_weight": [0, 0, "1.0"],
                "wavelet_weight": None,
            },
            {
                "azimuth": None,
                "component": None,
                "data_error": ["2.9627585", 0, 0],
                "derivative": False,
                "distance": None,
                "dt": None,
                "duration": 1,
                "file": "None",
                "location": [-30.838970184326172, -70.68913269042969],
                "name": "pedr",
                "observed": ["-3.6849258", 0, 0],
                "start_signal": 10,
                "synthetic": None,
                "trace_weight": ["0.5", 0, 0],
                "wavelet_weight": None,
            },
        ]

        for idx, d in enumerate(target_no_gnss_data):
            match1 = False
            match2 = False
            for idy, props in enumerate(static_json1):
                if props["name"] == d["name"]:
                    assert static_json1[idy] == static_info1[idy] == d
                    match = True
            if not match:
                raise Exception(f"No data from name, {d['name']}, found!")

        shutil.copyfile(RESULTS_DIR / "NP1" / "gnss_data", tempdir / "gnss_data")

        static_info2 = static_data(
            CMT,
            directory=tempdir,
        )
        with open(tempdir / "static_data.json", "r") as f:
            static_json2 = json.load(f)
        for idx, d in enumerate(STATIC_DATA):
            for key in [
                "duration",
                "start_signal",
                "trace_weight",
                "observed",
                "name",
                "component",
                "azimuth",
                "distance",
                "location",
                "derivative",
                "data_error",
            ]:
                test_values = static_json2[idx][key]
                target_values = static_info2[idx][key]
                print(target_values)
                if isinstance(target_values, list):
                    if isinstance(target_values, list):
                        try:  # some are lists of floats but one is list of components (alpha)
                            target_values = [float(v) for v in target_values]
                            test_values = [float(v) for v in test_values]
                            np.testing.assert_almost_equal(
                                test_values, target_values, decimal=5
                            )
                        except ValueError:
                            assert test_values == target_values
                else:
                    assert test_values == target_values
    finally:
        shutil.rmtree(tempdir)


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
            for key, target_values in d.items():
                test_values = strong_info[idx][key]
                if isinstance(target_values, float):
                    np.testing.assert_almost_equal(
                        target_values, test_values, decimal=5
                    )
                else:
                    assert test_values == target_values
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
            for key, target_values in d.items():
                test_values = tele_info[idx][key]
                if isinstance(target_values, float):
                    np.testing.assert_almost_equal(
                        target_values, test_values, decimal=5
                    )
                else:
                    assert test_values == target_values
    finally:
        shutil.rmtree(tempdir)


def test_tele_surf_traces():
    tempdir = pathlib.Path(tempfile.mkdtemp())
    try:
        new_surf_waves = update_manager_file_locations(
            SURF_WAVES, tempdir, replace_dir=str(RESULTS_DIR / "data")
        )
        os.mkdir(tempdir / "LOVE")
        os.mkdir(tempdir / "RAYLEIGH")
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
            for key, target_values in d.items():
                test_values = surf_info[idx][key]
                if isinstance(target_values, float):
                    np.testing.assert_almost_equal(
                        target_values, test_values, decimal=5
                    )
                else:
                    assert test_values == target_values
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


def test_wavelets_body_waves():
    s1, s2 = SAMPLING_FILTER["wavelet_scales"]
    filter = SAMPLING_FILTER["tele_filter"]
    for wavelet in TELE_WAVES:
        comp = wavelet["component"]
        w0, w1 = wavelets_body_waves(wavelet["duration"], filter, wavelet["dt"], s1, s2)
        if comp == "BHZ":
            assert w0 == wavelet["wavelet_weight"]
        else:
            assert w1 == wavelet["wavelet_weight"]


def test_wavelets_dart():
    assert __wavelets_dart(1, 8) == "0 0 2 2 2 2 2 2\n"


def test_wavelets_surf_waves():
    s1, s2 = SAMPLING_FILTER["wavelet_scales"]
    filter = SAMPLING_FILTER["surf_filter"]
    for wavelet in SURF_WAVES:
        comp = wavelet["component"]
        w = wavelets_surf_tele(filter, s1, s2)
        assert w == wavelet["wavelet_weight"]


def test_wavelets_strong_motion():
    s1, s2 = SAMPLING_FILTER["wavelet_scales"]
    filter = SAMPLING_FILTER["strong_filter"]
    for wavelet in STRONG_WAVES:
        comp = wavelet["component"]
        w = wavelets_strong_motion(wavelet["duration"], filter, wavelet["dt"], s1, s2)
        assert w == wavelet["wavelet_weight"]
    assert (
        wavelets_strong_motion(wavelet["duration"], filter, wavelet["dt"], s1, s2, True)
        == "0 2 2 2 2 2 2 0\n"
    )
