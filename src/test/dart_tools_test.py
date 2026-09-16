import json
import os
import pathlib
import shutil
import tempfile
from itertools import islice
from operator import itemgetter
from unittest import mock

import numpy as np
import pytest
import requests
from obspy import UTCDateTime, read  # type: ignore

from ffm.dart_tools import (
    dart_gf_time_freq_domain,
    dart_high_res_locations,
    get_dart_observations,
    get_subfault_locations,
    make_DART_GFs_DS_SS,
    make_dart_json_no_observations,
    make_ffm_kinematic_dtopo,
    make_ffm_static_dtopo,
    make_topo,
    plot_dtopo_snapshots,
    query_dart_api,
    run_ffm_tsunami,
)
from ffm.management import default_dirs

from .testutils import (
    END_TO_END_DIR,
    HOME,
    RESULTS_DIR,
    get_dart_waves_json,
    get_segments_data,
    get_tensor_info,
)


def test_make_DART_GFs_DS_SS():
    tempdir = tempfile.mkdtemp()
    try:
        tempdir = pathlib.Path(tempdir)
        tensor_info = get_tensor_info()
        segments_data = get_segments_data()
        dart_data = get_dart_waves_json()
        ndip = segments_data["segments"][0]["dip_subfaults"]
        nstk = segments_data["segments"][0]["stk_subfaults"]
        num_subfaults = float(ndip) * float(nstk)
        subfault_start = int(num_subfaults - 1)
        make_DART_GFs_DS_SS(
            tempdir,
            tensor_info,
            segments_data,
            dart_data,
            dt=60,
            tfinal=10800,
            subfault_start=subfault_start,
            time_domain=True,
            plots=False,
        )
        num_dart = len(dart_data)
        for kdart in range(num_dart):
            id_dart = dart_data[kdart]["name"]
            gf_file_freq_dip = f"{id_dart}_dip_gf.txt"
            gf_file_time_dip = (
                f"dart.{id_dart}.{str(subfault_start).rjust(4, '0')}.GF.dip.txt"
            )
            gf_file_freq_stk = f"{id_dart}_stk_gf.txt"
            gf_file_time_stk = (
                f"dart.{id_dart}.{str(subfault_start).rjust(4, '0')}.GF.stk.txt"
            )
            content1_freq_dip = np.genfromtxt(
                tempdir / gf_file_freq_dip, skip_header=1
            ).reshape(-1, 2)
            content2_freq_dip = np.genfromtxt(
                RESULTS_DIR / "Tsunami_Output" / gf_file_freq_dip, skip_header=1
            ).reshape(-1, 2)
            np.testing.assert_almost_equal(
                content1_freq_dip[:, 0], content2_freq_dip[:, 0], decimal=6
            )
            np.testing.assert_almost_equal(
                content1_freq_dip[:, 1], content2_freq_dip[:, 1], decimal=6
            )

            content1_time_dip = np.genfromtxt(
                tempdir / "_GF" / gf_file_time_dip, skip_header=1
            ).reshape(-1, 2)
            content2_time_dip = np.genfromtxt(
                RESULTS_DIR / "Tsunami_Output" / gf_file_time_dip, skip_header=1
            ).reshape(-1, 2)
            np.testing.assert_almost_equal(
                content1_time_dip[:, 0], content2_time_dip[:, 0], decimal=6
            )
            np.testing.assert_almost_equal(
                content1_time_dip[:, 1], content2_time_dip[:, 1], decimal=6
            )

            content1_freq_strk = np.genfromtxt(
                tempdir / gf_file_freq_stk, skip_header=1
            ).reshape(-1, 2)
            content2_freq_strk = np.genfromtxt(
                RESULTS_DIR / "Tsunami_Output" / gf_file_freq_stk, skip_header=1
            ).reshape(-1, 2)
            np.testing.assert_almost_equal(
                content1_freq_strk[:, 0], content2_freq_strk[:, 0], decimal=6
            )
            np.testing.assert_almost_equal(
                content1_freq_strk[:, 1], content2_freq_strk[:, 1], decimal=6
            )

            content1_time_strk = np.genfromtxt(
                tempdir / "_GF" / gf_file_time_stk, skip_header=1
            ).reshape(-1, 2)
            content2_time_strk = np.genfromtxt(
                RESULTS_DIR / "Tsunami_Output" / gf_file_time_stk, skip_header=1
            ).reshape(-1, 2)
            np.testing.assert_almost_equal(
                content1_time_strk[:, 0], content2_time_strk[:, 0], decimal=6
            )
            np.testing.assert_almost_equal(
                content1_time_strk[:, 1], content2_time_strk[:, 1], decimal=6
            )
    finally:
        shutil.rmtree(tempdir)


def test_get_subfault_locations():
    segments_data = get_segments_data()
    tensor_info = get_tensor_info()
    subfault_lons_target = np.array(
        [
            -72.43570925,
            -72.10855538,
            -72.13523064,
            -71.80807677,
            -71.83475204,
            -71.86142731,
            -71.53427344,
            -71.56094871,
            -71.58762398,
            -71.26047011,
            -71.28714538,
            -71.31382065,
            -70.98666678,
            -71.01334204,
        ]
    )
    subfault_lats_target = np.array(
        [
            -32.79368777,
            -30.38975994,
            -31.68645655,
            -29.28252872,
            -30.57922533,
            -31.87592194,
            -29.47199411,
            -30.76869072,
            -32.06538733,
            -29.6614595,
            -30.95815611,
            -32.25485272,
            -29.85092489,
            -31.1476215,
        ]
    )
    subfault_depths_target = np.array(
        [
            2.688,
            2.688,
            7.616,
            7.616,
            12.544,
            17.472,
            17.472,
            22.4,
            27.328,
            27.328,
            32.256,
            37.184,
            37.184,
            42.112,
        ]
    )
    subfault_length_target = 17.92760869565217
    subfault_width_target = 14.92435321
    subfault_strike_target = 6.61391231
    subfault_dip_target = 19.28082797
    subfault_rake_target = 109.27817172
    (
        subfault_lons,
        subfault_lats,
        subfault_depths,
        subfault_length,
        subfault_width,
        subfault_strike,
        subfault_dip,
        subfault_rake,
    ) = get_subfault_locations(segments_data, tensor_info)
    np.testing.assert_allclose(subfault_lons[::15], subfault_lons_target, atol=1e-6)
    np.testing.assert_allclose(subfault_lats[::15], subfault_lats_target, atol=1e-6)
    np.testing.assert_allclose(subfault_depths[::15], subfault_depths_target, atol=1e-6)
    np.testing.assert_almost_equal(
        subfault_length[0], subfault_length_target, decimal=6
    )
    np.testing.assert_almost_equal(subfault_width[0], subfault_width_target, decimal=6)
    np.testing.assert_almost_equal(
        subfault_strike[0], subfault_strike_target, decimal=6
    )
    np.testing.assert_almost_equal(subfault_dip[0], subfault_dip_target, decimal=6)
    np.testing.assert_almost_equal(subfault_rake[0], subfault_rake_target, decimal=6)


def test_make_dart_json_no_observations():
    tempdir = tempfile.mkdtemp()
    try:
        tempdir = pathlib.Path(tempdir)
        shutil.copyfile(
            RESULTS_DIR / "NP1" / "tensor_info.json", tempdir / "tensor_info.json"
        )
        dart_station_file = HOME / "DART_station_locations.txt"
        make_dart_json_no_observations(tempdir, dart_station_file)
        with open(tempdir / "dart_waves.json", "r") as f:
            dartj = json.load(f)
        with open(RESULTS_DIR / "NP1" / "dart_waves_no_obs.json", "r") as g:
            target_dartj = json.load(g)
        sorted_dartj = sorted(dartj, key=itemgetter("name"))
        sorted_target_dartj = sorted(target_dartj, key=itemgetter("name"))
        for d1, d2 in zip(sorted_dartj, sorted_target_dartj):
            for key in [
                "azimuth",
                "distance",
                "dt",
                "duration",
                "location",
                "trace_weight",
            ]:
                np.testing.assert_almost_equal(d1[key], d2[key], decimal=7)
    finally:
        shutil.rmtree(tempdir)


def test_make_topo():
    tempdir = tempfile.mkdtemp()
    try:
        west, east, south, north = -75, -72, -33, -29
        topo_crop_out = make_topo(west, east, south, north, tempdir)
        assert (topo_crop_out).exists()
        topo_crop_target = RESULTS_DIR / "Tsunami_Output" / "topo.topo"
        with open(topo_crop_target) as f:
            topo_target = list(islice(f, 10))
        with open(topo_crop_out) as g:
            topo = list(islice(g, 10))
        assert topo == topo_target
    finally:
        shutil.rmtree(tempdir)


def test_make_ffm_kinematic_dtopo():
    tempdir = tempfile.mkdtemp()
    try:
        tempdir = pathlib.Path(tempdir)
        ffm_solution = RESULTS_DIR / "NP1" / "Solution.txt"
        dtopo_fname = make_ffm_kinematic_dtopo(ffm_solution, tempdir, makeplots=False)
        with open(dtopo_fname) as f:
            dtopo_10_lines = list(islice(f, 10))
        with open(RESULTS_DIR / "NP1" / "kinematic_dtopo.tt1") as g:
            dtopo_target_10_lines = list(islice(g, 10))
        assert dtopo_10_lines == dtopo_target_10_lines
    finally:
        shutil.rmtree(tempdir)


def test_make_ffm_static_dtopo():
    tempdir = tempfile.mkdtemp()
    try:
        tempdir = pathlib.Path(tempdir)
        ffm_solution = RESULTS_DIR / "NP1" / "Solution.txt"
        dtopo_fname = make_ffm_static_dtopo(ffm_solution, tempdir, makeplots=False)
        with open(dtopo_fname) as f:
            dtopo_10_lines = list(islice(f, 10))
        with open(RESULTS_DIR / "NP1" / "static_dtopo.tt1") as g:
            dtopo_target_10_lines = list(islice(g, 10))
        assert dtopo_10_lines == dtopo_target_10_lines
    finally:
        shutil.rmtree(tempdir)


@mock.patch(target="requests.get")
def test_query_dart_api(mock_request):
    tempdir = tempfile.mkdtemp()
    try:
        mock_request.return_value.status_code = 200
        mock_request.return_value.text = pathlib.Path(
            END_TO_END_DIR / "data" / "32402_dart_api.txt"
        ).read_text()
        tensor_info = get_tensor_info()
        origin_time = UTCDateTime(tensor_info["datetime"])
        bandpass_dart = query_dart_api(origin_time, "32402")
        print(bandpass_dart[0])
        bandpass_dart_target = np.genfromtxt(
            END_TO_END_DIR / "data" / "32402_dart_filt.txt"
        )
        print(bandpass_dart_target[:, 0])
        assert np.all(bandpass_dart[0] == bandpass_dart_target[:, 0])
        assert np.all(bandpass_dart[1] == bandpass_dart_target[:, 1])
    finally:
        shutil.rmtree(tempdir)


@mock.patch(target="requests.get")
def test_get_dart_observations(mock_request):
    tempdir = tempfile.mkdtemp()
    try:
        mock_request.return_value.status_code = 200
        mock_request.return_value.text = pathlib.Path(
            END_TO_END_DIR / "data" / "32402_dart_api.txt"
        ).read_text()
        tensor_info = get_tensor_info()
        get_dart_observations(tensor_info, tempdir, max_distance=10, plot=True)
        assert (pathlib.Path(tempdir) / "32402_dart.png").exists()
        assert (pathlib.Path(tempdir) / "processed_32402_dart.sac").exists()
        st = read(pathlib.Path(tempdir) / "processed_32402_dart.sac")
        st_target = read(END_TO_END_DIR / "data" / "processed_32402_dart.sac")
        assert np.all(st[0].data == st_target[0].data)
    finally:
        shutil.rmtree(tempdir)


def test_dart_high_res_locations():
    lon_dart_target = [-73.43278, -73.9775]
    lat_dart_target = [-20.46861, -26.73472]
    gauges_target = {
        "0": {
            "gauge_id": "32401",
            "lon": -73.43278,
            "lat": -20.46861,
            "t1": 0,
            "t2": 10000000000.0,
        },
        "1": {
            "gauge_id": "32402",
            "lon": -73.9775,
            "lat": -26.73472,
            "t1": 0,
            "t2": 10000000000.0,
        },
    }
    high_res_regions_target = {
        "0": {
            "minlevel": 3,
            "maxlevel": 3,
            "t1": 0.0,
            "t2": 10000.0,
            "x1": -73.93278,
            "x2": -72.93278,
            "y1": -20.96861,
            "y2": -19.96861,
        },
        "1": {
            "minlevel": 3,
            "maxlevel": 3,
            "t1": 0.0,
            "t2": 10000.0,
            "x1": -74.4775,
            "x2": -73.4775,
            "y1": -27.23472,
            "y2": -26.23472,
        },
        "2": {
            "minlevel": 3,
            "maxlevel": 3,
            "t1": 0.0,
            "t2": 300.0,
            "x1": -72.63570924588406,
            "x2": -70.57342920594063,
            "y1": -33.110501912886875,
            "y2": -29.067926954701772,
        },
    }

    dart_data = get_dart_waves_json()
    segments_data = get_segments_data()
    tensor_info = get_tensor_info()
    high_res_regions, gauges, lon_dart, lat_dart = dart_high_res_locations(
        dart_data, segments_data, tensor_info
    )
    assert high_res_regions == high_res_regions_target
    assert gauges == gauges_target
    assert lon_dart == lon_dart_target
    assert lat_dart == lat_dart_target


def test_plot_dtopo_snapshots():
    tempdir = tempfile.mkdtemp()
    try:
        tempdir = pathlib.Path(tempdir)
        dtopo_file = RESULTS_DIR / "NP1" / "kinematic_dtopo.tt1"
        hyp_lon, hyp_lat = -71.674, -31.573
        output_dir = tempdir / "dtopo_pngs"
        plot_dtopo_snapshots(dtopo_file, hyp_lon, hyp_lat, output_dir)
        assert (output_dir / "dtopo_t50.00.png").exists()
        assert (output_dir / "dtopo_animation.gif").exists()
    finally:
        shutil.rmtree(tempdir)


def test_run_ffm_tsunami():
    tempdir = tempfile.mkdtemp()

    model_region_target = {
        "west": -98.9775,
        "east": -45.77342920594063,
        "south": -57.91050191288687,
        "north": 4.531389999999998,
    }
    high_res_regions_target = {
        "0": {
            "minlevel": 3,
            "maxlevel": 3,
            "t1": 0.0,
            "t2": 10000.0,
            "x1": -73.93278,
            "x2": -72.93278,
            "y1": -20.96861,
            "y2": -19.96861,
        },
        "1": {
            "minlevel": 3,
            "maxlevel": 3,
            "t1": 0.0,
            "t2": 10000.0,
            "x1": -74.4775,
            "x2": -73.4775,
            "y1": -27.23472,
            "y2": -26.23472,
        },
        "2": {
            "minlevel": 3,
            "maxlevel": 3,
            "t1": 0.0,
            "t2": 300.0,
            "x1": -72.63570924588406,
            "x2": -70.57342920594063,
            "y1": -33.110501912886875,
            "y2": -29.067926954701772,
        },
    }
    gauges_target = {
        "0": {
            "gauge_id": "32401",
            "lon": -73.43278,
            "lat": -20.46861,
            "t1": 0,
            "t2": 10000000000.0,
        },
        "1": {
            "gauge_id": "32402",
            "lon": -73.9775,
            "lat": -26.73472,
            "t1": 0,
            "t2": 10000000000.0,
        },
    }
    topo_path_target = "topo.topo"
    dtopo_path_target = "kinematic_dtopo.tt1"
    dtopo_type_target = 1
    dt_max_dtopo_target = 5
    tfinal_target = 10800
    num_output_times_target = 36
    wave_tolerance_target = 0.0001
    try:
        tempdir = pathlib.Path(tempdir)
        ffm_solution = RESULTS_DIR / "NP1" / "Solution.txt"
        shutil.copyfile(
            RESULTS_DIR / "NP1" / "tensor_info.json", tempdir / "tensor_info.json"
        )
        shutil.copyfile(
            RESULTS_DIR / "NP1" / "segments_data.json", tempdir / "segments_data.json"
        )
        shutil.copyfile(
            RESULTS_DIR / "NP1" / "dart_waves_no_obs.json", tempdir / "dart_waves.json"
        )
        run_ffm_tsunami(tempdir, ffm_solution, run_geoclaw=False)
        geoclaw_parameter_file = tempdir / "geoclaw_parameter_file.json"
        with open(geoclaw_parameter_file) as g:
            geoclaw_parameter = json.load(g)
        assert geoclaw_parameter["model_region"] == model_region_target
        assert geoclaw_parameter["high_res_regions"] == high_res_regions_target
        assert geoclaw_parameter["gauges"] == gauges_target
        assert os.path.basename(geoclaw_parameter["topo_path"]) == topo_path_target
        assert os.path.basename(geoclaw_parameter["dtopo_path"]) == dtopo_path_target
        assert geoclaw_parameter["dtopo_type"] == dtopo_type_target
        assert geoclaw_parameter["dt_max_dtopo"] == dt_max_dtopo_target
        assert geoclaw_parameter["tfinal"] == tfinal_target
        assert geoclaw_parameter["num_output_times"] == num_output_times_target
        assert geoclaw_parameter["wave_tolerance"] == wave_tolerance_target
    finally:
        shutil.rmtree(tempdir)


def test_dart_gf_time_freq_domain():
    tempdir = tempfile.mkdtemp()
    try:
        tempdir = pathlib.Path(tempdir)
        id_dart = 32401
        subfault_start = 206
        dt = 60
        gf_file_time_dip = (
            f"dart.{id_dart}.{str(subfault_start).rjust(4, '0')}.GF.dip.txt"
        )
        timegf = np.genfromtxt(RESULTS_DIR / "Tsunami_Output" / gf_file_time_dip)
        time_out = timegf[:, 0]
        gf_fil = timegf[:, 1]
        gf_file_freq_dip = f"{id_dart}_dip_gf.txt"
        freq_dip_target = RESULTS_DIR / "Tsunami_Output" / gf_file_freq_dip
        dart_gf_time_freq_domain(
            id_dart,
            int(dt),
            subfault_start,
            time_out,
            gf_fil,
            tempdir / gf_file_freq_dip,
        )
        with (
            open(tempdir / gf_file_freq_dip, "r") as f1,
            open(freq_dip_target, "r") as f2,
        ):
            content1 = f1.read()
            content2 = f2.read()
        assert content1 == content2
    finally:
        shutil.rmtree(tempdir)
