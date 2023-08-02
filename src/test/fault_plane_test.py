import json
import pathlib
import shutil
import tempfile

import numpy as np
from obspy import UTCDateTime

import wasp.seismic_tensor as tensor
from wasp.fault_plane import (
    __default_vel_of_eq,
    __epicenter_location,
    __fault_plane_properties,
    __hypocenter_location2,
    __lat_lon,
    __plane_tensor_def,
    __point_sources_general,
    __rise_time_parameters,
    __save_plane_data,
    __source_layer,
    __subfaults_properties,
    __write_event_mult_in,
    _point_sources_def,
    create_finite_fault,
    event_mult_in_to_json,
    is_fault_correct,
    point_sources_param,
    shear_modulous,
)

from .testutils import (
    RESULTS_DIR,
    get_segments_data,
    get_tensor_info,
    get_velmodel_data,
)

TENSOR = get_tensor_info()
SEGMENTS = get_segments_data()


def test_create_finite_fault():
    info_np1, info_np2 = tensor.planes_from_tensor(TENSOR)
    data = create_finite_fault(
        TENSOR,
        info_np1["plane_info"],
        ["cgps", "gps", "insar", "strong_motion", "tele_body", "surf_waves"],
    )
    assert data == SEGMENTS


def test_epicenter_location():
    assert __epicenter_location(100, 50) == {
        "hyp_stk": 100,
        "hyp_dip": 50,
    }


def test_event_mult_in_to_json():
    tempdir = tempfile.mkdtemp()
    try:
        event_file = pathlib.Path(tempdir) / "Event_mult.in"
        with open(RESULTS_DIR / "NP1" / "Event_mult.in") as i:
            test_mult_in = i.read()
        with open(event_file, "w") as f:
            f.write(test_mult_in)

        event_mult_in_to_json(tempdir)
        with open(pathlib.Path(tempdir) / "segments_data.json", "r") as f:
            fault_data = json.load(f)
        assert fault_data == SEGMENTS
    finally:
        shutil.rmtree(tempdir)


def test__fault_plane_properties():
    info_np1, info_np2 = tensor.planes_from_tensor(TENSOR)
    pinfo = info_np1["plane_info"]
    plane_info = __plane_tensor_def(pinfo["strike"], pinfo["dip"], pinfo["rake"], 2.5)
    fault_dimensions = __fault_plane_properties(137.445, TENSOR, plane_info, 0)
    segments_target = {
        i: SEGMENTS["segments"][0][i]
        for i in SEGMENTS["segments"][0]
        if i
        not in [
            "hyp_dip",
            "hyp_stk",
            "neighbours",
            "rake",
            "rupture_vel",
            "strike",
            "dip",
        ]
    }
    assert fault_dimensions == segments_target


def test_hypocenter_location2():
    info_np1, info_np2 = tensor.planes_from_tensor(TENSOR)
    pinfo = info_np1["plane_info"]
    plane_info = __plane_tensor_def(pinfo["strike"], pinfo["dip"], pinfo["rake"], 2.5)
    assert __hypocenter_location2(
        plane_info,
        SEGMENTS["segments"][0],
        TENSOR,
        0,
        SEGMENTS["rise_time"],
    ) == {"hyp_dip": 5, "hyp_stk": 9}


def test_is_fault_correct():
    assert (
        is_fault_correct(TENSOR, {"dip": 40, "delta_dip": 10, "hyp_dip": 20}) == False
    )
    assert is_fault_correct(TENSOR, {"dip": 40, "delta_dip": 1, "hyp_dip": 2}) == True


def test__lat_lon():
    # pull an example from running point_sources_param
    assert __lat_lon(
        6.613912311529926,
        19.280827965117993,
        258.15756521739127,
        65.66715411624098,
        -31.57,
        -71.67,
    ) == (-29.326477058851623, -70.70558238456391)


def test_plane_tensor_def():
    assert __plane_tensor_def(20, 15, 30, 2) == {
        "strike": 20,
        "dip": 15,
        "rake": 30,
        "rupture_vel": 2,
    }


def test_point_sources_def():
    assert _point_sources_def(SEGMENTS["rise_time"], 2.5, SEGMENTS["segments"][0],) == {
        "dip_ps": 5,
        "dx": 3.585521739130434,
        "dy": 2.9848706416473174,
        "strike_ps": 5,
    }


def test_point_sources_general():
    assert __point_sources_general(10, 30, 1, 2) == {
        "strike_ps": 10,
        "dip_ps": 30,
        "dx": 1,
        "dy": 2,
    }

    def test_rise_time_parameters():
        assert (
            __rise_time_parameters(
                TENSOR,
                SEGMENTS["segments"][0],
                "tele_body",
            )
            == SEGMENTS["rise_time"]
        )

    assert __rise_time_parameters(
        {
            "date_origin": UTCDateTime("2022-5-19T10:00:00Z"),
            "moment_mag": 2.4021744462980653e26,
            "lat": -54.1325,
            "lon": 159.0268,
            "depth": 10,
            "time_shift": 12,
        },
        {
            "delay_segment": 0.0,
            "delta_dip": 2.6999999999999997,
            "delta_strike": 3.1271428571428572,
            "dip": 89.92771983437935,
            "dip_subfaults": 11,
            "hyp_dip": 4,
            "hyp_stk": 11,
            "neighbours": [],
            "rake": -173.31774652086895,
            "rupture_vel": 2.5,
            "stk_subfaults": 21,
            "strike": 197.57412152089023,
        },
        "regular",
    ) == {"delta_rise": 2.0, "min_rise": 2.0, "windows": 3}
    assert __rise_time_parameters(
        {
            "date_origin": UTCDateTime("2022-5-19T10:00:00Z"),
            "moment_mag": 2.4021744462980653e26,
            "lat": -54.1325,
            "lon": 159.0268,
            "depth": 10,
            "time_shift": 30,
        },
        {
            "delay_segment": 0.0,
            "delta_dip": 2.6999999999999997,
            "delta_strike": 3.1271428571428572,
            "dip": 89.92771983437935,
            "dip_subfaults": 11,
            "hyp_dip": 4,
            "hyp_stk": 11,
            "neighbours": [],
            "rake": -173.31774652086895,
            "rupture_vel": 2.5,
            "stk_subfaults": 21,
            "strike": 197.57412152089023,
        },
        "regular",
    ) == {"delta_rise": 2.5, "min_rise": 2.5, "windows": 2}
    assert __rise_time_parameters(
        {
            "date_origin": UTCDateTime("2022-5-19T10:00:00Z"),
            "moment_mag": 2.4021744462980653e26,
            "lat": -54.1325,
            "lon": 159.0268,
            "depth": 10,
            "time_shift": 50,
        },
        {
            "delay_segment": 0.0,
            "delta_dip": 2.6999999999999997,
            "delta_strike": 3.1271428571428572,
            "dip": 89.92771983437935,
            "dip_subfaults": 11,
            "hyp_dip": 4,
            "hyp_stk": 11,
            "neighbours": [],
            "rake": -173.31774652086895,
            "rupture_vel": 2.5,
            "stk_subfaults": 21,
            "strike": 197.57412152089023,
        },
        "regular",
    ) == {"delta_rise": 4.0, "min_rise": 4.0, "windows": 2}
    assert __rise_time_parameters(
        {
            "date_origin": UTCDateTime("2022-5-19T10:00:00Z"),
            "moment_mag": 2.4021744462980653e26,
            "lat": -54.1325,
            "lon": 159.0268,
            "depth": 100000,
            "time_shift": 5,
        },
        {
            "delay_segment": 0.0,
            "delta_dip": 2.6999999999999997,
            "delta_strike": 3.1271428571428572,
            "dip": 89.92771983437935,
            "dip_subfaults": 11,
            "hyp_dip": 4,
            "hyp_stk": 11,
            "neighbours": [],
            "rake": -173.31774652086895,
            "rupture_vel": 2.5,
            "stk_subfaults": 21,
            "strike": 197.57412152089023,
        },
        ["tele_body"],
    ) == {"delta_rise": 1.0, "min_rise": 1.0, "windows": 4}


def test_point_sources_param():
    point_sources = point_sources_param(
        SEGMENTS["segments"],
        TENSOR,
        SEGMENTS["rise_time"],
    )
    assert np.max([np.max(p.flatten()) for p in point_sources]) == 359.97359539607316


def test_save_plane_data():
    tempdir = tempfile.mkdtemp()
    try:
        test_dict = __save_plane_data(
            {"plane_test": "test_plane_value"},
            {"subfault_test": "test_subfault_value"},
            {"epicenter_test": "test_epicenter_value"},
            {"rise_test": "test_rise_value"},
            tempdir,
        )
        with open(pathlib.Path(tempdir) / "segments_data.json", "r") as f:
            sdata = json.load(f)
        assert sdata == {
            "rise_time": {"rise_test": "test_rise_value"},
            "segments": [
                {
                    "epicenter_test": "test_epicenter_value",
                    "neighbours": [],
                    "plane_test": "test_plane_value",
                    "subfault_test": "test_subfault_value",
                }
            ],
        }
    finally:
        shutil.rmtree(tempdir)


def test_shear_modulous():
    point_sources = point_sources_param(
        SEGMENTS["segments"],
        TENSOR,
        SEGMENTS["rise_time"],
    )
    velmodel = get_velmodel_data()
    shear = shear_modulous(point_sources, velmodel=velmodel)
    for idx, target in enumerate(
        [
            352833668787.89996,
            352833668787.89996,
            424030485870.0,
            424030485870.0,
            424030485870.0,
            564461557600.0,
            564461557600.0,
            564461557600.0,
            675340884666.0,
        ]
    ):
        assert shear[0][idx][0] == target


def test_source_layer():
    thicknesses = [10, 1, 50]
    source_depth = 12
    assert __source_layer(thicknesses, source_depth) == 2
    thicknesses = [10, 1, 50]
    source_depth = 10.1
    assert __source_layer(thicknesses, source_depth) == 1


def test_subfaults_properties():
    assert __subfaults_properties(2, 3, 4, 5) == {
        "delay_segment": 0,
        "delta_strike": 2,
        "delta_dip": 3,
        "stk_subfaults": 4,
        "dip_subfaults": 5,
    }


def test_write_event_mult_in():
    tempdir = tempfile.mkdtemp()
    try:
        info_np1, info_np2 = tensor.planes_from_tensor(TENSOR)
        pinfo = info_np1["plane_info"]
        plane_info = __plane_tensor_def(
            pinfo["strike"], pinfo["dip"], pinfo["rake"], 2.5
        )

        __write_event_mult_in(
            TENSOR,
            plane_info,
            SEGMENTS["segments"][0],
            {"hyp_dip": 5, "hyp_stk": 9},
            SEGMENTS["rise_time"],
            tempdir,
        )
        with open(pathlib.Path(tempdir) / "Event_mult.in", "r") as f:
            event_mult_in = f.readlines()
        with open(RESULTS_DIR / "NP1" / "Event_mult.in") as t:
            test_mult_in = t.read()
        target = test_mult_in.split("\n")
        for i in range(len(target)):
            np.testing.assert_allclose(
                np.array(event_mult_in[i].split(), dtype=np.float64),
                np.array(target[i].split(), dtype=np.float64),
            )
    finally:
        shutil.rmtree(tempdir)
