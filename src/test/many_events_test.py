import json
import pathlib
import shutil
import tempfile
from test.testutils import RESULTS_DIR

from wasp.many_events import (
    get_model_space_events,
    get_moment_events,
    get_segments_events,
    get_waveforms_events,
    get_waveforms_events2,
    save_dict,
    select_segments_event,
    select_waveforms_event,
)


def test_get_model_space_events():
    tempdir = pathlib.Path(tempfile.mkdtemp())
    try:
        get_model_space_events(
            [
                [
                    {
                        "regularization": {
                            "neighbour_down": {"segment": 1},
                            "neighbour_left": {"segment": 2},
                            "neighbour_right": {"segment": 3},
                            "neighbour_up": {"segment": 4},
                        },
                    }
                ]
            ],
            directory=tempdir,
        )
        with open(tempdir / "model_space.json") as f:
            d = json.load(f)
        assert d == [
            {
                "regularization": {
                    "neighbour_down": {"segment": 1},
                    "neighbour_left": {"segment": 2},
                    "neighbour_right": {"segment": 3},
                    "neighbour_up": {"segment": 4},
                }
            }
        ]
    finally:
        shutil.rmtree(tempdir)


def test_get_moment_events():
    tempdir = pathlib.Path(tempfile.mkdtemp())
    try:
        get_moment_events(
            [{"seismic_moment": "moment1"}, {"seismic_moment": "moment2"}],
            directory=tempdir,
        )
        with open(tempdir / "moment_events.txt") as f:
            d = f.read()
        assert d == "moment1\nmoment2\n"
    finally:
        shutil.rmtree(tempdir)


def test_get_segments_events():
    tempdir = pathlib.Path(tempfile.mkdtemp())
    try:
        get_segments_events(
            [
                {
                    "rise_time": {
                        "name": "rise_time1",
                    },
                    "connections": [
                        {"segment1": 2, "segment2": 1},
                        {"segment1": 3, "segment2": 4},
                    ],
                    "segments": [{"name": "segment1"}],
                },
                {
                    "rise_time": {
                        "name": "rise_time2",
                    },
                    "segments": [{"name": "segment2"}],
                },
            ],
            [{"lat": 1, "lon": 2, "depth": 3}, {"lat": 4, "lon": 5, "depth": 6}],
            directory=tempdir,
        )
        with open(tempdir / "segments_events.txt") as f:
            txt = f.read()
        assert txt == "segment event\n1 1\n2 2\n"
        with open(tempdir / "segments_data.json") as f:
            dict = json.load(f)

        assert dict == {
            "connections": [
                {"segment1": 2, "segment2": 1},
                {"segment1": 3, "segment2": 4},
            ],
            "rise_time": {"name": "rise_time1"},
            "segments": [
                {
                    "event": 1,
                    "hypocenter": {"depth": 3, "lat": 1, "lon": 2},
                    "name": "segment1",
                },
                {
                    "event": 2,
                    "hypocenter": {"depth": 6, "lat": 4, "lon": 5},
                    "name": "segment2",
                },
            ],
        }
    finally:
        shutil.rmtree(tempdir)


def test_get_waveforms_events():
    tempdir = pathlib.Path(tempfile.mkdtemp())
    try:
        for datatype, txt_file, dict_file in zip(
            ["cgnss", "strong", "body", "surf", "gnss"],
            [
                "cgnss_events.txt",
                "strong_motion_events.txt",
                "tele_events.txt",
                "surf_events.txt",
                "static_events.txt",
            ],
            [
                "cgnss_waves.json",
                "strong_motion_waves.json",
                "tele_waves.json",
                "surf_waves.json",
                "static_data.json",
            ],
        ):
            get_waveforms_events(
                [
                    [
                        {
                            "component": "component1",
                            "name": "waveform1",
                        }
                    ],
                    [
                        {
                            "component": "component2",
                            "name": "waveform2",
                        }
                    ],
                ],
                [datatype],
                directory=tempdir,
            )
            with open(tempdir / txt_file) as f:
                txt = f.read()
            assert txt == "waveform1 component1 1\nwaveform2 component2 2\n"
            with open(tempdir / dict_file) as f:
                dict = json.load(f)
            assert dict == [
                {"event": 1, "component": "component1", "name": "waveform1"},
                {"event": 2, "component": "component2", "name": "waveform2"},
            ]
    finally:
        shutil.rmtree(tempdir)


def test_get_waveforms_events2():
    tempdir = pathlib.Path(tempfile.mkdtemp())
    try:
        dict_file = tempdir / "test_file.txt"
        get_waveforms_events2(
            [
                [{"name": "waveform1", "component": "component1"}],
                [{"name": "waveform2", "component": "component2"}],
            ],
            dict_file,
        )
        with open(dict_file) as f:
            d = f.read()
        assert d == "waveform1 component1 1\nwaveform2 component2 2\n"
    finally:
        shutil.rmtree(tempdir)


def test_save_dict():
    tempdir = pathlib.Path(tempfile.mkdtemp())
    try:
        dict_file = tempdir / "test_dict.json"
        save_dict([[{"name": "event1"}], [{"name": "event2"}]], dict_file)
        with open(dict_file) as f:
            d = json.load(f)
        assert d == [{"event": 1, "name": "event1"}, {"event": 2, "name": "event2"}]
    finally:
        shutil.rmtree(tempdir)


def test_select_segments_event():
    tempdir = pathlib.Path(tempfile.mkdtemp())
    try:
        d = select_segments_event(
            {
                "rise_time": "rise_time1",
                "segments": [
                    {"event": 1, "segment": "segment1"},
                    {"event": 1, "segment": "segment2"},
                    {"event": 2, "segment": "segment3"},
                    {"event": 2, "segment": "segment4"},
                ],
            },
            2,
        )
        assert d == {
            "rise_time": "rise_time1",
            "segments": [
                {"event": 2, "segment": "segment3"},
                {"event": 2, "segment": "segment4"},
            ],
        }
    finally:
        shutil.rmtree(tempdir)


def test_select_waveforms_event():
    tempdir = pathlib.Path(tempfile.mkdtemp())
    try:
        d = select_waveforms_event(
            [
                {"name": "name1", "event": 1},
                {"name": "name2", "event": 2},
                {"name": "name3", "event": 2},
                {"name": "name4", "event": 3},
            ],
            2,
        )
        assert d == [{"name": "name2", "event": 2}, {"name": "name3", "event": 2}]
    finally:
        shutil.rmtree(tempdir)
