import json
import pathlib
import shutil
import tempfile
from unittest import mock

from wasp.modify_jsons import modify_channels

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


def test_modify_channels_no_input():
    tempdir = tempfile.mkdtemp()
    channels = pathlib.Path(tempdir) / "channels.json"
    try:
        with open(channels, "w") as f:
            json.dump(CHANNELS, f)

        # test downweight
        modify_channels(
            json_file=channels,
            method="downweight",
            input=False,
            station_dict={"ABCD": ["SH"], "EFGH": ["BH1", "BHZ"]},
        )
        with open(channels, "r") as f:
            downweighted = json.load(f)
        assert downweighted == [
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

        # test delete
        modify_channels(
            json_file=channels,
            method="delete",
            input=False,
            station_dict={"ABCD": ["SH"], "EFGH": ["BH2", "BHZ"]},
        )
        with open(channels, "r") as f:
            deleted = json.load(f)
        assert deleted == [
            {
                "component": "BH1",
                "name": "EFGH",
                "trace_weight": 0,
            }
        ]
    finally:
        print("Cleaning up test directory.")
        shutil.rmtree(tempdir)


@mock.patch(
    "builtins.input",
    side_effect=["ABCD", "SH", "exit", "EFGH", "BH1", "EFGH", "BHZ", "exit"],
)
def test_modify_channels_input(mock_input):
    tempdir = tempfile.mkdtemp()
    channels = pathlib.Path(tempdir) / "channels.json"
    try:
        with open(channels, "w") as f:
            json.dump(CHANNELS, f)

        # test downweight
        modify_channels(
            json_file=channels,
            method="downweight",
        )
        with open(channels, "r") as f:
            downweighted = json.load(f)
        assert downweighted == [
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
                "trace_weight": 0,
            },
        ]

        # test delete
        modify_channels(
            json_file=channels,
            method="delete",
        )
        with open(channels, "r") as f:
            deleted = json.load(f)
        assert deleted == [
            {
                "component": "BH2",
                "name": "EFGH",
                "trace_weight": 1.0,
            },
            {
                "component": "SH",
                "name": "ABCD",
                "trace_weight": 0,
            },
        ]
    finally:
        print("Cleaning up test directory.")
        shutil.rmtree(tempdir)
