import importlib
import json
import pathlib
from datetime import datetime, timezone
from shutil import rmtree
from tempfile import mkdtemp
from test.testutils import DATA_DIR, MockResponse
from unittest import mock

from ffm.acquisition.cmt import Cmt

DETAIL = {
    "properties": {
        "place": "here, there",
        "products": {
            "origin": [
                {
                    "source": "not",
                    "properties": {
                        "magnitude": 5,
                        "eventtime": 1771036067951,
                        "latitude": -15,
                        "longitude": 150,
                        "depth": 6,
                    },
                },
                {
                    "source": "us",
                    "properties": {
                        "magnitude": 6.4,
                        "eventtime": "2026-02-14T02:27:37.951Z",
                        "latitude": -14.8934,
                        "longitude": 166.6013,
                        "depth": 10,
                        "review-status": "reviewed",
                    },
                },
                {
                    "source": "us",
                    "properties": {
                        "magnitude": 5,
                        "eventtime": 1771036057951,
                        "latitude": -14.8934,
                        "longitude": 170.6013,
                        "depth": 4,
                        "review-status": "not-reviewed",
                    },
                },
                {
                    "source": "us2",
                    "properties": {
                        "magnitude": 2,
                        "eventtime": 1771036057951,
                        "latitude": -12.8934,
                        "longitude": 171.6013,
                        "depth": 12,
                        "review-status": "not-reviewed",
                    },
                },
            ],
            "moment-tensor": [
                {
                    "source": "not",
                    "properties": {
                        "tensor-mpp": "-3E+17",
                        "tensor-mrp": "-3E+18",
                        "tensor-mrr": "1E+18",
                        "tensor-mrt": "1E+18",
                        "tensor-mtp": "1E+18",
                        "tensor-mtt": "-1E+18",
                        "sourcetime-duration": "8",
                    },
                },
                {
                    "source": "us",
                    "properties": {
                        "tensor-mpp": "-2E+17",
                        "tensor-mrp": "-2E+18",
                        "tensor-mrr": "3E+18",
                        "tensor-mrt": "3E+18",
                        "tensor-mtp": "3E+18",
                        "tensor-mtt": "-2E+18",
                        "sourcetime-duration": "7",
                        "review-status": "reviewed",
                        "derived-magnitude-type": "Mww",
                        "derived-depth": 10.1,
                        "derived-latitude": -14.89341,
                        "derived-longitude": 166.60131,
                    },
                },
                {
                    "source": "us",
                    "properties": {
                        "derived-magnitude-type": "Mww",
                        "sourcetime-duration": "2",
                        "tensor-mpp": "-1E+17",
                        "tensor-mrp": "-2E+18",
                        "tensor-mrr": "3E+18",
                        "tensor-mrt": "5E+18",
                        "tensor-mtp": "3E+18",
                        "tensor-mtt": "-1E+18",
                    },
                },
                {
                    "source": "us2",
                    "properties": {
                        "tensor-mpp": "-4E+17",
                        "tensor-mrp": "-4E+18",
                        "tensor-mrr": "2E+18",
                        "tensor-mrt": "2E+18",
                        "tensor-mtp": "2E+18",
                        "sourcetime-duration": "5",
                    },
                },
            ],
        },
    },
}


def test_cmt_from_id():
    with mock.patch(target="requests.get") as mock_requests:
        mock_requests.return_value = MockResponse(json.dumps(DETAIL), 200)
        detail_cmt = Cmt.from_id("eventid", "us")
        assert detail_cmt.depth == 10
        assert detail_cmt.eventid == "eventid"
        assert detail_cmt.latitude == -14.8934
        assert detail_cmt.longitude == 166.6013
        assert detail_cmt.magnitude == 6.4
        assert detail_cmt.magnitude_type == "Mww"
        assert detail_cmt.mpp == -2e17
        assert detail_cmt.mrp == -2e18
        assert detail_cmt.mrr == 3e18
        assert detail_cmt.mrt == 3e18
        assert detail_cmt.mtp == 3e18
        assert detail_cmt.mtt == -2e18
        assert detail_cmt.duration == 7
        assert detail_cmt.place == "there"
        assert detail_cmt.time == datetime(
            2026, 2, 14, 2, 27, 37, 951000, tzinfo=timezone.utc
        )


def test_cmt_from_detail():
    detail_cmt = Cmt.from_detail(DETAIL, "eventid", "us")
    assert detail_cmt.depth == 10
    assert detail_cmt.eventid == "eventid"
    assert detail_cmt.latitude == -14.8934
    assert detail_cmt.longitude == 166.6013
    assert detail_cmt.magnitude == 6.4
    assert detail_cmt.magnitude_type == "Mww"
    assert detail_cmt.mpp == -2e17
    assert detail_cmt.mrp == -2e18
    assert detail_cmt.mrr == 3e18
    assert detail_cmt.mrt == 3e18
    assert detail_cmt.mtp == 3e18
    assert detail_cmt.mtt == -2e18
    assert detail_cmt.duration == 7
    assert detail_cmt.place == "there"
    assert detail_cmt.time == datetime(
        2026, 2, 14, 2, 27, 37, 951000, tzinfo=timezone.utc
    )


def test_cmt_from_moment_tensor():
    detail_cmt = Cmt.from_moment_tensor(
        depth=10,
        eventid="eventid",
        latitude=-14.8934,
        longitude=166.6013,
        magnitude=6.4,
        moment_tensor=DETAIL["properties"]["products"]["moment-tensor"][1],
        place="there",
        time=datetime(2026, 2, 14, 2, 27, 37, 951000, tzinfo=timezone.utc),
    )
    assert detail_cmt.depth == 10
    assert detail_cmt.eventid == "eventid"
    assert detail_cmt.latitude == -14.8934
    assert detail_cmt.longitude == 166.6013
    assert detail_cmt.magnitude == 6.4
    assert detail_cmt.magnitude_type == "Mww"
    assert detail_cmt.mpp == -2e17
    assert detail_cmt.mrp == -2e18
    assert detail_cmt.mrr == 3e18
    assert detail_cmt.mrt == 3e18
    assert detail_cmt.mtp == 3e18
    assert detail_cmt.mtt == -2e18
    assert detail_cmt.duration == 7
    assert detail_cmt.place == "there"
    assert detail_cmt.time == datetime(
        2026, 2, 14, 2, 27, 37, 951000, tzinfo=timezone.utc
    )


def test_cmt_write():
    with mock.patch(target="requests.get") as mock_requests:
        mock_requests.return_value = MockResponse(json.dumps(DETAIL), 200)
        detail_cmt = Cmt.from_id("eventid", "us")

    tempdir = pathlib.Path(mkdtemp())
    try:
        detail_cmt.write(tempdir / "cmt")
        with open(tempdir / "cmt") as f:
            data = f.read()
        assert (
            data
            == """ US 2026 02 14 02 27 37.951000 -14.8934 166.6013 10.0  0.0 6.4 there
event name: eventid
time shift: 3.5
half duration: 3.5
latitude: -14.89341
longitude: 166.60131
depth: 10.1
Mrr: 3e+25
Mtt: -2e+25
Mpp: -2e+24
Mrt: 3e+25
Mrp: -2e+25
Mtp: 3e+25
"""
        )
    finally:
        rmtree(tempdir)
