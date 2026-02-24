import json
import pathlib
from datetime import datetime, timezone
from shutil import rmtree
from tempfile import mkdtemp
from test.testutils import MockResponse
from unittest import mock

import numpy as np
from obspy.core import Stream, Trace

from ffm.acquisition.irisquery import IrisQuery

DETAIL = {
    "properties": {
        "products": {
            "origin": [
                {
                    "source": "us",
                    "properties": {
                        "depth": "10",
                        "eventtime": "2026-02-14T02:27:37.951Z",
                        "latitude": "-14.8934",
                        "longitude": "166.6013",
                    },
                }
            ]
        },
    }
}


def test_irisquery_from_id():
    with mock.patch(target="requests.get") as mock_requests:
        mock_requests.return_value = MockResponse(json.dumps(DETAIL), 200)
        detail_irisquery = IrisQuery.from_id("eventid", "us")
        assert detail_irisquery.depth == 10
        assert detail_irisquery.latitude == -14.8934
        assert detail_irisquery.longitude == 166.6013
        assert detail_irisquery.duration == 3000
        assert detail_irisquery.event_time == datetime(
            2026, 2, 14, 2, 27, 37, 951000, tzinfo=timezone.utc
        )


def test_irisquery_from_detail():
    detail_irisquery = IrisQuery.from_detail(DETAIL, "us")
    assert detail_irisquery.depth == 10
    assert detail_irisquery.latitude == -14.8934
    assert detail_irisquery.longitude == 166.6013
    assert detail_irisquery.duration == 3000
    assert detail_irisquery.event_time == datetime(
        2026, 2, 14, 2, 27, 37, 951000, tzinfo=timezone.utc
    )


def test_irisquery_from_moment_tensor():
    detail_irisquery = IrisQuery.from_origin(
        origin=DETAIL["properties"]["products"]["origin"][0]
    )
    assert detail_irisquery.depth == 10
    assert detail_irisquery.latitude == -14.8934
    assert detail_irisquery.longitude == 166.6013
    assert detail_irisquery.duration == 3000
    assert detail_irisquery.event_time == datetime(
        2026, 2, 14, 2, 27, 37, 951000, tzinfo=timezone.utc
    )


def test_irisquery_end_to_end():
    detail_irisquery = IrisQuery.from_detail(DETAIL, "us")
    assert detail_irisquery.depth == 10
    assert detail_irisquery.latitude == -14.8934
    assert detail_irisquery.longitude == 166.6013
    assert detail_irisquery.duration == 3000
    assert detail_irisquery.event_time == datetime(
        2026, 2, 14, 2, 27, 37, 951000, tzinfo=timezone.utc
    )
    with mock.patch(target="ffm.acquisition.irisquery.Client") as mock_client:
        mocked_client = mock.MagicMock()
        mock_client.return_value = mocked_client

        # test getting data
        mocked_client.get_waveforms.return_value = Stream(
            [
                Trace(
                    np.array([1, 2, 3]),
                    header={
                        "station": "STAT1",
                        "channel": "BHZ",
                        "location": "00",
                        "network": "US",
                    },
                ),
                Trace(
                    np.array([3, 4, 5]),
                    header={
                        "station": "STAT1",
                        "channel": "BHZ",
                        "location": "10",
                        "network": "US",
                    },
                ),
            ]
        )

        data = detail_irisquery.get_data(
            networks=["US"], stations={"US": {"STAT1": ["BH*"]}}
        )
        assert len(data) == 2

    with mock.patch(target="ffm.acquisition.irisquery.Iris") as mock_iris:
        mocked_iris = mock.MagicMock()
        mock_iris.return_value = mocked_iris

        # test getting response
        mocked_iris.sacpz.side_effect = [b"response 1", b"response 2"]
        data_with_response = detail_irisquery.get_responses(data)
        assert data_with_response[0].stats.response == b"response 1"
        assert data_with_response[1].stats.response == b"response 2"

    # test writing the data
    tempdir = pathlib.Path(mkdtemp())
    try:
        detail_irisquery.write_data(tempdir, data_with_response)
        assert (tempdir / "US_STAT1_BHZ_00.sac").exists()
        assert (tempdir / "US_STAT1_BHZ_10.sac").exists()
        assert (tempdir / "SAC_PZs_US_STAT1_BHZ_00.sac").exists()
        assert (tempdir / "SAC_PZs_US_STAT1_BHZ_10.sac").exists()
    finally:
        rmtree(tempdir)
