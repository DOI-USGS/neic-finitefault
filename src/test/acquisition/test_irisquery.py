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
        detail_irisquery = IrisQuery("teleseismic", "eventid", "us")
        assert detail_irisquery.query.depth == 10
        assert detail_irisquery.query.latitude == -14.8934
        assert detail_irisquery.query.longitude == 166.6013
        assert detail_irisquery.query.duration == 3000
        assert detail_irisquery.query.event_time == datetime(
            2026, 2, 14, 2, 27, 37, 951000, tzinfo=timezone.utc
        )


def test_irisquery_end_to_end():
    with mock.patch(target="requests.get") as mock_requests:
        mock_requests.return_value = MockResponse(json.dumps(DETAIL), 200)
        detail_irisquery = IrisQuery("teleseismic", "eventid", "us")
    assert detail_irisquery.query.depth == 10
    assert detail_irisquery.query.latitude == -14.8934
    assert detail_irisquery.query.longitude == 166.6013
    assert detail_irisquery.query.duration == 3000
    assert detail_irisquery.query.event_time == datetime(
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

        data = detail_irisquery.get_teleseismic_data(
            networks=["US"], stations={"US": {"STAT1": ["BH*"]}}
        )
        assert len(data) == 2

        data = detail_irisquery.get_strongmotion_data(
            networks=["US"],
        )
