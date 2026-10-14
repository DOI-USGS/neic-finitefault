import pathlib
import shutil
from tempfile import mkdtemp
from unittest import TestCase, mock
from datetime import datetime
from ffm.acquisition.dartquery import DartQuery
from obspy import read
import numpy as np
import json

DETAIL = {
    "properties": {
        "products": {
            "origin": [
                {
                    "source": "us",
                    "properties": {
                        "depth": "10",
                        "eventtime": "2015-09-16T22:54:32.00Z",
                        "latitude": "-31.57",
                        "longitude": "-71.67",
                    },
                }
            ],
        },
    }
}


class TestDartQuery(TestCase):
    def setUp(self):
        self.tempdir = pathlib.Path(mkdtemp())

    def tearDown(self):
        shutil.rmtree(self.tempdir)

    @mock.patch(target="requests.get")
    def test_query_and_write(self, mock_request):
        mock_request.return_value.status_code = 200
        mock_request.return_value.text = (
            pathlib.Path(__file__).parent / "data" / "32402_dart_api.txt"
        ).read_text()

        # setup query
        query = DartQuery(
            latitude=-31.57,
            longitude=-71.67,
            event_time=datetime(2015, 9, 16, 22, 54, 32),
            hours_after=7,
        )

        # get stations and data
        stations = query.get_stations()
        stations = query.get_data(stations)
        query.write_data(stations, directory=self.tempdir, plot=True)

        # compare processed data to target
        assert (self.tempdir / "32402_dart.png").exists()
        assert (self.tempdir / "DART" / "processed_32402_dart.sac").exists()
        st = read(self.tempdir / "DART" / "processed_32402_dart.sac")
        st_target = read(
            pathlib.Path(__file__).parent / "data" / "processed_32402_dart.sac"
        )
        np.testing.assert_allclose(st[0].data, st_target[0].data)

    @mock.patch(target="requests.get")
    def test_from_id(self, mock_request):

        mock_detail = mock.MagicMock()
        mock_detail.status_code = 200
        mock_detail.json.return_value = DETAIL
        mock_api_query = mock.MagicMock()
        mock_api_query.status_code = 200
        mock_api_query.text = (
            pathlib.Path(__file__).parent / "data" / "32402_dart_api.txt"
        ).read_text()
        mock_request.side_effect = [mock_detail, mock_api_query, mock_api_query]

        # setup query
        query = DartQuery.from_id(eventid="eventid")
        query.hours_after = 7

        # get stations and data
        stations = query.get_stations()
        stations = query.get_data(stations)
        query.write_data(stations, directory=self.tempdir, plot=True)

        # compare processed data to target
        assert (self.tempdir / "32402_dart.png").exists()
        assert (self.tempdir / "DART" / "processed_32402_dart.sac").exists()
        st = read(self.tempdir / "DART" / "processed_32402_dart.sac")
        st_target = read(
            pathlib.Path(__file__).parent / "data" / "processed_32402_dart.sac"
        )
        np.testing.assert_allclose(st[0].data, st_target[0].data)

    @mock.patch(target="requests.get")
    def test_geonet(self, mock_request):
        mock_request.return_value.status_code = 200
        mock_request.return_value.json.return_value = [
            {
                "data": [
                    {
                        "val": 5779.116,
                        "err": 0.0,
                        "qc": "",
                        "ts": "2021-03-05T00:00:00Z",
                    },
                    {
                        "val": 5779.122,
                        "err": 0.0,
                        "qc": "",
                        "ts": "2021-03-05T00:00:30Z",
                    },
                    {
                        "val": 5779.114,
                        "err": 0.0,
                        "qc": "",
                        "ts": "2021-03-05T00:00:15Z",
                    },
                ],
                "series": {
                    "domain": "dart",
                    "station": "NZE",
                    "name": "water-height",
                    "sensorCode": "40",
                    "method": "15s",
                    "aspect": "nil",
                },
                "latitude": -36.0493,
                "longitude": -177.708,
                "elevationM": 0,
                "datum": "WGS84",
                "relativeHeightM": -5779,
                "valueUnit": "m",
                "errorUnit": "m",
            },
            {
                "data": [],
                "series": {
                    "sensorCode": "41",
                },
            },
            {
                "data": [],
                "series": {
                    "sensorCode": "42",
                },
            },
            {
                "data": [],
                "series": {
                    "sensorCode": "43",
                },
            },
        ]

        # setup query
        query = DartQuery(
            latitude=-39.10518,
            longitude=176.35803,
            event_time=datetime(2021, 3, 5, 17, 46, 58),
            hours_after=5,
        )

        # get stations and data
        stations = query.get_stations()
        stations = query.get_data(stations)
        query.write_data(stations, directory=self.tempdir, raw=True)

        # compare raw data to target
        raw_data_file = self.tempdir / "NZE_dart.json"
        self.assertTrue(raw_data_file.exists())
        with open(raw_data_file) as f:
            dart_data = json.load(f)

        self.assertDictContainsSubset(
            {
                "data": {
                    "time": [
                        "2021-03-05T00:00:00+00:00",
                        "2021-03-05T00:00:15+00:00",
                        "2021-03-05T00:00:30+00:00",
                    ],
                    "height": [5779.116, 5779.114, 5779.122],
                    "unit": "m",
                    "interpolated_time": None,
                    "interpolated_height": None,
                    "datum": "WGS84",
                    "elevation": 0.0,
                    "error_unit": "m",
                    "error": [0.0, 0.0, 0.0],
                    "relative_height": -5779.0,
                    "sensor_code": "40",
                    "value_unit": "m",
                },
            },
            dart_data,
        )
        self.assertDictContainsSubset(
            {
                "id": "NZE",
                "latitude": -36.0493,
                "longitude": -177.708,
                "region": "NewZealand",
                "type": "geonet",
            },
            dart_data,
        )
