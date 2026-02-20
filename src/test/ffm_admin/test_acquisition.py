import json
import pathlib
from shutil import rmtree
from tempfile import mkdtemp
from unittest import mock

from typer.testing import CliRunner

from ffm.ffm_admin.acquisition import app
from test.testutils import MockResponse

runner = CliRunner()


DETAIL = {
    "properties": {
        "mag": 6.4,
        "time": 1771036057951,
        "products": {
            "moment-tensor": [
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
                    },
                },
            ]
        },
    },
    "geometry": {
        "type": "Point",
        "coordinates": [166.6013, -14.8934, 10],
    },
}


def test_app():
    with mock.patch(target="requests.get") as mock_requests:
        mock_requests.return_value = MockResponse(json.dumps(DETAIL), 200)
        tempdir = pathlib.Path(mkdtemp())
        try:
            result = runner.invoke(app, ["eventid", "-d", str(tempdir)])
            assert result.exit_code == 0
            with open(tempdir / "eventid_cmt_CMT") as f:
                data = f.read()
            assert (
                data
                == """ US 2026 02 14 02 27 37.951000 -14.8934 166.6013 10.0  0.0 6.4 PlanetEarth
event name: eventid
time shift: 3.5
half duration: 3.5
latitude: -14.8934
longitude: 166.6013
depth: 10.0
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
