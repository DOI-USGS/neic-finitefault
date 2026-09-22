import json
import pathlib
from shutil import rmtree
from tempfile import mkdtemp
from test.testutils import MockResponse
from unittest import mock

from typer.testing import CliRunner

runner = CliRunner()


DETAIL = {
    "properties": {
        "mag": 6.4,
        "place": "here, there",
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
            ],
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
            ],
        },
    },
    "geometry": {
        "type": "Point",
        "coordinates": [166.6013, -14.8934, 10],
    },
}


def test_teleseismic():
    with mock.patch(target="requests.get") as mock_requests:
        mock_requests.return_value = MockResponse(json.dumps(DETAIL), 200)
        from ffm.ffm_admin.acquisition import app

        tempdir = pathlib.Path(mkdtemp())
        try:
            result = runner.invoke(
                app, ["teleseismic", "eventid", "-d", str(tempdir), "-n", "US"]
            )
            assert result.exit_code == 0
        finally:
            rmtree(tempdir)
