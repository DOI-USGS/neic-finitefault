import pathlib
import shutil
from tempfile import mkdtemp
from unittest import TestCase, mock
from typer.testing import CliRunner

from ffm.ffm_admin.acquisition import app

runner = CliRunner()

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


class TestDart(TestCase):
    def setUp(self):
        self.tempdir = pathlib.Path(mkdtemp())

    def tearDown(self):
        shutil.rmtree(self.tempdir)

    def test_dart(self):
        with mock.patch(target="requests.get") as mock_requests:
            mock_requests.return_value.json.return_value = DETAIL
            mock_requests.return_value.status_code = 200

            from ffm.ffm_admin.acquisition import app

            result = runner.invoke(
                app,
                [
                    "dart",
                    "--eventid",
                    "eventid",
                    "-d",
                    str(self.tempdir),
                ],
            )
            assert result.exit_code == 0
