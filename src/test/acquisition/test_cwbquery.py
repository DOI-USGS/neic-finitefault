from datetime import datetime, timezone
import pathlib
from unittest import mock

from ffm.acquisition.cwbquery import CwbQuery

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


def test_irisquery_end_to_end():
    detail_cwbquery = CwbQuery.from_detail(DETAIL, "us")
    assert detail_cwbquery.depth == 10
    assert detail_cwbquery.latitude == -14.8934
    assert detail_cwbquery.longitude == 166.6013
    assert detail_cwbquery.duration == 3000
    assert detail_cwbquery.event_time == datetime(
        2026, 2, 14, 2, 27, 37, 951000, tzinfo=timezone.utc
    )
    with mock.patch(target="os.chdir"):
        with mock.patch(
            target="ffm.acquisition.cwbquery.subprocess"
        ) as mock_subprocess:
            mocked_process = mock.MagicMock()
            mock_subprocess.Popen.return_value = mocked_process
            mocked_process.communicate.return_value = (None, None)
            mocked_process.returncode = 0

            commands = detail_cwbquery.get_data(
                directory="here",
                edge_cwb_jar_path="CWBQuery.jar",
                java="/java",
                networks=["US"],
                stations={"US": {"STAT1": ["BH*"]}},
            )
            assert len(commands) == 2
            assert (
                commands[0]
                == '/java -jar CWBQuery.jar -s US.....BH. -delazc 30:90:-14.8934:166.6013 -b "2026/02/14 02:27:37" -d 3000.0 -nogaps -sacpz nm'
            )
            assert (
                commands[1]
                == '/java -jar CWBQuery.jar -s USSTAT1BH* -delazc 30:90:-14.8934:166.6013 -b "2026/02/14 02:27:37" -d 3000.0 -nogaps -sacpz nm'
            )
