import pathlib
import shutil
import tempfile
from glob import glob
from unittest import mock

from obspy import UTCDateTime  # type:ignore
from obspy.core.inventory.channel import Channel  # type:ignore

from .data import MockObspyClient


@mock.patch("obspy.clients.fdsn.Client", return_value=MockObspyClient)
@mock.patch("obspy.clients.iris.Client", return_value=MockObspyClient)
def test_acquisition(mock_iris, mock_fdsn):
    from wasp.data_acquisition import acquisition

    tempdir = tempfile.mkdtemp()
    try:
        event_time = UTCDateTime("2023-07-02T09:29:49Z")
        lat_ep = 33.826
        lon_ep = -118.881
        depth = 10.7
        acquisition(
            event_time,
            lat_ep,
            lon_ep,
            depth,
            ["strong", "body"],
            waveform_directory=pathlib.Path(tempdir),
        )
        num_files = glob(f"{tempdir}/*")
        assert len(num_files) == 5
    finally:
        print("Cleaning up test directory.")
        shutil.rmtree(tempdir)


def test_get_channel_information_manual():
    from wasp.data_acquisition import __get_channel_information_manual

    channel = Channel("code", "loc", 56, 78, 10, 12)
    channel.azimuth = 20
    channel.dip = 30

    metadata = __get_channel_information_manual(
        channel=channel, lat_ep=12, lon_ep=34, depth=10000
    )
    assert metadata == {
        "stla": 56.0,
        "stlo": 78.0,
        "khole": "loc",
        "cmpaz": 20.0,
        "cmpinc": 120.0,
        "evla": 12,
        "evlo": 34,
        "evdp": 10000,
    }
