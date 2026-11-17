import unittest
from unittest import mock

import pathlib
import numpy as np
from ffm.acquisition.parsing import parse_noaa_dart
from ffm.data_processing import process_dart
from ffm.data_types.dartstation import DartStation


class TestDataProcessing(unittest.TestCase):
    def test_dart_processing(self):
        dart_data_station = parse_noaa_dart(
            DartStation(
                id="station",
                azimuth=0,
                distance=10,
                latitude=20,
                longitude=100,
                region="test region",
            ),
            pathlib.Path(__file__).parent
            / "acquisition"
            / "data"
            / "32402_dart_api.txt",
        )

        processed_dart_data = process_dart(
            dart_data_station.data,
            frequency_corners=[5, 120],
            dt=60,
            filter_peak_noise=True,
        )

        bandpass_dart_target = np.genfromtxt(
            pathlib.Path(__file__).parent
            / "acquisition"
            / "data"
            / "32402_dart_filt.txt"
        )

        np.testing.assert_array_equal(
            processed_dart_data.interpolated_time, bandpass_dart_target[:, 0]
        )
        np.testing.assert_array_almost_equal(
            processed_dart_data.interpolated_height, bandpass_dart_target[:, 1]
        )
