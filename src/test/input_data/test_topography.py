from unittest import TestCase
from ffm.input_data.topography import Topography
import pandas as pd

import numpy as np


class TestTopography(TestCase):
    @classmethod
    def setUpClass(cls):
        cls.default_data = pd.DataFrame(
            {"x": [0, 1, 2, 3], "y": [4, 5, 6, 7], "z": [8, 9, 10, 11]}
        )

    def test_bad_dataframe(self):
        with self.assertRaises(KeyError):
            Topography(pd.DataFrame())

    def test_bad_z_unit(self):
        with self.assertRaises(ValueError):
            Topography(self.default_data, z_unit="unknown")

    def test_bad_zero(self):
        with self.assertRaises(ValueError):
            Topography(self.default_data, zero="unknown")

    def test_cm_to_m(self):
        topography = Topography(self.default_data, z_unit="cm")
        np.testing.assert_array_almost_equal(
            topography.data.z.values, [0.08, 0.09, 0.10, 0.11], decimal=16
        )

    def test_km_to_m(self):
        topography = Topography(self.default_data, z_unit="km")
        np.testing.assert_array_almost_equal(
            topography.data.z.values, [8000, 9000, 10000, 11000], decimal=16
        )

    def test_mm_to_m(self):
        topography = Topography(self.default_data, z_unit="mm")
        np.testing.assert_array_almost_equal(
            topography.data.z.values, [0.008, 0.009, 0.010, 0.011], decimal=16
        )
