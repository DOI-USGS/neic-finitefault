import pathlib
from shutil import rmtree
from tempfile import mkdtemp
from unittest import TestCase
from ffm.input_data.topographynetcdf import TopographyNetCDF
import pandas as pd
import numpy as np


class TestTopographyNetCDF(TestCase):

    def setUp(self):
        self.tempdir = pathlib.Path(mkdtemp())
        self.data_dir = pathlib.Path(__file__).parent.parent / "data" / "input_data"

    def tearDown(self):
        rmtree(self.tempdir)

    def test_bad_data_type(self):
        with self.assertRaises(TypeError):
            topography = TopographyNetCDF(b"")

    def test_bad_filepath(self):
        with self.assertRaises(FileNotFoundError):
            topography = TopographyNetCDF(self.tempdir / "topo.nc")

    def test_from_file(self):
        topography = TopographyNetCDF(self.data_dir / "topo.nc")

        # check against xyz
        target_z = pd.read_csv(
            self.data_dir / "topo.xyz", sep="\t", names=["x", "y", "z"]
        ).z.values
        np.testing.assert_array_almost_equal(
            topography.data.z.values, target_z, decimal=15
        )
