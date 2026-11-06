import pathlib
from shutil import rmtree
from tempfile import mkdtemp
from unittest import TestCase
from ffm.input_data.topographyesriascii import TopographyEsriAscii
import pandas as pd
import numpy as np


class TestTopographyEsriAscii(TestCase):

    def setUp(self):
        self.tempdir = pathlib.Path(mkdtemp())
        self.data_dir = pathlib.Path(__file__).parent.parent / "data" / "input_data"

    def tearDown(self):
        rmtree(self.tempdir)

    def test_bad_data_type(self):
        with self.assertRaises(TypeError):
            topography = TopographyEsriAscii(b"")

    def test_bad_filepath(self):
        with self.assertRaises(FileNotFoundError):
            topography = TopographyEsriAscii(self.tempdir / "topo.asc")

    def test_from_file(self):
        topography = TopographyEsriAscii(self.data_dir / "topo.asc")

        # check against xyz
        target_z = pd.read_csv(
            self.data_dir / "topo.xyz", sep="\t", names=["x", "y", "z"]
        ).z.values
        np.testing.assert_array_almost_equal(
            topography.data.z.values, target_z, decimal=15
        )
