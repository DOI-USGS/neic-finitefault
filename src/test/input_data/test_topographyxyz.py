import pathlib
from shutil import rmtree
from tempfile import mkdtemp
from unittest import TestCase

import xarray
from ffm.input_data.topographyesriascii import TopographyEsriAscii
from ffm.input_data.topographynetcdf import TopographyNetCDF
from ffm.input_data.topographyxyz import TopographyXYZ
import pandas as pd

import numpy as np


class TestTopographyXYZ(TestCase):

    def setUp(self):
        self.tempdir = pathlib.Path(mkdtemp())
        self.default_data = {"x": [0, 1, 2, 3], "y": [4, 5, 6, 7], "z": [8, 9, 10, 11]}
        self.data_dir = pathlib.Path(__file__).parent.parent / "data" / "input_data"

    def tearDown(self):
        rmtree(self.tempdir)

    def test_bad_data_type(self):
        with self.assertRaises(TypeError):
            topography = TopographyXYZ(b"")

    def test_bad_filepath(self):
        with self.assertRaises(FileNotFoundError):
            topography = TopographyXYZ(self.tempdir / "topo.xyz")

    def test_from_content(self):
        content = ""
        for x, y, z in zip(
            self.default_data["x"], self.default_data["y"], self.default_data["z"]
        ):
            content += f"{x}\t{y}\t{z}\n"
        topography = TopographyXYZ(content)
        np.testing.assert_array_equal(topography.data.z.values, [8, 9, 10, 11])

    def test_from_dataframe(self):
        topography = TopographyXYZ(pd.DataFrame(self.default_data))
        np.testing.assert_array_equal(topography.data.z.values, [8, 9, 10, 11])

    def test_from_file(self):
        content = ""
        for x, y, z in zip(
            self.default_data["x"], self.default_data["y"], self.default_data["z"]
        ):
            content += f"{x}\t{y}\t{z}\n"

        filepath = self.tempdir / "file.xyz"
        with open(filepath, "w") as f:
            f.write(content)

        # filepath as path
        topography = TopographyXYZ(filepath)
        np.testing.assert_array_equal(topography.data.z.values, [8, 9, 10, 11])

        # filepath as string
        topography = TopographyXYZ(str(filepath))
        np.testing.assert_array_equal(topography.data.z.values, [8, 9, 10, 11])

    def test_to_esri_ascii(self):
        topography = TopographyXYZ(self.data_dir / "topo.xyz")
        filepath = self.tempdir / "file.asc"
        esri = topography.write_esri_ascii(filepath)

        self.assertTrue(filepath.exists())
        with open(self.data_dir / "topo.asc") as f:
            target_lines = f.readlines()
        with open(filepath) as f:
            lines = f.readlines()
        # ncols
        self.assertEqual(lines.pop(0), target_lines.pop(0))
        # nrows
        self.assertEqual(lines.pop(0), target_lines.pop(0))
        # xlower
        self.assertEqual(
            round(float(lines.pop(0).split(" ")[-1])),
            float(target_lines.pop(0).split(" ")[-1]),
        )
        # ylower
        self.assertEqual(
            round(float(lines.pop(0).split(" ")[-1])),
            float(target_lines.pop(0).split(" ")[-1]),
        )
        # cellsize
        np.testing.assert_almost_equal(
            float(lines.pop(0).split(" ")[-1]),
            float(target_lines.pop(0).split(" ")[-1]),
            decimal=10,
        )
        # nodata_value
        self.assertEqual(lines.pop(0), target_lines.pop(0))
        # compare rest of data values
        for line, target_line in zip(lines, target_lines):
            data = [float(v.strip()) for v in line.split("\t")]
            target_data = [float(v.strip()) for v in target_line.split("\t")]
            np.testing.assert_array_almost_equal(data, target_data, decimal=10)

    def test_to_geoclaw_ascii(self):
        topography = TopographyXYZ(self.data_dir / "topo.xyz")
        filepath = self.tempdir / "file.topo"
        return_value = topography.write_geoclaw_ascii(filepath)

        self.assertIsNone(return_value)
        self.assertTrue(filepath.exists())
        with open(self.data_dir / "topo.topo") as f:
            target_lines = f.readlines()
        with open(filepath) as f:
            lines = f.readlines()
        # ncols
        self.assertEqual(lines.pop(0), target_lines.pop(0))
        # nrows
        self.assertEqual(lines.pop(0), target_lines.pop(0))
        # xlower
        np.testing.assert_almost_equal(
            float(lines.pop(0).split("xlower")[0]),
            float(target_lines.pop(0).split("xlower")[0]),
            decimal=8,
        )
        # ylower
        np.testing.assert_almost_equal(
            float(lines.pop(0).split("ylower")[0]),
            float(target_lines.pop(0).split("ylower")[0]),
            decimal=10,
        )
        # cellsize
        np.testing.assert_almost_equal(
            float(lines.pop(0).split("cellsize")[0]),
            float(target_lines.pop(0).split("cellsize")[0]),
            decimal=10,
        )
        # nodata_value
        np.testing.assert_almost_equal(
            float(lines.pop(0).split("nodata_value")[0]),
            float(target_lines.pop(0).split("nodata_value")[0]),
            decimal=10,
        )
        # compare rest of data values
        data = [float(v.strip()) for v in lines]
        target_data = [float(v.strip()) for v in target_lines]
        np.testing.assert_array_almost_equal(data, target_data, decimal=10)

    def test_to_netcdf(self):
        topography = TopographyXYZ(self.data_dir / "topo.xyz")
        filepath = self.tempdir / "file.nc"
        topography.write_netcdf(filepath)

        self.assertTrue(filepath.exists())
        array = xarray.open_dataset(self.data_dir / "topo.nc")
        target_array = xarray.open_dataset(filepath)
        np.testing.assert_array_almost_equal(
            array["lon"], target_array["lon"], decimal=9
        )
        np.testing.assert_array_almost_equal(
            array["lat"], target_array["lat"], decimal=10
        )
        np.testing.assert_array_almost_equal(array["z"], target_array["z"], decimal=10)

    def test_write(self):
        topography = TopographyXYZ(pd.DataFrame(self.default_data))
        filepath = self.tempdir / "file.xyz"
        topography.write_xyz(filepath)

        # check file
        with open(filepath) as f:
            lines = f.readlines()
        self.assertListEqual(
            lines, ["0\t4\t8\n", "1\t5\t9\n", "2\t6\t10\n", "3\t7\t11\n"]
        )
