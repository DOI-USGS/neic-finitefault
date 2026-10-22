import pathlib
from shutil import rmtree
from tempfile import mkdtemp

import numpy as np
import pandas as pd
import pygmt
from geoclaw.src.python.geoclaw import topotools


class Topography:
    def __init__(
        self,
        data: pd.DataFrame,
        spacing: str = "30s",
        z_unit: str = "m",
        zero: str = "sea-level",
    ):
        """
        Topography data. Metadata defaults to assuming use of SRTM30_PLUS
        (https://topex.ucsd.edu/WWW_html/srtm30_plus.html) with 30s spacing

        :param self: XYZ instnace
        :param data: The xyz topography data
        :type data: pd.DataFrame
        :param spacing: The spacing of the data. Use GMT/PyGMT equivalent annotations.
        :type spacing: str
        :param z_unit: The unit of the z value (for SRTM30_PLUS topography, depth in meters).
                    All z values are converted to m for consistency.
        :type z_unit: str
        :param zero: The zero reference point (for SRTM30_PLUS topography, sea level)
        :type zero: str
        """
        # set  data
        for column in ["x", "y", "z"]:
            if column not in data.columns:
                raise KeyError(
                    f"Column name '{column}', must be available in the dataframe"
                )
        self.data = data.copy(deep=True)

        # convert z to meters if necessary
        supported_units = ["mm", "cm", "m", "km"]
        if z_unit not in supported_units:
            raise ValueError(
                f"Specified z_unit '{z_unit}' is currently not supported. Supported units include:{supported_units}"
            )
        match z_unit:
            case "mm":
                self.data.z = self.data.z.values * 0.001
            case "cm":
                self.data.z = self.data.z.values * 0.01
            case "m":
                pass
            case "km":
                self.data.z = self.data.z.values * 1000
            case _:
                raise ValueError(
                    f"Specified z_unit '{z_unit}' is currently not supported. Supported units include:{supported_units}"
                )
        self.z_unit = "m"
        self.spacing = spacing

        # convert 0/360 longitude data to -180/180
        self.data.x = [x - 360 if x > 180 else x for x in self.data.x.values]

        # set zero
        if zero != "sea-level":
            raise ValueError(
                f"Specified zero '{zero}' is currently not supported. Supported zeros include: 'sea-level'"
            )
        self.zero = zero

    @property
    def max_x(self) -> float:
        return np.max(self.data.x)

    @property
    def max_y(self) -> float:
        return np.max(self.data.y)

    @property
    def min_x(self) -> float:
        return np.min(self.data.x)

    @property
    def min_y(self) -> float:
        return np.min(self.data.y)

    def write_esri_ascii(self, outfile: pathlib.Path | str):
        """
        Write Esri Ascii topography data to the specified outfile

        :param self: Topography instance
        :param outfile: Path where the file format will be written
        :type outfile: pathlib.Path | str
        """
        outfile = pathlib.Path(outfile)
        pygmt.xyz2grd(
            x=self.data.x,
            y=self.data.y,
            z=self.data.z,
            spacing=self.spacing,
            outgrid=str(outfile.resolve()) + "=ei",
            region=[self.min_x, self.max_x, self.min_y, self.max_y],
        )

    def write_geoclaw_ascii(self, outfile: pathlib.Path):
        """
        Write geoclaw ascii file

        :param self: Topography instance
        :param outfile: Path where the file format will be written
        :type outfile: pathlib.Path | str
        """
        outfile = pathlib.Path(outfile)
        try:
            tempdir = pathlib.Path(mkdtemp())
            tempfile = tempdir / "topo.asc"
            self.write_esri_ascii(tempfile)
            topo = topotools.Topography(tempfile, topo_type=2)
            topo.write(outfile, no_data_value=-999999)
        finally:
            rmtree(tempdir)

    def write_netcdf(self, outfile: pathlib.Path | str):
        """
        Write NetCDF topography data to the specified outfile

        :param self: Topographyinstance
        :param outfile: Path where the file format will be written
        :type outfile: pathlib.Path | str
        """
        outfile = pathlib.Path(outfile)
        pygmt.xyz2grd(
            x=self.data.x,
            y=self.data.y,
            z=self.data.z,
            spacing=self.spacing,
            outgrid=str(outfile.resolve()) + "=nf",
            region=[self.min_x, self.max_x, self.min_y, self.max_y],
        )

    def write_xyz(self, outfile: pathlib.Path | str):
        """
        Write XYZ topography data to the specified outfile

        :param self: Topography instance
        :param outfile: Path where the file format will be written
        :type outfile: pathlib.Path | str
        """
        outfile = pathlib.Path(outfile)
        self.data.to_csv(
            outfile,
            columns=["x", "y", "z"],
            index=False,
            header=False,
            sep="\t",
        )
