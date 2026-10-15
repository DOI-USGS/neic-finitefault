import io
import pathlib

import pandas as pd
import requests
from pydantic import BaseModel, Field

from ffm.input_data.topographyesriascii import TopographyEsriAscii
from ffm.input_data.topographynetcdf import TopographyNetCDF
from ffm.input_data.topographyxyz import TopographyXYZ


class TopographyQuery(BaseModel):
    """
    Class to query for bathymetry data
    """

    north: float = Field(..., description="The north (maximum latitude) extent")
    south: float = Field(..., description="The north (minimum latitude) extent")
    east: float = Field(..., description="The east (maximum longitude) extent")
    west: float = Field(..., description="The west (minimum longitude) extent")

    def get_data(
        self,
        directory: pathlib.Path,
        filename: str | None = None,
        format: str = "geoclaw_ascii",
    ) -> pathlib.Path | TopographyEsriAscii | TopographyNetCDF:
        url = "https://topex.ucsd.edu/cgi-bin/get_srtm30.cgi"
        params = {
            "north": self.north,
            "south": self.south,
            "west": self.west,
            "east": self.east,
        }

        response = requests.post(url, data=params)
        response.raise_for_status()

        # convert 360 longitudes to 180
        content = io.StringIO(response.text)
        data = pd.read_csv(content, sep="\t", names=["x", "y", "z"])
        bathymetry = TopographyXYZ(data=data)

        format = format.lower()
        filepath = None
        match format:
            case "esri_ascii":
                filepath = directory / (filename or "topo.asc")
                bathymetry.write_esri_ascii(filepath)
            case "geoclaw_ascii":
                filepath = directory / (filename or "topo.topo")
                bathymetry.write_geoclaw_ascii(filepath)
            case "netcdf":
                filepath = directory / (filename or "topo.nc")
                bathymetry.write_netcdf(filepath)
            case "xyz":
                filepath = directory / (filename or "topo.xyz")
                bathymetry.write_xyz(filepath)
            case _:
                raise ValueError(
                    "Invalid format specified (format). Currently supported formats are 'esri_ascii', 'geoclaw_ascii, 'netcdf', and 'xyz'"
                )
        return filepath
