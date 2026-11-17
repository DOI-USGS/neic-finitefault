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

    west: float = Field(..., description="The west (minimum longitude) extent")
    east: float = Field(..., description="The east (maximum longitude) extent")
    south: float = Field(..., description="The south (minimum latitude) extent")
    north: float = Field(..., description="The north (maximum latitude) extent")

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
        print(params)

        # check that lon/lat extent is reasonable
        if not (-90 <= params["north"] <= 90 and -90 <= params["south"] <= 90):
            raise ValueError("Latitude values must be between -90 and 90")
        if not (-180 <= params["west"] <= 360 and -180 <= params["east"] <= 360):
            raise ValueError(
                "Longitude values must be between 0 and 360 or -180 to 180"
            )
        if params["south"] > params["north"]:
            raise ValueError("South latitude must be lower than North latitude")
        if params["west"] > params["east"]:
            raise ValueError("West longitude must be lower than East longitude")

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
