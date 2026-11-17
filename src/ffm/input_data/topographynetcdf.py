import pathlib
from shutil import rmtree
from tempfile import mkdtemp

import numpy as np
import pandas as pd
import pygmt

from ffm.input_data.topography import Topography


class TopographyNetCDF(Topography):
    def __init__(self, data: str | pathlib.Path | pd.DataFrame, **kwargs):
        """
        Topography NetCDF grid data. Metadata defaults to assuming use of SRTM30_PLUS
        (https://topex.ucsd.edu/WWW_html/srtm30_plus.html) with 30s spacing

        :param self: TopographyNetCDF instnace
        :param data: The filepath or xarray
            (dataframe's must have the columns: 'x', 'y', and 'z')
        :type data: str | pathlib.Path | pd.DataFrame
        :param kwargs: Keyword arguments of Topography
        """

        # check that a filepath hasn't been passed as a string
        if isinstance(data, str) and pathlib.Path(data).exists():
            data = pathlib.Path(data)
        if isinstance(data, pathlib.Path) and not data.exists():
            raise FileNotFoundError(f"Specified file '{data.resolve()}' does not exist")

        # get data
        if isinstance(data, pathlib.Path):
            data = pygmt.load_dataarray(data)
            data = pygmt.grd2xyz(data)
            data.rename(columns={"lon": "x", "lat": "y"}, inplace=True)
        elif isinstance(data, pd.DataFrame):
            pass
        else:
            raise TypeError(f"Invalid data type: {type(data)}")

        super().__init__(data=data, **kwargs)
