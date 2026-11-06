import io
import pathlib

import pandas as pd

from ffm.input_data.topography import Topography


class TopographyXYZ(Topography):
    def __init__(self, data: str | pathlib.Path | pd.DataFrame, **kwargs):
        """
        Topography XYZ data. Metadata defaults to assuming use of SRTM30_PLUS
        (https://topex.ucsd.edu/WWW_html/srtm30_plus.html) with 30s spacing

        :param self: XYZ instnace
        :param data: The data as a string (content), filepath, or pandas dataframe
            (dataframe's must have the columns: 'x', 'y', and 'z')
        :type data: str | pathlib.Path | pd.DataFrame
        :param kwargs: Keyword arguments of Topography
        """

        # check that a filepath hasn't been passed as a string
        if isinstance(data, str) and pathlib.Path(data).exists():
            data = pathlib.Path(data)

        # get data
        if isinstance(data, str):
            content = io.StringIO(data)
            content.seek(0)
            data = pd.read_csv(content, names=["x", "y", "z"], sep="\t")
            print(data)
        elif isinstance(data, pd.DataFrame):
            data = data
        elif isinstance(data, pathlib.Path):
            if not data.exists():
                raise FileNotFoundError(
                    f"Specified file '{data.resolve()}' does not exist"
                )
            data = pd.read_csv(data, names=["x", "y", "z"], sep="\t")
        else:
            raise TypeError(f"Invalid data type '{type(data)}'")

        super().__init__(data=data, **kwargs)
