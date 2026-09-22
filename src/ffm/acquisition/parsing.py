import io
import json
import pathlib
from datetime import datetime, timezone
from typing import List, Union

import pandas as pd  # type: ignore [import-untyped]

from ffm.data_types.dartdata import GeoNetDart, NoaaDart
from ffm.data_types.dartstation import DartStation, GeoNetDartStation, NoaaDartStation


class EmptyDataError(Exception):
    """
    Class for reporting empty data when parsing
    """

    pass


def parse_geonet_dart(
    station: DartStation,
    dart_data: Union[pathlib.Path, io.StringIO, str, list],
) -> List[GeoNetDartStation]:
    """
    Parse response from GeoNet dart data api
        example: https://tilde.geonet.org.nz/v4/data/dart/NZE/water-height/-/15s/-/2021-03-05/2021-03-05

    :param station: The basic station information
    :param dart_data: The response data
    :return: Description
    :rtype: The GeoNet dart data station with data
    """
    if isinstance(dart_data, str):
        geonet_data = json.loads(dart_data)
    elif isinstance(dart_data, list):
        geonet_data = dart_data
    elif isinstance(dart_data, pathlib.Path):
        with open(dart_data) as f:
            geonet_data = json.load(f)
    elif isinstance(dart_data, io.StringIO):
        geonet_data = json.load(dart_data)
    else:
        raise ValueError(f"Invalid dart data type: {type(dart_data)}")

    data_stations = []
    for sensor in geonet_data:
        if len(sensor["data"]) == 0:
            continue
        data_station = GeoNetDartStation(**station.model_dump())
        data_station.latitude = sensor["latitude"]
        data_station.longitude = sensor["longitude"]
        data = pd.DataFrame(sensor["data"])
        # sort just in case the values come in unordered
        data = data.sort_values(by="ts")
        data_station.data = GeoNetDart(  # type: ignore [call-arg]
            datum=sensor["datum"],
            elevation=sensor["elevationM"],
            error=data.err.values,
            error_unit=sensor["errorUnit"],
            height=data.val.values,
            relative_height=sensor["relativeHeightM"],
            sensor_code=sensor["series"]["sensorCode"],
            time=[
                pd.Timestamp(t).to_pydatetime().replace(tzinfo=timezone.utc)
                for t in pd.to_datetime(data.ts, utc=True).values
            ],
            unit=sensor["valueUnit"],
            value_unit=sensor["valueUnit"],
        )
        data_stations += [data_station]
    return data_stations


def parse_noaa_dart(
    station: DartStation,
    dart_data: Union[pathlib.Path, io.StringIO, str, list],
) -> NoaaDartStation:
    """
    Parse response from NOAA dart data api
        example: https://www.ndbc.noaa.gov/dart_data.php?station=32402&startmonth=9&startday=16&startyear=2015

    :param station: The basic station information
    :param dart_data: The response data
    :return: The NOAA dart data station with data
    """
    if isinstance(dart_data, str):
        lines = dart_data.split("\n")
    elif isinstance(dart_data, list):
        lines = dart_data
    elif isinstance(dart_data, pathlib.Path):
        with open(dart_data) as f:
            lines = f.readlines()
    elif isinstance(dart_data, io.StringIO):
        dart_data.seek(0)
        lines = dart_data.readlines()
    else:
        raise ValueError(f"Invalid dart data type: {type(dart_data)}")

    noa_station = NoaaDartStation(**station.model_dump())

    if len(lines) < 10:
        raise EmptyDataError(f"No dart data available in NOAA response")
    data = pd.read_fwf(
        io.StringIO("\n".join(lines[10:])),
        names=["year", "month", "day", "hour", "minute", "second", "T", "height"],
    )
    data = data.drop(data[data.height == 9999.000].index)
    data = data.reset_index(drop=True)

    # Noaa provides data for multiple resolutions in the same response
    # this means that there may be duplicate timestamps with different values
    # in this process, we have chosen to prioritize the higher resolution value
    #  For example:
    # 2015 09 16 21 15 00 2 4053.827
    # 2015 09 16 21 15 00 1 4053.829
    #
    # the first line with T=2 would be used and the second line with T=1 would be dropped
    # first sort data by resolution (T column) with the highest resolution first
    data = data.sort_values(by="T", ascending=False)
    # drop duplicate timestamps and keep highest precision time
    data = data.drop_duplicates(
        subset=["year", "month", "day", "hour", "minute", "second"], keep="first"
    )

    # get datetime column
    data["datetime"] = [
        datetime.strptime(f"{y}-{m}-{d} {H}:{M}:{S}", "%Y-%m-%d %H:%M:%S").replace(
            tzinfo=timezone.utc
        )
        for y, m, d, H, M, S in zip(
            data["year"].values,
            data["month"].values,
            data["day"].values,
            data["hour"].values,
            data["minute"].values,
            data["second"].values,
        )
    ]
    # sort by timestamp
    data = data.sort_values(by="datetime")
    noa_station.data = NoaaDart(  # type: ignore [call-arg]
        time=[
            pd.Timestamp(t).to_pydatetime().replace(tzinfo=timezone.utc)
            for t in data["datetime"].values
        ],
        height=data["height"].values,
        unit="m",
    )
    return noa_station
