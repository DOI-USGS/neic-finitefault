import io
import json
import logging
import pathlib
from datetime import datetime, timedelta, timezone, tzinfo
from typing import List, Tuple

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd  # type: ignore [import-untyped]
import requests
from obspy import Stream, Trace, UTCDateTime  # type: ignore [import-untyped]
from pydantic import BaseModel, Field

import ffm.management as mng
from ffm.acquisition.event import get_event_detail, get_product
from ffm.acquisition.parsing import (
    EmptyDataError,
    parse_geonet_dart,
    parse_noaa_dart,
)
from ffm.data_management import _dict_trace
from ffm.data_processing import process_dart
from ffm.data_types.dartdata import Dart, GeoNetDart
from ffm.data_types.dartstation import (
    DartStation,
    DartStationUnion,
    GeoNetDartStation,
    NoaaDartStation,
)
from ffm.seismic_tensor import get_tensor


class DartQuery(BaseModel):
    """
    Class to query for DART data
    """

    dt: float = Field(
        60.0,
        description="Desired data sampling rate, in seconds, default is 60s",
    )
    event_time: datetime
    hours_after: int = Field(
        8,
        description="How many hours after the event to get data for",
    )
    hours_before: int = Field(
        1,
        description="How many hours before the event to get data for",
    )
    latitude: float = Field(..., description="The hypocenter latitude", ge=-90, le=90)
    longitude: float = Field(
        ..., description="The hypocenter longitude", ge=-180, le=180
    )
    max_distance: float = Field(
        20,
        description="Maximum station hypocentral distance, in kilometers, default is 20",
    )

    @property
    def endtime(self) -> datetime:
        """
        Calculate query endtime

        :param self: DartQuery instance
        :return: A datetime representing hours_after after the event_time
        """
        return self.event_time + timedelta(seconds=self.hours_after * 60 * 60)

    @classmethod
    def from_detail(
        cls,
        detail: dict,
        source: str = "us",
    ):
        """
        Get DartQuery from USGS event detail

        :param cls: DartQuery
        :param detail: USGS event detail dictionary
        :param source: The source of the origin product to get for event information
        """
        # get origin
        origin = get_product(detail=detail, product_type="origin", source=source)
        if origin is None:
            raise Exception(
                f"No origin with source '{source}' found. Try specifying a different source."
            )
        return cls.from_origin(origin=origin)

    @classmethod
    def from_gcmt_tensor_file(cls, filepath: pathlib.Path):
        """
        Get DartQuery from gcmt file

        :param cls: DartQuery
        :param filepath: Path to the gcmt tensor file
        """
        tensor_info = get_tensor(cmt_file=filepath)
        return cls(  # type: ignore [call-arg]
            event_time=datetime.strptime(
                tensor_info["datetime"], "%Y-%m-%dT%H:%M:%S"
            ).replace(tzinfo=timezone.utc),
            latitude=tensor_info["lat"],
            longitude=tensor_info["lon"],
        )

    @classmethod
    def from_id(
        cls,
        eventid: str,
        source: str = "us",
    ):
        """
        Get DartQuery from USGS feeds or ComCat

        :param cls: DartQuery
        :param eventid: The event ID to search the catalog for
        :param source: The source of the origin product to get for event information
        """
        # get event detail
        detail = get_event_detail(eventid)
        if detail is None:
            raise Exception(f"Error getting event with id '{eventid}'")
        return cls.from_detail(detail=detail, source=source)

    @classmethod
    def from_origin(cls, origin: dict):
        """
        Get DartQuery from origin product

        :param cls: DartQuery
        :param origin: An origin product dictionary
        """
        props = origin["properties"]
        event_time = datetime.strptime(
            props["eventtime"], "%Y-%m-%dT%H:%M:%S.%fZ"
        ).replace(tzinfo=timezone.utc)
        return cls(
            event_time=event_time,
            latitude=float(props["latitude"]),
            longitude=float(props["longitude"]),
        )  # type: ignore [call-arg]

    @classmethod
    def from_qcmt_tensor_file(cls, filepath: pathlib.Path):
        """
        Get DartQuery from qcmt file

        :param cls: DartQuery
        :param filepath: Path to the qcmt tensor file
        """
        tensor_info = get_tensor(quake_file=filepath)
        return cls(  # type: ignore [call-arg]
            event_time=datetime.strptime(
                tensor_info["datetime"], "%Y-%m-%dT%H:%M:%S"
            ).replace(tzinfo=timezone.utc),
            latitude=tensor_info["lat"],
            longitude=tensor_info["lon"],
        )

    def get_data(self, stations: List[DartStation]) -> List[DartStationUnion]:
        """
        Update stations with their specific station type and get data from the relevant source

        :param self: DartQuery instance
        :param stations: Stations with basic information
        :return: Stations specific to the data source with data
        """
        stations_with_data: List[DartStationUnion] = []
        for station in stations:
            try:
                if station.id.isnumeric():
                    stations_with_data += self.query_noaa(station)
                else:
                    stations_with_data += self.query_geonet(station)
            except Exception as e:
                logging.warning(f"Unable to get data for {station.id}: {str(e)}")

        return stations_with_data

    def get_stations(
        self,
        dart_station_file: pathlib.Path = pathlib.Path(__file__).parent
        / "DART_station_locations.txt",
    ) -> List[DartStation]:
        """
        Read a DART stations file

        :param self: DartQuery instance
        :param dart_station_file: Path to a DART stations file in format
                Station_ID Latitude Longitude Region
                <id> <latitude> <longitude> <region description>
        :return: DartStations with basic information
        """
        df_dart = pd.read_csv(
            dart_station_file,
            header=1,
            index_col=False,
            names=["id", "latitude", "longitude", "region"],
            sep=r"\s+",
        )

        stations = []
        for station in df_dart.to_dict(orient="records"):
            distance_km, azimuth, _ = mng._distazbaz(
                station_lat=station["latitude"],
                station_lon=station["longitude"],
                event_lat=self.latitude,
                event_lon=self.longitude,
            )
            dart_station = DartStation(
                **{**{"azimuth": azimuth, "distance": distance_km}, **station}
            )
            dist_degrees = distance_km / 111.11
            if dist_degrees < self.max_distance:
                stations += [dart_station]
        return stations

    def query_geonet(self, station: DartStation) -> List[GeoNetDartStation]:
        """
        Query the GeoNet DART data API

        :param self: DartQuery instance
        :param station: DartStation with basic station information
        :return: GeoNetDartStations with data
        """
        # query noaa for the station data
        response = requests.get(
            (
                "https://tilde.geonet.org.nz/v4/data/"
                f"dart/{station.id}/water-height/-/15s/-/"
                f"/{self.starttime.year}-{self.starttime.month:02d}-{self.starttime.day:02d}/"
                f"{self.endtime.year}-{self.endtime.month:02d}-{self.endtime.day:02d}"
            ),
        )
        response.raise_for_status()
        dart_data_stations = parse_geonet_dart(station, response.json())
        return dart_data_stations

    def query_noaa(self, station: DartStation) -> List[NoaaDartStation]:
        """
        Query the NOAA DART data API

        :param self: DartQuery instance
        :param station: DartStation with basic station information
        :return: NoaaDartStations with data
        """
        # query noaa for the station data
        response = requests.get(
            "http://www.ndbc.noaa.gov/dart_data.php",
            params={
                "station": station.id,
                "startmonth": str(self.starttime.month),
                "startday": str(self.starttime.day),
                "startyear": str(self.starttime.year),
                "endmonth": str(self.endtime.month),
                "endday": str(self.endtime.day),
                "endyear": str(self.endtime.year),
            },
            headers={"Accept": "text/csv"},
        )
        response.raise_for_status()
        # check that we have data
        try:
            station = parse_noaa_dart(station, response.text)
        except EmptyDataError:
            raise Exception(f"No data available in response from: {response.url}")

        return [station]

    @property
    def starttime(self) -> datetime:
        """
        Calculate query starttime

        :param self: DartQuery instance
        :return: A datetime representing hours_before before the event_time
        """
        return self.event_time - timedelta(seconds=self.hours_before * 60 * 60)

    def write_data(
        self,
        stations: List[DartStationUnion],
        directory: pathlib.Path,
        frequency_corners: Tuple[float, float] = (5, 120),
        dt: float = 60,
        raw: bool = False,
        plot: bool = False,
    ):
        """
        Write dart data to a processed file or to a raw data file

        :param self: DartQuery instance
        :param stations: GeoNetDartStations and/or NoaaDartStations with data
        :param directory: Path to the location where files will be written (processed data is put in the DART subfolder)
        :param frequency_corners: List of corners to use for processing bandpass filter (not used if raw = True)
        :param dt: The data time delta (not used if raw = True)
        :param raw: Dump the raw data as json
        :param plot: Plot the data (not used if raw = True)
        """
        metadata = []
        for station in stations:
            filename = f"{station.id}_dart"
            if not raw:
                if station.data is None:
                    raise ValueError("Cannot write data for a station with empty data")
                filename = f"processed_{filename}.sac"
                station.data = process_dart(
                    station.data,
                    frequency_corners=frequency_corners,
                    dt=dt,
                    filter_peak_noise=(
                        True if isinstance(station, NoaaDartStation) else False
                    ),
                )
                # setup trace and stats
                trace = Trace(data=station.data.interpolated_height)
                trace.stats.delta = dt
                trace.stats.station = station.id
                trace.stats.channel = "dart"
                trace.stats.starttime = UTCDateTime(station.data.time[0])
                trace.stats.sac = {"stla": station.latitude, "stlo": station.longitude}
                # in place, trim from one hour before origin time to 7 hours later
                trace.trim(
                    starttime=UTCDateTime(self.starttime),
                    endtime=UTCDateTime(self.endtime),
                )
                # demean data with first time before the origin time
                pre_event = trace.copy().trim(
                    starttime=UTCDateTime(self.starttime),
                    endtime=UTCDateTime(self.event_time),
                )
                if len(trace.data) == 0 or len(trace.data[~np.isnan(trace.data)]) == 0:
                    print(
                        f"No data left after processing and filtering for {station.id}"
                    )
                    continue
                trace.data = trace.data - np.mean(pre_event.data)
                # create and write stream
                dart_directory = directory / "DART"
                if not dart_directory.exists():
                    dart_directory.mkdir(parents=True)
                filepath = dart_directory / filename

                print(f"Writing data for {station.id}")
                Stream([trace]).write(str(filepath), format="SAC", byteorder=0)
                # get metadata for waveform to add to dart_waves.json
                metadata += [
                    _dict_trace(
                        file=str(filepath.resolve()),
                        name_station=str(station.id),
                        azimuth=station.azimuth,
                        channel="dart",
                        distance=station.distance / 111.11,
                        dt=dt,
                        duration=len(trace.data),
                        n_start_obs=len(pre_event),
                        trace_weight=1.0,
                        wavelet_weight="0 0 2 2 2 2 2 2\n",
                        observed_trace=None,
                        synthetic_trace=None,
                        location=[station.latitude, station.longitude],
                    )
                ]
                if plot:
                    # trace.plot(filename=str(directory / f"{station.id}_dart.png"))
                    plt.figure(figsize=(15, 5))
                    plt.plot(
                        trace.times(reftime=UTCDateTime(self.event_time)) / (3600),
                        trace.data * 100.0,
                    )
                    plt.xlim((-1, 7))
                    plt.xlabel("Time After Origin (hrs)")
                    plt.ylabel("Tsunami Amplitude (cm)")
                    plt.savefig(directory / f"{station.id}_dart.png")
                    plt.close()
            else:
                filename = f"{filename}.json"
                with open(directory / filename, "w") as f:
                    f.write(station.model_dump_json())

        if len(metadata) > 0:
            with open(directory / "dart_waves.json", "w") as f:
                json.dump(
                    metadata,
                    f,
                    sort_keys=True,
                    indent=4,
                    separators=(",", ": "),
                    ensure_ascii=False,
                )
