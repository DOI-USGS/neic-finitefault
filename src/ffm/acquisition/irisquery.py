import json
import multiprocessing
import pathlib
import time
from datetime import datetime, timedelta, timezone
from typing import Dict, List, Optional

from obspy import Stream, Trace, UTCDateTime  # type: ignore [import-untyped]
from obspy.clients.fdsn import Client  # type: ignore [import-untyped]
from obspy.clients.fdsn.header import (  # type: ignore [import-untyped]
    FDSNNoDataException,
)
from obspy.clients.iris.client import Client as Iris  # type: ignore [import-untyped]
from obspy.core.inventory import Inventory  # type: ignore [import-untyped]
from pydantic import BaseModel, Field
from ratelimiter import RateLimiter  # type: ignore [import-untyped]

from ffm.acquisition.event import get_event_detail, get_product


def rate_limit_iris(until):
    duration = int(round(until - time.time()))
    print("Rate limited, sleeping for {:d} seconds".format(duration))


class IrisQuery(BaseModel):
    depth: float = Field(..., description="The hypocenter depth")
    event_time: datetime
    latitude: float = Field(..., description="The hypocenter latitude", ge=-90, le=90)
    longitude: float = Field(
        ..., description="The hypocenter longitude", ge=-180, le=180
    )
    min_distance: float = Field(
        30, description="The minimum distance from the hypocenter"
    )
    max_distance: float = Field(
        90, description="The maximum distance from the hypocenter"
    )
    seconds_before: float = Field(
        0,
        description="The seconds before the event time to include in the search",
    )
    seconds_after: float = Field(
        3000,
        description="The seconds after the event time to include in the search",
    )

    @property
    def default_stations(self) -> Dict[str, Dict[str, List[str]]]:
        """Get the default stations for each preferred network from default_stations.json"""
        with open(
            pathlib.Path(__file__).parent / "default_stations.json"
        ) as stations_file:
            default_stations: Dict[str, Dict[str, List[str]]] = json.load(stations_file)
        return default_stations

    @property
    def duration(self):
        """Calculate query duration"""
        return (self.endtime - self.starttime).total_seconds()

    @property
    def endtime(self):
        """Calculate query endtime"""
        return self.event_time + timedelta(seconds=self.seconds_after)

    @classmethod
    def from_detail(
        cls,
        detail: dict,
        source: str = "us",
    ) -> "IrisQuery":
        """Get IrisQuery from USGS event detail"""
        # get origin
        origin = get_product(detail=detail, product_type="origin", source=source)
        if origin is None:
            raise Exception(
                f"No origin with source '{source}' found. Try specifying a different source."
            )
        return cls.from_origin(origin=origin)

    @classmethod
    def from_id(
        cls,
        eventid: str,
        source: str = "us",
    ) -> "IrisQuery":
        """Get IrisQuery from id in USGS feeds or ComCat"""
        # get event detail
        detail = get_event_detail(eventid)
        if detail is None:
            raise Exception(f"Error getting event with id '{eventid}'")
        return cls.from_detail(detail=detail, source=source)

    @classmethod
    def from_origin(cls, origin: dict) -> "IrisQuery":
        """IrisQuery from moment-tensor product"""
        props = origin["properties"]
        event_time = datetime.strptime(
            props["eventtime"], "%Y-%m-%dT%H:%M:%S.%fZ"
        ).replace(tzinfo=timezone.utc)
        return cls(
            depth=float(props["depth"]),
            event_time=event_time,
            latitude=float(props["latitude"]),
            longitude=float(props["longitude"]),
        )  # type: ignore [call-arg]

    def get_data(
        self,
        networks: List[str] = [],
        stations: Optional[Dict[str, Dict[str, List[str]]]] = None,
        debug: bool = False,
    ) -> Stream:
        """Get data for specified networks (and optionally stations) from IRIS with Obspy"""
        client = Client("IRIS", debug=debug)
        # list of network, station
        parallelized_params = []
        for network in networks:
            print(f"Querying network: {network}")
            inventory: Inventory = client.get_stations(
                channel="BH*",
                endtime=UTCDateTime(self.endtime),
                latitude=self.latitude,
                level="response",
                longitude=self.longitude,
                maxradius=self.max_distance,
                minradius=self.min_distance,
                network=network,
                starttime=UTCDateTime(self.starttime),
            )
            for station in inventory.networks[0].stations:
                parallelized_params += [(network, station.code, "BH*", debug)]
        if stations is not None:
            for network_key, network_stations in stations.items():
                for station_key, channels in network_stations.items():
                    for channel in channels:
                        print(
                            f"Querying station: {network_key} {station_key} {channel}"
                        )
                        parallelized_params += [
                            (
                                network_key,
                                station_key,
                                channel,
                                debug,
                            )
                        ]

        with multiprocessing.Pool() as pool:
            traces = pool.starmap(self._get_data, parallelized_params)
        flat_traces = []
        for trace in traces:
            flat_traces += trace
        stream = Stream(flat_traces)
        stream.merge(-1)
        return stream

    def _get_data(
        self, network: str, station: str, channel: str, debug: bool = False
    ) -> list:
        """Get data for a given network, station, channel configuration"""
        client = Client("IRIS", debug=debug)
        traces = []
        try:
            stream: Stream = client.get_waveforms(
                network=network,
                station=station,
                channel=channel,
                starttime=UTCDateTime(self.starttime),
                location="*",
                endtime=UTCDateTime(self.endtime),
            )
            traces += stream.traces
        except FDSNNoDataException as e:
            print(f"No data available for {network} {station} {channel}")
        return traces

    def get_responses(self, stream: Stream, debug: bool = False) -> Stream:
        """Get station responses corresponding to the data in the provided stream"""
        response_traces = []
        with RateLimiter(max_calls=10, period=1, callback=rate_limit_iris):
            for trace in stream:
                response_traces += [self._get_response(trace, debug)]
        return Stream(response_traces)

    @RateLimiter(max_calls=10, period=1, callback=rate_limit_iris)
    def _get_response(self, trace: Trace, debug: bool = False) -> Trace:
        """Get the trace station's response from IRIS and set the traces
        stats with the response poles and zeros"""
        client: Iris = Iris(debug=debug)
        response = client.sacpz(
            network=trace.stats.network,
            station=trace.stats.station,
            location=trace.stats.location,
            channel=trace.stats.channel,
            starttime=UTCDateTime(self.starttime),
            endtime=UTCDateTime(self.endtime),
        )

        trace.stats.response = response
        return trace

    @property
    def starttime(self):
        """Calculate query starttime"""
        return self.event_time - timedelta(seconds=self.seconds_before)

    def write_data(self, directory: pathlib.Path, stream: Stream):
        """Write all the data in the stream to the provided directory"""
        with multiprocessing.Pool() as pool:
            pool.starmap(
                self._write_data, [(directory, trace) for trace in stream.traces]
            )

    def _write_data(self, directory: pathlib.Path, trace: Trace):
        """Write the data (sac) and station response (poles and zeros) for the provided trace"""
        filename = f"{trace.stats.network}_{trace.stats.station}_{trace.stats.channel}_{trace.stats.location}.sac"
        trace.write(
            filename=str(directory / filename),
            format="sac",
        )
        with open(directory / ("SAC_PZs_" + filename), "w") as f:
            f.write(trace.stats.response.decode())
