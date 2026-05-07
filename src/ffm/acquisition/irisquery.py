import json
import math
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

from ffm.acquisition.teleseismicquery import TeleseismicQuery


def rate_limit_iris(until):
    duration = int(math.ceil(until - time.time()))
    print("Rate limited, sleeping for {:d} seconds".format(duration))


class IrisQuery(TeleseismicQuery):

    def get_data(
        self,
        networks: List[str] = [],
        stations: Optional[Dict[str, Dict[str, List[str]]]] = None,
        debug: bool = False,
    ) -> Stream:
        """Get data for specified networks (and optionally stations) from IRIS with Obspy"""
        client = Client("IRIS", debug=debug)
        # list of network, station
        traces = []
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
                traces += self._get_data(network, station.code, "BH*", debug)
        if stations is not None:
            for network_key, network_stations in stations.items():
                for station_key, channels in network_stations.items():
                    for channel in channels:
                        print(
                            f"Querying station: {network_key} {station_key} {channel}"
                        )
                    traces += self._get_data(
                        network_key,
                        station_key,
                        channel,
                        debug,
                    )
        stream = Stream(traces)
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
