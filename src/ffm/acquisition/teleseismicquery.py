import json
import math
import multiprocessing
import pathlib
import time
from datetime import datetime, timedelta, timezone
from typing import Dict, List, Optional

from pydantic import BaseModel, Field

from ffm.acquisition.event import get_event_detail, get_product


class TeleseismicQuery(BaseModel):
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
    ):
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
    ):
        """Get IrisQuery from id in USGS feeds or ComCat"""
        # get event detail
        detail = get_event_detail(eventid)
        if detail is None:
            raise Exception(f"Error getting event with id '{eventid}'")
        return cls.from_detail(detail=detail, source=source)

    @classmethod
    def from_origin(cls, origin: dict):
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

    @property
    def starttime(self):
        """Calculate query starttime"""
        return self.event_time - timedelta(seconds=self.seconds_before)
