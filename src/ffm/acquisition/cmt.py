import pathlib
from datetime import datetime, timezone

from pydantic import BaseModel

from ffm.acquisition.event import get_event_detail, get_moment_tensor


class Cmt(BaseModel):
    """Adapted from usID2CMT.py by wyeck created Tue Aug 27 08:16:36 2024"""

    depth: float
    duration: float
    eventid: str
    latitude: float
    longitude: float
    magnitude: float
    magnitude_type: str
    mpp: float
    mrp: float
    mrr: float
    mrt: float
    mtp: float
    mtt: float
    place: str
    time: datetime

    @classmethod
    def from_detail(cls, detail: dict, eventid: str, source: str = "us"):
        """Get CMT from USGS event detail"""
        # get moment_tensor
        moment_tensor = get_moment_tensor(detail, source=source)
        if moment_tensor is None:
            raise Exception(
                f"No moment tensor with source '{source}' found. Try specifying a different source."
            )
        props = detail["properties"]
        coordinates = detail["geometry"]["coordinates"]
        depth = float(coordinates[2])
        latitude = float(coordinates[1])
        longitude = float(coordinates[0])
        magnitude = float(props["mag"])
        place = props["place"].split(",")[-1].strip()
        time = datetime.fromtimestamp(float(props["time"]) / 1000, tz=timezone.utc)
        return cls.from_moment_tensor(
            depth=depth,
            eventid=eventid,
            latitude=latitude,
            longitude=longitude,
            magnitude=magnitude,
            moment_tensor=moment_tensor,
            place=place,
            time=time,
        )

    @classmethod
    def from_id(cls, eventid: str, source: str = "us"):
        """Get CMT from id in USGS feeds or ComCat"""
        # get event detail
        detail = get_event_detail(eventid)
        if detail is None:
            raise Exception(f"Error getting event with id '{eventid}'")
        return cls.from_detail(detail=detail, eventid=eventid, source=source)

    @classmethod
    def from_moment_tensor(
        cls,
        depth: float,
        eventid: str,
        latitude: float,
        longitude: float,
        magnitude: float,
        moment_tensor: dict,
        place: str,
        time: datetime,
    ):
        """CMT from moment-tensor product"""
        if time.tzinfo != timezone.utc:
            raise Exception("Datetimes should be in UTC timezone")
        props = moment_tensor["properties"]
        magnitude_type = props.get("derived-magnitude-type", "??")
        mpp = float(props["tensor-mpp"])
        mrr = float(props["tensor-mrr"])
        mtt = float(props["tensor-mtt"])
        mrt = float(props["tensor-mrt"])
        mtp = float(props["tensor-mtp"])
        mrp = float(props["tensor-mrp"])
        duration = float(props["sourcetime-duration"])
        return cls(
            depth=depth,
            duration=duration,
            eventid=eventid,
            latitude=latitude,
            longitude=longitude,
            magnitude=magnitude,
            magnitude_type=magnitude_type,
            mpp=mpp,
            mrp=mrp,
            mrr=mrr,
            mrt=mrt,
            mtp=mtp,
            mtt=mtt,
            place=place,
            time=time,
        )

    def write(self, filepath: pathlib.Path):
        """Write CMT to file"""
        with open(filepath, "w") as f:
            time = self.time.strftime("%Y %m %d %H %M %S.%f")
            f.write(
                f" US {time} {self.latitude} {self.longitude} {self.depth}  0.0 {self.magnitude} {self.place}\n"
            )
            f.write(f"event name: {self.eventid}\n")
            f.write(f"time shift: {self.duration / 2.0}\n")
            f.write(f"half duration: {self.duration / 2.0}\n")
            f.write(f"latitude: {self.latitude}\n")
            f.write(f"longitude: {self.longitude}\n")
            f.write(f"depth: {self.depth}\n")
            f.write(f"Mrr: {self.mrr * 10**7}\n")
            f.write(f"Mtt: {self.mtt * 10**7}\n")
            f.write(f"Mpp: {self.mpp * 10**7}\n")
            f.write(f"Mrt: {self.mrt*10**7}\n")
            f.write(f"Mrp: {self.mrp*10**7}\n")
            f.write(f"Mtp: {self.mtp*10**7}\n")
