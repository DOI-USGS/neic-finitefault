import pathlib
from datetime import datetime, timezone

from pydantic import BaseModel

from ffm.acquisition.event import get_event_detail, get_product


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
    tensor_depth: float
    tensor_latitude: float
    tensor_longitude: float
    time: datetime

    @classmethod
    def from_detail(
        cls,
        detail: dict,
        eventid: str,
        source: str = "us",
        allow_unreviewed: bool = False,
    ):
        """Get CMT from USGS event detail"""
        additional_properties = []
        if not allow_unreviewed:
            additional_properties += [("review-status", "reviewed")]
        # Get origin
        origin = get_product(
            detail=detail,
            product_type="origin",
            source=source,
            additional_properties=additional_properties,
        )
        if origin is None:
            origin_error = f"Unable to get origin with source ({source})"
            if len(additional_properties) > 0:
                origin_error += (
                    f" and additional property filters: {additional_properties}"
                )
            raise Exception(origin_error)
        # get moment_tensor
        tensor_props = additional_properties + [
            ("derived-magnitude-type", "Mww"),
        ]
        moment_tensor = get_product(
            detail=detail,
            product_type="moment-tensor",
            source=source,
            additional_properties=tensor_props,
            additional_keys=["sourcetime-duration"],
        )
        if moment_tensor is None:
            tensor_error = f"Unable to get moment tensor with source ({source})"
            if len(tensor_props) > 0:
                tensor_error += f" and additional property filters: {tensor_props}"
            raise Exception(tensor_error)
        detail_props = detail["properties"]
        depth = float(origin["properties"]["depth"])
        latitude = float(origin["properties"]["latitude"])
        longitude = float(origin["properties"]["longitude"])
        magnitude = float(origin["properties"]["magnitude"])
        place = detail_props["place"].split(",")[-1].strip()
        time = datetime.strptime(
            origin["properties"]["eventtime"], "%Y-%m-%dT%H:%M:%S.%fZ"
        ).replace(tzinfo=timezone.utc)
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
    def from_id(cls, eventid: str, source: str = "us", allow_unreviewed: bool = False):
        """Get CMT from id in USGS feeds or ComCat"""
        # get event detail
        detail = get_event_detail(eventid)
        if detail is None:
            raise Exception(f"Error getting event with id '{eventid}'")
        return cls.from_detail(
            detail=detail,
            eventid=eventid,
            source=source,
            allow_unreviewed=allow_unreviewed,
        )

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
        tensor_depth = float(props["derived-depth"])
        tensor_latitude = float(props["derived-latitude"])
        tensor_longitude = float(props["derived-longitude"])
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
            tensor_depth=tensor_depth,
            tensor_latitude=tensor_latitude,
            tensor_longitude=tensor_longitude,
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
            f.write(f"latitude: {self.tensor_latitude}\n")
            f.write(f"longitude: {self.tensor_longitude}\n")
            f.write(f"depth: {self.tensor_depth}\n")
            f.write(f"Mrr: {self.mrr * 10**7}\n")
            f.write(f"Mtt: {self.mtt * 10**7}\n")
            f.write(f"Mpp: {self.mpp * 10**7}\n")
            f.write(f"Mrt: {self.mrt*10**7}\n")
            f.write(f"Mrp: {self.mrp*10**7}\n")
            f.write(f"Mtp: {self.mtp*10**7}\n")
