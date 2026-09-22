from typing import Annotated, Optional, Union

from pydantic import BaseModel, Field

from ffm.data_types.dartdata import Dart, GeoNetDart, NoaaDart


class DartStation(BaseModel):
    id: str = Field(..., description="The id of the dart station")
    azimuth: float = Field(..., description="The azimuth of the station")
    distance: float = Field(
        ..., description="The distance between the station and the event (km)"
    )
    latitude: float = Field(
        ..., description="The latitude of the dart station", ge=-90, le=90
    )
    longitude: float = Field(
        ..., description="The longitude of the dart station", ge=-180, le=180
    )
    region: Optional[str] = Field(None, description="The region of the dart station")


class GeoNetDartStation(DartStation):
    type: str = "geonet"
    data: Optional[GeoNetDart] = Field(
        None, description="Data and metadata at the station"
    )


class NoaaDartStation(DartStation):
    type: str = "noaa"
    data: Optional[NoaaDart] = Field(
        None, description="Data and metadata at the station"
    )


DartStationUnion = Annotated[
    Union[NoaaDartStation, GeoNetDartStation], Field(discriminator="type")
]
