from datetime import datetime
from typing import Annotated, List, Optional, Union

import numpy as np
from pydantic import BaseModel, ConfigDict, Field, field_serializer


class Dart(BaseModel):
    model_config = ConfigDict(arbitrary_types_allowed=True)
    time: List[datetime] = Field(..., description="Time of the observation")
    height: np.ndarray = Field(..., description="Height observation")
    unit: str = Field(..., description="The unit of the height observation")

    # data that is created after processing
    interpolated_time: Optional[np.ndarray] = Field(
        None, description="Time (seconds) for processed interpolated data"
    )
    interpolated_height: Optional[np.ndarray] = Field(
        None, description="Processed interpolated height data"
    )

    @field_serializer("time", when_used="json")
    def serialize_times(self, time: List[datetime]):
        return [t.isoformat() for t in time]

    @field_serializer(
        "height", "interpolated_time", "interpolated_height", when_used="json"
    )
    def serialize_numpy_data(self, data: np.ndarray):
        if data is None:
            return None
        return [d for d in data]


class GeoNetDart(Dart):
    datum: str = Field(..., description="The datum used for the location measurement")
    elevation: float = Field(..., description="The station elevation (m)")
    error_unit: str = Field(..., description="The error unit")
    error: np.ndarray = Field(..., description="The error recoreded")
    relative_height: float = Field(
        ..., description="The relative height of the station (m)"
    )
    sensor_code: str = Field(
        ..., description="The code of the sensor, for multiple sensors at one location"
    )
    value_unit: str = Field(..., description="The data unit")

    @field_serializer("error", when_used="json")
    def serialize_error(self, data: np.ndarray):
        return [d for d in data]


class NoaaDart(Dart):
    pass


DartDataUnion = Annotated[Union[NoaaDart, GeoNetDart], Field(discriminator="type")]
