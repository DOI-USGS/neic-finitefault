from typing import List, Tuple

from pydantic import BaseModel, Field

from ffm.model.result.segmentdata import SegmentData
from ffm.model.setup.segmentmetadata import SegmentMetadata


class Segment(BaseModel):
    """
    Class recording a segment in a model solution
    """

    data: SegmentData
    dx: float = Field(..., description="The delta (km) of subfaults in the x direction")
    dy: float = Field(..., description="The delta (km) of subfaults in the y direction")
    metadata: SegmentMetadata
    nx: float = Field(..., description="The number of subfaults in the x direction")
    ny: float = Field(..., description="The number of subfaults in the y direction")
    polygon: List[Tuple[float, float, float]] = Field(
        ..., description="The polygon defining the extents of the segment"
    )
