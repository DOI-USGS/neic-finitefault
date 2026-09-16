from typing import List, Optional

from pydantic import BaseModel, Field

from ffm.model.setup.connection import Connection
from ffm.model.setup.risetime import RiseTime
from ffm.model.setup.segmentmetadata import SegmentMetadata


class SegmentsData(BaseModel):
    rise_time: RiseTime
    segments: List[SegmentMetadata]
    connections: Optional[List[Connection]] = Field(
        None, description="Description of connection between two segments"
    )
