from typing import Optional

from pydantic import BaseModel, Field


class Hypocenter(BaseModel):
    depth: float
    lat: float
    lon: float


class SegmentMetadata(BaseModel):
    # TODO: get descriptions for each of these fields
    delay_segment: float = Field(0, description="")
    delta_dip: float = Field(..., description="")
    delta_strike: float = Field(..., description="")
    dip_subfaults: int = Field(..., description="")
    dip: float = Field(..., description="")
    hypocenter: Optional[Hypocenter] = Field(
        None, description="The hypocenter location"
    )
    hyp_dip: int = Field(..., description="")
    hyp_stk: int = Field(..., description="")
    max_vel: Optional[float] = Field(None, description="")
    min_vel: Optional[float] = Field(None, description="")
    neighbours: list = Field(..., description="")
    rake: float = Field(..., description="")
    rupture_vel: float = Field(..., description="")
    stk_subfaults: int = Field(..., description="")
    strike: float = Field(..., description="")
