from pydantic import BaseModel, Field


class RiseTime(BaseModel):
    delta_rise: float = Field(..., description="The delta of the risetime function")
    min_rise: float = Field(..., description="The minimuum of the risetime function")
    windows: int = Field(
        ..., description="The number of windows in the risetime function"
    )
