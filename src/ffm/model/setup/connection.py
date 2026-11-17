from typing import List

from pydantic import BaseModel


class Connection(BaseModel):
    segment1: int
    segment2: int
    segment1_subfault: List[int]
    segment2_subfault: List[int]
