import json
import math
import multiprocessing
import pathlib
import time
from datetime import datetime, timedelta, timezone
from typing import Dict, List, Optional

from pydantic import BaseModel, Field

from ffm.acquisition.event import get_event_detail, get_product
from ffm.acquisition.waveformquery import WaveformQuery


class StrongMotionQuery(WaveformQuery):
    min_distance: float = Field(
        0, description="The minimum distance from the hypocenter"
    )
    max_distance: float = Field(
        5, description="The maximum distance from the hypocenter"
    )
    seconds_before: float = Field(
        60,
        description="The seconds before the event time to include in the search",
    )
    seconds_after: float = Field(
        300,
        description="The seconds after the event time to include in the search",
    )
