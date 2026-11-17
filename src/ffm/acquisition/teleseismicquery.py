import json
import pathlib
from typing import Dict, List

from pydantic import Field

from ffm.acquisition.waveformquery import WaveformQuery


class TeleseismicQuery(WaveformQuery):
    min_distance: float = Field(
        30, description="The minimum distance from the hypocenter"
    )
    max_distance: float = Field(
        90, description="The maximum distance from the hypocenter"
    )
    seconds_before: float = Field(
        0,
        description="The seconds before the event time to include in the search",
    )
    seconds_after: float = Field(
        3000,
        description="The seconds after the event time to include in the search",
    )

    @property
    def default_stations(self) -> Dict[str, Dict[str, List[str]]]:
        """Get the default stations for each preferred network from default_stations.json"""
        with open(
            pathlib.Path(__file__).parent / "default_stations.json"
        ) as stations_file:
            default_stations: Dict[str, Dict[str, List[str]]] = json.load(stations_file)
        return default_stations
