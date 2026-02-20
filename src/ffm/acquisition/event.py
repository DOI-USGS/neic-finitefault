import logging
from typing import Optional

import requests


def get_event_detail(eventid: str, retries: int = 3) -> Optional[dict]:
    """Get an event detail from USGS Feeds (reroutes to ComCat if event is unavailable)"""
    url = f"https://earthquake.usgs.gov/earthquakes/feed/v1.0/detail/{eventid}.geojson"
    for i in range(retries):
        try:
            response = requests.get(url)
            response.raise_for_status()
        except:
            logging.warning(f"Failed to get detail for {eventid}")
    if response.status_code != 200:
        return None
    return response.json()


def get_moment_tensor(detail: dict, source: str = "us") -> Optional[dict]:
    """Get an moment-tensor product from an event detail. Default is preferred 'us' source moment tensor"""
    products = detail["properties"]["products"]
    for moment_tensor_product in products.get("moment-tensor", []):
        if moment_tensor_product["source"] == source:
            return moment_tensor_product
    return None
