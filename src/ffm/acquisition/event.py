import logging
from typing import List, Optional, Tuple

import requests


def get_event_detail(eventid: str, retries: int = 3) -> Optional[dict]:
    """Get an event detail from USGS Feeds (reroutes to ComCat if event is unavailable)"""
    url = f"https://earthquake.usgs.gov/earthquakes/feed/v1.0/detail/{eventid}.geojson"
    for i in range(retries):
        try:
            response = requests.get(url)
            response.raise_for_status()
            break
        except:
            logging.warning(f"Failed to get detail for {eventid}")
    if response.status_code != 200:
        return None
    return response.json()


def get_product(
    detail: dict,
    product_type: str,
    source: str = "us",
    additional_properties: List[Tuple[str, str]] = [],
    additional_keys: List[str] = [],
):
    """Get product from event detail"""
    products = detail["properties"]["products"]
    for product in products.get(product_type, []):
        if product["source"] == source:
            invalid = False
            for additional_prop in additional_properties:
                key = additional_prop[0]
                value = additional_prop[1]
                if str(product.get("properties", {}).get(key)) != value:
                    invalid = True
                    break
            for key in additional_keys:
                if key not in list(product.get("properties", {}).keys()):
                    invalid = True
                    break
            if invalid:
                continue
            return product
    return None
