import pathlib

from obspy import UTCDateTime
from obspy.clients.fdsn import Client
from obspy.clients.iris import Client as IRIS_Client

DATA_DIR = pathlib.Path(__file__).parent


def get_test_event():
    """M 1.3 - 8 km WNW of Port Hadlock-Irondale, Washington

    Data source: https://earthquake.usgs.gov/earthquakes/eventpage/uw61938401/executive
    Data retrieved: 7/6/2023
    """
    return {
        "eventid": "uw61938401",
        "lat": 48.062,
        "lon": -122.885,
        "time": UTCDateTime("2023-07-05T09:09:27Z"),
    }


def generate_strongmotion_iris_inventory():
    """Get a strong motion inventory for a small event to use for testing"""
    event = get_test_event()

    # Get the inventory
    client = Client("IRIS")
    inventory = client.get_stations(
        starttime=event["time"] - 1 * 60,
        endtime=event["time"] + 1 * 60,
        network="C,C1,II,IU",
        channel="HN*",
        level="response",
        maxradius=10,
        latitude=event["lat"],
        longitude=event["lon"],
    )
    inventory.write(
        str(DATA_DIR / "strongmotion_iris_inventory.xml"),
        format="STATIONXML",
    )


def generate_teleseismic_iris_inventory():
    """Get a teleseismic inventory for a small event to use for testing"""
    event = get_test_event()

    # Get the inventory
    client = Client("IRIS")
    inventory = client.get_stations(
        starttime=event["time"] - 1 * 60,
        endtime=event["time"] + 1 * 60,
        network="II,G,IU,GE",
        channel="BH*",
        level="response",
        minradius=30,
        maxradius=35,
        latitude=event["lat"],
        longitude=event["lon"],
    )
    inventory.write(
        str(DATA_DIR / "teleseismic_iris_inventory.xml"),
        format="STATIONXML",
    )


if __name__ == "__main__":
    generate_strongmotion_iris_inventory()
    generate_teleseismic_iris_inventory()
