import pathlib
import trace
from typing import Literal

from obspy import Stream, Trace
from obspy.core.inventory import Inventory, read_inventory

DATA_DIR = pathlib.Path(__file__).parent


def get_test_inventory(data: Literal["strong", "body"]) -> Inventory:
    """Get a test inventory"""
    if data == "strong":
        return read_inventory(DATA_DIR / "strongmotion_iris_inventory.xml")
    elif data == "body":
        return read_inventory(DATA_DIR / "strongmotion_iris_inventory.xml")
    else:
        raise ValueError(f"'{data}' is not a valid option for getting a test inventory")


def get_test_waveforms(data: Literal["strong", "body"]) -> Stream:
    trace1 = Trace()
    if data == "strong":
        trace2 = Trace()
        return Stream(traces=[trace1, trace2])
    elif data == "body":
        trace1.stats["channel"] = "BHZ"
        return Stream(traces=[trace1])
    else:
        raise ValueError(f"'{data}' is not a valid option for getting a test inventory")


class MockObspyClient(object):
    def __init__(*args, **kwargs) -> None:
        pass

    def get_stations(*args, **kwargs):
        if kwargs["channel"].startswith("BH"):
            data = "body"
        else:
            data = "strong"
        return get_test_inventory(data=data)

    def get_waveforms(*args, **kwargs):
        chan = args[3]
        if chan.endswith("2"):
            data = "strong"
        else:
            data = "body"
        return get_test_waveforms(data)

    def sacpz(netwk, statn, loc_code, channel, t1, t2, filename):
        with open(filename, "w"):
            pass
