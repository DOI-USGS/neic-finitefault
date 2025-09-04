from enum import Enum
from typing import List

DEFAULT_MANAGEMENT_FILES = {
    "cgnss": "cgnss_waves.json",
    "gnss": "static_data.json",
    "insar": "insar_data.json",
    "strong": "strong_motion_waves.json",
    "surf": "surf_waves.json",
    "body": "tele_waves.json",
}


class AcquireDataTypes(str, Enum):
    strong = "strong"
    body = "body"


class ManagedDataTypes(str, Enum):
    cgnss = "cgnss"
    gnss = "gnss"
    insar = "insar"
    strong = "strong"
    surf = "surf"
    body = "body"


class ModifiableDataTypes(str, Enum):
    cgnss = "cgnss"
    gnss = "gnss"
    strong_motion = "strong"
    surf = "surf"
    body = "body"


ProcessDataTypes = ManagedDataTypes


class ShiftMatchDataTypes(str, Enum):
    cgnss = "cgnss"
    strong = "strong"
    surf = "surf"
    body = "body"
