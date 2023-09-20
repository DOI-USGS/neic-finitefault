from enum import Enum
from typing import List

DEFAULT_MANAGEMENT_FILES = {
    "cgps": "cgps_waves.json",
    "strong": "strong_motion_waves.json",
    "surf": "surf_waves.json",
    "body": "tele_waves.json",
}


class AcquireDataTypes(str, Enum):
    strong = "strong"
    body = "body"


class ManagedDataTypes(str, Enum):
    cgps = "cgps"
    gps = "gps"
    insar = "insar"
    strong = "strong"
    surf = "surf"
    body = "body"


class ModifiableDataTypes(str, Enum):
    cgps = "cgps"
    gps = "gps"
    strong_motion = "strong"
    surf = "surf"
    body = "body"


ProcessDataTypes = ManagedDataTypes


class ShiftMatchDataTypes(str, Enum):
    cgps = "cgps"
    strong = "strong"
    surf = "surf"
    body = "body"
