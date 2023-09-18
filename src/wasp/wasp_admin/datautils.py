from enum import Enum
from typing import List

DEFAULT_MANAGEMENT_FILES = {
    "cgps": "cgps_waves.json",
    "gps": "gps_waves.json",
    "strong": "strong_motion_waves.json",
    "surf": "surf_waves.json",
    "tele": "tele_waves.json",
}


class AcquireDataTypes(str, Enum):
    strong = "strong"
    tele = "tele"


class ManagedDataTypes(str, Enum):
    cgps = "cgps"
    gps = "gps"
    insar = "insar"
    strong_motion = "strong"
    surf_tele = "surf"
    tele_body = "tele"


class ModifiableDataTypes(str, Enum):
    cgps = "cgps"
    gps = "gps"
    strong_motion = "strong"
    surf_tele = "surf"
    tele_body = "tele"


ProcessDataTypes = ManagedDataTypes


class ShiftMatchDataTypes(str, Enum):
    cgps = "cgps"
    strong = "strong"
    surf = "surf"
    tele = "tele"
