import numpy as np
from pydantic import BaseModel, ConfigDict


class SegmentData(BaseModel):
    """
    Class recording teh segment data in the model solution
    """

    model_config = ConfigDict(arbitrary_types_allowed=True)
    depth: np.ndarray
    dip: np.ndarray
    lat: np.ndarray
    lon: np.ndarray
    moment: np.ndarray
    rake: np.ndarray
    slip: np.ndarray
    strike: np.ndarray
    tfall: np.ndarray
    trise: np.ndarray
    trup: np.ndarray
