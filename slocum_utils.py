#
# Shared utility functions for the Slocum Data Harvester
#
# June-2023, Pat Welch, pat@mousebrains.com

import math
import numpy as np


def mkDegrees_scalar(degmin: float) -> float:
    """Convert DDMM.MM format to decimal degrees (scalar version)."""
    sgn = -1 if degmin < 0 else +1
    degmin = abs(degmin)
    deg = math.floor(degmin / 100)
    minutes = degmin % 100
    return sgn * (deg + minutes / 60)


def mkDegrees(degmin: np.ndarray) -> np.ndarray:
    """Convert DDMM.MM format to decimal degrees (vectorized version)."""
    sgn = np.where(degmin < 0, -1.0, 1.0)
    degmin = np.abs(degmin)
    deg = np.floor(degmin / 100)
    minutes = np.mod(degmin, 100)
    deg = sgn * (deg + minutes / 60)
    deg[np.abs(deg) > 180] = np.nan
    return deg
