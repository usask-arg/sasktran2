# ruff: noqa: F401, E402
from __future__ import annotations

import os

os.environ.setdefault("OPENBLAS_NUM_THREADS", "1")
os.environ.setdefault("RUST_BACKTRACE", "full")

# from ._core import Geodetic

from . import (
    appconfig,
    climatology,
    constituent,
    database,
    mie,
    optical,
    solar,
    spectroscopy,
    test_util,
    util,
    viewinggeo,
)
from ._core_rust import (
    EmissionSource,
    GeometryType,
    InputValidationMode,
    InterpolationMethod,
    MultipleScatterSource,
    OccultationSource,
    SingleScatterSource,
    StokesBasis,
    ThreadingLib,
    ThreadingModel,
)
from .atmosphere import Atmosphere
from .config import Config
from .engine import Engine
from .geodetic import WGS84, Geodetic, SphericalGeoid
from .geometry import Geometry1D
from .viewinggeo.wrappers import (
    GroundViewingSolar,
    SolarAnglesObserverLocation,
    TangentAltitudeSolar,
    ViewingGeometry,
)
