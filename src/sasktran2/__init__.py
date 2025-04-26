# ruff: noqa: F401, E402
from __future__ import annotations

import os

os.environ.setdefault("OPENBLAS_NUM_THREADS", "1")

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

from ._core_rust import PyConfig as Config
from ._core_rust import PyGroundViewingSolar as GroundViewingSolar
from ._core_rust import PyViewingGeometry as ViewingGeometry
from ._core_rust import PyTangentAltitudeSolar as TangentAltitudeSolar
from ._core_rust import PySolarAnglesObserverLocation as SolarAnglesObserverLocation
from .rsk import Engine

from ._core_rust import (
    EmissionSource,
    GeometryType,
    InputValidationMode,
    InterpolationMethod,
    MultipleScatterSource,
    OccultationSource,
    SingleScatterSource,
    StokesBasis,
    ThreadingModel
)

from ._core import (
    Geodetic,
)
from ._version import __version__
from .geometry import Geometry1D
from .atmosphere import Atmosphere
from .geodetic import WGS84, SphericalGeoid
from .output import Output, OutputDerivMapped, OutputIdeal
