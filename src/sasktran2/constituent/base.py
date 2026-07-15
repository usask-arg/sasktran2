from __future__ import annotations

import abc
from typing import Literal

from sasktran2.atmosphere import Atmosphere


class Constituent(abc.ABC):
    @property
    def volume_spatial_mode(self) -> Literal["altitude_profile", "native_2d"]:
        """Spatial input mode used when adding the constituent to an atmosphere.

        Existing constituents are altitude profiles and are broadcast across the
        horizontal dimension of a 2D atmosphere. New constituents that consume
        native 2D fields can override this with ``"native_2d"``.
        """
        return "altitude_profile"

    @abc.abstractmethod
    def add_to_atmosphere(self, atmo: Atmosphere):
        pass

    @abc.abstractmethod
    def register_derivative(self, atmo: Atmosphere, name: str):
        pass
