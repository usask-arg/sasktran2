from __future__ import annotations

import numpy as np

from sasktran2._core_rust import PyMie, PyMieOutput


class MieValueAccessor:
    def __init__(self, owner):
        self._owner = owner

    def __getattr__(self, name):
        # Called when accessing "a.value.X"
        return getattr(self._owner, name)


class MieOutput:
    _mie_output: PyMieOutput

    def __init__(self, mie_output: PyMieOutput):
        self._mie_output = mie_output

    @property
    def values(self) -> MieValueAccessor:
        return MieValueAccessor(self)

    @property
    def cos_angles(self) -> np.ndarray:
        return self._mie_output.cos_angles

    @property
    def size_parameter(self) -> np.ndarray:
        return self._mie_output.size_param

    @property
    def Qsca(self) -> np.ndarray:
        return self._mie_output.Qsca

    @property
    def Qext(self) -> np.ndarray:
        return self._mie_output.Qext

    @property
    def S1(self) -> np.ndarray:
        return self._mie_output.S1

    @property
    def S2(self) -> np.ndarray:
        return self._mie_output.S2


class LinearizedMie:
    _mie: PyMie

    def __init__(self, num_threads: int = 1):  # noqa: ARG002
        self._mie = PyMie()

    def calculate(
        self,
        size_param: np.ndarray,
        refractive_index,
        cos_angles,
        calculate_derivatives=False,
    ):
        res = self._mie.calculate(
            np.atleast_1d(size_param).astype(np.float64),
            complex(refractive_index),
            np.atleast_1d(cos_angles).astype(np.float64),
            calculate_derivatives,
        )

        return MieOutput(res)
