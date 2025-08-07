from __future__ import annotations

from sasktran2._core_rust import LineDatabaseType
from sasktran2.database.aer_line import AERLineDatabase
from sasktran2.optical.hitran import LineAbsorber


class AERLineAbsorber(LineAbsorber):
    def __init__(self, molecule: str, **kwargs):
        super().__init__(LineDatabaseType.AER, AERLineDatabase(), molecule, **kwargs)
