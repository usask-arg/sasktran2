from pathlib import Path
import sasktran2 as sk
import numpy as np
from sasktran2.constants import K_BOLTZMANN

from . import base
from . import database


class O3DBM(database.OpticalDatabaseGenericAbsorber):
    def __init__(self, version: str = None) -> None:
        if version is not None:
            dbm_file = sk.appconfig.database_root().joinpath('cross_sections/o3/dbm_v{}'.format(version))
        else:
            # Take the last version
            dbm_file = list(sk.appconfig.database_root().joinpath('cross_sections/o3').glob('dbm_*'))[-1]

        if dbm_file.exists():
            super().__init__(dbm_file)
        else:
            msg = 'Could not find DBM file'
            raise OSError(msg)



def pressure_temperature_to_numberdensity(pressure_pa: np.array, temperature_k: np.array):
    """_summary_

    Parameters
    ----------
    pressure_pa : np.array
        _description_
    temperature_k : np.array
        _description_
    """
    # Ideal gas law is PV = N k T
    # Or (N/V) = P / (k T)

    return pressure_pa / (K_BOLTZMANN * temperature_k)
