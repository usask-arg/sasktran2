import numpy as np

from sasktran2.constants import K_BOLTZMANN


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
