import numpy as np
from pathlib import Path
import sasktran2 as sk

from sasktran2.optical.base import OpticalProperty, OpticalQuantities



class OpticalDatabaseGenericAbsorber(OpticalProperty):
    def __init__(self, db_filepath: Path) -> None:
        self._file = db_filepath

        try:
            import xarray as xr
        except ImportError:
            msg = 'xarray must be installed to use OpticalDatabaseGenericAbsorber'
            raise msg from OSError

        try:
            import scipy
        except ImportError:
            msg = 'scipy must be installed to use OpticalDatabaseGenericAbsorber'
            raise msg from OSError

        self._database = xr.open_dataset(self._file)
        self._validate_db()

    def _validate_db(self):
        pass

    def atmosphere_quantities(self, atmo: sk.Atmosphere) -> OpticalQuantities:
        quants = OpticalQuantities(ssa=np.zeros_like(atmo.storage.ssa))

        coords = self._database['xs'].coords

        interp_handler = {}

        if 'temperature' in coords:
            if atmo.temperature_k is None:
                msg = 'temperature_k must be specified in Atmosphere to use OpticalDatabaseGenericAbsorber'
                raise ValueError(msg)

            interp_handler['temperature'] = atmo.temperature_k

        if 'pressure' in coords:
            if atmo.pressure_pa is None:
                msg = 'pressure_pa must be specified in Atmosphere to use OpticalDatabaseGenericAbsorber'
                raise ValueError(msg)

            interp_handler['pressure'] = atmo.pressure_pa

        if 'wavelength_nm' in coords:
            if atmo.wavelengths_nm is None:
                msg = 'wavelengths_nm must be specified in Atmosphere to use OpticalDatabaseGenericAbsorber'
            interp_handler['wavelength_nm'] = atmo.wavelengths_nm

        if 'wavenumber_cminv' in coords:
            if atmo.wavenumber_cminv is None:
                msg = 'wavenumber_cminv must be specified in Atmosphere to use OpticalDatabaseGenericAbsorber'
            interp_handler['wavenumber_cminv'] = atmo.wavenumber_cminv

        quants.extinction = self._database['xs'].interp(**interp_handler).to_numpy()

        # Out of bounds, set to 0
        quants.extinction[np.isnan(quants.extinction)] = 0

        return quants