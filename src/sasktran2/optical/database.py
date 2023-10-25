from pathlib import Path

import numpy as np
import xarray as xr

import sasktran2 as sk
from sasktran2.atmosphere import Atmosphere, NativeGridDerivative
from sasktran2.optical.base import OpticalProperty, OpticalQuantities


class OpticalDatabase(OpticalProperty):
    def __init__(self, db_filepath: Path) -> None:
        """
        An optical property that is defined by a database file.  This is just a base class to handle file loading,
        derived classes must be used as actual optical properties.

        Parameters
        ----------
        db_filepath : Path
            Path to the optical database file

        Raises
        ------
        OSError
            If xarray is not installed
        """
        super().__init__()
        self._file = db_filepath

        try:
            import xarray as xr
        except ImportError:
            msg = "xarray must be installed to use OpticalDatabaseGenericAbsorber"
            raise msg from OSError

        self._database = xr.open_dataset(self._file)
        self._validate_db()


class OpticalDatabaseGenericAbsorber(OpticalDatabase):
    def __init__(self, db_filepath: Path) -> None:
        """
        A purely absorbing optical property defined by a database file.  The database must contain the following

        - xs : The absorption cross section in [m^2]

        xs must be a function of either wavelength_nm or wavenumber_cminv, and optionally temperature and pressure.

        Parameters
        ----------
        db_filepath : Path
            Path to the database file
        """
        super().__init__(db_filepath)

        self._validate_db()

    def _validate_db(self):
        # Old file format used "temperature" and "pressure" instead of the standard keys
        # "temperature_k" and "pressure_pa", if so we rename them

        if "temperature" in self._database:
            self._database = self._database.rename({"temperature": "temperature_k"})
        if "pressure" in self._database:
            self._database = self._database.rename({"pressure": "pressure_pa"})

        if "xs" not in self._database:
            msg = "xs must be defined in the optical database"
            raise ValueError(msg)

    def optical_derivatives(self, atmo: Atmosphere, **kwargs) -> dict:
        derivs = {}

        interp_handler = self._construct_interp_handler(atmo, **kwargs)

        # Split the interpolators into ones that are 'z' dependent and ones that are not

        interp_handler_z = {}
        interp_handler_noz = {}

        for key, val in interp_handler.items():
            if val[0] == "z":
                interp_handler_z[key] = val
            else:
                interp_handler_noz[key] = val

        # Get the derivatives of the cross section with respect to the z dependent variables
        partial_interp = self._database["xs"].interp(**interp_handler_noz)

        for key, val in interp_handler_z.items():
            # Interpolate over the other variables
            new_interpolants = {k: v for k, v in interp_handler_z.items() if k != key}

            partial_interp2 = partial_interp.interp(**new_interpolants)

            dT = partial_interp2.diff(key) / partial_interp2[key].diff(key)

            interp_index = np.argmax(dT[key].to_numpy() > val[1][:, np.newaxis], axis=1)

            if "z" in dT.dims:
                dT = dT.isel(
                    {
                        "z": xr.DataArray(list(range(len(interp_index))), dims="z"),
                        key: xr.DataArray(interp_index, dims="z"),
                    }
                )
            else:
                dT = dT.isel({key: interp_index})

            derivs[key] = NativeGridDerivative(d_extinction=dT.to_numpy())
            derivs[key].d_extinction[np.isnan(derivs[key].d_extinction)] = 0

        return derivs

    def _construct_interp_handler(self, atmo: Atmosphere, **kwargs) -> dict:
        coords = self._database["xs"].coords

        interp_handler = {}

        # TODO: this could probably be refactored to iterate over atmosphere properties?
        if "temperature_k" in coords:
            if atmo.temperature_k is None:
                msg = "temperature_k must be specified in Atmosphere to use OpticalDatabaseGenericAbsorber"
                raise ValueError(msg)

            interp_handler["temperature_k"] = ("z", atmo.temperature_k)

        if "pressure_pa" in coords:
            if atmo.pressure_pa is None:
                msg = "pressure_pa must be specified in Atmosphere to use OpticalDatabaseGenericAbsorber"
                raise ValueError(msg)

            interp_handler["pressure_pa"] = ("z", atmo.pressure_pa)

        if "wavelength_nm" in coords:
            if atmo.wavelengths_nm is None:
                msg = "wavelengths_nm must be specified in Atmosphere to use OpticalDatabaseGenericAbsorber"
            interp_handler["wavelength_nm"] = atmo.wavelengths_nm

        if "wavenumber_cminv" in coords:
            if atmo.wavenumber_cminv is None:
                msg = "wavenumber_cminv must be specified in Atmosphere to use OpticalDatabaseGenericAbsorber"
            interp_handler["wavenumber_cminv"] = atmo.wavenumber_cminv

        return interp_handler

    def atmosphere_quantities(self, atmo: sk.Atmosphere, **kwargs) -> OpticalQuantities:
        quants = OpticalQuantities(ssa=np.zeros_like(atmo.storage.ssa))

        interp_handler = self._construct_interp_handler(atmo, **kwargs)

        quants.extinction = self._database["xs"].interp(**interp_handler).to_numpy()

        # Out of bounds, set to 0
        quants.extinction[np.isnan(quants.extinction)] = 0

        return quants


class OpticalDatabaseGenericScatterer(OpticalDatabase):
    def __init__(self, db_filepath: Path) -> None:
        """
        A purely scattering optical property defined by a database file.  The database must contain the following

        - xs_total : The total cross section in [m^2]
        - xs_scattering : The scattering cross section in [m^2]
        - lm_a1 : the legendre coefficients for the phase function

        All variables must be a function of either wavelength_nm or wavenumber_cminv, and optionally any other dimension such
        as particle size.

        Parameters
        ----------
        db_filepath : Path
            Path to the database file
        """
        super().__init__(db_filepath)

        self._validate_db()

    def _validate_db(self):
        self._database["lm_a1"] /= self._database["lm_a1"].isel(legendre=0)

    def _construct_interp_handler(self, atmo: Atmosphere, **kwargs) -> dict:
        coords = self._database["xs_total"].coords
        interp_handler = {}

        if "wavelength_nm" in coords:
            if atmo.wavelengths_nm is None:
                msg = "wavelengths_nm must be specified in Atmosphere to use OpticalDatabaseGenericScatterer"
                raise ValueError(msg)
            interp_handler["wavelength_nm"] = atmo.wavelengths_nm

        if "wavenumber_cminv" in coords:
            if atmo.wavenumber_cminv is None:
                msg = "wavenumber_cminv must be specified in Atmosphere to use OpticalDatabaseGenericScatterer"
                raise ValueError(msg)
            interp_handler["wavenumber_cminv"] = atmo.wavenumber_cminv

        for name, vals in kwargs.items():
            if name in coords:
                interp_handler[name] = ("z", vals)

        return interp_handler

    def cross_sections(
        self, wavelengths_nm: np.array, altitudes_m: np.array, **kwargs  # noqa: ARG002
    ) -> OpticalQuantities:
        quants = OpticalQuantities()

        coords = self._database["xs_total"].coords
        interp_handler = {}

        interp_handler["wavelength_nm"] = wavelengths_nm

        for name, vals in kwargs.items():
            if name in coords:
                interp_handler[name] = ("z", vals)

        ds_interp = self._database.interp(**interp_handler)

        quants.extinction = ds_interp["xs_total"].to_numpy()
        quants.ssa = ds_interp["xs_scattering"].to_numpy() / quants.extinction

        quants.extinction[np.isnan(quants.extinction)] = 0
        quants.ssa[np.isnan(quants.ssa)] = 0

        return quants

    def atmosphere_quantities(self, atmo: Atmosphere, **kwargs) -> OpticalQuantities:
        quants = OpticalQuantities()

        interp_handler = self._construct_interp_handler(atmo, **kwargs)

        ds_interp = self._database.interp(**interp_handler)

        quants.extinction = ds_interp["xs_total"].to_numpy()
        quants.ssa = ds_interp["xs_scattering"].to_numpy()

        num_assign_legendre = min(
            atmo.storage.leg_coeff.shape[0], len(ds_interp["legendre"])
        )

        quants.leg_coeff = (
            (ds_interp["lm_a1"])
            .transpose("legendre", "z", "wavelength_nm")
            .to_numpy()[:num_assign_legendre]
        )

        # Renormalize leg_coeffs so a1 is always 1
        quants.leg_coeff /= quants.leg_coeff[0, :, :][np.newaxis, :, :]

        quants.extinction[np.isnan(quants.extinction)] = 0
        quants.ssa[np.isnan(quants.ssa)] = 0
        quants.leg_coeff[np.isnan(quants.leg_coeff)] = 0

        if atmo.nstokes > 1:
            # TODO: add pol properties
            pass

        return quants

    def optical_derivatives(self, atmo: Atmosphere, **kwargs) -> dict:
        derivs = {}

        interp_handler = self._construct_interp_handler(atmo, **kwargs)

        # Split the interpolators into ones that are 'z' dependent and ones that are not

        interp_handler_z = {}
        interp_handler_noz = {}

        for key, val in interp_handler.items():
            if val[0] == "z":
                interp_handler_z[key] = val
            else:
                interp_handler_noz[key] = val

        # Get the derivatives of the cross section with respect to the z dependent variables
        partial_interp = self._database.interp(**interp_handler_noz)

        num_assign_legendre = min(
            atmo.storage.leg_coeff.shape[0], len(partial_interp["legendre"])
        )

        for key, val in interp_handler_z.items():
            # Interpolate over the other variables
            new_interpolants = {k: v for k, v in interp_handler_z.items() if k != key}

            partial_interp2 = partial_interp.interp(**new_interpolants)

            dT = partial_interp2.diff(key) / partial_interp2[key].diff(key)

            interp_index = np.argmax(dT[key].to_numpy() > val[1][:, np.newaxis], axis=1)

            if "z" in dT.dims:
                dT = dT.isel(
                    {
                        "z": xr.DataArray(list(range(len(interp_index))), dims="z"),
                        key: xr.DataArray(interp_index, dims="z"),
                    }
                )
            else:
                dT = dT.isel({key: xr.DataArray(interp_index, dims="z")})

            derivs[key] = NativeGridDerivative(
                d_extinction=dT["xs_total"].to_numpy(),
                d_ssa=dT["xs_scattering"].to_numpy(),
                d_leg_coeff=dT["lm_a1"]
                .transpose("legendre", "z", "wavelength_nm")
                .to_numpy()[:num_assign_legendre],
            )

            derivs[key].d_extinction[np.isnan(derivs[key].d_extinction)] = 0
            derivs[key].d_ssa[np.isnan(derivs[key].d_ssa)] = 0
            derivs[key].d_leg_coeff[np.isnan(derivs[key].d_leg_coeff)] = 0

        return derivs
