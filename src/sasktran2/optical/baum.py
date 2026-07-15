from __future__ import annotations

import warnings
from pathlib import Path
from typing import Literal

import numpy as np
import xarray as xr

from sasktran2.atmosphere import Atmosphere
from sasktran2.database.web import StandardDatabase
from sasktran2.optical.base import OpticalQuantities
from sasktran2.optical.database import OpticalDatabaseGenericScattererRust

BaumParticleModel = Literal[
    "general_habit_mixture",
    "aggregate_solid_columns",
    "solid_columns",
]

_PARTICLE_MODELS = (
    "general_habit_mixture",
    "aggregate_solid_columns",
    "solid_columns",
)
_DEFAULT_DATABASE_MOMENTS = 256
_STANDARD_DATABASE_KEY = "cross_sections/ice/baum_ice_crystals_v3_6.nc"
_FULL_STANDARD_DATABASE_KEY = "cross_sections/ice/baum_ice_crystals_v3_6_full.nc"
_RUNTIME_VARIABLES = (
    "xs_total",
    "xs_scattering",
    "lm_a1",
    "lm_a2",
    "lm_a3",
    "lm_a4",
    "lm_b1",
    "lm_b2",
)


class BaumIceCrystal(OpticalDatabaseGenericScattererRust):
    """Baum V3.6 severely rough ice-crystal optical properties.

    This database provides extinction, single-scattering albedo, and the stored
    polarized phase-matrix expansion for three ice-crystal habit models. The
    tabulation contains 445 wavelengths from 199 nm to 99,000 nm (0.199 to
    99 microns) and 23 effective diameters from 10 to 120 microns in 5-micron
    increments. Wavelengths supplied to SASKTRAN2 remain in nanometres.

    The available ``particle_model`` values are:

    - ``"general_habit_mixture"``: the Baum general habit mixture.
    - ``"aggregate_solid_columns"``: aggregates of solid columns.
    - ``"solid_columns"``: solid columns.

    All three models use the severely rough particle treatment in the Baum V3.6
    tables. The large database is capped at 16,384 moments. Some sharply peaked
    short-wavelength phase matrices do not meet the database-generation
    reconstruction tolerances at that cap; the NetCDF generation diagnostics mark
    those cells explicitly.

    This optical property is intended for use with
    :class:`sasktran2.constituent.ExtinctionScatterer`. The effective diameter is
    selected by passing ``effective_diameter_um`` to the constituent; supply one
    value per point on the constituent altitude grid. Values between tabulated
    diameters are interpolated by the Rust database backend. The stored extinction
    cross section is extinction per ice-water content in ``m^2 g^-1``, so an
    extinction-profile constituent is preferred over treating the table as a
    per-particle number-density cross section.

    Parameters
    ----------
    particle_model : {"general_habit_mixture", "aggregate_solid_columns", "solid_columns"}, optional
        Ice-crystal habit model to select. By default ``"general_habit_mixture"``.
    max_moments : int or None, optional
        Maximum number of Greek-coefficient moments to load. Limits up to 256 use
        the smaller standard-database artifact. Values above 256 select the full
        database, which is several gigabytes. ``None`` selects the full database and
        loads all 16,384 stored moments, which can also require multiple gigabytes of
        memory. This value must be at least
        :attr:`sasktran2.Config.num_singlescatter_moments` for any atmosphere using
        the property. By default 256.
    db_filepath : pathlib.Path or None, optional
        Optional local database override. By default the file is obtained from the
        SASKTRAN2 standard database and cached locally.

    Raises
    ------
    TypeError
        If ``max_moments`` is not an integer or ``None``.
    ValueError
        If the particle model or moment limit is invalid, or when a requested
        effective diameter is outside 10--120 microns.
    OSError
        If the database cannot be resolved or opened.

    Examples
    --------
    Construct an extinction-profile ice cloud with a 40-micron effective diameter:

    >>> import numpy as np
    >>> import sasktran2 as sk
    >>> altitudes_m = np.arange(0.0, 20001.0, 1000.0)
    >>> extinction_per_m = np.full_like(altitudes_m, 1.0e-5)
    >>> effective_diameter_um = np.full_like(altitudes_m, 40.0)
    >>> ice_optical_property = sk.optical.BaumIceCrystal(
    ...     particle_model="general_habit_mixture",
    ...     max_moments=256,
    ... )
    >>> ice_cloud = sk.constituent.ExtinctionScatterer(
    ...     ice_optical_property,
    ...     altitudes_m,
    ...     extinction_per_m,
    ...     extinction_wavelength_nm=550.0,
    ...     effective_diameter_um=effective_diameter_um,
    ... )
    """

    def __init__(
        self,
        particle_model: BaumParticleModel = "general_habit_mixture",
        max_moments: int | None = 256,
        db_filepath: Path | None = None,
    ) -> None:
        if particle_model not in _PARTICLE_MODELS:
            valid = ", ".join(_PARTICLE_MODELS)
            msg = f"Unknown Baum particle model {particle_model!r}; expected one of {valid}"
            raise ValueError(msg)

        if max_moments is not None:
            if isinstance(max_moments, (bool, np.bool_)) or not isinstance(
                max_moments, (int, np.integer)
            ):
                msg = "max_moments must be an integer or None"
                raise TypeError(msg)
            if max_moments <= 0:
                msg = "max_moments must be positive or None"
                raise ValueError(msg)
            max_moments = int(max_moments)

        if db_filepath is None:
            use_full_database = (
                max_moments is None or max_moments > _DEFAULT_DATABASE_MOMENTS
            )
            database_key = (
                _FULL_STANDARD_DATABASE_KEY
                if use_full_database
                else _STANDARD_DATABASE_KEY
            )
            if use_full_database:
                warnings.warn(
                    "The requested max_moments selects the full Baum ice-crystal "
                    "database, whose download is several gigabytes. Runtime memory "
                    "use scales with max_moments and can also reach several gigabytes.",
                    UserWarning,
                    stacklevel=2,
                )
            db_filepath = StandardDatabase().path(database_key)
        if db_filepath is None:
            msg = "Could not resolve the Baum ice-crystal database"
            raise OSError(msg)

        db_filepath = Path(db_filepath)
        if not db_filepath.exists():
            msg = f"Baum ice-crystal database does not exist: {db_filepath}"
            raise OSError(msg)

        with xr.open_dataset(db_filepath) as source:
            missing = [name for name in _RUNTIME_VARIABLES if name not in source]
            if missing:
                msg = f"Baum database is missing required variables: {missing}"
                raise ValueError(msg)
            if "particle_model" not in source.coords:
                msg = "Baum database is missing the particle_model coordinate"
                raise ValueError(msg)
            available_models = tuple(
                str(value) for value in source.particle_model.values
            )
            if particle_model not in available_models:
                msg = (
                    f"Particle model {particle_model!r} is not present in {db_filepath}; "
                    f"available models are {available_models}"
                )
                raise ValueError(msg)

            available_moments = source.sizes.get("legendre", 0)
            if available_moments == 0:
                msg = "Baum database contains no Legendre moments"
                raise ValueError(msg)
            if max_moments is None:
                loaded_moments = available_moments
            else:
                if max_moments > available_moments:
                    msg = (
                        f"Requested {max_moments} moments, but the Baum database only "
                        f"contains {available_moments}"
                    )
                    raise ValueError(msg)
                loaded_moments = int(max_moments)

            selected = (
                source[list(_RUNTIME_VARIABLES)]
                .sel(particle_model=particle_model, drop=True)
                .isel(legendre=slice(0, loaded_moments))
                .load()
            )
            diameter = np.asarray(selected["effective_diameter_um"], dtype=np.float64)

        self._particle_model = particle_model
        self._available_moments = available_moments
        self._max_moments = loaded_moments
        self._diameter_range_um = (float(diameter.min()), float(diameter.max()))
        super().__init__(db=selected)
        # The Rust interpolator owns its packed coefficient array after construction;
        # retaining the selected xarray data would duplicate the largest allocation.
        self._database = None

    @property
    def particle_model(self) -> BaumParticleModel:
        """Selected ice-crystal habit model."""
        return self._particle_model

    @property
    def max_moments(self) -> int:
        """Number of Greek-coefficient moments loaded into the Rust backend."""
        return self._max_moments

    @property
    def available_moments(self) -> int:
        """Moments stored in the selected standard or local database file."""
        return self._available_moments

    def _validate_effective_diameter(self, kwargs: dict) -> None:
        if "effective_diameter_um" not in kwargs:
            msg = "effective_diameter_um must be supplied for BaumIceCrystal"
            raise ValueError(msg)
        diameter = np.asarray(kwargs["effective_diameter_um"], dtype=np.float64)
        lower, upper = self._diameter_range_um
        if not np.all(np.isfinite(diameter)):
            msg = "effective_diameter_um must contain only finite values"
            raise ValueError(msg)
        if np.any((diameter < lower) | (diameter > upper)):
            msg = f"effective_diameter_um must be within [{lower}, {upper}] microns"
            raise ValueError(msg)

    def _validate_atmosphere_moments(self, atmo: Atmosphere) -> None:
        requested = atmo.leg_coeff.a1.shape[0]
        if requested > self._max_moments:
            msg = (
                f"Atmosphere requests {requested} single-scatter moments, but this "
                f"BaumIceCrystal loaded only {self._max_moments}; increase max_moments"
            )
            raise ValueError(msg)

    def cross_sections(
        self, wavelengths_nm: np.ndarray, altitudes_m: np.ndarray, **kwargs
    ) -> OpticalQuantities:
        self._validate_effective_diameter(kwargs)
        return super().cross_sections(wavelengths_nm, altitudes_m, **kwargs)

    def atmosphere_quantities(self, atmo: Atmosphere, **kwargs) -> OpticalQuantities:
        self._validate_effective_diameter(kwargs)
        self._validate_atmosphere_moments(atmo)
        return super().atmosphere_quantities(atmo, **kwargs)

    def optical_derivatives(self, atmo: Atmosphere, **kwargs) -> dict:
        self._validate_effective_diameter(kwargs)
        self._validate_atmosphere_moments(atmo)
        return super().optical_derivatives(atmo, **kwargs)

    def cross_section_derivatives(
        self, wavelengths_nm: np.ndarray, altitudes_m: np.ndarray, **kwargs
    ) -> dict:
        self._validate_effective_diameter(kwargs)
        return super().cross_section_derivatives(wavelengths_nm, altitudes_m, **kwargs)

    def _into_rust_object(self):
        # Keep the thin Python validation layer so a moment mismatch is reported
        # before the Rust interpolator fills an atmosphere-sized output buffer.
        return None
