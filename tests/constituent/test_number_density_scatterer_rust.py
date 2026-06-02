from __future__ import annotations

import numpy as np
import pytest
import sasktran2 as sk
from sasktran2.test_util.reference_numdenscatterer import (
    ExtinctionScatterer as ReferenceExtinctionScatterer,
)
from sasktran2.test_util.reference_numdenscatterer import (
    NumberDensityScatterer as ReferenceNumberDensityScatterer,
)


def _atmosphere():
    config = sk.Config()
    config.num_streams = 4
    config.delta_m_scaling = False

    altitude_grid = np.arange(0, 65001, 1000.0)
    geometry = sk.Geometry1D(
        0.6,
        0,
        6327000,
        altitude_grid,
        sk.InterpolationMethod.LinearInterpolation,
        sk.GeometryType.Spherical,
    )

    atmosphere = sk.Atmosphere(
        geometry, config, wavelengths_nm=np.array([310, 525, 750])
    )
    sk.climatology.us76.add_us76_standard_atmosphere(atmosphere)
    atmosphere["rayleigh"] = sk.constituent.Rayleigh()

    return atmosphere


def _scatterer_inputs():
    altitudes_m = np.array([0, 10000, 30000, 70000], dtype=float)
    number_density = np.array([1.0e7, 2.0e7, 8.0e6, 2.0e6])
    radius = np.array([80.0, 105.0, 140.0, 180.0])

    return altitudes_m, number_density, radius


def _optical_property(cls):
    return cls(
        sk.database.StandardDatabase().path("cross_sections/mie/sulfate_test.nc")
    )


def _run_atmosphere(constituent):
    atmosphere = _atmosphere()
    atmosphere["aerosol"] = constituent
    atmosphere.internal_object()
    return atmosphere


def _assert_storage_parity(rust_atmo, reference_atmo):
    np.testing.assert_allclose(
        rust_atmo.storage.total_extinction,
        reference_atmo.storage.total_extinction,
        rtol=1e-14,
        atol=0,
    )
    np.testing.assert_allclose(
        rust_atmo.storage.ssa, reference_atmo.storage.ssa, rtol=1e-14, atol=0
    )
    np.testing.assert_allclose(
        rust_atmo.storage.leg_coeff,
        reference_atmo.storage.leg_coeff,
        rtol=1e-14,
        atol=0,
    )


def _assert_mapping_parity(rust_atmo, reference_atmo, name):
    rust_mapping = rust_atmo.storage.get_derivative_mapping(name)
    reference_mapping = reference_atmo.storage.get_derivative_mapping(name)

    for field in [
        "d_extinction",
        "d_ssa",
        "d_leg_coeff",
        "scat_factor",
        "interpolator",
    ]:
        np.testing.assert_allclose(
            getattr(rust_mapping, field),
            getattr(reference_mapping, field),
            rtol=1e-12,
            atol=1e-15,
        )

    assert rust_mapping.interp_dim == reference_mapping.interp_dim


@pytest.mark.parametrize(
    "optical_cls",
    [
        sk.optical.database.OpticalDatabaseGenericScattererRust,
        sk.optical.database.OpticalDatabaseGenericScatterer,
    ],
)
def test_number_density_scatterer_matches_python_reference(optical_cls):
    altitudes_m, number_density, radius = _scatterer_inputs()

    rust_atmo = _run_atmosphere(
        sk.constituent.NumberDensityScatterer(
            _optical_property(optical_cls),
            altitudes_m,
            number_density,
            lognormal_median_radius=radius,
        )
    )
    reference_atmo = _run_atmosphere(
        ReferenceNumberDensityScatterer(
            _optical_property(optical_cls),
            altitudes_m,
            number_density,
            lognormal_median_radius=radius,
        )
    )

    _assert_storage_parity(rust_atmo, reference_atmo)
    _assert_mapping_parity(rust_atmo, reference_atmo, "wf_aerosol_number_density")
    _assert_mapping_parity(
        rust_atmo, reference_atmo, "wf_aerosol_lognormal_median_radius"
    )


@pytest.mark.parametrize(
    "optical_cls",
    [
        sk.optical.database.OpticalDatabaseGenericScattererRust,
        sk.optical.database.OpticalDatabaseGenericScatterer,
    ],
)
def test_extinction_scatterer_matches_python_reference(optical_cls):
    altitudes_m, _, radius = _scatterer_inputs()
    extinction_per_m = np.array([2.0e-7, 3.5e-7, 1.2e-7, 4.0e-8])
    extinction_wavelength_nm = 525.0

    rust_atmo = _run_atmosphere(
        sk.constituent.ExtinctionScatterer(
            _optical_property(optical_cls),
            altitudes_m,
            extinction_per_m,
            extinction_wavelength_nm,
            lognormal_median_radius=radius,
        )
    )
    reference_atmo = _run_atmosphere(
        ReferenceExtinctionScatterer(
            _optical_property(optical_cls),
            altitudes_m,
            extinction_per_m,
            extinction_wavelength_nm,
            lognormal_median_radius=radius,
        )
    )

    _assert_storage_parity(rust_atmo, reference_atmo)
    _assert_mapping_parity(rust_atmo, reference_atmo, "wf_aerosol_extinction")
    _assert_mapping_parity(
        rust_atmo, reference_atmo, "wf_aerosol_lognormal_median_radius"
    )
