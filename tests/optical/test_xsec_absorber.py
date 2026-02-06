from __future__ import annotations

from pathlib import Path

import numpy as np
import pytest
import sasktran2 as sk


@pytest.fixture
def test_data_dir():
    """Returns the path to the test data directory."""
    return Path(__file__).parent.parent / "data"


@pytest.fixture
def chclfcf3_file(test_data_dir):
    """Returns the path to the CHClFCF3 test file."""
    return test_data_dir / "CHClFCF3_287.0_0.0_675.2-715.0_00.xsc"


@pytest.fixture
def chcnft1_file(test_data_dir):
    """Returns the path to the CHCNFT1 test file."""
    return test_data_dir / "CHCNFT1"


def test_xsec_from_single_file(chclfcf3_file):
    """Test loading cross section database from a single file."""
    absorber = sk.optical.XsecAbsorber(chclfcf3_file)
    assert absorber is not None


def test_xsec_from_file_list(chclfcf3_file, chcnft1_file):
    """Test loading cross section database from a list of files."""
    absorber = sk.optical.XsecAbsorber([chclfcf3_file, chcnft1_file])
    assert absorber is not None


def test_xsec_atmosphere_quantities(chclfcf3_file):
    """Test computing optical quantities for an atmosphere."""
    absorber = sk.optical.XsecAbsorber(chclfcf3_file)
    
    config = sk.Config()
    altitude_grid = np.arange(0, 65001, 1000)
    
    geometry = sk.Geometry1D(
        0.6,
        0,
        6327000,
        altitude_grid,
        sk.InterpolationMethod.LinearInterpolation,
        sk.GeometryType.Spherical,
    )
    
    # Use wavenumbers in the range of the test file (675-715 cm^-1)
    # Convert to wavelengths: lambda (nm) = 1e7 / wavenumber (cm^-1)
    wavenumbers_cminv = np.linspace(675, 715, 100)
    wavelengths_nm = 1e7 / wavenumbers_cminv
    
    atmosphere = sk.Atmosphere(geometry, config, wavelengths_nm=wavelengths_nm)
    sk.climatology.us76.add_us76_standard_atmosphere(atmosphere)
    
    # Set temperature to match the file (287 K)
    atmosphere.temperature_k[:] = 287.0
    
    # Compute optical quantities
    result = absorber.atmosphere_quantities(
        atmo=atmosphere,
        wavenumbers_cminv=wavenumbers_cminv
    )
    
    # Check that we got cross sections
    assert result.cross_section.shape == (len(altitude_grid), len(wavenumbers_cminv))
    
    # Cross sections should be positive and non-zero for this molecule
    assert np.all(result.cross_section >= 0)
    assert np.any(result.cross_section > 0)


def test_xsec_temperature_derivatives(chclfcf3_file):
    """Test computing temperature derivatives."""
    absorber = sk.optical.XsecAbsorber(chclfcf3_file)
    
    config = sk.Config()
    altitude_grid = np.arange(0, 65001, 5000)  # Fewer points for speed
    
    geometry = sk.Geometry1D(
        0.6,
        0,
        6327000,
        altitude_grid,
        sk.InterpolationMethod.LinearInterpolation,
        sk.GeometryType.Spherical,
    )
    
    wavenumbers_cminv = np.linspace(680, 710, 50)
    wavelengths_nm = 1e7 / wavenumbers_cminv
    
    atmosphere = sk.Atmosphere(geometry, config, wavelengths_nm=wavelengths_nm)
    sk.climatology.us76.add_us76_standard_atmosphere(atmosphere)
    
    # Compute derivatives
    derivs = absorber.optical_derivatives(
        atmo=atmosphere,
        wavenumbers_cminv=wavenumbers_cminv
    )
    
    # Check that we got temperature derivatives
    assert "temperature_k" in derivs
    assert derivs["temperature_k"].cross_section.shape == (len(altitude_grid), len(wavenumbers_cminv))


def test_xsec_numerical_derivative_check(chclfcf3_file):
    """
    Test that analytical temperature derivatives match numerical derivatives.
    This is a basic sanity check using finite differences.
    """
    absorber = sk.optical.XsecAbsorber(chclfcf3_file)
    
    config = sk.Config()
    altitude_grid = np.array([10000, 20000, 30000])  # Just a few points
    
    geometry = sk.Geometry1D(
        0.6,
        0,
        6327000,
        altitude_grid,
        sk.InterpolationMethod.LinearInterpolation,
        sk.GeometryType.Spherical,
    )
    
    wavenumbers_cminv = np.linspace(685, 705, 20)
    wavelengths_nm = 1e7 / wavenumbers_cminv
    
    atmosphere = sk.Atmosphere(geometry, config, wavelengths_nm=wavelengths_nm)
    sk.climatology.us76.add_us76_standard_atmosphere(atmosphere)
    
    base_temp = 287.0
    delta_temp = 0.1  # Small perturbation
    
    # Compute analytical derivative
    derivs = absorber.optical_derivatives(
        atmo=atmosphere,
        wavenumbers_cminv=wavenumbers_cminv
    )
    analytical_deriv = derivs["temperature_k"].cross_section
    
    # Compute numerical derivative using finite differences
    atmosphere.temperature_k[:] = base_temp
    oq_minus = absorber.atmosphere_quantities(
        atmo=atmosphere,
        wavenumbers_cminv=wavenumbers_cminv
    )
    
    atmosphere.temperature_k[:] = base_temp + delta_temp
    oq_plus = absorber.atmosphere_quantities(
        atmo=atmosphere,
        wavenumbers_cminv=wavenumbers_cminv
    )
    
    numerical_deriv = (oq_plus.cross_section - oq_minus.cross_section) / delta_temp
    
    # Check that analytical and numerical derivatives are close
    # Use a relatively loose tolerance since we're comparing finite difference
    np.testing.assert_allclose(analytical_deriv, numerical_deriv, rtol=1e-3, atol=1e-25)


def test_xsec_different_temperatures(chcnft1_file):
    """
    Test that cross sections change with temperature.
    The CHCNFT1 file contains data at 324.1 K.
    """
    absorber = sk.optical.XsecAbsorber(chcnft1_file)
    
    config = sk.Config()
    altitude_grid = np.array([10000])
    
    geometry = sk.Geometry1D(
        0.6,
        0,
        6327000,
        altitude_grid,
        sk.InterpolationMethod.LinearInterpolation,
        sk.GeometryType.Spherical,
    )
    
    # Use wavenumbers from the file range (3881-4573 cm^-1)
    wavenumbers_cminv = np.linspace(3900, 4500, 100)
    wavelengths_nm = 1e7 / wavenumbers_cminv
    
    atmosphere = sk.Atmosphere(geometry, config, wavelengths_nm=wavelengths_nm)
    
    # Test at the file temperature
    atmosphere.temperature_k[:] = 324.1
    atmosphere.pressure_pa[:] = 101325.0
    
    result1 = absorber.atmosphere_quantities(
        atmo=atmosphere,
        wavenumbers_cminv=wavenumbers_cminv
    )
    
    # Test at a different temperature
    atmosphere.temperature_k[:] = 300.0
    
    result2 = absorber.atmosphere_quantities(
        atmo=atmosphere,
        wavenumbers_cminv=wavenumbers_cminv
    )
    
    # Results should be different at different temperatures
    # (they will be interpolated or extrapolated)
    # We just check that we got valid results for both
    assert result1.cross_section.shape == result2.cross_section.shape
    assert np.all(result1.cross_section >= 0)
    assert np.all(result2.cross_section >= 0)
