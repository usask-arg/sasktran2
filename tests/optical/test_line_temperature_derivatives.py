from __future__ import annotations

import numpy as np
import pytest
import sasktran2 as sk
import sasktran2.atmosphere as atmosphere_module
from sasktran2._core_rust import LineDatabaseType, PyLineAbsorber
from sasktran2.constituent.base import Constituent
from sasktran2.optical.hitran import LineAbsorber


def _absorber(tmp_path, partition_function=None):
    """A local synthetic O2 line list, without downloads or a HAPI dependency."""
    records = []
    for center, energy in [(13100.0, 30.0), (13100.8, 450.0)]:
        record = (
            f"{7:2d}{1:1d}{center:12.6f}{1e-24:10.3e}{0.1:10.3e}"
            f"{0.06:5.3f}{0.10:5.3f}{energy:10.4f}{0.7:4.2f}{0.003:8.6f}"
            + " " * 79
            + f"{1.0:7.1f}{1.0:7.1f}"
        )
        assert len(record) == 160
        records.append(record)
    (tmp_path / "O2.data").write_text("\n".join(records) + "\n")

    partition_temperatures = []

    def partition(_mol, _iso, temperature):
        partition_temperatures.append(temperature)
        if partition_function is not None:
            return partition_function(temperature)
        return temperature**1.5

    absorber = LineAbsorber.__new__(LineAbsorber)
    absorber._internal = PyLineAbsorber(
        LineDatabaseType.HITRAN,
        "O2",
        str(tmp_path),
        py_tips=partition,
        py_molmass=lambda _mol, _iso: 31.9988,
    )
    return absorber, partition_temperatures


class _FixedEmission(Constituent):
    def add_to_atmosphere(self, atmo):
        profile = np.exp(-(((atmo.model_geometry.altitudes() - 60_000) / 20_000) ** 2))
        atmo.storage.emission_source[:] += profile[:, np.newaxis]

    def register_derivative(self, atmo, name):
        pass


def _scenario(
    spectral_mode=sk.SpectralGridMode.Monochromatic,
    wavenumbers=None,
    **atmosphere_options,
):
    config = sk.Config()
    config.spectral_grid_mode = spectral_mode
    config.single_scatter_source = sk.SingleScatterSource.NoSource
    config.multiple_scatter_source = sk.MultipleScatterSource.NoSource
    config.emission_source = sk.EmissionSource.VolumeEmissionRate
    altitudes = np.array([0, 10_000, 20_000, 40_000, 80_000, 120_000], dtype=float)
    geometry = sk.Geometry1D(
        0.6,
        0,
        6_372_000,
        altitudes,
        sk.InterpolationMethod.LinearInterpolation,
        sk.GeometryType.Spherical,
    )
    atmosphere = sk.Atmosphere(
        geometry,
        config,
        wavenumber_cminv=(
            np.linspace(13099.8, 13101.0, 151) if wavenumbers is None else wavenumbers
        ),
        **atmosphere_options,
    )
    atmosphere.temperature_k = np.array([288, 240, 220, 250, 210, 200], dtype=float)
    atmosphere.pressure_pa = np.array([1e5, 2.6e4, 6e3, 300, 1, 0.002])
    viewing = sk.ViewingGeometry()
    for tangent in [10_000, 40_000, 70_000]:
        viewing.add_ray(sk.TangentAltitudeSolar(tangent, 0, 200_000, 0.6))
    return atmosphere, sk.Engine(config, geometry, viewing)


@pytest.mark.parametrize(
    "spectral_mode",
    [
        sk.SpectralGridMode.Monochromatic,
        sk.SpectralGridMode.AtmosphereIntegratedLineShape,
    ],
)
def test_cross_section_temperature_derivative(tmp_path, spectral_mode):
    absorber, _ = _absorber(tmp_path)
    atmo, _ = _scenario(spectral_mode=spectral_mode)
    vmr = np.full_like(atmo.temperature_k, 0.21)
    analytic = absorber.optical_derivatives(atmo, vmr=vmr)[
        "temperature_k"
    ].cross_section
    step = 0.001
    atmo.temperature_k += step
    above = absorber.atmosphere_quantities(atmo, vmr=vmr).cross_section.copy()
    atmo.temperature_k -= 2 * step
    below = absorber.atmosphere_quantities(atmo, vmr=vmr).cross_section.copy()
    atmo.temperature_k += step
    numeric = (above - below) / (2 * step)
    assert np.max(np.abs(numeric)) > 0
    np.testing.assert_allclose(analytic, numeric, rtol=2e-5, atol=1e-35)


@pytest.mark.parametrize(
    "spectral_mode",
    [
        sk.SpectralGridMode.Monochromatic,
        sk.SpectralGridMode.AtmosphereIntegratedLineShape,
    ],
)
@pytest.mark.parametrize("temperature_derivative", [False, True])
@pytest.mark.parametrize("pressure_derivative", [False, True])
def test_combined_optical_evaluation(
    tmp_path, spectral_mode, temperature_derivative, pressure_derivative
):
    absorber, sampled = _absorber(tmp_path)
    atmo, _ = _scenario(
        spectral_mode=spectral_mode,
        temperature_derivative=temperature_derivative,
        pressure_derivative=pressure_derivative,
    )
    vmr = np.full_like(atmo.temperature_k, 0.21)
    quantities, derivatives = absorber.atmosphere_quantities_and_derivatives(
        atmo, vmr=vmr
    )
    # One reference Q plus one Q(T) per level, and two extra partition samples
    # per level only when derivatives are enabled. No repeated spectrum pass.
    assert len(sampled) == 1 + len(vmr) * (3 if temperature_derivative else 1)
    expected = absorber.atmosphere_quantities(atmo, vmr=vmr)
    np.testing.assert_allclose(
        quantities.cross_section, expected.cross_section, rtol=1e-10, atol=1e-35
    )
    np.testing.assert_array_equal(quantities.ssa, expected.ssa)
    expected_keys = {
        key
        for key, enabled in [
            ("temperature_k", temperature_derivative),
            ("pressure_pa", pressure_derivative),
        ]
        if enabled
    }
    assert derivatives.keys() == expected_keys
    if expected_keys:
        separate = absorber.optical_derivatives(atmo, vmr=vmr)
        for key in expected_keys:
            np.testing.assert_array_equal(
                derivatives[key].cross_section, separate[key].cross_section
            )
    else:
        assert derivatives == {}
        np.testing.assert_array_equal(quantities.cross_section, expected.cross_section)


@pytest.mark.parametrize(
    ("calculate_derivatives", "temperature_derivative"), [(False, True), (True, False)]
)
def test_disabled_temperature_derivatives_do_not_sample_partition_derivatives(
    tmp_path, calculate_derivatives, temperature_derivative
):
    absorber, sampled = _absorber(tmp_path)
    atmo, engine = _scenario(
        calculate_derivatives=calculate_derivatives,
        temperature_derivative=temperature_derivative,
        pressure_derivative=False,
    )
    atmo["o2"] = sk.constituent.VMRAltitudeAbsorber(
        absorber,
        atmo.model_geometry.altitudes(),
        np.full(6, 0.21),
    )
    atmo["emission"] = _FixedEmission()
    result = engine.calculate_radiance(atmo)
    assert "wf_temperature_k" not in result
    assert sampled
    assert all(t in [296.0, *atmo.temperature_k] for t in sampled)
    if not temperature_derivative:
        assert absorber.optical_derivatives(atmo) == {}


def test_temperature_jacobian_includes_self_absorption_by_default(tmp_path):
    absorber, sampled = _absorber(tmp_path)
    atmo, engine = _scenario()
    atmo["o2"] = sk.constituent.VMRAltitudeAbsorber(
        absorber,
        atmo.model_geometry.altitudes(),
        np.full(6, 0.21),
    )
    atmo["emission"] = _FixedEmission()
    result = engine.calculate_radiance(atmo)
    # One forward evaluation to assemble extinction, then one combined
    # evaluation for all derivative mappings, without a third line-list pass.
    assert len(sampled) == 2 + 4 * len(atmo.temperature_k)
    assert "wf_o2_temperature_k_xs" in atmo.storage.derivative_mapping_names()
    analytic = result.wf_temperature_k
    step = 0.001
    for level in range(len(atmo.temperature_k)):
        original = atmo.temperature_k[level]
        atmo.temperature_k[level] = original + step
        above = engine.calculate_radiance(atmo).radiance
        atmo.temperature_k[level] = original - step
        below = engine.calculate_radiance(atmo).radiance
        atmo.temperature_k[level] = original
        numeric = (above - below) / (2 * step)
        np.testing.assert_allclose(
            analytic.isel(altitude=level),
            numeric,
            rtol=5e-5,
            atol=2e-7,
        )


@pytest.mark.parametrize(
    "spectral_mode",
    [
        sk.SpectralGridMode.Monochromatic,
        sk.SpectralGridMode.AtmosphereIntegratedLineShape,
    ],
)
@pytest.mark.parametrize("vmr_value", [None, 0.0, 0.21, 1.0])
def test_pressure_cross_section_derivative(tmp_path, spectral_mode, vmr_value):
    absorber, sampled = _absorber(tmp_path)
    atmo, _ = _scenario(spectral_mode=spectral_mode, temperature_derivative=False)
    # Keep this finite-difference stencil above the clipped low-pressure wings.
    # Zero pressure and clipping are checked separately in the Rust kernel tests.
    atmo.pressure_pa = np.maximum(atmo.pressure_pa, 1.0)
    kwargs = (
        {} if vmr_value is None else {"vmr": np.full_like(atmo.pressure_pa, vmr_value)}
    )
    derivatives = absorber.optical_derivatives(atmo, **kwargs)
    assert derivatives.keys() == {"pressure_pa"}
    assert len(sampled) == 1 + len(atmo.temperature_k)
    analytic = derivatives["pressure_pa"].cross_section
    original = atmo.pressure_pa.copy()
    step = np.maximum(original * 1e-3, 100.0)
    # A fourth-order forward difference keeps even the smallest pressure positive
    # and resolves tiny pressure-induced shifts of the large line-center values.
    samples = []
    for i in range(5):
        atmo.pressure_pa = original + i * step
        samples.append(
            absorber.atmosphere_quantities(atmo, **kwargs).cross_section.copy()
        )
    atmo.pressure_pa = original
    numeric = sum(
        weight * value
        for weight, value in zip([-25, 48, -36, 16, -3], samples, strict=True)
    ) / (12 * step[:, None])
    scale = np.max(np.abs(numeric), axis=1, keepdims=True)
    assert np.all(scale > 0)
    np.testing.assert_allclose(analytic / scale, numeric / scale, rtol=2e-4, atol=2e-5)


@pytest.mark.parametrize("temperature_derivative", [False, True])
def test_pressure_radiance_jacobian(tmp_path, temperature_derivative):
    absorber, sampled = _absorber(tmp_path)
    atmo, engine = _scenario(temperature_derivative=temperature_derivative)
    atmo["o2"] = sk.constituent.VMRAltitudeAbsorber(
        absorber, atmo.model_geometry.altitudes(), np.full(6, 0.21)
    )
    atmo["emission"] = _FixedEmission()
    result = engine.calculate_radiance(atmo)
    assert "wf_o2_pressure_pa_xs" in atmo.storage.derivative_mapping_names()
    assert len(sampled) == 2 + len(atmo.temperature_k) * (
        4 if temperature_derivative else 2
    )
    if not temperature_derivative:
        assert "wf_temperature_k" not in result
    for level in range(len(atmo.pressure_pa)):
        original = atmo.pressure_pa[level]
        # Resolve the very small radiance change at the highest altitude while
        # keeping the stencil on the same side of the clipped line-wing boundary.
        step = min(original * 0.025, max(original * 1e-3, 1e-3))
        atmo.pressure_pa[level] = original + step
        above = engine.calculate_radiance(atmo).radiance
        atmo.pressure_pa[level] = original - step
        below = engine.calculate_radiance(atmo).radiance
        atmo.pressure_pa[level] = original
        numeric = (above - below) / (2 * step)
        np.testing.assert_allclose(
            result.wf_pressure_pa.isel(altitude=level), numeric, rtol=2e-4, atol=2e-6
        )


def test_pressure_derivative_switches(tmp_path):
    absorber, sampled = _absorber(tmp_path)
    for calculate_derivatives, pressure_derivative in [(False, True), (True, False)]:
        atmo, engine = _scenario(
            calculate_derivatives=calculate_derivatives,
            pressure_derivative=pressure_derivative,
            temperature_derivative=False,
        )
        atmo["o2"] = sk.constituent.VMRAltitudeAbsorber(
            absorber, atmo.model_geometry.altitudes(), np.full(6, 0.21)
        )
        atmo["emission"] = _FixedEmission()
        result = engine.calculate_radiance(atmo)
        assert "wf_pressure_pa" not in result
        assert "wf_o2_pressure_pa_xs" not in atmo.storage.derivative_mapping_names()
        if not pressure_derivative:
            assert absorber.optical_derivatives(atmo) == {}
        assert all(t in [296.0, *atmo.temperature_k] for t in sampled)


def test_line_derivatives_compose_with_python_optical_properties(tmp_path):
    absorber, _ = _absorber(tmp_path)
    atmo, engine = _scenario()
    vmr = np.full(6, 0.21)
    derivatives = absorber.optical_derivatives(atmo, vmr=vmr)
    combined = absorber + absorber
    summed = combined.optical_derivatives(atmo, vmr=vmr)
    for key, derivative in derivatives.items():
        np.testing.assert_array_equal(derivative.d_extinction, derivative.cross_section)
        np.testing.assert_array_equal(derivative.d_ssa, derivative.ssa)
        assert derivative.d_leg_coeff is None
        np.testing.assert_array_equal(
            summed[key].d_extinction, 2 * derivative.cross_section
        )

    atmo["o2"] = sk.constituent.VMRAltitudeAbsorber(
        combined, atmo.model_geometry.altitudes(), vmr
    )
    atmo["emission"] = _FixedEmission()
    result = engine.calculate_radiance(atmo)
    reference_atmo, reference_engine = _scenario()
    for name in ["o2a", "o2b"]:
        reference_atmo[name] = sk.constituent.VMRAltitudeAbsorber(
            absorber, reference_atmo.model_geometry.altitudes(), vmr
        )
    reference_atmo["emission"] = _FixedEmission()
    reference = reference_engine.calculate_radiance(reference_atmo)
    for name in ["radiance", "wf_temperature_k", "wf_pressure_pa"]:
        # Separate mappings sum the density and cross-section terms in a
        # different floating-point order, including where they nearly cancel.
        np.testing.assert_allclose(
            result[name], reference[name], rtol=1e-10, atol=1e-10
        )


def test_temperature_derivative_at_partition_table_boundaries(tmp_path):
    def partition(temperature):
        if not 150.0 <= temperature <= 350.0:
            msg = "Temperature outside the partition table"
            raise ValueError(msg)
        return temperature**1.5

    absorber, _ = _absorber(tmp_path, partition)
    atmo, _ = _scenario(pressure_derivative=False)
    original = np.array([150.0, 150.00001, 230.0, 296.0, 349.99999, 350.0])
    atmo.temperature_k = original.copy()
    analytic = absorber.optical_derivatives(atmo)["temperature_k"].cross_section
    # An independent fourth-order one-sided stencil stays within the table.
    step = np.where(original < 250.0, 0.01, -0.01)
    samples = []
    for i in range(5):
        atmo.temperature_k = original + i * step
        samples.append(absorber.atmosphere_quantities(atmo).cross_section.copy())
    numeric = sum(
        weight * value
        for weight, value in zip([-25, 48, -36, 16, -3], samples, strict=True)
    ) / (12 * step[:, None])
    np.testing.assert_allclose(analytic, numeric, rtol=2e-5, atol=1e-35)


def test_partition_derivative_callback_errors_are_propagated(tmp_path):
    def partition(temperature):
        if temperature != 296.0:
            msg = "Partition samples unavailable"
            raise ValueError(msg)
        return temperature**1.5

    absorber, _ = _absorber(tmp_path, partition)
    atmo, _ = _scenario(pressure_derivative=False)
    atmo.temperature_k[:] = 296.0
    absorber.atmosphere_quantities(atmo)
    with pytest.raises(RuntimeError, match="Partition samples unavailable"):
        absorber.optical_derivatives(atmo)
    atmo["o2"] = sk.constituent.VMRAltitudeAbsorber(
        absorber, atmo.model_geometry.altitudes(), np.full(6, 0.21)
    )
    with pytest.raises(RuntimeError, match="Partition samples unavailable"):
        atmo.internal_object()


def test_derivatives_reduce_from_distinct_fine_spectral_grid(tmp_path, monkeypatch):
    absorber, _ = _absorber(tmp_path)
    fine_atmo, _ = _scenario()
    fine_grid = sk.basis.Grid.from_triangles(fine_atmo.wavenumbers_cminv)
    coarse_wavenumbers = np.linspace(13099.8, 13101.0, 11)
    coarse_grid = sk.basis.Grid.from_triangles(coarse_wavenumbers)
    native_atmosphere = atmosphere_module.PyAtmosphere

    def with_fine_grid(*args):
        return native_atmosphere(*args[:7], fine_grid._internal_object(), *args[8:])

    # The public constructor currently uses the same fine and output grids.
    # Supply distinct native grids to check actual spectral reduction.
    with monkeypatch.context() as patch:
        patch.setattr(atmosphere_module, "PyAtmosphere", with_fine_grid)
        coarse_atmo, _ = _scenario(
            spectral_mode=sk.SpectralGridMode.AtmosphereIntegratedLineShape,
            wavenumbers=coarse_wavenumbers,
        )
    vmr = np.full(6, 0.21)
    fine, fine_derivatives = absorber.atmosphere_quantities_and_derivatives(
        fine_atmo, vmr=vmr
    )
    coarse, coarse_derivatives = absorber.atmosphere_quantities_and_derivatives(
        coarse_atmo, vmr=vmr
    )
    separate = absorber.optical_derivatives(coarse_atmo, vmr=vmr)
    mapping = fine_grid.mapping_to(coarse_grid, normalize=False)
    assert coarse.cross_section.shape == (6, 11)
    np.testing.assert_allclose(
        coarse.cross_section, fine.cross_section @ mapping.T, rtol=1e-12
    )
    for key in ["temperature_k", "pressure_pa"]:
        np.testing.assert_allclose(
            coarse_derivatives[key].cross_section,
            fine_derivatives[key].cross_section @ mapping.T,
            rtol=1e-12,
            atol=1e-35,
        )
        np.testing.assert_array_equal(
            separate[key].cross_section, coarse_derivatives[key].cross_section
        )
