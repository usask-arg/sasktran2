from __future__ import annotations

import numpy as np
import pytest
import sasktran2 as sk
import xarray as xr
from sasktran2.linearization import Linearization, _ParameterSpec


def _synthetic_linearization() -> Linearization:
    wavelength = xr.DataArray([300.0, 310.0], dims="wavelength")
    value = xr.DataArray(
        np.arange(4.0).reshape(2, 2, 1),
        dims=("wavelength", "los", "stokes"),
        coords={"wavelength": wavelength, "stokes": ["I"]},
    )
    profile = xr.DataArray(
        np.arange(12.0).reshape(3, 2, 2, 1),
        dims=("altitude", *value.dims),
        coords={"altitude": [0.0, 10_000.0, 20_000.0], **value.coords},
    )
    diagonal_surface = xr.DataArray(
        np.arange(4.0, 8.0).reshape(2, 2, 1),
        dims=value.dims,
        coords=value.coords,
    )
    scalar_surface = xr.DataArray(
        np.arange(8.0, 12.0).reshape(2, 2, 1),
        dims=value.dims,
        coords=value.coords,
    )
    structured = xr.DataArray(
        np.arange(24.0).reshape(2, 3, 2, 2, 1),
        dims=("horizontal_angle", "altitude", *value.dims),
        coords={
            "horizontal_angle": [-0.1, 0.1],
            "altitude": [0.0, 10_000.0, 20_000.0],
            **value.coords,
        },
    )
    result = xr.Dataset(
        {
            "radiance": value,
            "wf_profile": profile,
            "wf_surface": diagonal_surface,
            "wf_scalar": scalar_surface,
            "wf_structured": structured,
        }
    )
    specs = {
        "wf_profile": _ParameterSpec(("altitude",)),
        "wf_surface": _ParameterSpec(("wavelength",), ("wavelength",)),
        "wf_scalar": _ParameterSpec(()),
        "wf_structured": _ParameterSpec(("horizontal_angle", "altitude")),
    }
    return Linearization._from_result(result, specs)


def test_structured_jvp_vjp_and_adjoint_identity():
    lin = _synthetic_linearization()
    assert lin.parameters == ("profile", "surface", "scalar", "structured")
    assert lin.parameter_dims["surface"] == ("wavelength",)
    assert lin.backends == {
        "jvp": sk.LinearizationBackend.MaterializedJacobian,
        "vjp": sk.LinearizationBackend.MaterializedJacobian,
    }

    tangent = xr.Dataset(
        {
            "profile": xr.DataArray(
                [1.0, 2.0, 3.0],
                dims="altitude",
                coords={"altitude": lin.jacobian.altitude},
            ),
            "surface": xr.DataArray(
                [4.0, 5.0],
                dims="wavelength",
                coords={"wavelength": lin.value.wavelength},
            ),
            "scalar": xr.DataArray(6.0),
            "structured": xr.DataArray(
                np.arange(6.0).reshape(2, 3),
                dims=("horizontal_angle", "altitude"),
                coords={
                    "horizontal_angle": lin.jacobian.horizontal_angle,
                    "altitude": lin.jacobian.altitude,
                },
            ),
        }
    )

    expected = (
        (lin.jacobian["profile"] * tangent["profile"]).sum("altitude")
        + lin.jacobian["surface"] * tangent["surface"]
        + lin.jacobian["scalar"] * tangent["scalar"]
        + (lin.jacobian["structured"] * tangent["structured"]).sum(
            ("horizontal_angle", "altitude")
        )
    )
    xr.testing.assert_allclose(lin.jvp(tangent), expected)

    cotangent = xr.DataArray(
        np.arange(1.0, 5.0).reshape(2, 2, 1),
        dims=lin.value.dims,
        coords=lin.value.coords,
    )
    gradient = lin.vjp(cotangent)
    xr.testing.assert_allclose(
        gradient["profile"],
        (lin.jacobian["profile"] * cotangent).sum(lin.value.dims),
    )
    xr.testing.assert_allclose(
        gradient["surface"],
        (lin.jacobian["surface"] * cotangent).sum(("los", "stokes")),
    )
    xr.testing.assert_allclose(
        gradient["scalar"],
        (lin.jacobian["scalar"] * cotangent).sum(lin.value.dims),
    )

    lhs = float((lin.jvp(tangent) * cotangent).sum())
    rhs = sum(float((tangent[name] * gradient[name]).sum()) for name in tangent)
    np.testing.assert_allclose(lhs, rhs)


def test_vjp_can_select_parameters():
    lin = _synthetic_linearization()
    cotangent = xr.ones_like(lin.value)

    selected = lin.vjp(cotangent, parameters=("structured", "scalar"))

    assert tuple(selected.data_vars) == ("structured", "scalar")
    xr.testing.assert_allclose(
        selected["structured"],
        (lin.jacobian["structured"] * cotangent).sum(lin.value.dims),
    )
    xr.testing.assert_allclose(
        selected["scalar"],
        (lin.jacobian["scalar"] * cotangent).sum(lin.value.dims),
    )
    assert len(lin.vjp(cotangent, parameters=()).data_vars) == 0


def test_vjp_rejects_invalid_parameter_selection():
    lin = _synthetic_linearization()
    cotangent = xr.ones_like(lin.value)

    with pytest.raises(ValueError, match="Unknown VJP parameter"):
        lin.vjp(cotangent, parameters=("missing",))
    with pytest.raises(ValueError, match="must not contain duplicates"):
        lin.vjp(cotangent, parameters=("profile", "profile"))
    with pytest.raises(TypeError, match="must be strings"):
        lin.vjp(cotangent, parameters=("profile", 1))


def test_tangent_template_describes_parameter_domain():
    lin = _synthetic_linearization()
    template = lin.tangent_template

    assert tuple(template.data_vars) == lin.parameters
    assert template["profile"].dims == ("altitude",)
    assert template["profile"].shape == (3,)
    assert template["structured"].dims == ("horizontal_angle", "altitude")
    assert template["structured"].shape == (2, 3)
    assert template["surface"].dims == ("wavelength",)
    assert template["scalar"].dims == ()
    xr.testing.assert_identical(template["profile"].altitude, lin.jacobian.altitude)
    xr.testing.assert_identical(template["surface"].wavelength, lin.value.wavelength)
    assert all(np.count_nonzero(template[name]) == 0 for name in lin.parameters)

    template["profile"].data[:] = 1
    assert np.count_nonzero(lin.tangent_template["profile"]) == 0


def test_jvp_validation_and_omitted_parameters():
    lin = _synthetic_linearization()
    xr.testing.assert_identical(lin.jvp(xr.Dataset()), xr.zeros_like(lin.value))

    with pytest.raises(ValueError, match="Unknown tangent"):
        lin.jvp(xr.Dataset({"missing": xr.DataArray(1.0)}))

    with pytest.raises(ValueError, match="must have dimensions"):
        lin.jvp(
            xr.Dataset(
                {
                    "profile": xr.DataArray(
                        np.ones(3),
                        dims="wrong",
                    )
                }
            )
        )

    with pytest.raises(ValueError, match="coordinate"):
        lin.jvp(
            xr.Dataset(
                {
                    "surface": xr.DataArray(
                        np.ones(2),
                        dims="wavelength",
                        coords={"wavelength": [301.0, 311.0]},
                    )
                }
            )
        )


@pytest.mark.parametrize("dtype", [np.int64, np.float32])
def test_engine_products_accept_real_numeric_dtypes(dtype):
    engine, atmosphere = _raw_engine_scenario()
    lin = engine.linearize(atmosphere)

    tangent = xr.Dataset(
        {
            "extinction": xr.DataArray(
                np.ones(atmosphere.num_locations, dtype=dtype),
                dims="altitude",
                coords={"altitude": lin.tangent_template.altitude},
            )
        }
    )
    expected_jvp = lin.jacobian["extinction"].sum("altitude")
    xr.testing.assert_allclose(lin.jvp(tangent), expected_jvp)

    cotangent = xr.ones_like(lin.value).astype(dtype)
    expected_vjp = (lin.jacobian["extinction"] * cotangent).sum(lin.value.dims)
    xr.testing.assert_allclose(lin.vjp(cotangent)["extinction"], expected_vjp)


def test_engine_products_reject_complex_inputs():
    engine, atmosphere = _raw_engine_scenario()
    lin = engine.linearize(atmosphere)
    tangent = lin.tangent_template[["extinction"]].astype(np.complex128)

    with pytest.raises(ValueError, match="real numeric"):
        lin.jvp(tangent)
    with pytest.raises(ValueError, match="real numeric"):
        lin.vjp(xr.ones_like(lin.value).astype(np.complex128))


def _raw_engine_scenario(
    *,
    calculate_derivatives: bool = True,
    interpolation: sk.InterpolationMethod = sk.InterpolationMethod.LinearInterpolation,
):
    config = sk.Config()
    config.single_scatter_source = sk.SingleScatterSource.Exact
    geometry = sk.Geometry1D(
        0.6,
        0.0,
        6_372_000.0,
        np.arange(0.0, 30_001.0, 5_000.0),
        interpolation,
        sk.GeometryType.Spherical,
    )
    viewing = sk.ViewingGeometry()
    viewing.add_ray(sk.GroundViewingSolar(0.6, 0.0, 0.8, 100_000.0))
    if calculate_derivatives:
        atmosphere = sk.test_util.scenarios.default_pure_scattering_atmosphere(
            config, geometry, 0.8, albedo=0.3
        )
    else:
        atmosphere = sk.Atmosphere(
            geometry,
            config,
            numwavel=1,
            calculate_derivatives=False,
        )
    return sk.Engine(config, geometry, viewing), atmosphere


def _constituent_engine_scenario(surface, *, num_stokes: int = 1):
    config = sk.Config()
    config.single_scatter_source = sk.SingleScatterSource.Exact
    config.multiple_scatter_source = sk.MultipleScatterSource.NoSource
    config.num_stokes = num_stokes
    geometry = sk.Geometry1D(
        0.6,
        0.0,
        6_372_000.0,
        np.arange(0.0, 30_001.0, 5_000.0),
        sk.InterpolationMethod.LinearInterpolation,
        sk.GeometryType.Spherical,
    )
    viewing = sk.ViewingGeometry()
    viewing.add_ray(sk.GroundViewingSolar(0.6, 0.2, 0.8, 100_000.0))
    atmosphere = sk.Atmosphere(
        geometry, config, wavelengths_nm=np.array([300.0, 350.0, 400.0])
    )
    sk.climatology.us76.add_us76_standard_atmosphere(atmosphere)
    atmosphere["rayleigh"] = sk.constituent.Rayleigh()
    atmosphere["surface"] = surface
    return sk.Engine(config, geometry, viewing), atmosphere


def test_engine_linearize_matches_weighting_function_output():
    engine, atmosphere = _raw_engine_scenario()
    expected = engine.calculate_radiance(atmosphere)
    lin = engine.linearize(atmosphere)

    assert engine._engine._supports_linearization(0)
    assert engine._engine._supports_linearization(1)
    assert engine._engine._supports_linearization(2)
    assert engine._engine._linearization_backend(0) == 2
    assert engine._engine._linearization_backend(1) == 2
    assert engine._engine._linearization_backend(2) == 2
    assert lin.backends == {
        "jvp": sk.LinearizationBackend.Native,
        "vjp": sk.LinearizationBackend.Native,
    }
    xr.testing.assert_allclose(lin.value, expected["radiance"])
    xr.testing.assert_allclose(lin.jacobian["extinction"], expected["wf_extinction"])
    xr.testing.assert_allclose(lin.jacobian["ssa"], expected["wf_ssa"])
    xr.testing.assert_allclose(lin.jacobian["albedo"], expected["wf_albedo"])


@pytest.mark.parametrize(
    "surface",
    [
        sk.constituent.LambertianSurface(0.3),
        sk.constituent.LambertianSurface(
            [0.2, 0.4], wavelengths_nm=np.array([290.0, 410.0])
        ),
        sk.constituent.LambertianSurface([0.2, 0.3, 0.4]),
    ],
)
def test_native_products_cover_surface_parameterizations(surface):
    engine, atmosphere = _constituent_engine_scenario(surface)
    lin = engine.linearize(atmosphere)
    name = "surface_albedo"

    tangent = lin.tangent_template[[name]]
    tangent[name].data[...] = np.linspace(0.2, 0.6, tangent[name].size).reshape(
        tangent[name].shape
    )
    expected_jvp = (lin.jacobian[name] * tangent[name]).sum(lin.parameter_dims[name])
    xr.testing.assert_allclose(lin.jvp(tangent), expected_jvp)

    cotangent = xr.ones_like(lin.value) * 0.7
    expected_vjp = (lin.jacobian[name] * cotangent).sum(lin.value.dims)
    gradient = lin.vjp(cotangent, parameters=(name,))
    assert tuple(gradient.data_vars) == (name,)
    xr.testing.assert_allclose(gradient[name], expected_vjp)


def test_constant_surface_parameter_has_scalar_domain():
    engine, atmosphere = _constituent_engine_scenario(
        sk.constituent.LambertianSurface(0.3)
    )
    lin = engine.linearize(atmosphere)

    assert lin.parameter_dims["surface_albedo"] == ()
    assert lin.tangent_template["surface_albedo"].dims == ()
    assert lin.jacobian["surface_albedo"].dims == lin.value.dims


def test_native_products_apply_polarized_output_rotation():
    engine, atmosphere = _constituent_engine_scenario(
        sk.constituent.LambertianSurface([0.2, 0.3, 0.4]), num_stokes=3
    )
    lin = engine.linearize(atmosphere)
    tangent = lin.tangent_template[["surface_albedo"]]
    tangent["surface_albedo"].data[:] = [0.3, -0.2, 0.5]
    expected = (lin.jacobian["surface_albedo"] * tangent["surface_albedo"]).sum(
        lin.parameter_dims["surface_albedo"]
    )
    xr.testing.assert_allclose(lin.jvp(tangent), expected)

    cotangent = xr.ones_like(lin.value)
    xr.testing.assert_allclose(
        lin.vjp(cotangent)["surface_albedo"],
        (lin.jacobian["surface_albedo"] * cotangent).sum(lin.value.dims),
    )


def test_constituent_rebuild_invalidates_linearization():
    engine, atmosphere = _constituent_engine_scenario(
        sk.constituent.LambertianSurface(0.3)
    )
    lin = engine.linearize(atmosphere)
    revision = atmosphere.revision

    engine.calculate_radiance(atmosphere)

    assert atmosphere.revision == revision + 1
    with pytest.raises(sk.StaleLinearizationError):
        lin.jvp(lin.tangent_template[["surface_albedo"]])


def test_native_products_apply_assign_name_parameter_mapping():
    engine, atmosphere = _constituent_engine_scenario(
        sk.constituent.LambertianSurface(0.3)
    )
    lin = engine.linearize(atmosphere)
    assert "pressure_pa" in lin.parameters

    tangent = lin.tangent_template[["pressure_pa"]]
    tangent["pressure_pa"].data[:] = np.linspace(0.1, 0.7, 7)
    xr.testing.assert_allclose(
        lin.jvp(tangent),
        (lin.jacobian["pressure_pa"] * tangent["pressure_pa"]).sum(
            lin.parameter_dims["pressure_pa"]
        ),
    )

    cotangent = xr.ones_like(lin.value) * 0.2
    xr.testing.assert_allclose(
        lin.vjp(cotangent)["pressure_pa"],
        (lin.jacobian["pressure_pa"] * cotangent).sum(lin.value.dims),
    )


def test_scalar_surface_emission_products():
    config = sk.Config()
    config.single_scatter_source = sk.SingleScatterSource.NoSource
    config.multiple_scatter_source = sk.MultipleScatterSource.NoSource
    config.emission_source = sk.EmissionSource.Standard
    geometry = sk.Geometry1D(
        0.6,
        0.0,
        6_372_000.0,
        np.arange(0.0, 30_001.0, 5_000.0),
        sk.InterpolationMethod.LinearInterpolation,
        sk.GeometryType.Spherical,
    )
    viewing = sk.ViewingGeometry()
    viewing.add_ray(sk.GroundViewingSolar(0.6, 0.0, 0.8, 100_000.0))
    wavelengths = np.array([8_000.0, 10_000.0])
    atmosphere = sk.Atmosphere(geometry, config, wavelengths_nm=wavelengths)
    atmosphere["manual"] = sk.constituent.Manual(
        np.full((7, 2), 1.0e-5), np.zeros((7, 2))
    )
    atmosphere["surface"] = sk.constituent.SurfaceThermalEmission(290.0, 0.8)
    engine = sk.Engine(config, geometry, viewing)
    lin = engine.linearize(atmosphere)

    assert not engine._engine._supports_linearization(1)
    assert not engine._engine._supports_linearization(2)
    assert engine._engine._linearization_backend(1) == 1
    assert engine._engine._linearization_backend(2) == 1
    assert lin.backends == {
        "jvp": sk.LinearizationBackend.StreamingJacobian,
        "vjp": sk.LinearizationBackend.StreamingJacobian,
    }

    assert lin.parameter_dims["surface_temperature_k"] == ()
    assert lin.parameter_dims["surface_emissivity"] == ()
    tangent = xr.Dataset(
        {
            "surface_temperature_k": xr.DataArray(1.3),
            "surface_emissivity": xr.DataArray(-0.2),
        }
    )
    expected = (
        lin.jacobian["surface_temperature_k"] * tangent["surface_temperature_k"]
        + lin.jacobian["surface_emissivity"] * tangent["surface_emissivity"]
    )
    xr.testing.assert_allclose(lin.jvp(tangent), expected)

    cotangent = xr.ones_like(lin.value) * 0.6
    gradient = lin.vjp(cotangent)
    for name in tangent:
        xr.testing.assert_allclose(
            gradient[name], (lin.jacobian[name] * cotangent).sum(lin.value.dims)
        )
    selected_gradient = lin.vjp(cotangent, parameters=("surface_temperature_k",))
    assert tuple(selected_gradient.data_vars) == ("surface_temperature_k",)
    xr.testing.assert_allclose(
        selected_gradient["surface_temperature_k"],
        gradient["surface_temperature_k"],
    )


def test_log_radiance_mapping_is_converted_for_linearization():
    engine, atmosphere = _constituent_engine_scenario(
        sk.constituent.LambertianSurface(0.3)
    )
    atmosphere["amf"] = sk.constituent.AirMassFactor()
    legacy = engine.calculate_radiance(atmosphere)
    lin = engine.linearize(atmosphere)

    xr.testing.assert_allclose(
        lin.jacobian["air_mass_factor"],
        legacy["air_mass_factor"] * lin.value,
    )
    tangent = lin.tangent_template[["air_mass_factor"]]
    tangent["air_mass_factor"].data[:] = 1
    xr.testing.assert_allclose(
        lin.jvp(tangent),
        lin.jacobian["air_mass_factor"].sum("altitude"),
    )


def test_engine_native_products_match_materialized_jacobian():
    engine, atmosphere = _raw_engine_scenario()
    lin = engine.linearize(atmosphere)

    tangent = lin.tangent_template[["extinction", "ssa", "leg_coeff_2", "albedo"]]
    tangent["extinction"].data[:] = np.linspace(0.2, 0.8, 7)
    tangent["ssa"].data[:] = np.linspace(-0.3, 0.1, 7)
    tangent["leg_coeff_2"].data[:] = np.linspace(0.05, 0.2, 7)
    tangent["albedo"].data[:] = 0.4
    expected_jvp = (
        (lin.jacobian["extinction"] * tangent["extinction"]).sum("altitude")
        + (lin.jacobian["ssa"] * tangent["ssa"]).sum("altitude")
        + (lin.jacobian["leg_coeff_2"] * tangent["leg_coeff_2"]).sum("altitude")
        + lin.jacobian["albedo"] * tangent["albedo"]
    )
    xr.testing.assert_allclose(lin.jvp(tangent), expected_jvp)

    cotangent = xr.ones_like(lin.value) * 1.7
    gradient = lin.vjp(cotangent)
    xr.testing.assert_allclose(
        gradient["extinction"],
        (lin.jacobian["extinction"] * cotangent).sum(lin.value.dims),
    )
    xr.testing.assert_allclose(
        gradient["ssa"],
        (lin.jacobian["ssa"] * cotangent).sum(lin.value.dims),
    )
    xr.testing.assert_allclose(
        gradient["leg_coeff_2"],
        (lin.jacobian["leg_coeff_2"] * cotangent).sum(lin.value.dims),
    )
    xr.testing.assert_allclose(
        gradient["albedo"],
        (lin.jacobian["albedo"] * cotangent).sum(("los", "stokes")),
    )

    lhs = float((lin.jvp(tangent) * cotangent).sum())
    rhs = sum(float((tangent[name] * gradient[name]).sum()) for name in tangent)
    np.testing.assert_allclose(lhs, rhs)


def test_distinct_atmosphere_linearisations_coexist_on_one_engine():
    engine, atmosphere_a = _constituent_engine_scenario(
        sk.constituent.LambertianSurface(0.2)
    )
    _, atmosphere_b = _constituent_engine_scenario(
        sk.constituent.LambertianSurface(0.6)
    )

    linearization_a = engine.linearize(atmosphere_a)
    tangent_a = linearization_a.tangent_template[["surface_albedo"]]
    tangent_a["surface_albedo"].data[...] = 0.4
    expected_a = (
        linearization_a.jacobian["surface_albedo"] * tangent_a["surface_albedo"]
    )

    linearization_b = engine.linearize(atmosphere_b)
    tangent_b = linearization_b.tangent_template[["surface_albedo"]]
    tangent_b["surface_albedo"].data[...] = -0.3
    expected_b = (
        linearization_b.jacobian["surface_albedo"] * tangent_b["surface_albedo"]
    )

    assert not np.allclose(linearization_a.value, linearization_b.value)
    xr.testing.assert_allclose(linearization_a.jvp(tangent_a), expected_a)
    xr.testing.assert_allclose(linearization_b.jvp(tangent_b), expected_b)

    cotangent_a = xr.ones_like(linearization_a.value) * 0.7
    expected_gradient_a = (
        linearization_a.jacobian["surface_albedo"] * cotangent_a
    ).sum(linearization_a.value.dims)
    gradient_a = linearization_a.vjp(cotangent_a, parameters=("surface_albedo",))
    xr.testing.assert_allclose(gradient_a["surface_albedo"], expected_gradient_a)


def test_native_products_match_lower_interpolation_jacobian():
    engine, atmosphere = _raw_engine_scenario(
        interpolation=sk.InterpolationMethod.LowerInterpolation
    )
    expected = engine.calculate_radiance(atmosphere)
    lin = engine.linearize(atmosphere)
    xr.testing.assert_allclose(lin.value, expected["radiance"])

    tangent = lin.tangent_template[["extinction", "ssa"]]
    tangent["extinction"].data[:] = np.linspace(0.2, 0.8, 7)
    tangent["ssa"].data[:] = np.linspace(-0.3, 0.1, 7)
    expected_jvp = (lin.jacobian["extinction"] * tangent["extinction"]).sum(
        "altitude"
    ) + (lin.jacobian["ssa"] * tangent["ssa"]).sum("altitude")
    xr.testing.assert_allclose(lin.jvp(tangent), expected_jvp)

    cotangent = xr.ones_like(lin.value) * 1.7
    gradient = lin.vjp(cotangent)
    for name in tangent:
        xr.testing.assert_allclose(
            gradient[name],
            (lin.jacobian[name] * cotangent).sum(lin.value.dims),
        )


def test_native_products_do_not_materialize_jacobian():
    engine, atmosphere = _raw_engine_scenario()
    lin = engine.linearize(atmosphere)
    assert lin._backend._jacobian is None

    tangent = lin.tangent_template[["extinction"]]
    tangent["extinction"].data[:] = 1
    lin.jvp(tangent)
    lin.vjp(xr.ones_like(lin.value))

    assert lin._backend._jacobian is None


def test_engine_linearization_detects_changed_atmosphere():
    engine, atmosphere = _raw_engine_scenario()
    lin = engine.linearize(atmosphere)
    value_at_point = lin.value.copy(deep=True)
    jacobian_at_point = lin.jacobian.copy(deep=True)

    tangent = lin.tangent_template[["extinction"]]
    tangent["extinction"].data[:] = 1

    atmosphere.storage.total_extinction[:] = 0
    atmosphere.mark_changed()
    updated_value = engine.calculate_radiance(atmosphere)["radiance"]

    assert float(np.abs(updated_value - value_at_point).max()) > 0
    xr.testing.assert_identical(lin.value, value_at_point)
    xr.testing.assert_identical(lin.jacobian, jacobian_at_point)
    with pytest.raises(sk.StaleLinearizationError):
        lin.jvp(tangent)
    with pytest.raises(sk.StaleLinearizationError):
        lin.vjp(xr.ones_like(lin.value))


def test_unmaterialized_jacobian_rejects_changed_atmosphere():
    engine, atmosphere = _raw_engine_scenario()
    lin = engine.linearize(atmosphere)

    atmosphere.mark_changed()
    with pytest.raises(sk.StaleLinearizationError):
        _ = lin.jacobian


def test_engine_linearize_requires_derivatives():
    engine, atmosphere = _raw_engine_scenario(calculate_derivatives=False)
    with pytest.raises(ValueError, match="calculate_derivatives=True"):
        engine.linearize(atmosphere)
