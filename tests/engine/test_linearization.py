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


def _raw_engine_scenario(*, calculate_derivatives: bool = True):
    config = sk.Config()
    config.single_scatter_source = sk.SingleScatterSource.Exact
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


def test_engine_linearize_matches_weighting_function_output():
    engine, atmosphere = _raw_engine_scenario()
    expected = engine.calculate_radiance(atmosphere)
    lin = engine.linearize(atmosphere)

    assert engine._engine._supports_linearization(0)
    assert not engine._engine._supports_linearization(1)
    assert not engine._engine._supports_linearization(2)
    xr.testing.assert_allclose(lin.value, expected["radiance"])
    xr.testing.assert_allclose(lin.jacobian["extinction"], expected["wf_extinction"])
    xr.testing.assert_allclose(lin.jacobian["ssa"], expected["wf_ssa"])
    xr.testing.assert_allclose(lin.jacobian["albedo"], expected["wf_albedo"])


def test_engine_linearization_is_fixed_at_atmosphere_state():
    engine, atmosphere = _raw_engine_scenario()
    lin = engine.linearize(atmosphere)
    value_at_point = lin.value.copy(deep=True)

    tangent = lin.tangent_template[["extinction"]]
    tangent["extinction"].data[:] = 1
    product_at_point = lin.jvp(tangent).copy(deep=True)

    atmosphere.storage.total_extinction[:] = 0
    updated_value = engine.calculate_radiance(atmosphere)["radiance"]

    assert float(np.abs(updated_value - value_at_point).max()) > 0
    xr.testing.assert_identical(lin.value, value_at_point)
    xr.testing.assert_identical(lin.jvp(tangent), product_at_point)


def test_engine_linearize_requires_derivatives():
    engine, atmosphere = _raw_engine_scenario(calculate_derivatives=False)
    with pytest.raises(ValueError, match="calculate_derivatives=True"):
        engine.linearize(atmosphere)
