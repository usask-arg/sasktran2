from __future__ import annotations

import numpy as np
import pytest
import sasktran2 as sk

EARTH_RADIUS = 6372000.0
TANGENT_ALTITUDE = 20000.0
OBSERVER_ALTITUDE = 200000.0
TANGENT_COS_SZA = 0.6
RELATIVE_AZIMUTHS = np.deg2rad([-150, -120, -90, -60, -30, 30, 60, 90, 120, 150])

BASES = [sk.StokesBasis.Standard, sk.StokesBasis.Solar, sk.StokesBasis.Observer]


def _rotation(angle, axis):
    axis = axis / np.linalg.norm(axis)
    k = np.array(
        [[0, -axis[2], axis[1]], [axis[2], 0, -axis[0]], [-axis[1], axis[0], 0]]
    )
    return np.eye(3) + np.sin(angle) * k + (1 - np.cos(angle)) * k @ k


def _perpendicular(v, look):
    p = v - v.dot(look) * look
    return p / np.linalg.norm(p)


def _tangent_altitude_solar_vectors(cos_sza, saa, relative_azimuth):
    """z, sun, look and observer used internally for a TangentAltitudeSolar ray"""
    z = np.array([0.0, 0.0, 1.0])
    sun = np.array(
        [
            np.sqrt(1 - cos_sza**2) * np.cos(saa),
            np.sqrt(1 - cos_sza**2) * np.sin(saa),
            cos_sza,
        ]
    )
    normal = np.cross(sun, z)
    if np.linalg.norm(normal) < 1e-12:
        normal = np.array([0.0, 1.0, 0.0])
    tangent = _rotation(np.arccos(TANGENT_COS_SZA), normal) @ sun
    tangent *= EARTH_RADIUS + TANGENT_ALTITUDE

    up = tangent / np.linalg.norm(tangent)
    look = _rotation(-relative_azimuth, up) @ _perpendicular(sun, up)

    s = np.sqrt(
        (EARTH_RADIUS + OBSERVER_ALTITUDE) ** 2 - (EARTH_RADIUS + TANGENT_ALTITUDE) ** 2
    )
    return z, sun, look, tangent - s * look


def _rayleigh_single_scatter_qu(reference, sun, look):
    """
    Normalized (Q, U) of single scattered Rayleigh light, which is polarized
    perpendicular to the scattering plane.  The reference is projected
    perpendicular to the line of sight, and the second axis is n x reference
    with n = -look the propagation direction.
    """
    n = -look
    p = _perpendicular(reference, look)
    q = np.cross(n, p)
    e = np.cross(n, sun)
    chi = np.arctan2(e.dot(q), e.dot(p))
    return np.array([np.cos(2 * chi), np.sin(2 * chi)])


def _calculate(cos_sza, saa, rays, basis, multiple_scatter):
    config = sk.Config()
    config.num_stokes = 3
    config.stokes_basis = basis
    if multiple_scatter:
        config.multiple_scatter_source = sk.MultipleScatterSource.DiscreteOrdinates
        config.num_streams = 4
    else:
        config.multiple_scatter_source = sk.MultipleScatterSource.NoSource

    geometry = sk.Geometry1D(
        cos_sza=cos_sza,
        solar_azimuth=saa,
        earth_radius_m=EARTH_RADIUS,
        altitude_grid_m=np.arange(0, 100001, 1000.0),
        interpolation_method=sk.InterpolationMethod.LinearInterpolation,
        geometry_type=sk.GeometryType.Spherical,
    )

    viewing_geo = sk.ViewingGeometry()
    for ray in rays:
        viewing_geo.add_ray(ray)

    atmosphere = sk.Atmosphere(geometry, config, wavelengths_nm=np.array([600.0]))
    sk.climatology.us76.add_us76_standard_atmosphere(atmosphere)
    atmosphere["rayleigh"] = sk.constituent.Rayleigh()
    atmosphere["surface"] = sk.constituent.LambertianSurface(0.3)

    radiance = sk.Engine(config, geometry, viewing_geo).calculate_radiance(atmosphere)

    return radiance["radiance"].isel(wavelength=0).transpose("los", "stokes").to_numpy()


@pytest.mark.parametrize(
    ("cos_sza", "saa"),
    [
        # Tangent point at the reference point, the observer is in the plane
        # of z and the look vector
        (TANGENT_COS_SZA, 0.0),
        # Tangent point away from the reference point
        (0.85, 0.7),
    ],
)
@pytest.mark.parametrize(
    ("basis", "reference_index"),
    # Index of the reference vector in _tangent_altitude_solar_vectors
    [
        (sk.StokesBasis.Standard, 0),
        (sk.StokesBasis.Solar, 1),
        (sk.StokesBasis.Observer, 3),
    ],
    ids=["standard", "solar", "observer"],
)
def test_single_scatter_stokes_basis(cos_sza, saa, basis, reference_index):
    """
    Single scattered Rayleigh light has a known polarization direction, check
    that every basis returns it for suns on both sides of the line of sight
    """
    rays = [
        sk.TangentAltitudeSolar(
            TANGENT_ALTITUDE, raz, OBSERVER_ALTITUDE, TANGENT_COS_SZA
        )
        for raz in RELATIVE_AZIMUTHS
    ]
    radiance = _calculate(cos_sza, saa, rays, basis, multiple_scatter=False)

    for i, raz in enumerate(RELATIVE_AZIMUTHS):
        vectors = _tangent_altitude_solar_vectors(cos_sza, saa, raz)
        _, sun, look, _ = vectors

        qu = radiance[i, 1:] / np.hypot(radiance[i, 1], radiance[i, 2])

        np.testing.assert_allclose(
            qu,
            _rayleigh_single_scatter_qu(vectors[reference_index], sun, look),
            atol=1e-10,
        )


def test_stokes_basis_rotated_scene():
    """
    Rotating the whole scene about the sun leaves it physically unchanged but
    moves the observer out of the plane of z and the look vector.  The solar
    and observer bases must not change, and the standard basis must change
    only by the rotation from the observer reference to z, including the
    multiple scatter contribution.
    """
    psi = np.deg2rad(35.0)

    rays = []
    rotated_rays = []
    standard_angles = []
    for raz in RELATIVE_AZIMUTHS:
        rays.append(
            sk.TangentAltitudeSolar(
                TANGENT_ALTITUDE, raz, OBSERVER_ALTITUDE, TANGENT_COS_SZA
            )
        )

        z, sun, look, observer = _tangent_altitude_solar_vectors(
            TANGENT_COS_SZA, 0.0, raz
        )
        rotation = _rotation(psi, sun)
        tangent = rotation @ z
        look = rotation @ look
        observer = rotation @ observer

        # Angles of the rotated tangent point, and the viewing azimuth in its
        # local frame, as defined by sk.TangentAltitude
        x_unit = np.array([1.0, 0.0, 0.0])
        y_unit = np.array([0.0, 1.0, 0.0])
        theta = np.arctan2(tangent[0], tangent[2])
        phi = np.arcsin(-tangent[1])
        primary = _rotation(theta, y_unit)
        to_tangent = _rotation(phi, primary @ x_unit) @ primary
        viewing_azimuth = np.arctan2(
            look.dot(to_tangent @ y_unit), look.dot(to_tangent @ x_unit)
        )
        rotated_rays.append(
            sk.TangentAltitude(
                TANGENT_ALTITUDE, OBSERVER_ALTITUDE, theta, viewing_azimuth, phi
            )
        )

        p_obs = _perpendicular(observer, look)
        p_z = _perpendicular(z, look)
        standard_angles.append(
            np.arctan2(-look.dot(np.cross(p_obs, p_z)), p_obs.dot(p_z))
        )

    for basis in BASES:
        original = _calculate(TANGENT_COS_SZA, 0.0, rays, basis, True)
        rotated = _calculate(TANGENT_COS_SZA, 0.0, rotated_rays, basis, True)

        if basis == sk.StokesBasis.Standard:
            # In the original scene the standard and observer references
            # coincide, turn them to z in the rotated scene
            c = np.cos(2 * np.array(standard_angles))
            s = np.sin(2 * np.array(standard_angles))
            original = np.stack(
                [
                    original[:, 0],
                    c * original[:, 1] + s * original[:, 2],
                    -s * original[:, 1] + c * original[:, 2],
                ],
                axis=1,
            )

        np.testing.assert_allclose(
            rotated, original, atol=1e-5 * np.abs(original[:, 0]).max()
        )
