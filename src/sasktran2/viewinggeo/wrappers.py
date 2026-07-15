from __future__ import annotations

from sasktran2._core_rust import (
    PyFluxObserverSolar,
    PyGroundViewingSolar,
    PySolarAnglesObserverLocation,
    PyTangentAltitude,
    PyTangentAltitudeSolar,
    PyViewingGeometry,
)


class ViewingGeometry:
    _viewing_geometry: PyViewingGeometry

    def __init__(self):
        self._viewing_geometry = PyViewingGeometry()
        self._rays = []
        self._flux_observers = []

    def add_ray(self, ray):
        self._rays.append(ray)
        self._viewing_geometry.add_ray(ray._internal)

    def add_flux_observer(self, observer):
        self._flux_observers.append(observer)
        self._viewing_geometry.add_flux_observer(observer._internal)

    @property
    def observer_rays(self):
        return self._rays

    @property
    def flux_observers(self):
        return self._flux_observers


class TangentAltitude:
    """A geometry-relative limb viewing ray defined at its tangent point.

    This policy specifies the line of sight relative to the coordinate system
    of the model geometry, rather than relative to the Sun. It is therefore the
    recommended tangent-ray policy for :class:`sasktran2.Geometry2D`.

    ``Geometry2D`` defines a reference plane spanned by its reference-z and
    reference-x directions. The horizontal grid angle is zero at reference-z
    and increases toward reference-x. ``horizontal_angle_radians`` places the
    tangent point on that same angular coordinate, so a value of zero centers
    the tangent point at the Geometry2D reference point.

    ``viewing_azimuth_radians`` sets the propagation direction of the line of
    sight at the tangent point. The propagation direction is the direction from
    the observer, through the tangent point, and toward the far side of the
    atmosphere. Its convention is:

    - ``0``: in the Geometry2D plane toward increasing horizontal angle;
    - ``pi``: in the Geometry2D plane toward decreasing horizontal angle;
    - ``+pi/2``: along the positive invariant/reference-y direction;
    - ``-pi/2``: along the negative invariant/reference-y direction.

    A viewing azimuth of ``+/-pi/2`` consequently follows a constant
    Geometry2D horizontal coordinate for a straight ray. Intermediate values
    combine along-plane and out-of-plane motion.

    Solar geometry is deliberately not specified here. The Sun is fixed by
    ``cos_sza`` and ``solar_azimuth`` on the model geometry. Changing those
    values changes illumination and scattering angles without moving this
    observing ray relative to the atmospheric grid.

    Parameters
    ----------
    tangent_altitude_m : float
        Geometric, unrefracted tangent altitude above the spherical surface in
        metres. It must be finite and non-negative.
    observer_altitude_m : float
        Observer altitude above the spherical surface in metres. It must be
        finite and no lower than ``tangent_altitude_m``.
    horizontal_angle_radians : float
        Angular position of the tangent point in the Geometry2D reference
        plane, in radians. Zero is the reference point and positive angles
        point toward reference-x.
    viewing_azimuth_radians : float
        Line-of-sight azimuth in the local tangent plane, in radians, following
        the convention above. Angles are periodic and are not restricted to a
        particular wrapped interval.

    Notes
    -----
    The policy is evaluated against the model coordinate system when an
    :class:`sasktran2.Engine` is constructed, so the Geometry2D object does not
    need to be passed to this constructor. The policy is also valid for a
    spherical :class:`sasktran2.Geometry1D`; its horizontal placement is then
    physically irrelevant because the atmosphere is horizontally uniform.

    The tangent altitude describes the straight-line launch geometry. If a
    refracting ray tracer is used, the actual refracted path need not retain
    the same tangent altitude.

    Examples
    --------
    Construct an in-plane ray centered 0.1 radians from the reference point::

        import numpy as np
        import sasktran2 as sk

        ray = sk.TangentAltitude(
            tangent_altitude_m=20_000.0,
            observer_altitude_m=200_000.0,
            horizontal_angle_radians=0.1,
            viewing_azimuth_radians=0.0,
        )

    Construct a ray at the reference point that travels through the invariant
    direction and therefore stays at horizontal angle zero::

        ray = sk.TangentAltitude(
            tangent_altitude_m=20_000.0,
            observer_altitude_m=200_000.0,
            horizontal_angle_radians=0.0,
            viewing_azimuth_radians=np.pi / 2,
        )
    """

    _internal: PyTangentAltitude

    def __init__(
        self,
        tangent_altitude_m: float,
        observer_altitude_m: float,
        horizontal_angle_radians: float,
        viewing_azimuth_radians: float,
    ):
        self._internal = PyTangentAltitude(
            tangent_altitude_m,
            observer_altitude_m,
            horizontal_angle_radians,
            viewing_azimuth_radians,
        )
        self._tangent_altitude_m = tangent_altitude_m
        self._observer_altitude_m = observer_altitude_m
        self._horizontal_angle_radians = horizontal_angle_radians
        self._viewing_azimuth_radians = viewing_azimuth_radians

    @property
    def tangent_altitude_m(self) -> float:
        """Geometric tangent altitude in metres."""
        return self._tangent_altitude_m

    @property
    def observer_altitude_m(self) -> float:
        """Observer altitude in metres."""
        return self._observer_altitude_m

    @property
    def horizontal_angle_radians(self) -> float:
        """Tangent-point angle in the Geometry2D reference plane."""
        return self._horizontal_angle_radians

    @property
    def viewing_azimuth_radians(self) -> float:
        """Line-of-sight azimuth in the local tangent plane."""
        return self._viewing_azimuth_radians

    def __repr__(self):
        return (
            "Tangent Viewing Ray: "
            f"tangent_altitude_m: {self._tangent_altitude_m}, "
            f"observer_altitude_m: {self._observer_altitude_m}, "
            f"horizontal_angle_radians: {self._horizontal_angle_radians}, "
            f"viewing_azimuth_radians: {self._viewing_azimuth_radians}"
        )


class TangentAltitudeSolar:
    _internal: PyTangentAltitudeSolar

    def __init__(
        self,
        tangent_altitude_m: float,
        relative_azimuth: float,
        observer_altitude_m: float,
        cos_sza: float,
    ):
        """
        Defines a viewing ray from the observer altitude, and tangent point parameters. Note that all of
        these parameters assume straight line paths (i.e. no atmospheric refraction)

        Parameters
        ----------
        tangent_altitude_m: float
            Tangent altitude in [m]
        relative_azimuth: float
            Relative azimuth angle to the sun. An angle of 0 degrees corresponds to the forward scattering plane. [rad]
        observer_altitude_m: float
            Observer altitude relative to the earth [m]
        cos_sza: float
            Cosine of the solar zenith angle at the tangent point [unitless]
        """
        self._internal = PyTangentAltitudeSolar(
            tangent_altitude_m, relative_azimuth, observer_altitude_m, cos_sza
        )

        self._tangent_altitude_m = tangent_altitude_m
        self._relative_azimuth = relative_azimuth
        self._observer_altitude_m = observer_altitude_m
        self._cos_sza = cos_sza

    def __repr__(self):
        # "Tangent Viewing Ray: tangentaltitude: {}, relative_azimuth_angle: "
        #    "{}, observeraltitude: {}, theta: {}, phi: {}",

        return f"Tangent Viewing Ray: tangentaltitude: {self._tangent_altitude_m}, relative_azimuth_angle: {self._relative_azimuth}, observeraltitude: {self._observer_altitude_m}, cos_sza: {self._cos_sza}"


class GroundViewingSolar:
    _internal: PyGroundViewingSolar

    def __init__(
        self,
        cos_sza: float,
        relative_azimuth: float,
        cos_viewing_zenith: float,
        observer_altitude_m: float,
    ):
        """
                Defines a viewing ray that is looking at the ground from angles defined at the ground location. Note that
        all of these parameters assumes straight line paths (i.e. no atmospheric refraction)

        Parameters
        ----------
        cos_sza: float
            Cosine of solar zenith angle at the ground point [unitless]
        relative_azimuth: float
            Relative azimuth angle to the sun [rad] at the ground point. An angle of 0 degrees corresponds to the forward scattering plane.
        observer_altitude_m: float
            Observer altitude relative to the earth [m]
        cos_viewing_zenith: float
            Cosine of the viewing zenith angle at the ground point [unitless]
        """
        self._internal = PyGroundViewingSolar(
            cos_sza, relative_azimuth, cos_viewing_zenith, observer_altitude_m
        )

        self._cos_sza = cos_sza
        self._relative_azimuth = relative_azimuth
        self._cos_viewing_zenith = cos_viewing_zenith
        self._observer_altitude_m = observer_altitude_m

    def __repr__(self):
        #            "Ground Viewing Ray: cos_sza: {}, relative_azimuth_angle: {}, "
        #    "cos_viewing_zenith: {}, observer_altitude: {}",
        return f"Ground Viewing Ray: cos_sza: {self._cos_sza}, relative_azimuth_angle: {self._relative_azimuth}, cos_viewing_zenith: {self._cos_viewing_zenith}, observer_altitude_m: {self._observer_altitude_m}"


class SolarAnglesObserverLocation:
    _internal: PySolarAnglesObserverLocation

    def __init__(
        self,
        cos_sza: float,
        relative_azimuth: float,
        cos_viewing_zenith: float,
        observer_altitude_m: float,
    ):
        """
        Defines a viewing ray that is defined at a location defined from the solar angles. Note that
        all of these parameters assumes straight line paths (i.e. no atmospheric refraction).
        This differs from sk.GroundViewingSolar in that the angles are defined at the observer location, not the ground location.

        Parameters
        ----------
        cos_sza: float
            Cosine of solar zenith angle at the observer point [unitless]
        relative_azimuth: float
            Relative azimuth angle to the sun [rad] at the observer point. An angle of 0 degrees corresponds to the forward scattering plane.
        cos_viewing_zenith: float
            Cosine of the viewing zenith angle at the observer point.  Positive angles are viewing up,
            negative angles are viewing down. [unitless]
        observer_altitude_m: float
            Observer altitude relative to the earth [m]
        """
        self._internal = PySolarAnglesObserverLocation(
            cos_sza, relative_azimuth, cos_viewing_zenith, observer_altitude_m
        )

        self._cos_sza = cos_sza
        self._relative_azimuth = relative_azimuth
        self._cos_viewing_zenith = cos_viewing_zenith
        self._observer_altitude_m = observer_altitude_m

    def __repr__(self):
        #            "Up Viewing Ray: cos_sza: {}, relative_azimuth_angle: {}, "
        #    "cos_viewing_zenith: {}, observer_altitude: {}",
        return f"Up Viewing Ray: cos_sza: {self._cos_sza}, relative_azimuth_angle: {self._relative_azimuth}, cos_viewing_zenith: {self._cos_viewing_zenith}, observer_altitude: {self._observer_altitude_m}"


class FluxObserverSolar:
    _internal: PyFluxObserverSolar

    def __init__(
        self,
        cos_sza: float,
        observer_altitude_m: float,
    ):
        """
        Defines a flux observer that is defined at a location defined from the solar angles.

        Parameters
        ----------
        cos_sza: float
            Cosine of solar zenith angle at the observer point [unitless]
        observer_altitude_m: float
            Observer altitude relative to the earth [m]
        """
        self._internal = PyFluxObserverSolar(cos_sza, observer_altitude_m)

        self._cos_sza = cos_sza
        self._observer_altitude_m = observer_altitude_m

    def __repr__(self):
        return f"Flux Observer: cos_sza: {self._cos_sza}, observer_altitude: {self._observer_altitude_m}"
