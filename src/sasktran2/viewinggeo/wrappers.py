from __future__ import annotations

from sasktran2._core_rust import (
    PyGroundViewingSolar,
    PySolarAnglesObserverLocation,
    PyTangentAltitudeSolar,
    PyViewingGeometry,
)


class ViewingGeometry:
    _viewing_geometry: PyViewingGeometry

    def __init__(self):
        self._viewing_geometry = PyViewingGeometry()
        self._rays = []

    def add_ray(self, ray):
        self._rays.append(ray)
        self._viewing_geometry.add_ray(ray._internal)

    @property
    def observer_rays(self):
        return self._rays


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
