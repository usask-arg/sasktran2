from __future__ import annotations

from copy import copy
from typing import ClassVar

import numpy as np
import pandas as pd
import sasktran2 as sk
import xarray as xr

from .fer_format import FrontEndRadiance


class FERGenerator:
    _optical_property_fns: ClassVar[dict[str, callable]] = {}
    _climatology_fns: ClassVar[dict[str, callable]] = {}
    _viewing_policy_fns: ClassVar[dict[str, callable]] = {}
    _atmospheric_state_fns: ClassVar[dict[str, callable]] = {}

    @classmethod
    def register_optical_property(cls, species_name: str):
        def decorator(optical_property_fn: callable):
            cls._optical_property_fns[species_name] = optical_property_fn
            return optical_property_fn

        return decorator

    @classmethod
    def register_climatology(cls, type: str):
        def decorator(climatology_fn: callable):
            cls._climatology_fns[type] = climatology_fn
            return climatology_fn

        return decorator

    @classmethod
    def register_viewing_policy(cls, type: str):
        def decorator(viewing_policy_fn: callable):
            cls._viewing_policy_fns[type] = viewing_policy_fn
            return viewing_policy_fn

        return decorator

    @classmethod
    def register_atmospheric_state(cls, type: str):
        def decorator(atmospheric_state_fn: callable):
            cls._atmospheric_state_fns[type] = atmospheric_state_fn
            return atmospheric_state_fn

        return decorator

    def __init__(self, cfg: dict | None = None):
        self._cfg = cfg if cfg is not None else {}

    @staticmethod
    def _convert_timestamp(time_value: dict | str | pd.Timestamp) -> pd.Timestamp:
        if isinstance(time_value, str):
            # Convert string to Timestamp
            return pd.Timestamp(time_value)
        if isinstance(time_value, dict):
            # Convert dictionary to Timestamp by unpacking the keys
            return pd.Timestamp(**time_value)
        if isinstance(time_value, pd.Timestamp):
            # Already a Timestamp, no conversion needed
            return time_value
        msg = "Unsupported time format"
        raise ValueError(msg)

    def _construct_atmosphere(
        self, model_geometry: sk.Geometry1D, sk_config: sk.Config, atmosphere: dict
    ) -> sk.Atmosphere:
        atmo = sk.Atmosphere(
            model_geometry=model_geometry,
            config=sk_config,
            wavelengths_nm=self._construct_model_wavelength(),
            calculate_derivatives=False,
        )

        self._atmospheric_state_fns[atmosphere["atmospheric_state"]["type"]](
            self, atmo, **atmosphere["atmospheric_state"]
        )

        if atmosphere["atmospheric_state"]["include_rayleigh_scattering"]:
            atmo["rayleigh"] = sk.constituent.Rayleigh()

        for name, cfg in atmosphere.get("constituents", {}).items():
            if cfg["type"] in self._climatology_fns:
                atmo[name] = self._climatology_fns[cfg["type"]](self, name, **cfg)

        return atmo

    def _construct_solar_handler(
        self, acquisition: dict
    ) -> sk.solar.SolarGeometryHandlerBase:
        solar_options = acquisition["solar_geometry"]

        if solar_options["type"] == "forced":
            return sk.solar.SolarHandlerForced(
                solar_options["solar_zenith"], solar_options["solar_azimuth"]
            )

        msg = f"Unknown solar geometry type: {solar_options['type']}"
        raise ValueError(msg)

    def _construct_viewing_geometry(
        self, acquisition: dict, solar_handler: sk.solar.SolarGeometryHandlerBase
    ) -> sk.ViewingGeometry:
        cfg = copy(acquisition)
        cfg["time"] = self._convert_timestamp(cfg["time"])

        return self._viewing_policy_fns[acquisition["type"]](self, solar_handler, **cfg)

    def _construct_model_geometry(
        self, acquisition: dict  # noqa: ARG002
    ) -> sk.Geometry1D:
        return sk.Geometry1D(
            0.6,
            0,
            6372000,
            np.arange(0, 60000.0, 1000.0),
            sk.InterpolationMethod.LinearInterpolation,
            sk.GeometryType.Spherical,
        )

    def _construct_sk2_config(self, acquisition: dict) -> sk.Config:  # noqa: ARG002
        return sk.Config()

    def _construct_model_wavelength(self) -> np.array:
        return np.array([745])

    def generate(self, acquistions: dict, atmosphere: dict):
        results = {}
        for acquisition, cfg in acquistions.items():
            solar_handler = self._construct_solar_handler(cfg)

            sk2_view, geometry_ds = self._construct_viewing_geometry(cfg, solar_handler)
            model = self._construct_model_geometry(cfg)
            config = self._construct_sk2_config(cfg)
            atmo = self._construct_atmosphere(model, config, atmosphere)

            engine = sk.Engine(config, model, sk2_view)

            # Convert to FER format
            results[acquisition] = FrontEndRadiance(
                engine.calculate_radiance(atmo), geometry_ds
            )

        return results


# Register all the default optical properties
@FERGenerator.register_optical_property("o3")
def o3_optical_property(*args, **kwargs):
    return sk.optical.O3DBM()


@FERGenerator.register_optical_property("no2")
def no2_optical_property(*args, **kwargs):
    return sk.optical.NO2Vandaele()


@FERGenerator.register_optical_property("bro")
def bro_optical_property(*args, **kwargs):
    return sk.optical.HITRANUV("BrO")


@FERGenerator.register_optical_property("so2")
def so2_optical_property(*args, **kwargs):
    return sk.optical.HITRANUV("SO2")


@FERGenerator.register_climatology("mipas_climatology")
def mipas_climatology(cls, name: str, *args, **kwargs):
    return sk.climatology.mipas.constituent(name, cls._optical_property_fns[name]())


@FERGenerator.register_atmospheric_state("us76")
def us76_atmospheric_state(cls, atmo, *args, **kwargs):  # noqa: ARG001
    sk.climatology.us76.add_us76_standard_atmosphere(atmo)


@FERGenerator.register_viewing_policy("limb_vertical_image")
def limb_vertical_image(
    cls,  # noqa: ARG001
    solar_handler,
    tangent_altitudes,
    tangent_latitude,
    tangent_longitude,
    time,
    inst_loc,
    *args,
    **kwargs,
):
    viewing_geo = sk.ViewingGeometry()

    for alt in tangent_altitudes:
        solar_zenith, solar_azimuth = solar_handler.target_solar_angles(
            tangent_latitude, tangent_longitude, alt, time
        )

        viewing_geo.add_ray(
            sk.TangentAltitudeSolar(
                tangent_altitude_m=alt,
                relative_azimuth=np.deg2rad(solar_azimuth),
                observer_altitude_m=inst_loc["altitude"],
                cos_sza=np.cos(np.deg2rad(solar_zenith)),
            )
        )

    geometry_ds = xr.Dataset(
        {
            "tangent_altitude": (["los"], tangent_altitudes),
            "tangent_longitude": (
                ["los"],
                np.ones_like(tangent_altitudes) * tangent_longitude,
            ),
            "tangent_latitude": (
                ["los"],
                np.ones_like(tangent_altitudes) * tangent_latitude,
            ),
        },
    )
    return viewing_geo, geometry_ds
