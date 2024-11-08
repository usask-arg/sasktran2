import numpy as np
import xarray as xr

import sasktran2 as sk


def limb_vertical_image(
    solar_handler,
    tangent_altitudes,
    tangent_latitude,
    tangent_longitude,
    time,
    inst_loc,
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


class LimbVerticalImage(sk.ViewingGeometry, sk.viewinggeo.ViewingGeometryContainer):
    def __init__(
        self,
        solar_handler,
        tangent_altitudes,
        tangent_latitude,
        tangent_longitude,
        time,
        inst_loc,
    ):
        sk.ViewingGeometry.__init__(self)
        for alt in tangent_altitudes:
            solar_zenith, solar_azimuth = solar_handler.target_solar_angles(
                tangent_latitude, tangent_longitude, alt, time
            )

            self.add_ray(
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

        sk.viewinggeo.ViewingGeometryContainer.__init__(self, geometry_ds)
