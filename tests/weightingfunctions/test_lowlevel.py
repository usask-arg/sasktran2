import sasktran2 as sk
import numpy as np
from sasktran2.test_util.wf import validate_wf


def _raw_scenarios() -> dict:
    config = sk.Config()

    altitude_grid = np.arange(0, 65001, 1000.0)
    geometry = sk.Geometry1D(0.6, 0, 6327000, altitude_grid, sk.InterpolationMethod.LinearInterpolation, sk.GeometryType.Spherical)

    viewing_geo = sk.ViewingGeometry()

    for tan_alt in np.arange(10000, 60000, 2000):
        viewing_geo.add_ray(sk.TangentAltitudeSolar(tan_alt, 0, 600000, 0.6))

    scen = {}

    scen['default'] = {'config': config,
                        'geometry': geometry,
                        'viewing_geo': viewing_geo,
                        'atmosphere': sk.test_util.scenarios.default_pure_scattering_atmosphere(config, geometry, 0.8)}

    return scen


def test_wf_extinction():
    D_FRACTION = 1e-5
    test_scens = _raw_scenarios()

    for name, scen in test_scens.items():
        engine = sk.Engine(scen['config'], scen['geometry'], scen['viewing_geo'])


        radiance = engine.calculate_radiance(scen['atmosphere'])

        numeric_wf = np.zeros_like(radiance['wf_extinction'].values)

        for i in range(len(radiance['altitude'])):
            dk = D_FRACTION * scen['atmosphere'].storage.total_extinction[i]

            scen['atmosphere'].storage.total_extinction[i] += dk
            radiance_above = engine.calculate_radiance(scen['atmosphere'])

            scen['atmosphere'].storage.total_extinction[i] -= 2*dk
            radiance_below = engine.calculate_radiance(scen['atmosphere'])

            scen['atmosphere'].storage.total_extinction[i] += dk

            numeric_wf[:, :, i, :] = (radiance_above['radiance'].values - radiance_below['radiance'].values) / (2*dk)

        radiance['wf_extinction_numeric'] = (['stokes', 'los', 'altitude', 'wavelength'], numeric_wf)
        validate_wf(radiance['wf_extinction'], radiance['wf_extinction_numeric'], decimal=5)


def test_wf_ssa():
    D_FRACTION = 1e-5
    test_scens = _raw_scenarios()

    for name, scen in test_scens.items():
        engine = sk.Engine(scen['config'], scen['geometry'], scen['viewing_geo'])

        radiance = engine.calculate_radiance(scen['atmosphere'])

        numeric_wf = np.zeros_like(radiance['wf_ssa'].values)

        for i in range(len(radiance['altitude'])):
            d_ssa = D_FRACTION * scen['atmosphere'].storage.ssa[i]

            scen['atmosphere'].storage.ssa[i] += d_ssa
            radiance_above = engine.calculate_radiance(scen['atmosphere'])

            scen['atmosphere'].storage.ssa[i] -= 2*d_ssa
            radiance_below = engine.calculate_radiance(scen['atmosphere'])

            scen['atmosphere'].storage.ssa[i] += d_ssa

            numeric_wf[:, :, i, :] = (radiance_above['radiance'].values - radiance_below['radiance'].values) / (2*d_ssa)

        radiance['wf_ssa_numeric'] = (['stokes', 'los', 'altitude', 'wavelength'], numeric_wf)
        validate_wf(radiance['wf_ssa'], radiance['wf_ssa_numeric'], decimal=6)
