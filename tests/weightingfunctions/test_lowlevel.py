import numpy as np
import sasktran2 as sk
from sasktran2.test_util.wf import validate_wf


def _raw_scenarios() -> list:
    """
    Defines the set of scenarios that we later on test dI/dk, dI/dssa, dI/dleg_coeff, dI/dalbedo on

    Right now we test for

    Configurations

    - Single scatter only
    - DO Source

    Viewing geometries

    - Limb viewing
    - Nadir viewing

    Model Geometry

    - 1D spherical

    Atmosphere

    - ficticious atmosphere with only Rayleigh scattering

    """
    configs = []

    configs.append(sk.Config())
    configs[-1].multiple_scatter_source = sk.MultipleScatterSource.DiscreteOrdinates
    configs[-1].num_streams = 4

    configs.append(sk.Config())

    geometrys = []
    altitude_grid = np.arange(0, 65001, 1000.0)

    geometrys.append(sk.Geometry1D(0.6, 0, 6327000, altitude_grid, sk.InterpolationMethod.LinearInterpolation, sk.GeometryType.Spherical))

    viewing_geos = []
    viewing_geos.append(sk.ViewingGeometry())

    for tan_alt in np.arange(10000, 60000, 2000):
        viewing_geos[-1].add_ray(sk.TangentAltitudeSolar(tan_alt, 0, 600000, 0.6))


    viewing_geos.append(sk.ViewingGeometry())
    viewing_geos[-1].add_ray(sk.GroundViewingSolar(0.6, 0, -0.8, 200000))

    scen = []

    for config in configs:
        for geometry in geometrys:
            for viewing_geo in viewing_geos:
                scen.append({'config': config,
                             'geometry': geometry,
                             'viewing_geo': viewing_geo,
                             'atmosphere': sk.test_util.scenarios.default_pure_scattering_atmosphere(config, geometry, 0.8)})

    return scen


def test_wf_extinction():
    """
    Verifies that the extinction derivative is accurate to 5 decimal places
    """
    D_FRACTION = 1e-4
    test_scens = _raw_scenarios()

    for scen in test_scens:
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

            numeric_wf[:, :, :, i] = (radiance_above['radiance'].values - radiance_below['radiance'].values) / (2*dk)

        radiance['wf_extinction_numeric'] = (['stokes', 'wavelength', 'los', 'altitude'], numeric_wf)
        validate_wf(radiance['wf_extinction'], radiance['wf_extinction_numeric'], decimal=5)


def test_wf_ssa():
    """
    Verifies that the single scatter albedo derivative is accurate to 6 decimal places
    """
    D_FRACTION = 1e-5
    test_scens = _raw_scenarios()

    for scen in test_scens:
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

            numeric_wf[:, :, :, i] = (radiance_above['radiance'].values - radiance_below['radiance'].values) / (2*d_ssa)

        radiance['wf_ssa_numeric'] = (['stokes', 'wavelength', 'los', 'altitude'], numeric_wf)
        validate_wf(radiance['wf_ssa'], radiance['wf_ssa_numeric'], decimal=6)
