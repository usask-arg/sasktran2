use std::f64::consts::PI;

use anyhow::Result;
use ndarray::s;
use sasktran2_rs::bindings::prelude::*;

#[test]
fn test_coulson() -> Result<()> {
    let mut atmosphere = Atmosphere::new(1, 2, 40, true, false, Stokes::Stokes3);

    atmosphere.storage.ssa.fill(1.0);
    atmosphere.storage.total_extinction.fill(0.5);
    atmosphere.storage.solar_irradiance.fill(PI);

    // Set the rayleigh scattering legendre coefficients
    //    lephasef[0].a1 = 1;
    // lephasef[2].a2 = 6 * (1 - d) / (2 + d);
    // lephasef[2].a1 = (1 - d) / (2 + d);
    // lephasef[2].b1 = sqrt(6.0) * (1 - d) / (2 + d);
    atmosphere
        .storage
        .leg_coeff
        .slice_mut(s![0, .., ..])
        .fill(1.0);
    atmosphere
        .storage
        .leg_coeff
        .slice_mut(s![8, .., ..])
        .fill(1.0 / 2.0);
    atmosphere
        .storage
        .leg_coeff
        .slice_mut(s![9, .., ..])
        .fill(3.0);
    atmosphere
        .storage
        .leg_coeff
        .slice_mut(s![11, .., ..])
        .fill(-6.0_f64.sqrt() / 2.0);

    let mut binding = Config::new();
    let config = binding
        .with_num_streams(40)?
        .with_multiple_scatter_source(MultipleScatterSource::DiscreteOrdinates)?
        .with_single_scatter_source(SingleScatterSource::DiscreteOrdinates)?
        .with_num_stokes(3)
        .unwrap();

    let altitude_grid = vec![0.0, 1.0];

    let geometry = Geometry1D::new(
        0.2,
        0.0,
        6371000.0,
        altitude_grid,
        InterpolationMethod::Lower,
        GeometryType::PlaneParallel,
    );

    let mut viewing_geometry = ViewingGeometry::new();
    viewing_geometry.add_ground_viewing_solar(0.2, 0.0, 200000.0, 0.02);

    let engine = Engine::new(config, &geometry, &viewing_geometry)?;
    let output = engine.calculate_radiance(&atmosphere).unwrap();
    println!("Radiance: {:?}", output.radiance);

    Ok(())
}
