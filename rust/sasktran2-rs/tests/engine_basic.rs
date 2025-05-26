use ndarray::s;
use sasktran2_rs::bindings::prelude::*;

#[test]
fn test_engine_basic() -> Result<()> {
    let mut atmosphere = Atmosphere::new(10, 50, 16, true, true, Stokes::Stokes1);

    atmosphere.storage.ssa.fill(1.0);
    atmosphere.storage.total_extinction.fill(0.0001);
    atmosphere
        .storage
        .leg_coeff
        .slice_mut(s![0, .., ..])
        .fill(1.0);

    let mut altitude_grid = Vec::new();
    for i in 0..50 {
        altitude_grid.push(i as f64 * 1000.0);
    }

    let geometry = Geometry1D::new(
        0.6,
        0.0,
        6371000.0,
        altitude_grid,
        InterpolationMethod::Linear,
        GeometryType::Spherical,
    );

    let mut config = Config::new();
    let _ =
        config.with_emission_source(sasktran2_rs::bindings::config::EmissionSource::Standard)?;

    let mut viewing_geometry = ViewingGeometry::new();
    viewing_geometry.add_ground_viewing_solar(0.6, 0.0, 200000.0, 1.0);

    let engine = Engine::new(&config, &geometry, &viewing_geometry)?;

    let output = engine.calculate_radiance(&atmosphere).unwrap();
    println!("Radiance: {:?}", output.radiance);

    Ok(())
}
