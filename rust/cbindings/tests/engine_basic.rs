use sasktran2_bindings::prelude::*;
use ndarray::s;

#[test]
fn test_engine_basic() {
    let mut atmosphere = Atmosphere::new(10, 50, 16, 0, false);

    atmosphere.storage.ssa.fill(1.0);
    atmosphere.storage.total_extinction.fill(0.1);
    atmosphere.storage.leg_coeff.slice_mut(s![0, .., ..]).fill(1.0);

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

    let config = Config::new();

    let mut viewing_geometry = ViewingGeometry::new();
    viewing_geometry.add_ground_viewing_solar(
        0.5,
        0.0,
        200000.0,
        -0.99,
    );

    let engine = Engine::new(&config, &geometry, &viewing_geometry);

    let output = engine.calculate_radiance(&atmosphere).unwrap();
    println!("Radiance: {:?}", output.radiance);

}
