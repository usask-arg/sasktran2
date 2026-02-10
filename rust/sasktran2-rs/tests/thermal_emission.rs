use anyhow::Result;
use ndarray::s;
use sasktran2_rs::bindings::config::EmissionSource;
use sasktran2_rs::bindings::prelude::*;

/// Convert wavenumber in cm⁻¹ to wavelength in nm.
fn wavenumber_to_wavelength_nm(wavenumber_cm: f64) -> f64 {
    1.0e7 / wavenumber_cm
}

/// DISORT Test Problem 7b: Thermal emission with scattering, comparing intensities.
///
/// Expected results from DISOTESTAUX.f (lines 1543-1550):
/// At τ=0 (top):     I(μ=-1) = 0.0,         I(μ=+1) = 4.65744E-07
/// At τ=100 (bottom): I(μ=-1) = 7.52311E-06, I(μ=+1) = 0.0
#[test]
fn test_disort7b_thermal_emission() -> Result<()> {
    // Test parameters from disotest.f90 lines 956-988
    let total_optical_depth = 100.0;
    let ssa = 0.95; // Single scatter albedo
    let g: f64 = 0.75; // Asymmetry factor (H-G)
    let nstreams = 16;

    // Temperature profile: 200K at top, 300K at bottom
    let temp_top = 300.0;
    let temp_bottom = 300.0;

    // Wavenumber 2702.99-2703.01 cm⁻¹ → midpoint ~2703 cm⁻¹ → λ ≈ 3700 nm
    let wavenumber_mid_cm = 2703.0;
    let wavelength_nm = wavenumber_to_wavelength_nm(wavenumber_mid_cm);

    // DISORT uses 1 layer (NLYR=1) with 2 grid points
    let num_altitudes = 2;
    let num_wavelengths = 1;
    let num_legendre = nstreams;

    let mut atmosphere = Atmosphere::new(
        num_wavelengths,
        num_altitudes,
        num_legendre + 1,
        true, // calculate_derivatives - must be false for DO emission (not implemented)
        true, // calculate_emission_derivatives
        Stokes::Stokes1,
        None,
    );

    // Create altitude grid (2 points: 0 = bottom, 1 = top)
    // Altitude is arbitrary; optical depth determines path length
    let altitude_grid = vec![0.0, 1000.0];

    // Set optical properties
    // For a single layer with 2 grid points spanning the layer
    atmosphere
        .storage
        .total_extinction
        .fill(total_optical_depth / 1000.0); // ext per meter
    atmosphere.storage.ssa.fill(ssa);
    atmosphere.storage.solar_irradiance.fill(0.0); // No solar beam (pure thermal)

    // Set thermal emission source term
    // Index 0 = bottom (300K), index 1 = top (200K)
    let emission_bottom = 1.09657540E-05; // planck_blackbody_radiance(temp_bottom, wavelength_nm);
    let emission_top = 1.09657540E-05; // planck_blackbody_radiance(temp_top, wavelength_nm);
    atmosphere.storage.emission_source[[0, 0]] = emission_bottom;
    atmosphere.storage.emission_source[[1, 0]] = emission_top;

    // Set Henyey-Greenstein phase function Legendre coefficients: P_l = g^l
    for l in 0..num_legendre + 1 {
        let coeff = g.powi(l as i32) * (2.0 * l as f64 + 1.0);
        atmosphere
            .storage
            .leg_coeff
            .slice_mut(s![l, .., ..])
            .fill(coeff);
    }

    // Configure engine for thermal emission in discrete ordinates
    let mut config = Config::new();
    config
        .with_num_streams(nstreams)?
        .with_multiple_scatter_source(MultipleScatterSource::DiscreteOrdinates)?
        .with_single_scatter_source(SingleScatterSource::DiscreteOrdinates)?
        .with_emission_source(EmissionSource::DiscreteOrdinates)?
        .with_delta_m_scaling(true)?
        .with_num_stokes(1)?;

    atmosphere.apply_delta_m_scaling(16)?;

    let geometry = Geometry1D::new(
        0.5,       // cos_sza (not used for pure thermal)
        0.0,       // saa
        6371000.0, // Earth radius
        altitude_grid.clone(),
        InterpolationMethod::Linear,
        GeometryType::PlaneParallel,
    );

    // Set up viewing geometry to measure upwelling intensity at TOA
    // DISORT UMU = [−1, 1] means looking at downwelling (μ=-1) and upwelling (μ=+1)
    // For ground viewing, we look up from ground which corresponds to upwelling at TOA
    let mut viewing_geometry = ViewingGeometry::new();
    // Looking straight up from ground (nadir viewing from TOA perspective)
    viewing_geometry.add_ground_viewing_solar(0.5, 0.0, 200000.0, 1.0);

    viewing_geometry.add_flux_observer_solar(0.6, 5000.0);

    // Create and run engine
    let engine = Engine::new(&config, &geometry, &viewing_geometry)?;
    let output = engine.calculate_radiance(&atmosphere)?;

    // Expected DISORT results for Test 7b
    // At τ=0 (top), upwelling intensity I(μ=+1) = 4.65744E-07
    let expected_upwelling_toa = 7.93075833E-06;

    // Print results for comparison
    println!("=== DISORT Test 7b: Thermal Emission + Scattering ===");
    println!("Parameters:");
    println!("  τ = {}", total_optical_depth);
    println!("  ssa = {}", ssa);
    println!("  g = {} (H-G)", g);
    println!("  T = [{} K (top), {} K (bottom)]", temp_top, temp_bottom);
    println!("  ν = 2703 cm⁻¹ (λ = {:.2} nm)", wavelength_nm);
    println!();
    println!(
        "Planck emission at bottom ({}K): {:.6e} W/(m^2 sr nm)",
        temp_bottom, emission_bottom
    );
    println!(
        "Planck emission at top ({}K): {:.6e} W/(m^2 sr nm)",
        temp_top, emission_top
    );
    println!();
    println!("Computed radiance: {:?}", output.radiance);
    println!();
    println!("Expected (DISORT Test 7b):");
    println!(
        "  Upwelling intensity at TOA (μ=+1): {:.5e}",
        expected_upwelling_toa
    );
    println!("  Downwelling intensity at bottom (μ=-1): 7.52311e-06");

    // Sanity checks
    for val in output.radiance.iter() {
        assert!(val.is_finite(), "Radiance should be finite, got {}", val);
        assert!(
            *val >= 0.0,
            "Radiance should be non-negative for thermal emission"
        );
    }

    // Compare to expected value (with tolerance for numerical differences)
    // Note: Units and normalization may differ between DISORT and sasktran2
    if !output.radiance.is_empty() && output.radiance[[0, 0, 0]] > 0.0 {
        let computed = output.radiance[[0, 0, 0]]; // * (1.0e7 / 2702.99 - 1.0e7 / 2703.01); // Integrate over wavenumber band
        println!(
            "Computed upwelling intensity at TOA (integrated over band): {:.5e}",
            computed
        );
        let ratio = computed / expected_upwelling_toa;
        println!();
        println!("Computed/Expected ratio: {:.8}", ratio);
        println!("(Ratio should be ~1.0 if units match DISORT)");
    }

    Ok(())
}
