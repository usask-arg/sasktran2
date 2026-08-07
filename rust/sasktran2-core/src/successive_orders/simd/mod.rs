//! Fixed-width wavelength batching for the successive-orders source.
//!
//! The scalar solver does not depend on this module. Packed values are kept in
//! one workspace so their memory lifetime and CXX-facing lifecycle remain
//! isolated from the ordinary per-wavelength state.

pub(super) mod batch;

#[cfg(not(test))]
use super::SolverDiagnostics;
#[cfg(not(test))]
use batch::LANES;

#[cfg(not(test))]
#[derive(Default)]
pub(super) struct Workspace {
    pub(super) batch4_transport_values: Vec<f64>,
    pub(super) batch4_scattering_coefficients: Vec<f64>,
    pub(super) batch4_boundary_scattering_values: Vec<f64>,
    pub(super) batch4_first_order_forcing: Vec<f64>,
    pub(super) batch4_initial_guess: Vec<f64>,
    pub(super) batch4_solution: Vec<f64>,
    pub(super) batch4_diagnostics: [SolverDiagnostics; LANES],
    pub(super) batch4_has_initial_guess: bool,
    pub(super) batch4_staged_lanes: usize,
    pub(super) batch4_transport_value_tangent: Vec<f64>,
    pub(super) batch4_scattering_coefficient_tangent: Vec<f64>,
    pub(super) batch4_boundary_scattering_value_tangent: Vec<f64>,
    pub(super) batch4_first_order_forcing_tangent: Vec<f64>,
    pub(super) batch4_solution_jvp: Vec<f64>,
    pub(super) batch4_jvp_staged_lanes: usize,
    pub(super) batch4_extinction_tangent: Vec<f64>,
    pub(super) batch4_single_scatter_albedo_tangent: Vec<f64>,
    pub(super) batch4_legendre_coefficient_tangent: Vec<f64>,
    pub(super) batch4_end_of_ray_source: Vec<f64>,
    pub(super) batch4_end_of_ray_source_tangent: Vec<f64>,
    pub(super) batch4_solution_cotangent: Vec<f64>,
    pub(super) batch4_scattering_coefficient_gradient: Vec<f64>,
    pub(super) batch4_boundary_scattering_value_gradient: Vec<f64>,
    pub(super) batch4_first_order_forcing_gradient: Vec<f64>,
    pub(super) batch4_vjp_staged_lanes: usize,
    pub(super) batch4_ray_vjp_ready: bool,
    pub(super) batch4_extinction_gradient: Vec<f64>,
    pub(super) batch4_single_scatter_albedo_gradient: Vec<f64>,
    pub(super) batch4_legendre_coefficient_gradient: Vec<f64>,
    pub(super) batch4_end_of_ray_source_gradient: Vec<f64>,
    pub(super) batch4_lane_solution: Vec<f64>,
}
