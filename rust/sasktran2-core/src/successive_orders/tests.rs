use super::*;

fn scalar_problem(multiplier: f64, forcing: f64, config: SolverConfig) -> SuccessiveOrdersSolver {
    let transport = CsrMatrix::new(1, 1, vec![0, 1], vec![0], vec![multiplier]).unwrap();
    let scattering =
        BlockDiagonalMatrix::new(vec![0, 1], vec![0, 1], vec![0, 1], vec![1.0]).unwrap();
    let problem = FixedPointProblem::new(transport, scattering, vec![forcing]).unwrap();
    SuccessiveOrdersSolver::new(problem, config).unwrap()
}

#[test]
fn picard_converges_to_geometric_series() {
    let config = SolverConfig {
        max_iterations: 200,
        relative_tolerance: 1.0e-12,
        absolute_tolerance: 1.0e-14,
        ..SolverConfig::default()
    };
    let mut solver = scalar_problem(0.8, 2.0, config);
    let solution = solver.solve(None).unwrap()[0];
    assert!((solution - 10.0).abs() < 1.0e-10);
    assert!(solver.diagnostics().converged);
    assert_eq!(solver.diagnostics().reason, ConvergenceReason::Tolerance);
    assert!(solver.diagnostics().final_residual < solver.diagnostics().initial_residual);
}

#[test]
fn fixed_iteration_mode_reports_maximum_iterations() {
    let config = SolverConfig {
        max_iterations: 3,
        relative_tolerance: 0.0,
        absolute_tolerance: 0.0,
        ..SolverConfig::default()
    };
    let mut solver = scalar_problem(0.5, 1.0, config);
    let solution = solver.solve(None).unwrap()[0];
    assert_eq!(solution, 1.75);
    assert!(!solver.diagnostics().converged);
    assert_eq!(solver.diagnostics().iterations, 3);
    assert_eq!(
        solver.diagnostics().reason,
        ConvergenceReason::MaximumIterations
    );
}

#[test]
fn anderson_accelerates_nearly_conservative_problem() {
    let base = SolverConfig {
        max_iterations: 200,
        relative_tolerance: 1.0e-10,
        absolute_tolerance: 1.0e-12,
        ..SolverConfig::default()
    };
    let mut picard = scalar_problem(0.98, 1.0, base);
    picard.solve(None).unwrap();

    let mut anderson = scalar_problem(
        0.98,
        1.0,
        SolverConfig {
            anderson_depth: 2,
            ..base
        },
    );
    let solution = anderson.solve(None).unwrap()[0];
    assert!((solution - 50.0).abs() < 1.0e-7);
    assert!(anderson.diagnostics().converged);
    assert!(anderson.diagnostics().iterations < picard.diagnostics().iterations);
}

#[test]
fn block_diagonal_operator_supports_rectangular_blocks() {
    let matrix = BlockDiagonalMatrix::new(
        vec![0, 2, 3],
        vec![0, 1, 3],
        vec![0, 2, 4],
        vec![2.0, 3.0, 4.0, -1.0],
    )
    .unwrap();
    let mut output = vec![0.0; 3];
    matrix.apply(&[5.0, 2.0, 7.0], &mut output).unwrap();
    assert_eq!(output, vec![10.0, 15.0, 1.0]);
}

#[test]
fn compressed_stencils_accept_1d_and_2d_widths() {
    let stencils = CompressedStencils::new(
        5,
        vec![0, 2, 6],
        vec![0, 1, 1, 2, 3, 4],
        vec![0.25, 0.75, 0.1, 0.2, 0.3, 0.4],
    )
    .unwrap();
    let mut output = vec![0.0; 2];
    stencils
        .interpolate(&[2.0, 4.0, 6.0, 8.0, 10.0], &mut output)
        .unwrap();
    assert_eq!(output, vec![3.5, 8.0]);
}
