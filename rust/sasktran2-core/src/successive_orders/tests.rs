use super::*;

fn scalar_problem(config: SolverConfig) -> SuccessiveOrdersSolver {
    let transport = CsrMatrix::new(1, 1, &[0, 1], &[0]).unwrap();
    let scattering =
        BlockDiagonalMatrix::new(vec![0, 1], vec![0, 1], vec![0, 1], vec![1.0]).unwrap();
    let problem = FixedPointProblem::new(transport, scattering).unwrap();
    SuccessiveOrdersSolver::new(problem, config).unwrap()
}

fn solve_scalar(solver: &mut SuccessiveOrdersSolver, multiplier: f64, forcing: f64) -> f64 {
    solver
        .solve(&[0, 1], &[0], &[multiplier], &[forcing], None)
        .unwrap()[0]
}

#[test]
fn picard_converges_to_geometric_series() {
    let config = SolverConfig {
        max_iterations: 200,
        relative_tolerance: 1.0e-12,
        absolute_tolerance: 1.0e-14,
        ..SolverConfig::default()
    };
    let mut solver = scalar_problem(config);
    let solution = solve_scalar(&mut solver, 0.8, 2.0);
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
    let mut solver = scalar_problem(config);
    let solution = solve_scalar(&mut solver, 0.5, 1.0);
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
    let mut picard = scalar_problem(base);
    solve_scalar(&mut picard, 0.98, 1.0);

    let mut anderson = scalar_problem(SolverConfig {
        anderson_depth: 2,
        ..base
    });
    let solution = solve_scalar(&mut anderson, 0.98, 1.0);
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

#[test]
fn scalar_coefficient_scattering_matches_dense_legendre_kernel() {
    let incoming = [
        [1.0, 0.0, 0.0],
        [0.0, 1.0, 0.0],
        [0.0, 0.0, 1.0],
        normalize([1.0, -2.0, 0.5]),
    ];
    let incoming_weights = [0.1, 0.2, 0.3, 0.4];
    let outgoing = [
        normalize([0.2, 0.7, -0.4]),
        normalize([-0.5, 0.1, 0.8]),
        normalize([0.6, -0.3, 0.2]),
    ];
    let coefficients = [1.0, 0.3, -0.07, 0.02];
    let basis = ScalarCoefficientBasis::from_directions(
        &incoming,
        &incoming_weights,
        &outgoing,
        coefficients.len(),
    )
    .unwrap();
    let mut scattering = ScalarCoefficientScattering::new(
        vec![0, outgoing.len(), outgoing.len() + 1],
        vec![0, incoming.len(), incoming.len() + 2],
        1,
        basis,
    )
    .unwrap();
    scattering.set_coefficients(&coefficients).unwrap();
    scattering.set_dense_values(&[0.25, 0.75]).unwrap();

    let atmosphere_input = [2.0, -1.0, 0.5, 3.0];
    let input = [2.0, -1.0, 0.5, 3.0, 4.0, 8.0];
    let mut output = vec![0.0; outgoing.len() + 1];
    scattering.apply(&input, &mut output).unwrap();

    let mut expected = vec![0.0; outgoing.len() + 1];
    for (output_index, output_direction) in outgoing.iter().enumerate() {
        for (input_index, input_direction) in incoming.iter().enumerate() {
            let cosine = dot(*input_direction, *output_direction);
            let kernel = coefficients
                .iter()
                .enumerate()
                .map(|(degree, coefficient)| coefficient * legendre(degree, cosine))
                .sum::<f64>();
            expected[output_index] +=
                incoming_weights[input_index] * kernel * atmosphere_input[input_index];
        }
    }
    expected[outgoing.len()] = 0.25 * 4.0 + 0.75 * 8.0;

    for (actual, expected) in output.iter().zip(expected) {
        assert!(
            (actual - expected).abs() < 2.0e-13,
            "{actual} != {expected}"
        );
    }
}

fn normalize(direction: [f64; 3]) -> [f64; 3] {
    let norm = dot(direction, direction).sqrt();
    [
        direction[0] / norm,
        direction[1] / norm,
        direction[2] / norm,
    ]
}

fn dot(left: [f64; 3], right: [f64; 3]) -> f64 {
    left.iter()
        .zip(right)
        .map(|(left, right)| left * right)
        .sum()
}

fn legendre(degree: usize, cosine: f64) -> f64 {
    if degree == 0 {
        return 1.0;
    }
    if degree == 1 {
        return cosine;
    }
    let mut previous = 1.0;
    let mut current = cosine;
    for order in 2..=degree {
        let next = ((2 * order - 1) as f64 * cosine * current - (order - 1) as f64 * previous)
            / order as f64;
        previous = current;
        current = next;
    }
    current
}
