use super::*;
use crate::math::wigner::WignerDCalculator;

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
fn fixed_point_jvp_and_vjp_match_finite_difference_and_duality() {
    let config = SolverConfig {
        max_iterations: 4,
        relative_tolerance: 0.0,
        absolute_tolerance: 0.0,
        ..SolverConfig::default()
    };
    let multiplier = 0.42;
    let forcing = 1.3;
    let transport_tangent = [0.17];
    let scattering_tangent = [-0.08];
    let forcing_tangent = [0.23];
    let mut solver = scalar_problem(config);
    solve_scalar(&mut solver, multiplier, forcing);
    let jvp = solver
        .solve_jvp(
            &[0, 1],
            &[0],
            &[multiplier],
            &transport_tangent,
            &[forcing],
            &forcing_tangent,
            &[],
            &scattering_tangent,
        )
        .unwrap();

    let epsilon = 1.0e-6;
    let mut plus = scalar_problem(config);
    plus.problem_mut()
        .set_scattering_values(&[1.0 + epsilon * scattering_tangent[0]])
        .unwrap();
    let plus_value = solve_scalar(
        &mut plus,
        multiplier + epsilon * transport_tangent[0],
        forcing + epsilon * forcing_tangent[0],
    );
    let mut minus = scalar_problem(config);
    minus
        .problem_mut()
        .set_scattering_values(&[1.0 - epsilon * scattering_tangent[0]])
        .unwrap();
    let minus_value = solve_scalar(
        &mut minus,
        multiplier - epsilon * transport_tangent[0],
        forcing - epsilon * forcing_tangent[0],
    );
    let numeric = (plus_value - minus_value) / (2.0 * epsilon);
    assert!((jvp[0] - numeric).abs() < 2.0e-9, "{} != {numeric}", jvp[0]);

    let solution_cotangent = [1.7];
    let gradient = solver
        .solve_vjp(
            &[0, 1],
            &[0],
            &[multiplier],
            &[forcing],
            &solution_cotangent,
        )
        .unwrap();
    let tangent_product = solution_cotangent[0] * jvp[0];
    let gradient_product = gradient.transport_values[0] * transport_tangent[0]
        + gradient.dense_scattering_values[0] * scattering_tangent[0]
        + gradient.forcing[0] * forcing_tangent[0];
    assert!((tangent_product - gradient_product).abs() < 2.0e-12);
}

#[test]
fn implicit_products_match_converged_system_with_anderson() {
    let config = SolverConfig {
        max_iterations: 200,
        relative_tolerance: 1.0e-12,
        absolute_tolerance: 1.0e-14,
        anderson_depth: 2,
        ..SolverConfig::default()
    };
    let multiplier = 0.8;
    let forcing = 2.0;
    let transport_tangent = [0.17];
    let scattering_tangent = [-0.08];
    let forcing_tangent = [0.23];
    let mut solver = scalar_problem(config);
    let solution = solve_scalar(&mut solver, multiplier, forcing);
    assert!(solver.diagnostics().converged);

    let jvp = solver
        .solve_jvp(
            &[0, 1],
            &[0],
            &[multiplier],
            &transport_tangent,
            &[forcing],
            &forcing_tangent,
            &[],
            &scattering_tangent,
        )
        .unwrap();
    let right_hand_side =
        scattering_tangent[0] * solution + forcing_tangent[0] + transport_tangent[0] * solution;
    let expected_jvp = right_hand_side / (1.0 - multiplier);
    assert!((jvp[0] - expected_jvp).abs() < 2.0e-10);

    let solution_cotangent = [1.7];
    let gradient = solver
        .solve_vjp(
            &[0, 1],
            &[0],
            &[multiplier],
            &[forcing],
            &solution_cotangent,
        )
        .unwrap();
    let adjoint = solution_cotangent[0] / (1.0 - multiplier);
    assert!((gradient.transport_values[0] - adjoint * solution).abs() < 2.0e-10);
    assert!((gradient.dense_scattering_values[0] - adjoint * solution).abs() < 2.0e-10);
    assert!((gradient.forcing[0] - adjoint).abs() < 2.0e-10);
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

    let input_tangent = [0.2, -0.1, 0.04, 0.3, -0.5, 0.7];
    let coefficient_tangent = [0.03, -0.02, 0.01, 0.04];
    let dense_tangent = [-0.06, 0.09];
    let mut output_tangent = vec![0.0; output.len()];
    scattering
        .apply_jvp(
            &input,
            &input_tangent,
            &coefficient_tangent,
            &dense_tangent,
            &mut output_tangent,
        )
        .unwrap();
    let output_cotangent = [0.4, -0.7, 0.2, 0.9];
    let mut input_cotangent = vec![0.0; input.len()];
    let mut coefficient_gradient = vec![0.0; coefficients.len()];
    let mut dense_gradient = vec![0.0; 2];
    scattering
        .apply_vjp(
            &input,
            &output_cotangent,
            &mut input_cotangent,
            &mut coefficient_gradient,
            &mut dense_gradient,
        )
        .unwrap();
    let mut transpose_input = vec![0.0; input.len()];
    scattering
        .apply_transpose(&output_cotangent, &mut transpose_input)
        .unwrap();
    for (transpose, vjp) in transpose_input.iter().zip(&input_cotangent) {
        assert!((transpose - vjp).abs() < 2.0e-13);
    }
    let forward_product = inner_product(&output_tangent, &output_cotangent);
    let reverse_product = inner_product(&input_tangent, &input_cotangent)
        + inner_product(&coefficient_tangent, &coefficient_gradient)
        + inner_product(&dense_tangent, &dense_gradient);
    assert!((forward_product - reverse_product).abs() < 3.0e-13);
}

#[test]
fn vector_coefficient_scattering_matches_dense_greek_kernel() {
    let incoming = [
        normalize([1.0, 0.2, 0.4]),
        [0.0, 0.0, 1.0],
        [0.0, 0.0, -1.0],
        normalize([-0.8, -0.2, 0.3]),
    ];
    let incoming_weights = [0.17, 0.23, 0.29, 0.31];
    let outgoing = [incoming[0], [0.0, 0.0, 1.0], normalize([0.6, -0.3, 0.2])];
    let num_coefficients = 5;
    let coefficients = [
        1.0, 0.0, 0.0, 0.0, 0.21, -0.04, 0.03, 0.02, -0.13, 0.34, -0.08, 0.11, 0.07, -0.12, 0.19,
        -0.05, -0.02, 0.06, -0.09, 0.04,
    ];
    let input = [
        2.0, -0.3, 0.2, 0.7, 0.1, -0.4, 1.2, -0.2, -0.1, 0.4, 0.3, 0.15,
    ];

    for (output_index, &output_direction) in outgoing.iter().enumerate() {
        for (input_index, &input_direction) in incoming.iter().enumerate() {
            let pair_basis = VectorCoefficientBasis::from_directions(
                &[input_direction],
                &[incoming_weights[input_index]],
                &[output_direction],
                num_coefficients,
            )
            .unwrap();
            let mut pair_scattering =
                VectorCoefficientScattering::new(vec![0, 3], vec![0, 3], 1, pair_basis).unwrap();
            pair_scattering.set_coefficients(&coefficients).unwrap();
            let mut pair_output = [0.0; 3];
            pair_scattering
                .apply(
                    &input[input_index * 3..input_index * 3 + 3],
                    &mut pair_output,
                )
                .unwrap();
            let mut pair_expected = [0.0; 3];
            for degree in 0..num_coefficients {
                let matrix = polarized_greek_matrix(
                    input_direction,
                    output_direction,
                    degree,
                    &coefficients[degree * 4..degree * 4 + 4],
                );
                for row in 0..3 {
                    for column in 0..3 {
                        pair_expected[row] += incoming_weights[input_index]
                            * matrix[row][column]
                            * input[input_index * 3 + column];
                    }
                }
            }
            for row in 0..3 {
                assert!(
                    (pair_output[row] - pair_expected[row]).abs() < 3.0e-12,
                    "pair ({input_index}, {output_index}) row {row}: {} != {}",
                    pair_output[row],
                    pair_expected[row]
                );
            }
        }
    }

    let basis = VectorCoefficientBasis::from_directions(
        &incoming,
        &incoming_weights,
        &outgoing,
        num_coefficients,
    )
    .unwrap();
    let mut scattering = VectorCoefficientScattering::new(
        vec![0, outgoing.len() * 3],
        vec![0, incoming.len() * 3],
        1,
        basis,
    )
    .unwrap();
    scattering.set_coefficients(&coefficients).unwrap();
    let mut output = vec![0.0; outgoing.len() * 3];
    scattering.apply(&input, &mut output).unwrap();

    let mut expected = vec![0.0; outgoing.len() * 3];
    for (output_index, output_direction) in outgoing.iter().enumerate() {
        for (input_index, input_direction) in incoming.iter().enumerate() {
            for degree in 0..num_coefficients {
                let matrix = polarized_greek_matrix(
                    *input_direction,
                    *output_direction,
                    degree,
                    &coefficients[degree * 4..degree * 4 + 4],
                );
                for row in 0..3 {
                    for column in 0..3 {
                        expected[output_index * 3 + row] += incoming_weights[input_index]
                            * matrix[row][column]
                            * input[input_index * 3 + column];
                    }
                }
            }
        }
    }

    for (actual, expected) in output.iter().zip(expected) {
        assert!(
            (actual - expected).abs() < 3.0e-12,
            "{actual} != {expected}"
        );
    }

    let input_tangent = [
        0.2, -0.1, 0.04, -0.3, 0.08, 0.17, 0.14, -0.22, 0.09, 0.31, -0.07, 0.11,
    ];
    let coefficient_tangent = [
        0.01, -0.03, 0.02, 0.04, -0.02, 0.05, -0.01, 0.03, 0.04, -0.02, 0.06, -0.01, 0.03, 0.02,
        -0.04, 0.01, -0.05, 0.03, 0.02, -0.02,
    ];
    let mut output_tangent = vec![0.0; output.len()];
    scattering
        .apply_jvp(
            &input,
            &input_tangent,
            &coefficient_tangent,
            &[],
            &mut output_tangent,
        )
        .unwrap();
    let output_cotangent = [0.4, -0.7, 0.2, 0.9, -0.3, 0.5, -0.2, 0.6, 0.8];
    let mut input_cotangent = vec![0.0; input.len()];
    let mut coefficient_gradient = vec![0.0; coefficients.len()];
    scattering
        .apply_vjp(
            &input,
            &output_cotangent,
            &mut input_cotangent,
            &mut coefficient_gradient,
            &mut [],
        )
        .unwrap();
    let mut transpose_input = vec![0.0; input.len()];
    scattering
        .apply_transpose(&output_cotangent, &mut transpose_input)
        .unwrap();
    for (transpose, vjp) in transpose_input.iter().zip(&input_cotangent) {
        assert!((transpose - vjp).abs() < 2.0e-12);
    }
    let forward_product = inner_product(&output_tangent, &output_cotangent);
    let reverse_product = inner_product(&input_tangent, &input_cotangent)
        + inner_product(&coefficient_tangent, &coefficient_gradient);
    assert!((forward_product - reverse_product).abs() < 2.0e-11);
}

fn polarized_greek_matrix(
    incoming: [f64; 3],
    outgoing: [f64; 3],
    degree: usize,
    coefficients: &[f64],
) -> [[f64; 3]; 3] {
    let (theta, c1, c2, s1, s2) = stokes_scattering_factors(negate(incoming), negate(outgoing));
    let d00 = WignerDCalculator::new(0, 0).d(theta, degree as i32);
    let d22 = WignerDCalculator::new(2, 2).d(theta, degree as i32);
    let d02 = WignerDCalculator::new(0, 2).d(theta, degree as i32);
    let d2m2 = WignerDCalculator::new(2, -2).d(theta, degree as i32);
    let w_add = d22 + d2m2;
    let w_minus = d22 - d2m2;
    let [a1, a2, a3, b1] = coefficients.try_into().unwrap();
    let mut result = [[0.0; 3]; 3];
    result[0][0] = a1 * d00;
    result[1][0] = -b1 * c2 * d02;
    result[2][0] = -b1 * s2 * d02;
    result[0][1] = -b1 * c1 * d02;
    result[0][2] = b1 * s1 * d02;
    result[1][1] = 0.5
        * (a2 * (c1 * c2 * w_add - s1 * s2 * w_minus) + a3 * (c1 * c2 * w_minus - s1 * s2 * w_add));
    result[2][1] = 0.5
        * (a2 * (c1 * s2 * w_add + s1 * c2 * w_minus) + a3 * (c1 * s2 * w_minus + s1 * c2 * w_add));
    result[1][2] = -0.5
        * (a2 * (s1 * c2 * w_add + c1 * s2 * w_minus) + a3 * (s1 * c2 * w_minus + c1 * s2 * w_add));
    result[2][2] = 0.5
        * (a2 * (-s1 * s2 * w_add + c1 * c2 * w_minus)
            + a3 * (-s1 * s2 * w_minus + c1 * c2 * w_add));
    result
}

fn stokes_scattering_factors(incoming: [f64; 3], outgoing: [f64; 3]) -> (f64, f64, f64, f64, f64) {
    let cosine = dot(incoming, outgoing).clamp(-1.0, 1.0);
    let theta = cosine.acos();
    let sin_scatter = theta.sin();
    if sin_scatter.abs() < 1.0e-8 {
        return (theta, 1.0, 1.0, 0.0, 0.0);
    }
    let cos_incoming = incoming[2].clamp(-1.0, 1.0);
    let cos_outgoing = outgoing[2].clamp(-1.0, 1.0);
    let sin_incoming = cos_incoming.acos().sin();
    let sin_outgoing = cos_outgoing.acos().sin();
    if sin_incoming.abs() < 1.0e-8 || sin_outgoing.abs() < 1.0e-8 {
        return (theta, 1.0, 1.0, 0.0, 0.0);
    }
    let cosine1 =
        ((cos_outgoing - cos_incoming * cosine) / (sin_incoming * sin_scatter)).clamp(-1.0, 1.0);
    let cosine2 =
        ((cos_incoming - cos_outgoing * cosine) / (sin_outgoing * sin_scatter)).clamp(-1.0, 1.0);
    let mut sine1 = 2.0 * (1.0 - cosine1 * cosine1).sqrt() * cosine1;
    let mut sine2 = 2.0 * (1.0 - cosine2 * cosine2).sqrt() * cosine2;
    let phi_incoming = incoming[1].atan2(incoming[0]);
    let phi_outgoing = outgoing[1].atan2(outgoing[0]);
    if (phi_incoming - phi_outgoing).rem_euclid(2.0 * std::f64::consts::PI) < std::f64::consts::PI {
        sine1 *= -1.0;
        sine2 *= -1.0;
    }
    (
        theta,
        2.0 * cosine1 * cosine1 - 1.0,
        2.0 * cosine2 * cosine2 - 1.0,
        sine1,
        sine2,
    )
}

fn negate(direction: [f64; 3]) -> [f64; 3] {
    [-direction[0], -direction[1], -direction[2]]
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

fn inner_product(left: &[f64], right: &[f64]) -> f64 {
    left.iter()
        .zip(right)
        .map(|(&left, &right)| left * right)
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
