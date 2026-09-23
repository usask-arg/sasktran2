use super::*;
use crate::interpolation::linear::Interp1Weights;
use crate::optical::line::shape::{LineShape, LineShapeDirection, assign_with_derivatives};
use rayon::prelude::*;

#[derive(Clone, Copy)]
pub(super) enum LineDerivative<'a> {
    Temperature,
    /// d(self pressure)/d(total pressure), e.g. VMR for fixed composition.
    Pressure(ArrayView1<'a, f64>),
}

/// Follow exactly the same approximation regions as `split_and_assign`.
fn split_with_derivatives<const N: usize>(
    line: &AdjustedLineParameters,
    directions: &[LineShapeDirection; N],
    grid: &Grid1DView,
    values: &mut [f64],
    mut derivatives: [&mut [f64]; N],
) {
    const EPSILON: f64 = 1.0e-4;
    let n = grid.x.len();
    if n == 0 {
        return;
    }
    if 2.84 * line.y * line.y > 1.52 / EPSILON {
        assign_with_derivatives(
            grid.x,
            line,
            directions,
            LineShape::Lorentzian,
            values,
            derivatives,
        );
        return;
    }
    if grid.x[0].abs().max(grid.x[n - 1].abs()) < 2.15 - 2.53 * line.y / EPSILON {
        assign_with_derivatives(
            grid.x,
            line,
            directions,
            LineShape::Gaussian,
            values,
            derivatives,
        );
        return;
    }
    let split = (1.52 / EPSILON - 2.84 * line.y * line.y).sqrt() * line.doppler_width;
    let left = grid.lower_bound(line.line_center - split).min(n);
    let right = grid.lower_bound(line.line_center + split).min(n);
    for (range, shape) in [
        (0..left, LineShape::Lorentzian),
        (left..right, LineShape::Voigt),
        (right..n, LineShape::Lorentzian),
    ] {
        assign_with_derivatives(
            &grid.x[range.clone()],
            line,
            directions,
            shape,
            &mut values[range.clone()],
            derivatives.each_mut().map(|d| &mut d[range.clone()]),
        );
    }
}

impl LineAbsorber {
    /// Cross-section derivative at fixed total and self pressure.
    pub fn cross_section_temperature_derivative(
        &self,
        wavenumbers: ArrayView1<'_, f64>,
        temperature: ArrayView1<'_, f64>,
        pressure: ArrayView1<'_, f64>,
        pself: ArrayView1<'_, f64>,
    ) -> Result<Array2<f64>> {
        Ok(self
            .cross_section_with_temperature_derivative(wavenumbers, temperature, pressure, pself)?
            .1)
    }

    /// Cross sections and their temperature derivatives in one line-list pass.
    /// The value-only batch and SIMD kernels remain independent.
    pub fn cross_section_with_temperature_derivative(
        &self,
        wavenumbers: ArrayView1<'_, f64>,
        temperature: ArrayView1<'_, f64>,
        pressure: ArrayView1<'_, f64>,
        pself: ArrayView1<'_, f64>,
    ) -> Result<(Array2<f64>, Array2<f64>)> {
        let (values, [derivative]) = self.cross_section_with_derivatives(
            wavenumbers,
            temperature,
            pressure,
            pself,
            [LineDerivative::Temperature],
        )?;
        Ok((values, derivative))
    }

    /// Values and selected state derivatives, sharing all line-profile evaluations.
    pub(super) fn cross_section_with_derivatives<const N: usize>(
        &self,
        wavenumbers: ArrayView1<'_, f64>,
        temperature: ArrayView1<'_, f64>,
        pressure: ArrayView1<'_, f64>,
        pself: ArrayView1<'_, f64>,
        requests: [LineDerivative<'_>; N],
    ) -> Result<(Array2<f64>, [Array2<f64>; N])> {
        anyhow::ensure!(
            pressure.len() == temperature.len() && pself.len() == temperature.len(),
            "Temperature, pressure, and self pressure must have the same length"
        );
        anyhow::ensure!(
            temperature.iter().all(|t| t.is_finite() && *t > 0.0),
            "Temperature must be positive and finite"
        );
        for request in &requests {
            if let LineDerivative::Pressure(d_pself_dp) = request {
                anyhow::ensure!(
                    d_pself_dp.len() == temperature.len(),
                    "Self-pressure derivative and temperature must have the same length"
                );
            }
        }
        let mut values = Array2::zeros((temperature.len(), wavenumbers.len()));
        let mut derivatives = std::array::from_fn(|_| Array2::zeros(values.raw_dim()));
        if wavenumbers.is_empty() || temperature.is_empty() {
            return Ok((values, derivatives));
        }

        let wavenumbers = wavenumbers.as_slice().unwrap();
        let permutation = (!wavenumbers.is_sorted()).then(|| argsort_f64(wavenumbers));
        let sorted_storage;
        let sorted = if let Some(permutation) = &permutation {
            sorted_storage = permutation
                .iter()
                .map(|&i| wavenumbers[i])
                .collect::<Vec<_>>();
            sorted_storage.as_slice()
        } else {
            wavenumbers
        };
        let grid = Grid1DView::new(sorted);
        let lines = self.db.between_slice(
            sorted[0] - self.line_contribution_width,
            sorted[sorted.len() - 1] + self.line_contribution_width,
        );
        let params = self.gen_mol_param(lines, temperature.as_slice().unwrap())?;
        let temperature_derivative = requests
            .iter()
            .any(|r| matches!(r, LineDerivative::Temperature));
        let d_log_partition: HashMap<_, Vec<_>> = if temperature_derivative {
            let partition = self
                .partition_generator
                .as_ref()
                .ok_or_else(|| anyhow!("Partition generator not set"))?;
            params
                .keys()
                .map(|&(mol, iso)| {
                    (
                        (mol, iso),
                        temperature
                            .iter()
                            .map(|&t| partition.log_temperature_derivative(mol, iso, t))
                            .collect(),
                    )
                })
                .collect()
        } else {
            HashMap::new()
        };
        let max_pself = pself.iter().copied().fold(0.0, f64::max);
        let max_pself = if max_pself == 0.0 {
            101325.0
        } else {
            max_pself
        };
        let prepared: Vec<_> = lines
            .iter()
            .filter_map(|line| {
                if line.line_intensity * max_pself / (1.38064852e-16 * 1e-7 * 296.0)
                    < self.cull_factor
                {
                    return None;
                }
                let start = sorted
                    .partition_point(|&v| v < line.line_center - self.line_contribution_width);
                let end = sorted
                    .partition_point(|&v| v < line.line_center + self.line_contribution_width);
                if start >= end {
                    return None;
                }
                let coupling = if self.enable_line_coupling && !line.y_coupling.is_empty() {
                    Some((
                        Array1::from_vec(line.coupling_temperature.clone()),
                        Array1::from_vec(line.y_coupling.clone()),
                        Array1::from_vec(line.g_coupling.clone()),
                    ))
                } else {
                    None
                };
                Some((line, start, end, coupling))
            })
            .collect();
        crate::threading::thread_pool()?.install(|| {
            // Group disjoint row views for Rayon without copying spectra or
            // allocating derivative buffers inside the line loop.
            let mut rows = derivatives.each_mut().map(|d| d.axis_iter_mut(Axis(0)));
            let derivative_rows: Vec<_> = (0..temperature.len())
                .map(|_| rows.each_mut().map(|rows| rows.next().unwrap()))
                .collect();
            values
                .axis_iter_mut(Axis(0))
                .into_par_iter()
                .zip(derivative_rows.into_par_iter())
                .enumerate()
                .for_each(|(i, (mut value, mut derivative_rows))| {
                    for (line, start, end, coupling) in &prepared {
                        let key = (line.mol_id, line.iso_id);
                        let mol = &params[&key];
                        let mut adjusted = line
                            .adjusted_parameters(
                                temperature[i],
                                pressure[i],
                                pself[i],
                                mol.partition_factor[i],
                                mol.mol_mass,
                            )
                            .unwrap();
                        let mut directions = requests.map(|request| match request {
                            LineDerivative::Temperature => line.adjusted_temperature_derivative(
                                &adjusted,
                                temperature[i],
                                d_log_partition[&key][i],
                            ),
                            LineDerivative::Pressure(d_pself_dp) => line
                                .adjusted_pressure_derivative(
                                    &adjusted,
                                    temperature[i],
                                    d_pself_dp[i],
                                ),
                        });
                        let range = *start..*end;
                        let out = &mut value.as_slice_mut().unwrap()[range.clone()];
                        let deriv = derivative_rows
                            .each_mut()
                            .map(|d| &mut d.as_slice_mut().unwrap()[range.clone()]);
                        if let Some((temps, ys, gs)) = coupling {
                            let y = ys.interp1(temps, temperature[i], OutOfBoundsMode::Extend);
                            let g = gs.interp1(temps, temperature[i], OutOfBoundsMode::Extend);
                            let (dy, dg): (f64, f64) = if temperature_derivative {
                                let weights =
                                    temps.interp1_weights(temperature[i], OutOfBoundsMode::Extend);
                                (
                                    weights.iter().map(|&(j, _, dw)| ys[j] * dw).sum(),
                                    weights.iter().map(|&(j, _, dw)| gs[j] * dw).sum(),
                                )
                            } else {
                                (0.0, 0.0)
                            };
                            let p = pressure[i] / 101325.0;
                            let amplitude = adjusted.line_intensity_re;
                            adjusted.line_intensity_re = amplitude * (1.0 + p * p * g);
                            adjusted.line_intensity_im = -amplitude * p * y;
                            for (direction, request) in directions.iter_mut().zip(requests) {
                                let d_amplitude = direction.line_intensity_re;
                                match request {
                                    LineDerivative::Temperature => {
                                        direction.line_intensity_re = d_amplitude
                                            * (1.0 + p * p * g)
                                            + amplitude * p * p * dg;
                                        direction.line_intensity_im =
                                            -p * (d_amplitude * y + amplitude * dy);
                                    }
                                    LineDerivative::Pressure(_) => {
                                        direction.line_intensity_re =
                                            amplitude * 2.0 * p * g / 101325.0;
                                        direction.line_intensity_im = -amplitude * y / 101325.0;
                                    }
                                }
                            }
                            assign_with_derivatives(
                                &sorted[*start..*end],
                                &adjusted,
                                &directions,
                                LineShape::Voigt,
                                out,
                                deriv,
                            );
                        } else {
                            split_with_derivatives(
                                &adjusted,
                                &directions,
                                &grid.slice(*start, *end),
                                out,
                                deriv,
                            );
                        }
                    }
                });
        });
        // The forward calculation clips negative cross sections after summing lines.
        for derivative in &mut derivatives {
            Zip::from(derivative)
                .and(&values)
                .for_each(|derivative, &value| {
                    if value <= 0.0 {
                        *derivative = 0.0;
                    }
                });
        }
        values.mapv_inplace(|value| value.max(0.0));
        if let Some(permutation) = permutation {
            let mut output = Array2::zeros(values.raw_dim());
            let mut d_output = std::array::from_fn(|_| Array2::zeros(values.raw_dim()));
            for (sorted_idx, &original_idx) in permutation.iter().enumerate() {
                output
                    .column_mut(original_idx)
                    .assign(&values.column(sorted_idx));
                for (out, derivative) in d_output.iter_mut().zip(&derivatives) {
                    out.column_mut(original_idx)
                        .assign(&derivative.column(sorted_idx));
                }
            }
            Ok((output, d_output))
        } else {
            Ok((values, derivatives))
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::optical::line::shape::assign_with_derivative;
    use std::sync::{
        Arc,
        atomic::{AtomicUsize, Ordering},
    };

    fn forward_shape(grid: &[f64], line: &AdjustedLineParameters, shape: LineShape) -> Vec<f64> {
        let mut values = vec![0.0; grid.len()];
        match shape {
            LineShape::Gaussian => gaussian_assign(
                grid,
                line.line_center,
                line.doppler_width,
                line.y,
                line.line_intensity_re,
                &mut values,
            ),
            LineShape::Lorentzian => lorentzian_assign(
                grid,
                line.line_center,
                line.doppler_width,
                line.y,
                line.line_intensity_re,
                &mut values,
            ),
            LineShape::Voigt => w_jpole_assign(
                grid,
                line.line_center,
                line.doppler_width,
                line.y,
                line.line_intensity_re,
                line.line_intensity_im,
                &mut values,
            ),
        }
        values
    }

    #[test]
    fn line_shape_parameter_derivatives_match_forward_kernels() {
        // 19 samples exercise SIMD blocks and their scalar tails, including wings.
        let grid: Vec<_> = [
            -200.0, -130.0, -80.0, -20.0, -8.0, -4.0, -2.0, -1.0, -0.3, 0.0, 0.2, 0.7, 1.0, 2.0,
            4.0, 8.0, 30.0, 100.0, 200.0,
        ]
        .iter()
        .map(|x| 7.0 + 0.2 * x)
        .collect();
        for (shape, imaginary_amplitude) in [
            (LineShape::Gaussian, 0.0),
            (LineShape::Lorentzian, 0.0),
            (LineShape::Voigt, 0.0),
            (LineShape::Voigt, 0.03),
        ] {
            for y in [0.0001, 0.3, 100.0] {
                let base = AdjustedLineParameters {
                    line_center: 7.0,
                    doppler_width: 0.2,
                    y,
                    line_intensity_re: 2.0,
                    line_intensity_im: imaginary_amplitude,
                };
                for parameter in 0..5 {
                    let mut direction = LineShapeDirection::default();
                    let step = if parameter == 2 && matches!(shape, LineShape::Lorentzian) {
                        y * 1.0e-3
                    } else {
                        1.0e-4
                    };
                    *match parameter {
                        0 => &mut direction.line_center,
                        1 => &mut direction.doppler_width,
                        2 => &mut direction.y,
                        3 => &mut direction.line_intensity_re,
                        _ => &mut direction.line_intensity_im,
                    } = 1.0;
                    let shifted = |h: f64| AdjustedLineParameters {
                        line_center: base.line_center + h * direction.line_center,
                        doppler_width: base.doppler_width + h * direction.doppler_width,
                        y: base.y + h * direction.y,
                        line_intensity_re: base.line_intensity_re + h * direction.line_intensity_re,
                        line_intensity_im: base.line_intensity_im + h * direction.line_intensity_im,
                    };
                    let mut values = vec![0.0; grid.len()];
                    let mut derivatives = vec![0.0; grid.len()];
                    assign_with_derivative(
                        &grid,
                        &base,
                        &direction,
                        shape,
                        &mut values,
                        &mut derivatives,
                    );
                    let expected = forward_shape(&grid, &base, shape);
                    let above = forward_shape(&grid, &shifted(step), shape);
                    let below = forward_shape(&grid, &shifted(-step), shape);
                    let above_half = forward_shape(&grid, &shifted(step / 2.0), shape);
                    let below_half = forward_shape(&grid, &shifted(-step / 2.0), shape);
                    for i in 0..grid.len() {
                        assert!(
                            (values[i] - expected[i]).abs() < 1.0e-12 * expected[i].abs() + 1.0e-14
                        );
                        // Richardson extrapolation removes the central difference's
                        // leading truncation error, including at stationary points.
                        let numeric = (4.0 * (above_half[i] - below_half[i]) / step
                            - (above[i] - below[i]) / (2.0 * step))
                            / 3.0;
                        assert!(
                            (derivatives[i] - numeric).abs()
                                < 3.0e-6 * numeric.abs()
                                    + 2.0e-9
                                    + 16.0 * f64::EPSILON * expected[i].abs() / step,
                            "parameter {parameter}, y {y}, sample {i}: {} vs {numeric}",
                            derivatives[i]
                        );
                    }
                }
            }
        }
    }

    #[test]
    fn split_derivative_matches_forward_approximation_regions() {
        for (center, y) in [(0.2, 0.0), (13100.0, 0.1), (13100.0, 100.0)] {
            let grid = Array1::linspace(center - 0.18, center + 0.18, 127);
            let grid = Grid1DView::new(grid.as_slice().unwrap());
            let line = AdjustedLineParameters {
                line_center: center,
                doppler_width: 0.001,
                y,
                line_intensity_re: 1.0,
                line_intensity_im: 0.0,
            };
            let direction = LineShapeDirection {
                line_intensity_re: 0.7,
                ..Default::default()
            };
            let mut reference = vec![0.0; 127];
            split_and_assign(&line, &grid, &mut reference);
            let mut values = vec![0.0; 127];
            let mut derivatives = vec![0.0; 127];
            split_with_derivatives(&line, &[direction], &grid, &mut values, [&mut derivatives]);
            for i in 0..127 {
                assert!((values[i] - reference[i]).abs() < 1e-12);
                assert!((derivatives[i] - 0.7 * reference[i]).abs() < 1e-12);
            }
        }
    }

    struct Partition {
        derivative_calls: Arc<AtomicUsize>,
    }
    impl PartitionFactor for Partition {
        fn partition_factor(&self, _: i32, _: i32, t: f64) -> f64 {
            t.powf(1.5)
        }
        fn log_temperature_derivative(&self, _: i32, _: i32, t: f64) -> f64 {
            self.derivative_calls.fetch_add(1, Ordering::Relaxed);
            1.5 / t
        }
    }
    struct Mass;
    impl MolecularMass for Mass {
        fn molecular_mass(&self, _: i32, _: i32) -> f64 {
            31.9988
        }
    }
    fn synthetic_absorber(coupling: bool, calls: Arc<AtomicUsize>) -> LineAbsorber {
        let lines = [13100.0, 13100.8]
            .into_iter()
            .map(|center| OpticalLine {
                line_center: center,
                line_intensity: 1.0e-24,
                einstein_a: None,
                lower_energy: 450.0,
                gamma_air: 0.06,
                gamma_self: 0.1,
                delta_air: 0.003,
                n_air: 0.7,
                mol_id: 7,
                iso_id: 1,
                upper_quanta: String::new(),
                lower_quanta: String::new(),
                upper_local_quanta: String::new(),
                lower_local_quanta: String::new(),
                upper_statistical_weight: None,
                lower_statistical_weight: None,
                y_coupling: vec![0.001, 0.002, 0.003],
                g_coupling: vec![0.001, 0.004, 0.002],
                coupling_temperature: vec![200.0, 250.0, 300.0],
            })
            .collect();
        LineAbsorber::new(OpticalLineDB { lines })
            .with_partition_generator(Box::new(Partition {
                derivative_calls: calls,
            }))
            .with_molecular_mass_generator(Box::new(Mass))
            .with_line_coupling(coupling)
            .with_line_contribution_width(2.0)
    }

    #[test]
    fn temperature_derivative_includes_strength_widths_and_line_mixing() {
        for mixing in [false, true] {
            let calls = Arc::new(AtomicUsize::new(0));
            let absorber = synthetic_absorber(mixing, calls.clone());
            let wv = Array1::linspace(13097.0, 13104.0, 701);
            let temp = array![190.0, 230.0, 280.0, 330.0];
            let pressure = array![1.0, 1000.0, 101325.0, 1.0e8];
            let pself = &pressure * 0.21;
            let step = 0.001;
            let above = absorber
                .cross_section(
                    wv.view(),
                    (&temp + step).view(),
                    pressure.view(),
                    pself.view(),
                )
                .unwrap();
            let below = absorber
                .cross_section(
                    wv.view(),
                    (&temp - step).view(),
                    pressure.view(),
                    pself.view(),
                )
                .unwrap();
            assert_eq!(
                calls.load(Ordering::Relaxed),
                0,
                "forward-only spectra must not evaluate partition derivatives"
            );
            let (combined_xs, derivative) = absorber
                .cross_section_with_temperature_derivative(
                    wv.view(),
                    temp.view(),
                    pressure.view(),
                    pself.view(),
                )
                .unwrap();
            assert_eq!(calls.load(Ordering::Relaxed), temp.len());
            let numeric = (above - below) / (2.0 * step);
            assert!(numeric.iter().any(|v| v.abs() > 1e-32));
            for (actual, expected) in derivative.iter().zip(numeric.iter()) {
                assert!(
                    (actual - expected).abs() < expected.abs() * 2e-5 + 1e-35,
                    "mixing {mixing}: {actual:e} vs {expected:e}"
                );
            }
            // A permutation that is not its own inverse catches incorrect restoration.
            let indices = [503, 100, 230, 340, 510];
            let unsorted = Array1::from_iter(indices.map(|i| wv[i]));
            let xs = absorber
                .cross_section(wv.view(), temp.view(), pressure.view(), pself.view())
                .unwrap();
            for (actual, expected) in combined_xs.iter().zip(xs.iter()) {
                assert!((actual - expected).abs() < expected.abs() * 1e-10 + 1e-35);
            }
            let unsorted_xs = absorber
                .cross_section(unsorted.view(), temp.view(), pressure.view(), pself.view())
                .unwrap();
            let (unsorted_combined_xs, unsorted_d) = absorber
                .cross_section_with_temperature_derivative(
                    unsorted.view(),
                    temp.view(),
                    pressure.view(),
                    pself.view(),
                )
                .unwrap();
            for (j, idx) in indices.into_iter().enumerate() {
                for i in 0..temp.len() {
                    assert!(
                        (unsorted_xs[[i, j]] - xs[[i, idx]]).abs()
                            < xs[[i, idx]].abs() * 1e-10 + 1e-35
                    );
                    assert!(
                        (unsorted_combined_xs[[i, j]] - xs[[i, idx]]).abs()
                            < xs[[i, idx]].abs() * 1e-10 + 1e-35
                    );
                    assert!(
                        (unsorted_d[[i, j]] - derivative[[i, idx]]).abs()
                            < derivative[[i, idx]].abs() * 1e-10 + 1e-35
                    );
                }
            }
        }
    }

    #[test]
    fn pressure_derivative_includes_self_broadening_shifts_and_mixing() {
        for mixing in [false, true] {
            for vmr in [0.0, 0.21, 1.0] {
                let calls = Arc::new(AtomicUsize::new(0));
                let absorber = synthetic_absorber(mixing, calls.clone());
                let wv = Array1::linspace(13097.0, 13104.0, 701);
                let temp = array![190.0, 230.0, 280.0, 330.0];
                let pressure = array![10.0, 1000.0, 101325.0, 1.0e8];
                let pself = &pressure * vmr;
                let d_pself_dp = Array1::from_elem(temp.len(), vmr);
                let (xs, [dp]) = absorber
                    .cross_section_with_derivatives(
                        wv.view(),
                        temp.view(),
                        pressure.view(),
                        pself.view(),
                        [LineDerivative::Pressure(d_pself_dp.view())],
                    )
                    .unwrap();
                assert_eq!(
                    calls.load(Ordering::Relaxed),
                    0,
                    "pressure-only evaluation must not sample partition temperature slopes"
                );
                let (joint_xs, [joint_dt, joint_dp]) = absorber
                    .cross_section_with_derivatives(
                        wv.view(),
                        temp.view(),
                        pressure.view(),
                        pself.view(),
                        [
                            LineDerivative::Temperature,
                            LineDerivative::Pressure(d_pself_dp.view()),
                        ],
                    )
                    .unwrap();
                assert_eq!(calls.load(Ordering::Relaxed), temp.len());
                let dt = absorber
                    .cross_section_temperature_derivative(
                        wv.view(),
                        temp.view(),
                        pressure.view(),
                        pself.view(),
                    )
                    .unwrap();
                for (left, right) in [(&xs, &joint_xs), (&dp, &joint_dp), (&dt, &joint_dt)] {
                    for (&a, &b) in left.iter().zip(right) {
                        assert!((a - b).abs() < b.abs() * 1e-11 + 1e-35);
                    }
                }
                let step = pressure.mapv(|p| (p * 1e-3_f64).max(1.0));
                let shifted = |factor: f64| {
                    let shifted_pressure = &pressure + &step * factor;
                    absorber
                        .cross_section(
                            wv.view(),
                            temp.view(),
                            shifted_pressure.view(),
                            (&shifted_pressure * vmr).view(),
                        )
                        .unwrap()
                };
                let numeric = ((shifted(0.5) - shifted(-0.5)) * 4.0
                    - (shifted(1.0) - shifted(-1.0)) / 2.0)
                    / 3.0
                    / step.view().insert_axis(Axis(1));
                for (actual, expected) in dp.rows().into_iter().zip(numeric.rows()) {
                    let peak = expected.iter().copied().map(f64::abs).fold(0.0, f64::max);
                    assert!(peak > 0.0);
                    for (&a, &b) in actual.iter().zip(expected) {
                        assert!(
                            (a - b).abs() < b.abs() * 2e-4 + peak * 2e-5 + 1e-35,
                            "mixing {mixing}, VMR {vmr}: {a:e} vs {b:e}"
                        );
                    }
                }
                let indices = [503, 100, 230, 340, 510];
                let unsorted = Array1::from_iter(indices.map(|i| wv[i]));
                let (_, [unsorted_dp]) = absorber
                    .cross_section_with_derivatives(
                        unsorted.view(),
                        temp.view(),
                        pressure.view(),
                        pself.view(),
                        [LineDerivative::Pressure(d_pself_dp.view())],
                    )
                    .unwrap();
                for (j, idx) in indices.into_iter().enumerate() {
                    for i in 0..temp.len() {
                        assert!(
                            (unsorted_dp[[i, j]] - dp[[i, idx]]).abs()
                                < dp[[i, idx]].abs() * 1e-10 + 1e-35
                        );
                    }
                }
            }
        }
    }

    #[test]
    fn pressure_derivative_is_finite_at_zero_pressure() {
        let absorber = synthetic_absorber(true, Arc::new(AtomicUsize::new(0)));
        let wv = Array1::linspace(13099.98, 13100.02, 35);
        let temp = array![230.0];
        let pressure = array![0.0];
        let vmr = array![0.21];
        let (_, [dp]) = absorber
            .cross_section_with_derivatives(
                wv.view(),
                temp.view(),
                pressure.view(),
                pressure.view(),
                [LineDerivative::Pressure(vmr.view())],
            )
            .unwrap();
        let forward = |p: f64| {
            absorber
                .cross_section(
                    wv.view(),
                    temp.view(),
                    array![p].view(),
                    array![0.21 * p].view(),
                )
                .unwrap()
        };
        let h = 10.0;
        let numeric = (-25.0 * forward(0.0) + 48.0 * forward(h) - 36.0 * forward(2.0 * h)
            + 16.0 * forward(3.0 * h)
            - 3.0 * forward(4.0 * h))
            / (12.0 * h);
        let peak = numeric.iter().copied().map(f64::abs).fold(0.0, f64::max);
        assert!(peak > 0.0);
        for (&a, &b) in dp.iter().zip(&numeric) {
            assert!(a.is_finite() && (a - b).abs() < 2e-5 * peak);
        }
    }

    #[test]
    fn default_partition_derivative_matches_power_law() {
        struct PowerLaw;
        impl PartitionFactor for PowerLaw {
            fn partition_factor(&self, _: i32, _: i32, t: f64) -> f64 {
                t.powf(1.5)
            }
        }
        for t in [100.0, 200.0, 296.0, 500.0] {
            assert!((PowerLaw.log_temperature_derivative(7, 1, t) - 1.5 / t).abs() < 1e-11);
        }
    }

    #[test]
    fn adjusted_temperature_derivative_includes_stimulated_emission() {
        let absorber = synthetic_absorber(false, Arc::new(AtomicUsize::new(0)));
        let mut line = absorber.db.lines.into_iter().next().unwrap();
        // Far-infrared lines exercise the stimulated-emission term, which is
        // negligible for the visible O2 lines used by the radiance tests.
        for center in [5.0, 500.0, 13100.0] {
            line.line_center = center;
            for t in [190.0, 296.0, 500.0] {
                let parameters = |temperature: f64| {
                    line.adjusted_parameters(
                        temperature,
                        80000.0,
                        16000.0,
                        (temperature / 296.0).powf(1.5),
                        31.9988,
                    )
                    .unwrap()
                };
                let base = parameters(t);
                let direction = line.adjusted_temperature_derivative(&base, t, 1.5 / t);
                let above = parameters(t + 0.001);
                let below = parameters(t - 0.001);
                for (actual, a, b) in [
                    (
                        direction.line_intensity_re,
                        above.line_intensity_re,
                        below.line_intensity_re,
                    ),
                    (
                        direction.doppler_width,
                        above.doppler_width,
                        below.doppler_width,
                    ),
                    (direction.y, above.y, below.y),
                ] {
                    let expected = (a - b) / 0.002;
                    assert!((actual - expected).abs() < 1e-6 * expected.abs());
                }
            }
        }
    }

    #[test]
    fn temperature_derivative_respects_cross_section_clipping() {
        let mut absorber = synthetic_absorber(true, Arc::new(AtomicUsize::new(0)));
        for line in &mut absorber.db.lines {
            line.y_coupling.fill(0.0);
            line.g_coupling.fill(-2.0);
        }
        let wv = Array1::linspace(13099.9, 13100.9, 31);
        let t = array![280.0];
        let p = array![101325.0];
        let pself = &p * 0.21;
        let xs = absorber
            .cross_section(wv.view(), t.view(), p.view(), pself.view())
            .unwrap();
        let (combined_xs, derivatives) = absorber
            .cross_section_with_derivatives(
                wv.view(),
                t.view(),
                p.view(),
                pself.view(),
                [
                    LineDerivative::Temperature,
                    LineDerivative::Pressure(array![0.21].view()),
                ],
            )
            .unwrap();
        assert!(xs.iter().all(|&v| v == 0.0));
        assert!(combined_xs.iter().all(|&v| v == 0.0));
        assert!(derivatives.iter().all(|d| d.iter().all(|&v| v == 0.0)));
    }
}
