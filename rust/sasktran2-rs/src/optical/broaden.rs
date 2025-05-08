#[cfg(feature = "simd")]
use crate::math::simd::f64s;
use std::f64::consts::FRAC_1_SQRT_2;

use crate::math::errorfunctions::optimized::{
    w_jpole_assign, w_jpole_real_assign, w_jpole_real_assign_uniform,
};
use crate::prelude::*;

pub const SQRT_PI: f64 = 1.772_453_850_905_516;

trait LineShapes {
    fn lorentzian(self, y: f64) -> f64;
    fn gaussian(self, y: f64) -> f64;
}

impl LineShapes for f64 {
    #[inline(always)]
    fn lorentzian(self, y: f64) -> f64 {
        let x = self;

        y / (SQRT_PI * (x * x + y * y))
    }

    #[inline(always)]
    fn gaussian(self, _y: f64) -> f64 {
        let x = self;

        (-x * x).exp()
    }
}

#[cfg(not(feature = "simd"))]
fn lorentzian_assign_uniform(x_start: f64, x_delta: f64, y: f64, scale: f64, result: &mut [f64]) {
    for (i, result) in result.iter_mut().enumerate() {
        let x = x_start + (i as f64) * x_delta;
        *result += x.lorentzian(y) * scale;
    }
}

#[cfg(feature = "simd")]
fn lorentzian_assign_uniform(x_start: f64, x_delta: f64, y: f64, scale: f64, result: &mut [f64]) {
    let lanes = f64s::LEN;
    let mut x = f64s::splat(x_start);
    let delta = f64s::from_slice(&[
        0.0,
        x_delta,
        2.0 * x_delta,
        3.0 * x_delta,
        4.0 * x_delta,
        5.0 * x_delta,
        6.0 * x_delta,
        7.0 * x_delta,
    ]);

    x += delta;

    let ys = f64s::splat(y);
    let constant = f64s::splat(scale / SQRT_PI) * ys;

    let chunks = result.chunks_exact_mut(lanes);

    for result in chunks {
        let denom = x * x + ys * ys;
        let lorentzian = constant / denom;

        let r = f64s::from_slice(result);
        let updated = r + lorentzian;
        updated.copy_to_slice(result);

        x += f64s::splat(lanes as f64 * x_delta);
    }

    let remainder_len = result.len() % lanes;
    let n = result.len();

    #[allow(clippy::needless_range_loop)]
    for i in n - remainder_len..n {
        let x = x_start + (i as f64) * x_delta;
        result[i] += x.lorentzian(y) * scale;
    }
}

fn gaussian_assign_uniform(x_start: f64, x_delta: f64, y: f64, scale: f64, result: &mut [f64]) {
    for (i, result) in result.iter_mut().enumerate() {
        let x = x_start + (i as f64) * x_delta;
        *result += x.gaussian(y) * scale;
    }
}

#[allow(clippy::too_many_arguments)]
pub fn voigt_broaden_uniform(
    line_center: ArrayView1<f64>,
    line_intensity: ArrayView1<f64>,
    lower_energy: ArrayView1<f64>,
    gamma_air: ArrayView1<f64>,
    gamma_self: ArrayView1<f64>,
    delta_air: ArrayView1<f64>,
    n_air: ArrayView1<f64>,
    iso_id: ArrayView1<i32>,
    partitions: ArrayView2<f64>,
    mol_mass: ArrayView1<f64>,
    pressure: ArrayView1<f64>,
    pself: ArrayView1<f64>,
    temperature: ArrayView1<f64>,
    first_wavenumber: f64,
    wavenumber_spacing: f64,
    mut result: ArrayViewMut2<f64>, // [geometry, wavenumber]
    line_contribution_width: f64,
    cull_factor: f64,
    _num_threads: usize,
    subtract_pedastal: bool,
) -> Result<()> {
    const EPSILON: f64 = 1e-4;
    const C2: f64 = 1.4387769;
    const SPEED_OF_LIGHT: f64 = 2.99792458e10;
    const NA: f64 = 6.02214179e23;
    const K_B: f64 = 1.38064852e-16;

    let n_wavenumber = result.len_of(Axis(1));

    let last_wavenumber: f64 = first_wavenumber + (n_wavenumber as f64 - 1.0) * wavenumber_spacing;

    result.fill(0.0);

    let start = first_wavenumber - line_contribution_width;
    // Binary search for the first element >= to_x
    let start_line_idx = match line_center
        .as_slice()
        .unwrap()
        .binary_search_by(|x| x.partial_cmp(&start).unwrap())
    {
        Ok(i) => i,  // Exact match found
        Err(i) => i, // Not found, i is the index where `to_x` would be inserted
    };

    let end = last_wavenumber + line_contribution_width;
    let end_line_idx = match line_center
        .as_slice()
        .unwrap()
        .binary_search_by(|x| x.partial_cmp(&end).unwrap())
    {
        Ok(i) => i,  // Exact match found
        Err(i) => i, // Not found, i is the index where `to_x` would be inserted
    };

    let max_p_self = pself.iter().copied().fold(0.0, f64::max);
    let max_p_self = match max_p_self {
        0.0 => 1.0,
        _ => max_p_self,
    };

    let mut zero_pivot: i64 = 0;
    Zip::indexed(line_center.slice(ndarray::s![start_line_idx..end_line_idx])).for_each(
        |i, &lc| {
            let i = i + start_line_idx;
            let line_intensity = line_intensity[i];

            if line_intensity * 101325.0 * max_p_self / (K_B * 1e-7 * 296.0) < cull_factor {
                return;
            }

            let start_wavenumber_idx: i64 = ((lc - line_contribution_width - first_wavenumber)
                / wavenumber_spacing)
                .floor() as i64;
            let end_wavenumber_idx: i64 = ((lc + line_contribution_width - first_wavenumber)
                / wavenumber_spacing)
                .ceil() as i64;

            let start_wavenumber_idx = start_wavenumber_idx.max(0) as usize;
            let end_wavenumber_idx = end_wavenumber_idx.min(n_wavenumber as i64) as usize;

            while zero_pivot < end_wavenumber_idx as i64
                && zero_pivot as f64 * wavenumber_spacing + first_wavenumber < lc
            {
                zero_pivot += 1;
            }

            if start_wavenumber_idx == end_wavenumber_idx {
                return;
            }

            let le = lower_energy[i];

            let denominator = (1.0 - (-C2 * lc / 296.0).exp()) * (-C2 * le / 296.0).exp();

            Zip::indexed(temperature)
                .and(pressure)
                .and(pself)
                .and(result.axis_iter_mut(Axis(0)))
                .par_for_each(|g, &temperature, &pressure, &pself, mut result| {
                    let numerator =
                        (-C2 * le / temperature).exp() * (1.0 - (-C2 * lc / temperature).exp());

                    let common_factor = numerator / denominator;
                    let iso_index = (iso_id[i] - 1) as usize;
                    let partition_factor = 1.0 / partitions[[g, iso_index]];
                    let adjusted_line_intensity = line_intensity * common_factor * partition_factor;

                    let mol_mass = mol_mass[iso_index];
                    let doppler_width = lc / SPEED_OF_LIGHT
                        * (NA * K_B * temperature / mol_mass).sqrt()
                        / FRAC_1_SQRT_2;
                    let g_air = gamma_air[i];
                    let g_self = gamma_self[i];
                    let da = delta_air[i];
                    let na = n_air[i];

                    let gamma_val = (296.0 / temperature).powf(na)
                        * (g_air * (pressure - pself) + g_self * pself);

                    let shifted_center = lc + da * pressure;
                    let y = gamma_val / doppler_width;

                    let norm_factor = 1.0 / (SQRT_PI * doppler_width);

                    let normalized_intensity = adjusted_line_intensity * norm_factor;

                    let x_start = (start_wavenumber_idx as f64 * wavenumber_spacing
                        + first_wavenumber
                        - shifted_center)
                        / doppler_width;
                    let x_delta = wavenumber_spacing / doppler_width;

                    if 2.84 * y * y > 1.52 / EPSILON {
                        Zip::indexed(
                            result.slice_mut(ndarray::s![start_wavenumber_idx..end_wavenumber_idx]),
                        )
                        .for_each(|w, result| {
                            *result +=
                                (x_start + w as f64 * x_delta).lorentzian(y) * normalized_intensity;
                        });
                    } else {
                        let max_x = ((end_wavenumber_idx as f64) * wavenumber_spacing
                            + first_wavenumber
                            - shifted_center)
                            / doppler_width;
                        let min_x = ((start_wavenumber_idx as f64) * wavenumber_spacing
                            + first_wavenumber
                            - shifted_center)
                            / doppler_width;

                        let max_abs_x = max_x.abs().max(min_x.abs());

                        if max_abs_x < 2.15 - 2.53 * y / EPSILON {
                            gaussian_assign_uniform(
                                x_start,
                                x_delta,
                                y,
                                normalized_intensity,
                                result
                                    .slice_mut(ndarray::s![
                                        start_wavenumber_idx..end_wavenumber_idx
                                    ])
                                    .as_slice_mut()
                                    .unwrap(),
                            );
                        } else {
                            let split_x = (1.52 / EPSILON - 2.84 * y * y).sqrt();

                            let mut lorentzian_split =
                                ((split_x * doppler_width) / wavenumber_spacing).floor() as i64;

                            if zero_pivot + lorentzian_split > end_wavenumber_idx as i64 {
                                lorentzian_split = end_wavenumber_idx as i64 - zero_pivot;
                            }

                            if zero_pivot - lorentzian_split < start_wavenumber_idx as i64 {
                                lorentzian_split = zero_pivot - start_wavenumber_idx as i64;
                            }

                            let left = (zero_pivot - lorentzian_split) as usize;
                            let right = (zero_pivot + lorentzian_split) as usize;

                            lorentzian_assign_uniform(
                                x_start,
                                x_delta,
                                y,
                                normalized_intensity,
                                result
                                    .slice_mut(ndarray::s![start_wavenumber_idx..left])
                                    .as_slice_mut()
                                    .unwrap(),
                            );

                            /*
                            Zip::indexed(result.slice_mut(ndarray::s![start_wavenumber_idx..left]))
                                .for_each(|w, result| {
                                    *result += (x_start + w as f64 * x_delta).lorentzian(y)
                                        * normalized_intensity;
                                });
                            */

                            let x_start =
                                x_start + ((left - start_wavenumber_idx) as f64 * x_delta);

                            w_jpole_real_assign_uniform(
                                x_start,
                                x_delta,
                                y,
                                normalized_intensity,
                                result
                                    .slice_mut(ndarray::s![left..right])
                                    .as_slice_mut()
                                    .unwrap(),
                            );

                            let x_start = x_start + ((right - left) as f64 * x_delta);

                            lorentzian_assign_uniform(
                                x_start,
                                x_delta,
                                y,
                                normalized_intensity,
                                result
                                    .slice_mut(ndarray::s![right..end_wavenumber_idx])
                                    .as_slice_mut()
                                    .unwrap(),
                            );

                            /*
                            Zip::indexed(result.slice_mut(ndarray::s![right..end_wavenumber_idx]))
                                .for_each(|w, result| {
                                    *result += (x_start + w as f64 * x_delta).lorentzian(y)
                                        * normalized_intensity;
                                });
                                */
                        }
                    }

                    if subtract_pedastal {
                        let x = line_contribution_width / doppler_width;
                        let pedastal = x.lorentzian(y) * normalized_intensity;
                        Zip::from(
                            result.slice_mut(ndarray::s![start_wavenumber_idx..end_wavenumber_idx]),
                        )
                        .for_each(|result| {
                            *result -= pedastal;
                        });
                    }
                });
        },
    );

    Ok(())
}

#[allow(clippy::too_many_arguments)]
pub fn voigt_broaden(
    line_center: ArrayView1<f64>,
    line_intensity: ArrayView1<f64>,
    lower_energy: ArrayView1<f64>,
    gamma_air: ArrayView1<f64>,
    gamma_self: ArrayView1<f64>,
    delta_air: ArrayView1<f64>,
    n_air: ArrayView1<f64>,
    iso_id: ArrayView1<i32>,
    partitions: ArrayView2<f64>,
    mol_mass: ArrayView1<f64>,
    pressure: ArrayView1<f64>,
    pself: ArrayView1<f64>,
    temperature: ArrayView1<f64>,
    wavenumber_grid: ArrayView1<f64>,
    mut result: ArrayViewMut2<f64>, // [geometry, wavenumber]
    line_contribution_width: f64,
    cull_factor: f64,
    _num_threads: usize,
    subtract_pedastal: bool,
) -> Result<()> {
    const C2: f64 = 1.4387769;
    const SPEED_OF_LIGHT: f64 = 2.99792458e10;
    const NA: f64 = 6.02214179e23;
    const K_B: f64 = 1.38064852e-16;

    let n_wavenumber = result.len_of(Axis(1));

    let last_wavenumber: f64 = wavenumber_grid[n_wavenumber - 1];

    result.fill(0.0);

    let start = wavenumber_grid[0] - line_contribution_width;
    // Binary search for the first element >= to_x
    let start_line_idx = match line_center
        .as_slice()
        .unwrap()
        .binary_search_by(|x| x.partial_cmp(&start).unwrap())
    {
        Ok(i) => i,  // Exact match found
        Err(i) => i, // Not found, i is the index where `to_x` would be inserted
    };

    let end = last_wavenumber + line_contribution_width;
    let end_line_idx = match line_center
        .as_slice()
        .unwrap()
        .binary_search_by(|x| x.partial_cmp(&end).unwrap())
    {
        Ok(i) => i,  // Exact match found
        Err(i) => i, // Not found, i is the index where `to_x` would be inserted
    };

    let max_p_self = pself.iter().copied().fold(0.0, f64::max);
    let max_p_self = match max_p_self {
        0.0 => 1.0,
        _ => max_p_self,
    };

    let mut start_wavenumber_idx: i64 = 0;
    let mut end_wavenumber_idx: i64 = 0;

    Zip::indexed(line_center.slice(ndarray::s![start_line_idx..end_line_idx])).for_each(
        |i, &lc| {
            let i = i + start_line_idx;
            let line_intensity = line_intensity[i];

            if line_intensity * 101325.0 * max_p_self / (K_B * 1e-7 * 296.0) < cull_factor {
                return;
            }

            while start_wavenumber_idx < n_wavenumber as i64
                && wavenumber_grid[start_wavenumber_idx as usize] < lc - line_contribution_width
            {
                start_wavenumber_idx += 1;
            }

            while end_wavenumber_idx < n_wavenumber as i64
                && wavenumber_grid[end_wavenumber_idx as usize] < lc + line_contribution_width
            {
                end_wavenumber_idx += 1;
            }

            if start_wavenumber_idx == end_wavenumber_idx {
                return;
            }

            let le = lower_energy[i];

            let denominator = (1.0 - (-C2 * lc / 296.0).exp()) * (-C2 * le / 296.0).exp();

            Zip::indexed(temperature)
                .and(pressure)
                .and(pself)
                .and(result.axis_iter_mut(Axis(0)))
                .par_for_each(|g, &temperature, &pressure, &pself, mut result| {
                    let numerator =
                        (-C2 * le / temperature).exp() * (1.0 - (-C2 * lc / temperature).exp());

                    let common_factor = numerator / denominator;

                    let iso_index = (iso_id[i] - 1) as usize;

                    let partition_factor = 1.0 / partitions[[g, iso_index]];
                    let adjusted_line_intensity = line_intensity * common_factor * partition_factor;

                    let mol_mass = mol_mass[iso_index];
                    let doppler_width = lc / SPEED_OF_LIGHT
                        * (NA * K_B * temperature / mol_mass).sqrt()
                        / FRAC_1_SQRT_2;
                    let g_air = gamma_air[i];
                    let g_self = gamma_self[i];
                    let da = delta_air[i];
                    let na = n_air[i];

                    let gamma_val = (296.0 / temperature).powf(na)
                        * (g_air * (pressure - pself) + g_self * pself);

                    let shifted_center = lc + da * pressure;
                    let y = gamma_val / doppler_width;

                    let norm_factor = 1.0 / (SQRT_PI * doppler_width);

                    let normalized_intensity = adjusted_line_intensity * norm_factor;
                    w_jpole_real_assign(
                        wavenumber_grid
                            .slice(ndarray::s![
                                start_wavenumber_idx as usize..end_wavenumber_idx as usize
                            ])
                            .as_slice()
                            .unwrap(),
                        shifted_center,
                        doppler_width,
                        y,
                        normalized_intensity,
                        result
                            .slice_mut(ndarray::s![
                                start_wavenumber_idx as usize..end_wavenumber_idx as usize
                            ])
                            .as_slice_mut()
                            .unwrap(),
                    );

                    if subtract_pedastal {
                        let x = line_contribution_width / doppler_width;
                        let pedastal = x.lorentzian(y) * normalized_intensity;
                        Zip::from(result.slice_mut(ndarray::s![
                            start_wavenumber_idx as usize..end_wavenumber_idx as usize
                        ]))
                        .for_each(|result| {
                            *result -= pedastal;
                        });
                    }
                });
        },
    );

    Ok(())
}

#[allow(clippy::too_many_arguments)]
pub fn voigt_broaden_with_line_coupling(
    line_center: ArrayView1<f64>,
    line_intensity: ArrayView1<f64>,
    lower_energy: ArrayView1<f64>,
    gamma_air: ArrayView1<f64>,
    gamma_self: ArrayView1<f64>,
    delta_air: ArrayView1<f64>,
    n_air: ArrayView1<f64>,
    iso_id: ArrayView1<i32>,
    partitions: ArrayView2<f64>,
    y_coupling: ArrayView2<f64>, // [line, geometry]
    g_coupling: ArrayView2<f64>, // [line, geometry]
    mol_mass: ArrayView1<f64>,
    pressure: ArrayView1<f64>,
    pself: ArrayView1<f64>,
    temperature: ArrayView1<f64>,
    wavenumber_grid: ArrayView1<f64>,
    mut result: ArrayViewMut2<f64>, // [geometry, wavenumber]
    line_contribution_width: f64,
    cull_factor: f64,
    _num_threads: usize,
    subtract_pedastal: bool,
) -> Result<()> {
    const C2: f64 = 1.4387769;
    const SPEED_OF_LIGHT: f64 = 2.99792458e10;
    const NA: f64 = 6.02214179e23;
    const K_B: f64 = 1.38064852e-16;

    let n_wavenumber = result.len_of(Axis(1));

    let last_wavenumber: f64 = wavenumber_grid[n_wavenumber - 1];

    result.fill(0.0);

    let start = wavenumber_grid[0] - line_contribution_width;
    // Binary search for the first element >= to_x
    let start_line_idx = match line_center
        .as_slice()
        .unwrap()
        .binary_search_by(|x| x.partial_cmp(&start).unwrap())
    {
        Ok(i) => i,  // Exact match found
        Err(i) => i, // Not found, i is the index where `to_x` would be inserted
    };

    let end = last_wavenumber + line_contribution_width;
    let end_line_idx = match line_center
        .as_slice()
        .unwrap()
        .binary_search_by(|x| x.partial_cmp(&end).unwrap())
    {
        Ok(i) => i,  // Exact match found
        Err(i) => i, // Not found, i is the index where `to_x` would be inserted
    };

    let max_p_self = pself.iter().copied().fold(0.0, f64::max);
    let max_p_self = match max_p_self {
        0.0 => 1.0,
        _ => max_p_self,
    };

    let mut start_wavenumber_idx: i64 = 0;
    let mut end_wavenumber_idx: i64 = 0;

    Zip::indexed(line_center.slice(ndarray::s![start_line_idx..end_line_idx])).for_each(
        |i, &lc| {
            let i = i + start_line_idx;
            let line_intensity = line_intensity[i];

            if line_intensity * 101325.0 * max_p_self / (K_B * 1e-7 * 296.0) < cull_factor {
                return;
            }

            while start_wavenumber_idx < n_wavenumber as i64
                && wavenumber_grid[start_wavenumber_idx as usize] < lc - line_contribution_width
            {
                start_wavenumber_idx += 1;
            }

            while end_wavenumber_idx < n_wavenumber as i64
                && wavenumber_grid[end_wavenumber_idx as usize] < lc + line_contribution_width
            {
                end_wavenumber_idx += 1;
            }

            if start_wavenumber_idx == end_wavenumber_idx {
                return;
            }

            let le = lower_energy[i];

            let denominator = (1.0 - (-C2 * lc / 296.0).exp()) * (-C2 * le / 296.0).exp();

            Zip::indexed(temperature)
                .and(pressure)
                .and(pself)
                .and(result.axis_iter_mut(Axis(0)))
                .par_for_each(|g, &temperature, &pressure, &pself, mut result| {
                    let numerator =
                        (-C2 * le / temperature).exp() * (1.0 - (-C2 * lc / temperature).exp());

                    let common_factor = numerator / denominator;
                    let iso_index = (iso_id[i] - 1) as usize;
                    let partition_factor = 1.0 / partitions[[g, iso_index]];
                    let adjusted_line_intensity = line_intensity * common_factor * partition_factor;

                    let mol_mass = mol_mass[iso_index];
                    let doppler_width = lc / SPEED_OF_LIGHT
                        * (NA * K_B * temperature / mol_mass).sqrt()
                        / FRAC_1_SQRT_2;
                    let g_air = gamma_air[i];
                    let g_self = gamma_self[i];
                    let da = delta_air[i];
                    let na = n_air[i];

                    let gamma_val = (296.0 / temperature).powf(na)
                        * (g_air * (pressure - pself) + g_self * pself);

                    let shifted_center = lc + da * pressure;
                    let y = gamma_val / doppler_width;

                    let norm_factor = 1.0 / (SQRT_PI * doppler_width);

                    let normalized_intensity = adjusted_line_intensity * norm_factor;

                    let scale_re =
                        normalized_intensity * (1.0 + pressure * pressure * g_coupling[[i, g]]);
                    let scale_im = normalized_intensity * (-pressure * y_coupling[[i, g]]);

                    w_jpole_assign(
                        wavenumber_grid
                            .slice(ndarray::s![
                                start_wavenumber_idx as usize..end_wavenumber_idx as usize
                            ])
                            .as_slice()
                            .unwrap(),
                        shifted_center,
                        doppler_width,
                        y,
                        scale_re,
                        scale_im,
                        result
                            .slice_mut(ndarray::s![
                                start_wavenumber_idx as usize..end_wavenumber_idx as usize
                            ])
                            .as_slice_mut()
                            .unwrap(),
                    );

                    if subtract_pedastal {
                        let x = line_contribution_width / doppler_width;
                        let pedastal = x.lorentzian(y) * normalized_intensity;
                        Zip::from(result.slice_mut(ndarray::s![
                            start_wavenumber_idx as usize..end_wavenumber_idx as usize
                        ]))
                        .for_each(|result| {
                            *result -= pedastal;
                        });
                    }
                });
        },
    );

    Ok(())
}
