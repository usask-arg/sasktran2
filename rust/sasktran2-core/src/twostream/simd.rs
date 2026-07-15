//! Portable-SIMD preprocessing kernels.

use std::simd::{Simd, num::SimdFloat};

const LANES: usize = 4;
type F64x4 = Simd<f64, LANES>;

fn exp(value: F64x4) -> F64x4 {
    let input = value.to_array();
    let mut output = wide::f64x4::from(input).exp().to_array();
    if input
        .iter()
        .any(|value| !value.is_finite() || value.abs() >= 708.0)
    {
        for (output, input) in output.iter_mut().zip(input) {
            if !input.is_finite() || input.abs() >= 708.0 {
                *output = input.exp();
            }
        }
    }
    F64x4::from_array(output)
}

#[cfg(test)]
#[test]
fn vector_exp_tracks_scalar_libm() {
    for index in (0..=100_000).step_by(LANES) {
        let input = std::array::from_fn(|lane| -707.9 * (index + lane) as f64 / 100_000.0);
        let actual = exp(F64x4::from_array(input)).to_array();
        for (actual, input) in actual.into_iter().zip(input) {
            let expected = input.exp();
            let relative_error = ((actual - expected) / expected).abs();
            assert!(
                relative_error <= 3.0e-15,
                "exp({input}) relative error {relative_error}"
            );
        }
    }
    let edge = [f64::NEG_INFINITY, -740.0, 0.0, f64::INFINITY];
    assert_eq!(exp(F64x4::from_array(edge)).to_array(), edge.map(f64::exp));
}

pub(super) fn prepare_solar_boundary(
    factors: &[f64],
    num_wavelengths: usize,
    optical_depth: &[f64],
    irradiance: &[f64],
    slant: &mut [f64],
    transmission: &mut [f64],
) {
    let chunks = num_wavelengths / LANES;
    for (layer, &factor) in factors.iter().enumerate() {
        if factor == 0.0 {
            continue;
        }
        let optical_depth = &optical_depth[layer * num_wavelengths..(layer + 1) * num_wavelengths];
        let factor = F64x4::splat(factor);
        for chunk in 0..chunks {
            let offset = chunk * LANES;
            let value = F64x4::from_slice(&slant[offset..])
                + factor * F64x4::from_slice(&optical_depth[offset..]);
            value.copy_to_slice(&mut slant[offset..offset + LANES]);
        }
        for wave in chunks * LANES..num_wavelengths {
            slant[wave] += factor[0] * optical_depth[wave];
        }
    }

    for chunk in 0..chunks {
        let offset = chunk * LANES;
        let negative_slant = -F64x4::from_slice(&slant[offset..]);
        negative_slant.copy_to_slice(&mut slant[offset..offset + LANES]);
        let value = exp(negative_slant) * F64x4::from_slice(&irradiance[offset..]);
        value.copy_to_slice(&mut transmission[offset..offset + LANES]);
    }
    for wave in chunks * LANES..num_wavelengths {
        slant[wave] = -slant[wave];
        transmission[wave] = slant[wave].exp() * irradiance[wave];
    }
}

pub(super) fn prepare_average_secant(
    top: &[f64],
    bottom: &[f64],
    optical_depth: &[f64],
    output: &mut [f64],
) {
    let chunks = output.len() / LANES;
    for chunk in 0..chunks {
        let offset = chunk * LANES;
        let value = (F64x4::from_slice(&top[offset..]) - F64x4::from_slice(&bottom[offset..]))
            / F64x4::from_slice(&optical_depth[offset..]);
        value.copy_to_slice(&mut output[offset..offset + LANES]);
    }
    for wave in chunks * LANES..output.len() {
        output[wave] = (top[wave] - bottom[wave]) / optical_depth[wave];
    }
}

#[allow(clippy::too_many_arguments)]
pub(super) fn prepare_optical(
    parallel: bool,
    num_layers: usize,
    num_wavelengths: usize,
    thickness: &[f64],
    extinction: &[f64],
    level_ssa: &[f64],
    level_b1: &[f64],
    optical_depth: &mut [f64],
    layer_ssa: &mut [f64],
    layer_b1: &mut [f64],
) {
    if parallel {
        use rayon::prelude::*;
        optical_depth
            .par_chunks_mut(num_wavelengths)
            .zip(layer_ssa.par_chunks_mut(num_wavelengths))
            .zip(layer_b1.par_chunks_mut(num_wavelengths))
            .enumerate()
            .for_each(|(layer, ((optical_depth, layer_ssa), layer_b1))| {
                prepare_optical_layer(
                    layer,
                    num_wavelengths,
                    thickness[layer],
                    extinction,
                    level_ssa,
                    level_b1,
                    optical_depth,
                    layer_ssa,
                    layer_b1,
                );
            });
        return;
    }

    for (layer, ((optical_depth, layer_ssa), layer_b1)) in optical_depth
        .chunks_mut(num_wavelengths)
        .zip(layer_ssa.chunks_mut(num_wavelengths))
        .zip(layer_b1.chunks_mut(num_wavelengths))
        .enumerate()
        .take(num_layers)
    {
        prepare_optical_layer(
            layer,
            num_wavelengths,
            thickness[layer],
            extinction,
            level_ssa,
            level_b1,
            optical_depth,
            layer_ssa,
            layer_b1,
        );
    }
}

#[allow(clippy::too_many_arguments)]
fn prepare_optical_layer(
    layer: usize,
    num_wavelengths: usize,
    layer_thickness: f64,
    extinction: &[f64],
    level_ssa: &[f64],
    level_b1: &[f64],
    optical_depth: &mut [f64],
    layer_ssa: &mut [f64],
    layer_b1: &mut [f64],
) {
    let clamp = F64x4::splat(1.0 - 1.0e-9);
    let top = layer * num_wavelengths;
    let bottom = (layer + 1) * num_wavelengths;
    let chunks = num_wavelengths / LANES;
    for chunk in 0..chunks {
        let offset = chunk * LANES;
        let kt = F64x4::from_slice(&extinction[top + offset..]);
        let kb = F64x4::from_slice(&extinction[bottom + offset..]);
        let st = F64x4::from_slice(&level_ssa[top + offset..]);
        let sb = F64x4::from_slice(&level_ssa[bottom + offset..]);
        let bt = F64x4::from_slice(&level_b1[top + offset..]);
        let bb = F64x4::from_slice(&level_b1[bottom + offset..]);
        let avg_k = (kt + kb) * F64x4::splat(0.5);
        let kst = kt * st;
        let ksb = kb * sb;
        let avg_ks = (kst + ksb) * F64x4::splat(0.5);
        let od = avg_k * F64x4::splat(layer_thickness);
        let ssa = (avg_ks / avg_k).simd_min(clamp);
        let b1 = ((kst * bt + ksb * bb) * F64x4::splat(0.5)) / avg_ks;
        od.copy_to_slice(&mut optical_depth[offset..offset + LANES]);
        ssa.copy_to_slice(&mut layer_ssa[offset..offset + LANES]);
        b1.copy_to_slice(&mut layer_b1[offset..offset + LANES]);
    }
    for wave in chunks * LANES..num_wavelengths {
        let ti = top + wave;
        let bi = bottom + wave;
        let avg_k = 0.5 * (extinction[ti] + extinction[bi]);
        let kst = extinction[ti] * level_ssa[ti];
        let ksb = extinction[bi] * level_ssa[bi];
        let avg_ks = 0.5 * (kst + ksb);
        optical_depth[wave] = avg_k * layer_thickness;
        layer_ssa[wave] = (avg_ks / avg_k).min(1.0 - 1.0e-9);
        layer_b1[wave] = 0.5 * (kst * level_b1[ti] + ksb * level_b1[bi]) / avg_ks;
    }
}
