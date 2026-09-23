//! Directional derivatives of the line shapes used by the line absorber.
//!
//! Values and derivatives share intermediates. These kernels are only called
//! when derivatives are requested; the value-only kernels remain independent.
use super::AdjustedLineParameters;
use crate::math::errorfunctions::optimized::{BJ, CJ, SQRT_PI};
use std::ops::{Add, Div, Mul, Sub};

/// A direction in the adjusted line parameters (e.g. their temperature derivatives).
/// `line_intensity_re` includes the inverse Doppler width normalization, just as
/// in `AdjustedLineParameters`. `y` is the Lorentz/Doppler width ratio.
#[derive(Clone, Copy, Debug, Default)]
pub struct LineShapeDirection {
    pub line_center: f64,
    pub doppler_width: f64,
    pub y: f64,
    pub line_intensity_re: f64,
    pub line_intensity_im: f64,
}

#[derive(Clone, Copy)]
pub enum LineShape {
    Gaussian,
    Lorentzian,
    Voigt,
}

#[inline(always)]
fn evaluate<F>(
    wavenumber: F,
    line: &AdjustedLineParameters,
    direction: &LineShapeDirection,
    shape: LineShape,
    scalar: impl Fn(f64) -> F,
    exp: impl Fn(F) -> F,
) -> (F, F)
where
    F: Copy + Add<Output = F> + Sub<Output = F> + Mul<Output = F> + Div<Output = F>,
{
    let x = (wavenumber - scalar(line.line_center)) * scalar(1.0 / line.doppler_width);
    let dx = scalar(-direction.line_center / line.doppler_width)
        - x * scalar(direction.doppler_width / line.doppler_width);
    let (real, imag, d_real, d_imag) = match shape {
        LineShape::Gaussian => {
            let value = exp(scalar(0.0) - x * x);
            (
                value,
                scalar(0.0),
                scalar(-2.0) * x * value * dx,
                scalar(0.0),
            )
        }
        LineShape::Lorentzian => {
            let denominator = x * x + scalar(line.y * line.y);
            let inverse = scalar(1.0) / denominator;
            let value = scalar(line.y / SQRT_PI) * inverse;
            let derivative = scalar(direction.y / SQRT_PI) * inverse
                - value * inverse * (scalar(2.0) * x * dx + scalar(2.0 * line.y * direction.y));
            (value, scalar(0.0), derivative, scalar(0.0))
        }
        LineShape::Voigt => {
            let mut real = scalar(0.0);
            let mut imag = scalar(0.0);
            let mut derivative_real = scalar(0.0);
            let mut derivative_imag = scalar(0.0);
            for j in 0..8 {
                let delta_x = x - scalar(CJ[j].re);
                let delta_y = line.y - CJ[j].im;
                let inverse = scalar(1.0) / (delta_x * delta_x + scalar(delta_y * delta_y));
                // Contribution to -i b_j / (sqrt(pi) (z - c_j)).
                let re = (scalar(BJ[j].im) * delta_x - scalar(BJ[j].re * delta_y)) * inverse;
                let im = (scalar(-BJ[j].re) * delta_x - scalar(BJ[j].im * delta_y)) * inverse;
                real = real + re;
                imag = imag + im;
                // Differentiate the rational approximation itself, rather than
                // applying the exact Faddeeva identity to an approximate value.
                derivative_real = derivative_real - (re * delta_x + im * scalar(delta_y)) * inverse;
                derivative_imag = derivative_imag - (im * delta_x - re * scalar(delta_y)) * inverse;
            }
            let norm = scalar(1.0 / SQRT_PI);
            (
                real * norm,
                imag * norm,
                (derivative_real * dx - derivative_imag * scalar(direction.y)) * norm,
                (derivative_imag * dx + derivative_real * scalar(direction.y)) * norm,
            )
        }
    };
    (
        scalar(line.line_intensity_re) * real - scalar(line.line_intensity_im) * imag,
        scalar(direction.line_intensity_re) * real - scalar(direction.line_intensity_im) * imag
            + scalar(line.line_intensity_re) * d_real
            - scalar(line.line_intensity_im) * d_imag,
    )
}

/// Accumulate a line and a directional derivative, including complex amplitudes
/// for line mixing. No per-line spectral derivative arrays are allocated.
pub fn assign_with_derivative(
    wavenumbers: &[f64],
    line: &AdjustedLineParameters,
    direction: &LineShapeDirection,
    shape: LineShape,
    values: &mut [f64],
    derivatives: &mut [f64],
) {
    assert_eq!(wavenumbers.len(), values.len());
    assert_eq!(wavenumbers.len(), derivatives.len());
    #[cfg(not(feature = "simd"))]
    let start = 0;
    #[cfg(feature = "simd")]
    let start = {
        use crate::math::simd::f64s;
        use std::simd::StdFloat;
        let lanes = f64s::LEN;
        for ((wv, value), derivative) in wavenumbers
            .chunks_exact(lanes)
            .zip(values.chunks_exact_mut(lanes))
            .zip(derivatives.chunks_exact_mut(lanes))
        {
            let (v, d) = evaluate(
                f64s::from_slice(wv),
                line,
                direction,
                shape,
                f64s::splat,
                |x| x.exp(),
            );
            (f64s::from_slice(value) + v).copy_to_slice(value);
            (f64s::from_slice(derivative) + d).copy_to_slice(derivative);
        }
        wavenumbers.len() / lanes * lanes
    };
    for ((&wv, value), derivative) in wavenumbers[start..]
        .iter()
        .zip(&mut values[start..])
        .zip(&mut derivatives[start..])
    {
        let (v, d) = evaluate(wv, line, direction, shape, |x| x, f64::exp);
        *value += v;
        *derivative += d;
    }
}
