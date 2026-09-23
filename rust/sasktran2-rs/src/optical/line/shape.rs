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

#[derive(Clone, Copy, Default)]
struct PolePair {
    numerator: [f64; 4],
    d_numerator: [f64; 4],
    denominator: [f64; 2],
    d_denominator: [f64; 4],
}

/// Coefficients that are constant across a line's spectral points. Projecting
/// the complex amplitude onto each pole first avoids forming four separate
/// complex value/derivative sums in the hot loop, including for line mixing.
struct PreparedShape {
    center: f64,
    inverse_width: f64,
    dx_offset: f64,
    dx_slope: f64,
    scale: f64,
    d_scale: f64,
    gamma_squared: f64,
    d_gamma_squared: f64,
    d_center: f64,
    poles: [PolePair; 4],
}

const GAUSSIAN: u8 = 0;
const LORENTZIAN: u8 = 1;
const VOIGT: u8 = 2;

impl PreparedShape {
    fn new<const SHAPE: u8>(line: &AdjustedLineParameters, direction: &LineShapeDirection) -> Self {
        let mut shape = Self {
            center: line.line_center,
            inverse_width: 1.0 / line.doppler_width,
            dx_offset: -direction.line_center / line.doppler_width,
            dx_slope: -direction.doppler_width / line.doppler_width,
            scale: line.line_intensity_re,
            d_scale: direction.line_intensity_re,
            gamma_squared: 0.0,
            d_gamma_squared: 0.0,
            d_center: direction.line_center,
            poles: [PolePair::default(); 4],
        };
        if SHAPE == LORENTZIAN {
            let width = line.doppler_width;
            let gamma = line.y * width;
            shape.gamma_squared = gamma * gamma;
            shape.d_gamma_squared =
                2.0 * gamma * (direction.y * width + line.y * direction.doppler_width);
            shape.scale = line.line_intensity_re * line.y / SQRT_PI * width * width;
            shape.d_scale = ((direction.line_intensity_re * line.y
                + line.line_intensity_re * direction.y)
                * width
                * width
                + 2.0 * line.line_intensity_re * line.y * width * direction.doppler_width)
                / SQRT_PI;
        } else if SHAPE >= VOIGT {
            let width = line.doppler_width;
            let a = line.line_intensity_re * width / SQRT_PI;
            let b = line.line_intensity_im * width / SQRT_PI;
            let da = (direction.line_intensity_re * width
                + line.line_intensity_re * direction.doppler_width)
                / SQRT_PI;
            let db = (direction.line_intensity_im * width
                + line.line_intensity_im * direction.doppler_width)
                / SQRT_PI;
            for (j, pole) in shape.poles.iter_mut().enumerate() {
                // Work in wavenumber units so changing the Doppler width moves
                // a pole, rather than rescaling every spectral coordinate.
                let c = width * CJ[j].re;
                let dc = direction.doppler_width * CJ[j].re;
                let y = width * (line.y - CJ[j].im);
                let dy = direction.doppler_width * (line.y - CJ[j].im) + width * direction.y;
                let sum = c * c + y * y;
                let difference = y * y - c * c;
                let d_sum = 2.0 * (c * dc + y * dy);
                let d_difference = 2.0 * (y * dy - c * dc);
                let r = BJ[j].re;
                let m = BJ[j].im;
                let even = m * c - r * y;
                let offset = m * c + r * y;
                let odd = r * difference + 2.0 * m * y * c;
                let numerator = [
                    -2.0 * a * offset * sum,
                    2.0 * b * odd,
                    2.0 * a * even,
                    2.0 * b * r,
                ];
                let d_numerator = [
                    -2.0 * ((da * offset + a * (m * dc + r * dy)) * sum + a * offset * d_sum)
                        - direction.line_center * numerator[1],
                    2.0 * (db * odd + b * (r * d_difference + 2.0 * m * (dy * c + y * dc)))
                        - 2.0 * direction.line_center * numerator[2],
                    2.0 * (da * even + a * (m * dc - r * dy))
                        - 3.0 * direction.line_center * numerator[3],
                    2.0 * db * r,
                ];
                // Poles j and j+4 have opposite real parts and equal imaginary
                // parts. Their sum is P3(x) / (x^4 + q2*x^2 + q0), halving
                // divisions without changing the rational approximation.
                *pole = PolePair {
                    numerator,
                    d_numerator,
                    denominator: [sum * sum, 2.0 * difference],
                    d_denominator: [
                        2.0 * sum * d_sum,
                        -4.0 * difference * direction.line_center,
                        2.0 * d_difference,
                        -4.0 * direction.line_center,
                    ],
                };
            }
        }
        shape
    }
}

#[inline(always)]
fn evaluate<const SHAPE: u8, const N: usize, const ODD: usize, F>(
    wavenumber: F,
    shapes: &[PreparedShape; N],
    scalar: impl Fn(f64) -> F,
    exp: impl Fn(F) -> F,
) -> (F, [F; N])
where
    F: Copy + Add<Output = F> + Sub<Output = F> + Mul<Output = F> + Div<Output = F>,
{
    let shape = &shapes[0];
    let delta = wavenumber - scalar(shape.center);
    if SHAPE == LORENTZIAN {
        let inverse = scalar(1.0) / (delta * delta + scalar(shape.gamma_squared));
        let value = scalar(shape.scale) * inverse;
        let derivatives = shapes.each_ref().map(|shape| {
            (scalar(shape.d_scale)
                - value * (scalar(shape.d_gamma_squared) - scalar(2.0 * shape.d_center) * delta))
                * inverse
        });
        (value, derivatives)
    } else {
        let x = delta * scalar(shape.inverse_width);
        if SHAPE == GAUSSIAN {
            let value = exp(scalar(0.0) - x * x);
            (
                scalar(shape.scale) * value,
                shapes.each_ref().map(|shape| {
                    let dx = scalar(shape.dx_offset) + x * scalar(shape.dx_slope);
                    (scalar(shape.d_scale) - scalar(2.0 * shape.scale) * x * dx) * value
                }),
            )
        } else {
            let mut value = scalar(0.0);
            let mut derivatives = [scalar(0.0); N];
            let squared = delta * delta;
            let fourth = squared * squared;
            for (j, pole) in shape.poles.iter().enumerate() {
                let inverse = scalar(1.0)
                    / (fourth
                        + scalar(pole.denominator[1]) * squared
                        + scalar(pole.denominator[0]));
                let mut numerator = scalar(pole.numerator[2]) * squared + scalar(pole.numerator[0]);
                if ODD & 1 != 0 {
                    numerator = numerator
                        + delta * (scalar(pole.numerator[3]) * squared + scalar(pole.numerator[1]));
                }
                let term = numerator * inverse;
                value = value + term;
                // Differentiate the rational approximation itself, preserving
                // consistency with the forward kernel even in the line wings.
                for (k, (derivative, shape)) in derivatives.iter_mut().zip(shapes).enumerate() {
                    let pole = &shape.poles[j];
                    let mut d_numerator =
                        scalar(pole.d_numerator[2]) * squared + scalar(pole.d_numerator[0]);
                    let mut d_denominator =
                        scalar(pole.d_denominator[2]) * squared + scalar(pole.d_denominator[0]);
                    if ODD == usize::MAX || ODD & (1 << (k + 1)) != 0 {
                        d_numerator = d_numerator
                            + delta
                                * (scalar(pole.d_numerator[3]) * squared
                                    + scalar(pole.d_numerator[1]));
                        d_denominator = d_denominator
                            + delta
                                * (scalar(pole.d_denominator[3]) * squared
                                    + scalar(pole.d_denominator[1]));
                    }
                    *derivative = *derivative + (d_numerator - term * d_denominator) * inverse;
                }
            }
            (value, derivatives)
        }
    }
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
    assign_with_derivatives(
        wavenumbers,
        line,
        &[*direction],
        shape,
        values,
        [derivatives],
    );
}

/// Accumulate multiple directions together, sharing profile values and the
/// expensive exponentials/divisions. N is fixed outside the spectral loops.
pub fn assign_with_derivatives<const N: usize>(
    wavenumbers: &[f64],
    line: &AdjustedLineParameters,
    directions: &[LineShapeDirection; N],
    shape: LineShape,
    values: &mut [f64],
    derivatives: [&mut [f64]; N],
) {
    assert!(N > 0);
    assert_eq!(wavenumbers.len(), values.len());
    assert!(derivatives.iter().all(|d| d.len() == wavenumbers.len()));
    // Select the kernel once per region, outside both scalar and SIMD loops.
    match shape {
        LineShape::Gaussian => {
            assign::<GAUSSIAN, N, 0>(wavenumbers, line, directions, values, derivatives)
        }
        LineShape::Lorentzian => {
            assign::<LORENTZIAN, N, 0>(wavenumbers, line, directions, values, derivatives)
        }
        LineShape::Voigt => {
            // Bit zero selects odd value terms; subsequent bits select odd
            // terms in each direction. A pressure shift need not make the
            // temperature derivative or the unshifted profile asymmetric.
            // Specialize one/two directions; larger requests and complex
            // amplitudes use the fully general kernel.
            let odd = if N <= 2 && line.line_intensity_im == 0.0 {
                directions.iter().enumerate().fold(0, |mask, (k, d)| {
                    mask | (usize::from(d.line_intensity_im != 0.0 || d.line_center != 0.0)
                        << (k + 1))
                })
            } else {
                usize::MAX
            };
            match odd {
                0 => assign::<VOIGT, N, 0>(wavenumbers, line, directions, values, derivatives),
                2 => assign::<VOIGT, N, 2>(wavenumbers, line, directions, values, derivatives),
                4 => assign::<VOIGT, N, 4>(wavenumbers, line, directions, values, derivatives),
                6 => assign::<VOIGT, N, 6>(wavenumbers, line, directions, values, derivatives),
                _ => assign::<VOIGT, N, { usize::MAX }>(
                    wavenumbers,
                    line,
                    directions,
                    values,
                    derivatives,
                ),
            }
        }
    }
}

fn assign<const SHAPE: u8, const N: usize, const ODD: usize>(
    wavenumbers: &[f64],
    line: &AdjustedLineParameters,
    directions: &[LineShapeDirection; N],
    values: &mut [f64],
    mut derivatives: [&mut [f64]; N],
) {
    if wavenumbers.is_empty() {
        return;
    }
    let shapes = directions.map(|direction| PreparedShape::new::<SHAPE>(line, &direction));
    #[cfg(not(feature = "simd"))]
    let start = 0;
    #[cfg(feature = "simd")]
    let start = {
        use crate::math::simd::f64s;
        use std::simd::StdFloat;
        let lanes = f64s::LEN;
        let mut derivative_chunks = derivatives.each_mut().map(|d| d.chunks_exact_mut(lanes));
        for (wv, value) in wavenumbers
            .chunks_exact(lanes)
            .zip(values.chunks_exact_mut(lanes))
        {
            let (v, gradients) =
                evaluate::<SHAPE, N, ODD, _>(f64s::from_slice(wv), &shapes, f64s::splat, |x| {
                    x.exp()
                });
            (f64s::from_slice(value) + v).copy_to_slice(value);
            for (chunks, gradient) in derivative_chunks.iter_mut().zip(gradients) {
                let derivative = chunks.next().unwrap();
                (f64s::from_slice(derivative) + gradient).copy_to_slice(derivative);
            }
        }
        wavenumbers.len() / lanes * lanes
    };
    for (i, (&wv, value)) in wavenumbers[start..]
        .iter()
        .zip(&mut values[start..])
        .enumerate()
    {
        let (v, gradients) = evaluate::<SHAPE, N, ODD, _>(wv, &shapes, |x| x, f64::exp);
        *value += v;
        for (derivative, gradient) in derivatives.iter_mut().zip(gradients) {
            derivative[start + i] += gradient;
        }
    }
}
