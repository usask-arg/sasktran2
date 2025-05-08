#[cfg(feature = "simd")]
use super::super::simd::f64s;
use num::complex::Complex64;

pub const SQRT_PI: f64 = 1.772_453_850_905_516;

/// w(z) via an 8‐pole Padé sum (J‐pole approximation).
#[inline(always)]
pub fn w_jpole(z: Complex64) -> Complex64 {
    const BJ: [Complex64; 8] = [
        Complex64::new(0.00383968430671409, -0.0119854387180615),
        Complex64::new(-0.321597857664957, -0.218883985607935),
        Complex64::new(2.55515264319988, 0.613958600684469),
        Complex64::new(-2.73739446984183, 5.69007914897806),
        Complex64::new(0.00383968430671409, 0.0119854387180615),
        Complex64::new(-0.321597857664957, 0.218883985607935),
        Complex64::new(2.55515264319988, -0.613958600684469),
        Complex64::new(-2.73739446984183, -5.69007914897806),
    ];

    const CJ: [Complex64; 8] = [
        Complex64::new(2.51506776338386, -1.60713668042405),
        Complex64::new(-1.68985621846204, -1.66471695485661),
        Complex64::new(0.981465428659098, -1.70017951305004),
        Complex64::new(-0.322078795578047, -1.71891780447016),
        Complex64::new(-2.51506776338386, -1.60713668042405),
        Complex64::new(1.68985621846204, -1.66471695485661),
        Complex64::new(-0.981465428659098, -1.70017951305004),
        Complex64::new(0.322078795578047, -1.71891780447016),
    ];

    // manually fully unrolled
    let sum: Complex64 = BJ[0] / (z - CJ[0])
        + BJ[1] / (z - CJ[1])
        + BJ[2] / (z - CJ[2])
        + BJ[3] / (z - CJ[3])
        + BJ[4] / (z - CJ[4])
        + BJ[5] / (z - CJ[5])
        + BJ[6] / (z - CJ[6])
        + BJ[7] / (z - CJ[7]);

    sum / SQRT_PI * Complex64::new(0.0, -1.0)
}

// Assigns the real part of result of the w(z) function to the provided slice, assuming uniform x spacing
#[cfg(not(feature = "simd"))]
pub fn w_jpole_real_assign_uniform(
    x_start: f64,
    x_delta: f64,
    y: f64,
    scale: f64,
    result: &mut [f64],
) {
    // real parts of the BJ poles
    const BJ_RE: [f64; 8] = [
        0.00383968430671409,
        -0.321597857664957,
        2.55515264319988,
        -2.73739446984183,
        0.00383968430671409,
        -0.321597857664957,
        2.55515264319988,
        -2.73739446984183,
    ];

    // imaginary parts of the BJ poles
    const BJ_IM: [f64; 8] = [
        -0.0119854387180615,
        -0.218883985607935,
        0.613958600684469,
        5.69007914897806,
        0.0119854387180615,
        0.218883985607935,
        -0.613958600684469,
        -5.69007914897806,
    ];

    // real parts of the CJ poles
    const CJ_RE: [f64; 8] = [
        2.51506776338386,
        -1.68985621846204,
        0.981465428659098,
        -0.322078795578047,
        -2.51506776338386,
        1.68985621846204,
        -0.981465428659098,
        0.322078795578047,
    ];

    // imaginary parts of the CJ poles
    const CJ_IM: [f64; 8] = [
        -1.60713668042405,
        -1.66471695485661,
        -1.70017951305004,
        -1.71891780447016,
        -1.60713668042405,
        -1.66471695485661,
        -1.70017951305004,
        -1.71891780447016,
    ];

    // do this *once* for each distinct y, not per‐x
    let mut dy = [0.0; 8];
    let mut dy2 = [0.0; 8];
    for j in 0..8 {
        dy[j] = y - CJ_IM[j];
        dy2[j] = dy[j] * dy[j];
    }

    const INV_SQRT_PI: f64 = 1.0 / SQRT_PI;

    for (i, result) in result.iter_mut().enumerate() {
        let x = x_start + i as f64 * x_delta;
        let mut sum_im = 0.0;
        for j in 0..8 {
            let dx = x - CJ_RE[j];
            let den = dx * dx + dy2[j];
            // numerator of complex division
            let num_im = BJ_IM[j] * dx - BJ_RE[j] * dy[j];
            sum_im += num_im / den;
        }
        // multiply by –i/√π:  (sum_re + i sum_im) * (0 – i)/√π
        *result += sum_im * INV_SQRT_PI * scale;
    }
}

#[cfg(feature = "simd")]
pub fn w_jpole_real_assign_uniform(
    x_start: f64,
    x_delta: f64,
    y: f64,
    scale: f64,
    result: &mut [f64],
) {
    // real parts of the BJ poles
    const BJ_RE: [f64; 8] = [
        0.00383968430671409,
        -0.321597857664957,
        2.55515264319988,
        -2.73739446984183,
        0.00383968430671409,
        -0.321597857664957,
        2.55515264319988,
        -2.73739446984183,
    ];

    // imaginary parts of the BJ poles
    const BJ_IM: [f64; 8] = [
        -0.0119854387180615,
        -0.218883985607935,
        0.613958600684469,
        5.69007914897806,
        0.0119854387180615,
        0.218883985607935,
        -0.613958600684469,
        -5.69007914897806,
    ];

    // real parts of the CJ poles
    const CJ_RE: [f64; 8] = [
        2.51506776338386,
        -1.68985621846204,
        0.981465428659098,
        -0.322078795578047,
        -2.51506776338386,
        1.68985621846204,
        -0.981465428659098,
        0.322078795578047,
    ];

    // imaginary parts of the CJ poles
    const CJ_IM: [f64; 8] = [
        -1.60713668042405,
        -1.66471695485661,
        -1.70017951305004,
        -1.71891780447016,
        -1.60713668042405,
        -1.66471695485661,
        -1.70017951305004,
        -1.71891780447016,
    ];

    // do this *once* for each distinct y, not per‐x
    let mut dy = [0.0; 8];
    let mut dy2 = [0.0; 8];
    for j in 0..8 {
        dy[j] = y - CJ_IM[j];
        dy2[j] = dy[j] * dy[j];
    }

    const INV_SQRT_PI: f64 = 1.0 / SQRT_PI;

    let lanes = f64s::LEN;
    let remainder = result.len() % lanes;
    let chunks = result.chunks_exact_mut(lanes);

    for (i, result) in chunks.enumerate() {
        let i = i * lanes;
        let x_start = x_start + i as f64 * x_delta;
        let x = f64s::splat(x_start)
            + f64s::splat(x_delta) * f64s::from_slice(&[0.0, 1.0, 2.0, 3.0, 4.0, 5.0, 6.0, 7.0]);

        let mut sum_im = f64s::splat(0.0);
        for j in 0..8 {
            let dx = x - f64s::splat(CJ_RE[j]);
            let den = dx * dx + f64s::splat(dy2[j]);
            // numerator of complex division
            let num_im = f64s::splat(BJ_IM[j]) * dx - f64s::splat(BJ_RE[j]) * f64s::splat(dy[j]);
            sum_im += num_im / den;
        }
        // multiply by –i/√π:  (sum_re + i sum_im) * (0 – i)/√π
        sum_im *= f64s::splat(INV_SQRT_PI * scale);
        let arr = sum_im.to_array();
        for j in 0..lanes {
            result[j] += arr[j];
        }
    }

    // Handle the remainder
    let n = result.len();
    #[allow(clippy::needless_range_loop)]
    for i in n - remainder..n {
        let x = x_start + i as f64 * x_delta;
        let mut sum_im = 0.0;
        for j in 0..8 {
            let dx = x - CJ_RE[j];
            let den = dx * dx + dy2[j];
            // numerator of complex division
            let num_im = BJ_IM[j] * dx - BJ_RE[j] * dy[j];
            sum_im += num_im / den;
        }
        // multiply by –i/√π:  (sum_re + i sum_im) * (0 – i)/√π
        result[i] += sum_im * INV_SQRT_PI * scale;
    }
}

// Assigns the real part of result of the w(z) function to the provided slice.
#[cfg(not(feature = "simd"))]
pub fn w_jpole_real_assign(
    wvnum: &[f64],
    c: f64,
    width: f64,
    y: f64,
    scale: f64,
    result: &mut [f64],
) {
    // real parts of the BJ poles
    const BJ_RE: [f64; 8] = [
        0.00383968430671409,
        -0.321597857664957,
        2.55515264319988,
        -2.73739446984183,
        0.00383968430671409,
        -0.321597857664957,
        2.55515264319988,
        -2.73739446984183,
    ];

    // imaginary parts of the BJ poles
    const BJ_IM: [f64; 8] = [
        -0.0119854387180615,
        -0.218883985607935,
        0.613958600684469,
        5.69007914897806,
        0.0119854387180615,
        0.218883985607935,
        -0.613958600684469,
        -5.69007914897806,
    ];

    // real parts of the CJ poles
    const CJ_RE: [f64; 8] = [
        2.51506776338386,
        -1.68985621846204,
        0.981465428659098,
        -0.322078795578047,
        -2.51506776338386,
        1.68985621846204,
        -0.981465428659098,
        0.322078795578047,
    ];

    // imaginary parts of the CJ poles
    const CJ_IM: [f64; 8] = [
        -1.60713668042405,
        -1.66471695485661,
        -1.70017951305004,
        -1.71891780447016,
        -1.60713668042405,
        -1.66471695485661,
        -1.70017951305004,
        -1.71891780447016,
    ];

    // do this *once* for each distinct y, not per‐x
    let mut dy = [0.0; 8];
    let mut dy2 = [0.0; 8];
    for j in 0..8 {
        dy[j] = y - CJ_IM[j];
        dy2[j] = dy[j] * dy[j];
    }

    const INV_SQRT_PI: f64 = 1.0 / SQRT_PI;

    for (wvnum, result) in wvnum.iter().zip(result.iter_mut()) {
        let x = (wvnum - c) / width;
        let mut sum_im = 0.0;
        for j in 0..8 {
            let dx = x - CJ_RE[j];
            let den = dx * dx + dy2[j];
            // numerator of complex division
            let num_im = BJ_IM[j] * dx - BJ_RE[j] * dy[j];
            sum_im += num_im / den;
        }
        // multiply by –i/√π:  (sum_re + i sum_im) * (0 – i)/√π
        *result += sum_im * INV_SQRT_PI * scale;
    }
}

// Assigns the real part of result of the w(z) function to the provided slice.
#[cfg(feature = "simd")]
pub fn w_jpole_real_assign(
    wvnum: &[f64],
    c: f64,
    width: f64,
    y: f64,
    scale: f64,
    result: &mut [f64],
) {
    // real parts of the BJ poles
    const BJ_RE: [f64; 8] = [
        0.00383968430671409,
        -0.321597857664957,
        2.55515264319988,
        -2.73739446984183,
        0.00383968430671409,
        -0.321597857664957,
        2.55515264319988,
        -2.73739446984183,
    ];

    // imaginary parts of the BJ poles
    const BJ_IM: [f64; 8] = [
        -0.0119854387180615,
        -0.218883985607935,
        0.613958600684469,
        5.69007914897806,
        0.0119854387180615,
        0.218883985607935,
        -0.613958600684469,
        -5.69007914897806,
    ];

    // real parts of the CJ poles
    const CJ_RE: [f64; 8] = [
        2.51506776338386,
        -1.68985621846204,
        0.981465428659098,
        -0.322078795578047,
        -2.51506776338386,
        1.68985621846204,
        -0.981465428659098,
        0.322078795578047,
    ];

    // imaginary parts of the CJ poles
    const CJ_IM: [f64; 8] = [
        -1.60713668042405,
        -1.66471695485661,
        -1.70017951305004,
        -1.71891780447016,
        -1.60713668042405,
        -1.66471695485661,
        -1.70017951305004,
        -1.71891780447016,
    ];

    // do this *once* for each distinct y, not per‐x
    let mut dy = [0.0; 8];
    let mut dy2 = [0.0; 8];
    for j in 0..8 {
        dy[j] = y - CJ_IM[j];
        dy2[j] = dy[j] * dy[j];
    }

    const INV_SQRT_PI: f64 = 1.0 / SQRT_PI;

    let lanes = f64s::LEN;
    let chunks = wvnum.chunks_exact(lanes);
    let remainder = chunks.remainder();

    // Process the chunks of 4 elements at a time
    for (wvnum, result) in chunks.zip(result.chunks_exact_mut(lanes)) {
        let wv = f64s::from_slice(wvnum);
        let x = (wv - f64s::splat(c)) / f64s::splat(width);

        let mut sum_im = f64s::splat(0.0);
        for j in 0..8 {
            let dx = x - f64s::splat(CJ_RE[j]);
            let den = dx * dx + f64s::splat(dy2[j]);
            // numerator of complex division
            let num_im = f64s::splat(BJ_IM[j]) * dx - f64s::splat(BJ_RE[j]) * f64s::splat(dy[j]);
            sum_im += num_im / den;
        }
        sum_im *= f64s::splat(INV_SQRT_PI * scale);
        sum_im += f64s::from_slice(result);

        sum_im.copy_to_slice(result);
    }
    let n = wvnum.len();
    for i in n - remainder.len()..n {
        let x = (wvnum[i] - c) / width;
        let mut sum_im = 0.0;
        for j in 0..8 {
            let dx = x - CJ_RE[j];
            let den = dx * dx + dy2[j];
            // numerator of complex division
            let num_im = BJ_IM[j] * dx - BJ_RE[j] * dy[j];
            sum_im += num_im / den;
        }
        // multiply by –i/√π:  (sum_re + i sum_im) * (0 – i)/√π
        result[i] += sum_im * INV_SQRT_PI * scale;
    }
}

// Assigns the the real component of (scale + i scale_im) * w(z) to the provided slice. Used for line coupling
// calculations
#[inline(always)]
#[cfg(not(feature = "simd"))]
pub fn w_jpole_assign(
    wvnum: &[f64],
    c: f64,
    width: f64,
    y: f64,
    scale_re: f64,
    scale_im: f64,
    result: &mut [f64],
) {
    // real parts of the BJ poles
    const BJ_RE: [f64; 8] = [
        0.00383968430671409,
        -0.321597857664957,
        2.55515264319988,
        -2.73739446984183,
        0.00383968430671409,
        -0.321597857664957,
        2.55515264319988,
        -2.73739446984183,
    ];

    // imaginary parts of the BJ poles
    const BJ_IM: [f64; 8] = [
        -0.0119854387180615,
        -0.218883985607935,
        0.613958600684469,
        5.69007914897806,
        0.0119854387180615,
        0.218883985607935,
        -0.613958600684469,
        -5.69007914897806,
    ];

    // real parts of the CJ poles
    const CJ_RE: [f64; 8] = [
        2.51506776338386,
        -1.68985621846204,
        0.981465428659098,
        -0.322078795578047,
        -2.51506776338386,
        1.68985621846204,
        -0.981465428659098,
        0.322078795578047,
    ];

    // imaginary parts of the CJ poles
    const CJ_IM: [f64; 8] = [
        -1.60713668042405,
        -1.66471695485661,
        -1.70017951305004,
        -1.71891780447016,
        -1.60713668042405,
        -1.66471695485661,
        -1.70017951305004,
        -1.71891780447016,
    ];

    // do this *once* for each distinct y, not per‐x
    let mut dy = [0.0; 8];
    let mut dy2 = [0.0; 8];
    for j in 0..8 {
        dy[j] = y - CJ_IM[j];
        dy2[j] = dy[j] * dy[j];
    }

    const INV_SQRT_PI: f64 = 1.0 / SQRT_PI;

    for (wvnum, result) in wvnum.iter().zip(result.iter_mut()) {
        let x = (wvnum - c) / width;
        let mut sum_im = 0.0;
        let mut sum_re = 0.0;
        for j in 0..8 {
            let dx = x - CJ_RE[j];
            let den = dx * dx + dy2[j];
            // numerator of complex division
            let num_im = BJ_IM[j] * dx - BJ_RE[j] * dy[j];
            let num_re = BJ_RE[j] * dx + BJ_IM[j] * dy[j];
            sum_im += num_im / den;
            sum_re += num_re / den;
        }
        // multiply by –i/√π:  (sum_re + i sum_im) * -i / sqrt(pi) *  (scale_re + i scale_im)
        // = 1 / sqrt(pi) * (sum_re + i sum_im) * (scale_im - i scale_re)
        *result += sum_im * INV_SQRT_PI * scale_re;
        *result += sum_re * INV_SQRT_PI * scale_im;
    }
}

// Assigns the the real component of (scale + i scale_im) * w(z) to the provided slice. Used for line coupling
// calculations
#[inline(always)]
#[cfg(feature = "simd")]
pub fn w_jpole_assign(
    wvnum: &[f64],
    c: f64,
    width: f64,
    y: f64,
    scale_re: f64,
    scale_im: f64,
    result: &mut [f64],
) {
    // real parts of the BJ poles
    const BJ_RE: [f64; 8] = [
        0.00383968430671409,
        -0.321597857664957,
        2.55515264319988,
        -2.73739446984183,
        0.00383968430671409,
        -0.321597857664957,
        2.55515264319988,
        -2.73739446984183,
    ];

    // imaginary parts of the BJ poles
    const BJ_IM: [f64; 8] = [
        -0.0119854387180615,
        -0.218883985607935,
        0.613958600684469,
        5.69007914897806,
        0.0119854387180615,
        0.218883985607935,
        -0.613958600684469,
        -5.69007914897806,
    ];

    // real parts of the CJ poles
    const CJ_RE: [f64; 8] = [
        2.51506776338386,
        -1.68985621846204,
        0.981465428659098,
        -0.322078795578047,
        -2.51506776338386,
        1.68985621846204,
        -0.981465428659098,
        0.322078795578047,
    ];

    // imaginary parts of the CJ poles
    const CJ_IM: [f64; 8] = [
        -1.60713668042405,
        -1.66471695485661,
        -1.70017951305004,
        -1.71891780447016,
        -1.60713668042405,
        -1.66471695485661,
        -1.70017951305004,
        -1.71891780447016,
    ];

    // do this *once* for each distinct y, not per‐x
    let mut dy = [0.0; 8];
    let mut dy2 = [0.0; 8];
    for j in 0..8 {
        dy[j] = y - CJ_IM[j];
        dy2[j] = dy[j] * dy[j];
    }

    const INV_SQRT_PI: f64 = 1.0 / SQRT_PI;

    let lanes = f64s::LEN;
    let chunks = wvnum.chunks_exact(lanes);
    let remainder = chunks.remainder();

    // Process the chunks of 4 elements at a time
    for (wvnum, result) in chunks.zip(result.chunks_exact_mut(lanes)) {
        let wv = f64s::from_slice(wvnum);
        let x = (wv - f64s::splat(c)) / f64s::splat(width);

        let mut sum_im = f64s::splat(0.0);
        let mut sum_re = f64s::splat(0.0);
        for j in 0..8 {
            let dx = x - f64s::splat(CJ_RE[j]);
            let den = dx * dx + f64s::splat(dy2[j]);
            // numerator of complex division
            let num_im = f64s::splat(BJ_IM[j]) * dx - f64s::splat(BJ_RE[j]) * f64s::splat(dy[j]);
            let num_re = f64s::splat(BJ_RE[j]) * dx + f64s::splat(BJ_IM[j]) * f64s::splat(dy[j]);
            sum_im += num_im / den;
            sum_re += num_re / den;
        }
        sum_im *= f64s::splat(INV_SQRT_PI * scale_re);
        sum_re *= f64s::splat(INV_SQRT_PI * scale_im);

        sum_im += f64s::from_slice(result) + sum_re;

        sum_im.copy_to_slice(result);
    }
    let n = wvnum.len();
    for i in n - remainder.len()..n {
        let x = (wvnum[i] - c) / width;
        let mut sum_im = 0.0;
        let mut sum_re = 0.0;
        for j in 0..8 {
            let dx = x - CJ_RE[j];
            let den = dx * dx + dy2[j];
            // numerator of complex division
            let num_im = BJ_IM[j] * dx - BJ_RE[j] * dy[j];
            let num_re = BJ_RE[j] * dx + BJ_IM[j] * dy[j];
            sum_im += num_im / den;
            sum_re += num_re / den;
        }
        // multiply by –i/√π:  (sum_re + i sum_im) * (0 – i)/√π
        result[i] += sum_im * INV_SQRT_PI * scale_re + sum_re * INV_SQRT_PI * scale_im;
    }
}
