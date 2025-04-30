use num::complex::Complex64;

pub const SQRT_PI: f64 = 1.7724538509055160272981674833411;

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

pub fn w_jpole_assign(x_start: f64, x_delta: f64, y: f64, scale: f64, result: &mut [f64]) {
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
