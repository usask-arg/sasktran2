#[cfg(feature = "simd")]
use std::simd::f64x8;

#[cfg(feature = "simd")]
#[allow(non_camel_case_types)]
pub type f64s = f64x8;

#[cfg(not(feature = "simd"))]
pub fn axpy(a: f64, x: &[f64], y: &mut [f64]) {
    for i in 0..x.len() {
        y[i] += a * x[i];
    }
}

#[cfg(feature = "simd")]
pub fn axpy(a: f64, x: &[f64], y: &mut [f64]) {
    let lanes = f64s::LEN;

    let chunks = x.chunks_exact(lanes);
    let remainder = x.len() % lanes;
    let n = x.len();

    let a_simd = f64s::splat(a);

    for (x, y) in chunks.zip(y.chunks_exact_mut(lanes)) {
        let x_simd = f64s::from_slice(x);
        let y_simd = f64s::from_slice(y);
        let result = y_simd + (a_simd * x_simd);
        result.copy_to_slice(y);
    }

    for i in n - remainder..n {
        y[i] += a * x[i];
    }
}
