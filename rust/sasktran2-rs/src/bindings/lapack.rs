use anyhow::{Result, anyhow};
use ndarray::{Array1, Array2};
use sasktran2_sys::ffi;

fn to_col_major_vec(a: &Array2<f64>) -> Vec<f64> {
    let (nrows, ncols) = a.dim();
    let mut out = vec![0.0; nrows * ncols];

    for j in 0..ncols {
        for i in 0..nrows {
            out[j * nrows + i] = a[[i, j]];
        }
    }

    out
}

fn from_col_major_vec(buf: &[f64], nrows: usize, ncols: usize) -> Array2<f64> {
    let mut out = Array2::<f64>::zeros((nrows, ncols));

    for j in 0..ncols {
        for i in 0..nrows {
            out[[i, j]] = buf[j * nrows + i];
        }
    }

    out
}

/// Solves a banded linear system using LAPACK dgbsv through the C API.
///
/// The input matrix `a` must have shape `(n, n)`.
///
/// The RHS matrix `b` must have shape `(n, nrhs)`.
///
/// Returns `(solution, ipiv)` on success where `solution` has the same shape
/// as the input `b` and `ipiv` has length `n`.
pub fn dgesv(a: &Array2<f64>, b: &Array2<f64>) -> Result<(Array2<f64>)> {
    let (nrows_a, ncols_a) = a.dim();
    let n = nrows_a;
    let (n_b, nrhs) = b.dim();

    if n == 0 {
        return Err(anyhow!("dgesv requires n > 0"));
    }
    if ncols_a != n {
        return Err(anyhow!(
            "coefficient matrix must be square: got shape ({}, {})",
            nrows_a,
            ncols_a
        ));
    }
    if n_b != n {
        return Err(anyhow!(
            "rhs row count mismatch: b has {} rows but matrix has n={} columns",
            n_b,
            n
        ));
    }

    let mut a_col_major = to_col_major_vec(a);
    let mut b_col_major = to_col_major_vec(b);
    let mut ipiv = vec![0_i64; n];

    let info = unsafe {
        ffi::sk_lapack_dgesv(
            n as i64,
            nrhs as i64,
            a_col_major.as_mut_ptr(),
            n as i64,
            ipiv.as_mut_ptr(),
            b_col_major.as_mut_ptr(),
            n as i64,
        )
    };

    if info != 0 {
        return Err(anyhow!("sk_lapack_dgesv failed with info={}", info));
    }

    let x = from_col_major_vec(&b_col_major, n, nrhs);
    Ok(x)
}
