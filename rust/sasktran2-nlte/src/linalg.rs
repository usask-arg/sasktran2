use crate::prelude::*;
use nalgebra::DMatrix;

/// Solves the dense system `a x = b` by LU decomposition with partial pivoting.
///
/// `a` must have shape `(n, n)` and `b` shape `(n, nrhs)`. Returns `x` with
/// the shape of `b`, or an error if the matrix is exactly singular.
pub fn solve_linear_system(a: &Array2<f64>, b: &Array2<f64>) -> Result<Array2<f64>> {
    let (n, ncols_a) = a.dim();
    let (n_b, nrhs) = b.dim();

    if n == 0 {
        return Err(anyhow!("linear solve requires n > 0"));
    }
    if ncols_a != n {
        return Err(anyhow!(
            "coefficient matrix must be square: got shape ({}, {})",
            n,
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

    let a_mat = DMatrix::from_fn(n, n, |i, j| a[[i, j]]);
    let b_mat = DMatrix::from_fn(n, nrhs, |i, j| b[[i, j]]);

    let x = a_mat
        .lu()
        .solve(&b_mat)
        .ok_or_else(|| anyhow!("coefficient matrix is singular"))?;

    Ok(Array2::from_shape_fn((n, nrhs), |(i, j)| x[(i, j)]))
}

#[cfg(test)]
mod tests {
    use super::*;
    use ndarray::array;

    #[test]
    fn solves_system_requiring_pivoting() {
        let a = array![[0.0, 2.0, 1.0], [1.0, 1.0, 0.0], [3.0, 0.0, 1.0]];
        let x_expected = array![[1.0, -2.0], [2.0, 0.5], [-1.0, 4.0]];
        let b = a.dot(&x_expected);

        let x = solve_linear_system(&a, &b).unwrap();

        for (actual, expected) in x.iter().zip(x_expected.iter()) {
            assert!((actual - expected).abs() < 1e-14);
        }
    }

    #[test]
    fn singular_matrix_is_an_error() {
        let a = array![[1.0, 2.0], [2.0, 4.0]];
        let b = array![[1.0], [2.0]];

        let err = solve_linear_system(&a, &b).unwrap_err();

        assert!(err.to_string().contains("singular"));
    }

    #[test]
    fn shape_mismatch_is_an_error() {
        let a = array![[1.0, 0.0], [0.0, 1.0]];
        let b = array![[1.0], [2.0], [3.0]];

        assert!(solve_linear_system(&a, &b).is_err());
    }
}
