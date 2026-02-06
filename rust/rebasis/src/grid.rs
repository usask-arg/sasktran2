use std::iter::Sum;
use std::ops::{Add, Mul};

use ndarray::{Array1, Array2, ArrayView2, s};

use crate::basis::BasisType;

/// Wrapper struct so that maybe eventually we change this to be
/// vectors of individual basis types
#[derive(Debug, Clone)]
struct BasisSet {
    basis: Vec<BasisType>,
}

/// A Grid is a set of values defined with an associated Basis
#[derive(Debug, Clone)]
pub struct Grid {
    basis_set: BasisSet,
}

#[derive(Debug, Clone)]
pub struct MappingMatrix {
    pub matrix: Array2<f64>,
    pub start_index: Array1<usize>,
    pub end_index: Array1<usize>,
}

impl MappingMatrix {
    pub fn dot(&self, output: ArrayView2<f64>) -> Array2<f64> {
        let mut result = Array2::<f64>::zeros((output.dim().0, self.start_index.len()));

        for i in 0..self.start_index.len() {
            let start = self.start_index[i];
            let end = self.end_index[i];

            for j in 0..output.dim().0 {
                let sum: f64 = self
                    .matrix
                    .slice(s![i, start..end])
                    .iter()
                    .zip(output.slice(s![j, start..end]).iter())
                    .map(|(m, v)| m * v)
                    .sum();
                result[[j, i]] = sum;
            }
        }

        result
    }
}

impl Grid {
    pub fn new(basis: Vec<BasisType>) -> Self {
        Grid {
            basis_set: BasisSet { basis },
        }
    }

    pub fn vec_basis(&self) -> &Vec<BasisType> {
        &self.basis_set.basis
    }
}

pub fn mapping_matrix(from_grid: &Grid, to_grid: &Grid) -> MappingMatrix {
    let n_from = from_grid.basis_set.basis.len();
    let n_to = to_grid.basis_set.basis.len();

    let mut matrix = Array2::<f64>::zeros((n_to, n_from));
    let mut start_index = Array1::<usize>::zeros(n_to);
    let mut end_index = Array1::<usize>::zeros(n_to);

    let mut from_index_low: usize = 0;
    let mut from_index_high: usize = 0;

    // We basically loop through the output values, and keep track of two indices that determine which
    // elements of the from_grid contribute

    for (i, b) in to_grid.basis_set.basis.iter().enumerate() {
        // this will increment until we find one that contributes
        while from_index_low + 1 < n_from
            && from_grid.basis_set.basis[from_index_low].upper_limit() < b.lower_limit()
        {
            from_index_low += 1;
        }

        // And increment until we find the last one that contributes
        while from_index_high < n_from
            && from_grid.basis_set.basis[from_index_high].lower_limit() < b.upper_limit()
        {
            from_index_high += 1;
        }

        for j in from_index_low..from_index_high {
            matrix[[i, j]] = b.overlap_integral(&from_grid.basis_set.basis[j]);
        }

        start_index[i] = from_index_low;
        end_index[i] = from_index_high;
    }

    MappingMatrix {
        matrix,
        start_index,
        end_index,
    }
}

// todo!, probably delete the clone here and just add in place into the result
pub fn transform_sorted<T>(from_grid: &Grid, to_grid: &mut Grid, from_vals: &[T]) -> Vec<T>
where
    T: Add + Mul<f64, Output = T> + Sum + Clone,
{
    let mut from_index_low: usize = 0;
    let mut from_index_high: usize = 0;

    let max_n = from_vals.len();

    // We basically loop through the output values, and keep track of two indices that determine which
    // elements of the from_grid contribute

    to_grid
        .basis_set
        .basis
        .iter()
        .map(|b| {
            // this will increment until we find one that contributes
            while from_grid.basis_set.basis[from_index_low].upper_limit() < b.lower_limit()
                && from_index_low < max_n
            {
                from_index_low += 1;
            }

            // And increment until we find the last one that contributes
            while from_grid.basis_set.basis[from_index_high].lower_limit() < b.upper_limit()
                && from_index_high < max_n
            {
                from_index_high += 1;
            }

            from_grid.basis_set.basis[from_index_low..from_index_high]
                .iter()
                .zip(from_vals[from_index_low..from_index_high].iter())
                .map(|(b2, v)| v.clone() * b.overlap_integral(b2))
                .sum()
        })
        .collect()
}
