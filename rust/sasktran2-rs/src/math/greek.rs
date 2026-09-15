//! Quadrature transform shared by tabulated phase matrices and the Mie integrator.

use crate::prelude::*;

use super::wigner::WignerDCalculator;

pub struct GreekTransform {
    basis: [Array2<f64>; 4],
}

impl GreekTransform {
    pub fn new(cos_angles: ArrayView1<f64>, num_coefficients: usize) -> Self {
        let basis = [(0, 0), (2, 2), (2, -2), (0, 2)].map(|(m, n)| {
            let mut table = WignerDCalculator::new(m, n).matrix_d(cos_angles, num_coefficients);
            for (l, mut row) in table.outer_iter_mut().enumerate() {
                row *= (2 * l + 1) as f64 / 2.0;
            }
            table
        });
        Self { basis }
    }

    /// Returns [a1, a2, a3, a4, b1, b2] for [p11, p12, p22, p33, p34, p44].
    /// The angle axis is contiguous in the basis so every dot product reads it
    /// sequentially. Weights and phase combinations are computed once per row.
    pub fn project(
        &self,
        phase: [ArrayView1<f64>; 6],
        weights: ArrayView1<f64>,
    ) -> [Array1<f64>; 6] {
        let [p11, p12, p22, p33, p34, p44] = phase;
        let [d00, d22, d2m2, d02] = &self.basis;
        let a1 = d00.dot(&(&p11 * &weights));
        let a4 = d00.dot(&(&p44 * &weights));
        let b1 = d02.dot(&(&p12 * &weights));
        let b2 = -d02.dot(&(&p34 * &weights));
        let plus = d22.dot(&((&p22 + &p33) * weights));
        let minus = d2m2.dot(&((&p22 - &p33) * weights));
        let a2 = (&plus + &minus) * 0.5;
        let a3 = (&plus - &minus) * 0.5;
        [a1, a2, a3, a4, b1, b2]
    }
}
