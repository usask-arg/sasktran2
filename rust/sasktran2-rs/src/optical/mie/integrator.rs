use num::complex::Complex64;
use rayon::prelude::*;

use crate::math::greek::GreekTransform;
use crate::prelude::*;

use super::mie_f;

pub struct MieIntegrator {
    cos_angles: Array1<f64>,
    num_legendre: usize,
    pool: Option<rayon::ThreadPool>,
    transform: GreekTransform,
}

struct DistributionRow<'a> {
    xs_total: &'a mut f64,
    xs_scattering: &'a mut f64,
    pdf: ArrayView1<'a, f64>,
    phase: [ArrayViewMut1<'a, f64>; 4],
    coefficients: [ArrayViewMut1<'a, f64>; 6],
}

impl MieIntegrator {
    pub fn new(
        cos_angles: ArrayView1<f64>,
        num_legendre: usize,
        num_threads: usize,
    ) -> Result<Self> {
        anyhow::ensure!(
            cos_angles.iter().all(|x| x.is_finite() && x.abs() <= 1.0),
            "cos_angles must be finite and within [-1, 1]"
        );
        // A private pool honors this integrator's thread count without changing
        // the engine/global Rayon settings. One thread avoids dispatch overhead;
        // zero uses Rayon's automatic thread count.
        let pool = if num_threads == 1 {
            None
        } else {
            Some(crate::util::create_pool(num_threads)?)
        };
        Ok(Self {
            cos_angles: cos_angles.into_owned(),
            num_legendre,
            pool,
            transform: GreekTransform::new(cos_angles, num_legendre),
        })
    }

    #[allow(clippy::too_many_arguments)]
    pub fn integrate(
        &self,
        wavelength: f64,
        refractive_index: Complex64,
        size_param: ArrayView1<f64>,           // [size_param]
        pdf: ArrayView2<f64>,                  // [distribution, size_param]
        size_weights: ArrayView1<f64>,         // [size_param]
        angle_weights: ArrayView1<f64>,        // [angle]
        mut xs_total: ArrayViewMut1<f64>,      // [distribution]
        mut xs_scattering: ArrayViewMut1<f64>, // [distribution]
        mut p11: ArrayViewMut2<f64>,           // [distribution, angle]
        mut p12: ArrayViewMut2<f64>,           // [distribution, angle]
        mut p33: ArrayViewMut2<f64>,           // [distribution, angle]
        mut p34: ArrayViewMut2<f64>,           // [distribution, angle]
        mut lm_a1: ArrayViewMut2<f64>,         // [distribution, legendre]
        mut lm_a2: ArrayViewMut2<f64>,         // [distribution, legendre]
        mut lm_a3: ArrayViewMut2<f64>,         // [distribution, legendre]
        mut lm_a4: ArrayViewMut2<f64>,         // [distribution, legendre]
        mut lm_b1: ArrayViewMut2<f64>,         // [distribution, legendre]
        mut lm_b2: ArrayViewMut2<f64>,         // [distribution, legendre]
    ) -> Result<()> {
        let num_distributions = pdf.nrows();
        let num_angles = self.cos_angles.len();
        anyhow::ensure!(
            pdf.ncols() == size_param.len() && size_weights.len() == size_param.len(),
            "pdf must have shape (distribution, size_param), matching size_weights"
        );
        anyhow::ensure!(
            angle_weights.len() == num_angles
                && xs_total.len() == num_distributions
                && xs_scattering.len() == num_distributions,
            "angle weights or cross section output dimensions do not match"
        );
        for shape in [p11.dim(), p12.dim(), p33.dim(), p34.dim()] {
            anyhow::ensure!(
                shape == (num_distributions, num_angles),
                "phase outputs must have shape (distribution, angle)"
            );
        }
        for shape in [
            lm_a1.dim(),
            lm_a2.dim(),
            lm_a3.dim(),
            lm_a4.dim(),
            lm_b1.dim(),
            lm_b2.dim(),
        ] {
            anyhow::ensure!(
                shape == (num_distributions, self.num_legendre),
                "coefficient outputs must have shape (distribution, legendre)"
            );
        }
        let k = 2.0 * std::f64::consts::PI / wavelength;
        let c = 4.0 * std::f64::consts::PI / (2.0 * k * k);
        let calculate_mie = || {
            mie_f::mie_with_threads(
                size_param,
                refractive_index,
                self.cos_angles.view(),
                self.pool.is_some(),
            )
        };
        let output = if let Some(pool) = &self.pool {
            pool.install(calculate_mie)
        } else {
            calculate_mie()
        };

        // Split every output into disjoint distribution rows before dispatch.
        // Reductions within a row keep their serial order, independent of the
        // number of workers, and accept both C/F layouts and strided views.
        let phase_rows = p11
            .outer_iter_mut()
            .zip(p12.outer_iter_mut())
            .zip(p33.outer_iter_mut())
            .zip(p34.outer_iter_mut())
            .map(|(((a, b), c), d)| [a, b, c, d]);
        let coefficient_rows = lm_a1
            .outer_iter_mut()
            .zip(lm_a2.outer_iter_mut())
            .zip(lm_a3.outer_iter_mut())
            .zip(lm_a4.outer_iter_mut())
            .zip(lm_b1.outer_iter_mut())
            .zip(lm_b2.outer_iter_mut())
            .map(|(((((a, b), c), d), e), f)| [a, b, c, d, e, f]);
        let rows: Vec<_> = xs_total
            .iter_mut()
            .zip(xs_scattering.iter_mut())
            .zip(pdf.outer_iter())
            .zip(phase_rows)
            .zip(coefficient_rows)
            .map(
                |((((xs_total, xs_scattering), pdf), phase), coefficients)| DistributionRow {
                    xs_total,
                    xs_scattering,
                    pdf,
                    phase,
                    coefficients,
                },
            )
            .collect();
        let integrate_row = |row: DistributionRow<'_>| {
            let DistributionRow {
                xs_total,
                xs_scattering,
                pdf,
                phase: [mut p11, mut p12, mut p33, mut p34],
                coefficients,
            } = row;
            Zip::indexed(size_param)
                .and(size_weights)
                .and(pdf)
                .for_each(|j, &size_param, &size_weight, &pdf| {
                    let area_weight =
                        size_weight * pdf * std::f64::consts::PI * (size_param / k).powi(2);
                    *xs_total += area_weight * output.Qext[j];
                    *xs_scattering += area_weight * output.Qsca[j];
                });
            let phase_norm = c / *xs_scattering;
            Zip::indexed(size_weights)
                .and(pdf)
                .for_each(|j, &size_weight, &pdf| {
                    let weight = phase_norm * size_weight * pdf;
                    Zip::from(output.S1.row(j))
                        .and(output.S2.row(j))
                        .and(&mut p11)
                        .and(&mut p12)
                        .and(&mut p33)
                        .and(&mut p34)
                        .for_each(|s1, s2, p11, p12, p33, p34| {
                            let s1_norm = s1.norm_sqr();
                            let s2_norm = s2.norm_sqr();
                            let cross = s1 * s2.conj();
                            *p11 += weight * (s1_norm + s2_norm);
                            *p12 += weight * (s1_norm - s2_norm);
                            *p33 += weight * 2.0 * cross.re;
                            *p34 += weight * 2.0 * cross.im;
                        });
                });
            let projected = self.transform.project(
                [
                    p11.view(),
                    p12.view(),
                    p11.view(),
                    p33.view(),
                    p34.view(),
                    p33.view(),
                ],
                angle_weights,
            );
            for (mut output, values) in coefficients.into_iter().zip(projected) {
                output += &values;
            }
        };
        if let Some(pool) = &self.pool {
            pool.install(|| rows.into_par_iter().for_each(integrate_row));
        } else {
            rows.into_iter().for_each(integrate_row);
        }
        Ok(())
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use ndarray::{Array1, Array2, array};

    #[test]
    fn mie_integrator_accumulates_a4_from_p33() {
        let cos_angles = array![-0.8, 0.0, 0.8];
        let angle_weights = array![0.5, 1.0, 0.5];
        let integrator = MieIntegrator::new(cos_angles.view(), 1, 1).unwrap();

        let mut xs_total = Array1::zeros(1);
        let mut xs_scattering = Array1::zeros(1);
        let mut p11 = Array2::zeros((1, 3));
        let mut p12 = Array2::zeros((1, 3));
        let mut p33 = Array2::zeros((1, 3));
        let mut p34 = Array2::zeros((1, 3));
        let mut lm_a1 = Array2::zeros((1, 1));
        let mut lm_a2 = Array2::zeros((1, 1));
        let mut lm_a3 = Array2::zeros((1, 1));
        let mut lm_a4 = Array2::zeros((1, 1));
        let mut lm_b1 = Array2::zeros((1, 1));
        let mut lm_b2 = Array2::zeros((1, 1));

        integrator
            .integrate(
                500.0,
                Complex64::new(1.5, 0.01),
                array![2.0].view(),
                array![[1.0]].view(),
                array![1.0].view(),
                angle_weights.view(),
                xs_total.view_mut(),
                xs_scattering.view_mut(),
                p11.view_mut(),
                p12.view_mut(),
                p33.view_mut(),
                p34.view_mut(),
                lm_a1.view_mut(),
                lm_a2.view_mut(),
                lm_a3.view_mut(),
                lm_a4.view_mut(),
                lm_b1.view_mut(),
                lm_b2.view_mut(),
            )
            .unwrap();

        let expected_a4 = 0.5 * p33.row(0).dot(&angle_weights);
        assert!((lm_a4[[0, 0]] - expected_a4).abs() < 1e-13);
        assert!(lm_a4[[0, 0]].abs() > 0.0);
    }
}
