#![allow(non_snake_case)]

use num::complex::Complex64;

use crate::math::wigner::WignerDCalculator;
use crate::prelude::*;

use super::mie_f;

pub struct MieIntegrator {
    cos_angles: Array1<f64>,
    _num_legendre: usize,
    _num_threads: usize,

    lpoly_00: Array2<f64>,
    lpoly_22: Array2<f64>,
    lpoly_2m2: Array2<f64>,
    lpoly_02: Array2<f64>,
}

impl MieIntegrator {
    pub fn new(
        cos_angles: ArrayView1<f64>,
        num_legendre: usize,
        num_threads: usize,
    ) -> Result<Self> {
        // Construct wigner polynomials
        let wigner_00 = WignerDCalculator::new(0, 0);
        let wigner_22 = WignerDCalculator::new(2, 2);
        let wigner_2m2 = WignerDCalculator::new(2, -2);
        let wigner_02 = WignerDCalculator::new(0, 2);

        let mut lpoly_00 = Array2::zeros((cos_angles.len(), num_legendre));
        let mut lpoly_22 = Array2::zeros((cos_angles.len(), num_legendre));
        let mut lpoly_2m2 = Array2::zeros((cos_angles.len(), num_legendre));
        let mut lpoly_02 = Array2::zeros((cos_angles.len(), num_legendre));

        Zip::from(cos_angles)
            .and(lpoly_00.axis_iter_mut(Axis(0)))
            .and(lpoly_22.axis_iter_mut(Axis(0)))
            .and(lpoly_2m2.axis_iter_mut(Axis(0)))
            .and(lpoly_02.axis_iter_mut(Axis(0)))
            .for_each(
                |cos_theta, mut lpoly_00, mut lpoly_22, mut lpoly_2m2, mut lpoly_02| {
                    let theta = (*cos_theta).acos();
                    wigner_00.vector_d(theta, lpoly_00.as_slice_mut().unwrap());
                    wigner_22.vector_d(theta, lpoly_22.as_slice_mut().unwrap());
                    wigner_2m2.vector_d(theta, lpoly_2m2.as_slice_mut().unwrap());
                    wigner_02.vector_d(theta, lpoly_02.as_slice_mut().unwrap());
                },
            );

        Ok(Self {
            cos_angles: cos_angles.into_owned(),
            _num_legendre: num_legendre,
            _num_threads: num_threads,
            lpoly_00,
            lpoly_22,
            lpoly_2m2,
            lpoly_02,
        })
    }

    #[allow(clippy::too_many_arguments)]
    pub fn integrate(
        &self,
        wavelength: f64,
        refractive_index: Complex64,
        size_param: ArrayView1<f64>,       // [size_param]
        pdf: ArrayView2<f64>,              // [size_param, distribution]
        size_weights: ArrayView1<f64>,     // [size_param]
        angle_weights: ArrayView1<f64>,    // [angle]
        xs_total: ArrayViewMut1<f64>,      // [distribution]
        xs_scattering: ArrayViewMut1<f64>, // [distribution]
        mut p11: ArrayViewMut2<f64>,       // [distribution, angle]
        mut p12: ArrayViewMut2<f64>,       // [distribution, angle]
        mut p33: ArrayViewMut2<f64>,       // [distribution, angle]
        mut p34: ArrayViewMut2<f64>,       // [distribution, angle]
        mut lm_a1: ArrayViewMut2<f64>,     // [distribution, legendre]
        mut lm_a2: ArrayViewMut2<f64>,     // [distribution, legendre]
        mut lm_a3: ArrayViewMut2<f64>,     // [distribution, legendre]
        mut _lm_a4: ArrayViewMut2<f64>,    // [distribution, legendre]
        mut lm_b1: ArrayViewMut2<f64>,     // [distribution, legendre]
        mut lm_b2: ArrayViewMut2<f64>,     // [distribution, legendre]
    ) -> Result<()> {
        let k = 2.0 * std::f64::consts::PI / wavelength;
        let c = 4.0 * std::f64::consts::PI / (2.0 * k * k);

        // Alloc arrays for the Mie parameters and do the calculation
        let output = mie_f::mie(size_param, refractive_index, self.cos_angles.view());

        Zip::indexed(xs_total)
            .and(xs_scattering)
            .for_each(|i, xs_total, xs_scattering| {
                // Index everything that we need
                let pdf = pdf.index_axis(Axis(0), i);
                let mut p11 = p11.index_axis_mut(Axis(0), i);
                let mut p12 = p12.index_axis_mut(Axis(0), i);
                let mut p33 = p33.index_axis_mut(Axis(0), i);
                let mut p34 = p34.index_axis_mut(Axis(0), i);
                let mut lm_a1 = lm_a1.index_axis_mut(Axis(0), i);
                let mut lm_a2 = lm_a2.index_axis_mut(Axis(0), i);
                let mut lm_a3 = lm_a3.index_axis_mut(Axis(0), i);
                let mut lm_b1 = lm_b1.index_axis_mut(Axis(0), i);
                let mut lm_b2 = lm_b2.index_axis_mut(Axis(0), i);

                // First we go through the size parameters to calculate the cross sections
                Zip::indexed(size_param)
                    .and(size_weights)
                    .and(pdf)
                    .for_each(|j, &size_param, &size_weight, &pdf| {
                        *xs_total += size_weight
                            * pdf
                            * std::f64::consts::PI
                            * output.Qext[j]
                            * (size_param * wavelength / (2.0 * std::f64::consts::PI)).powf(2.0);

                        *xs_scattering += size_weight
                            * pdf
                            * std::f64::consts::PI
                            * output.Qsca[j]
                            * (size_param * wavelength / (2.0 * std::f64::consts::PI)).powf(2.0);
                    });

                // And the phase normalization factor
                let phase_norm = c / *xs_scattering;

                // Then we go through and calculate the phase functions
                Zip::indexed(size_weights)
                    .and(pdf)
                    .for_each(|j, &size_weight, &pdf| {
                        let S1 = output.S1.slice(ndarray::s![j, ..]);
                        let S2 = output.S2.slice(ndarray::s![j, ..]);

                        Zip::from(S1)
                            .and(S2)
                            .and(&mut p11)
                            .and(&mut p12)
                            .and(&mut p33)
                            .and(&mut p34)
                            .for_each(|s1, s2, p11, p12, p33, p34| {
                                *p11 += phase_norm
                                    * size_weight
                                    * pdf
                                    * (s1.norm_sqr() + s2.norm_sqr());
                                *p12 += phase_norm
                                    * size_weight
                                    * pdf
                                    * (s1.norm_sqr() - s2.norm_sqr());
                                *p33 += phase_norm
                                    * size_weight
                                    * pdf
                                    * (s1 * s2.conj() + s2 * s1.conj()).re;
                                *p34 += phase_norm
                                    * size_weight
                                    * pdf
                                    * (s1 * s2.conj() - s2 * s1.conj()).im;
                            });
                    });

                // Then with the phase function calculated we can calculate the legendre polynomials
                Zip::indexed(&mut lm_a1)
                    .and(&mut lm_a2)
                    .and(&mut lm_a3)
                    .and(&mut lm_b1)
                    .and(&mut lm_b2)
                    .for_each(|l, lm_a1, lm_a2, lm_a3, lm_b1, lm_b2| {
                        let l_weight = 1.0 / (2.0 / (2.0 * l as f64 + 1.0));

                        Zip::indexed(angle_weights).for_each(|j, w| {
                            let w = *w * l_weight;
                            let lpoly_00 = self.lpoly_00[[j, l]];
                            let lpoly_22 = self.lpoly_22[[j, l]];
                            let lpoly_2m2 = self.lpoly_2m2[[j, l]];
                            let lpoly_02 = self.lpoly_02[[j, l]];

                            let p11 = p11[[j]];
                            let p12 = p12[[j]];
                            let p33 = p33[[j]];
                            let p34 = p34[[j]];

                            let temp1 = w * lpoly_22 * (p11 + p33);
                            let temp2 = w * lpoly_2m2 * (p11 - p33);

                            *lm_a1 += w * lpoly_00 * p11;
                            *lm_a2 += (temp1 + temp2) / 2.0;
                            *lm_a3 += (temp1 - temp2) / 2.0;

                            *lm_b1 += w * lpoly_02 * p12;
                            *lm_b2 += -w * lpoly_02 * p34;
                        });
                    });
            });
        Ok(())
    }
}
