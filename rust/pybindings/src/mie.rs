use core::num;
use std::slice;

use sk_core::mie as sk_mie;
use sk_core::wigner::WignerDCalculator;

use numpy::{Complex64, PyArrayMethods};
use numpy::ndarray::{Array1, Array2, Axis, Zip};
use numpy::{PyArray1, PyReadonlyArray1, PyReadonlyArray2, PyReadwriteArray1, PyReadwriteArray2};
use pyo3::prelude::*;

#[pyclass]
pub struct Mie {
    mie: sk_mie::Mie,
    lpoly_00: Array2<f64>,
    lpoly_22: Array2<f64>,
    lpoly_2m2: Array2<f64>,
    lpoly_02: Array2<f64>,
}

#[pymethods]
impl Mie {
    #[new]
    pub fn new<'py>(
        cos_angles: PyReadonlyArray1<'py, f64>,
        num_coeff: usize,
        num_threads: usize,
    ) -> Self {
        let num_angles = cos_angles.as_array().len();

        println!("num_angles: {}", num_angles);

        let wigner00 = WignerDCalculator::new(0, 0);
        let wigner22 = WignerDCalculator::new(2, 2);
        let wigner2m2 = WignerDCalculator::new(2, -2);
        let wigner02 = WignerDCalculator::new(0, 2);

        let mut lpoly_00 = Array2::zeros((num_angles, num_coeff));
        let mut lpoly_22 = Array2::zeros((num_angles, num_coeff));
        let mut lpoly_2m2 = Array2::zeros((num_angles, num_coeff));
        let mut lpoly_02 = Array2::zeros((num_angles, num_coeff));

        Zip::from(lpoly_00.rows_mut())
            .and(lpoly_02.rows_mut())
            .and(lpoly_22.rows_mut())
            .and(lpoly_2m2.rows_mut())
            .and(cos_angles.as_array())
            .for_each(|mut lpoly00_row, mut lpoly02_row, mut lpoly22_row, mut lpoly2m2_row, cos_angle| {
                let theta = cos_angle.acos();
                wigner00.vector_d(theta, lpoly00_row.as_slice_mut().unwrap());
                wigner22.vector_d(theta, lpoly22_row.as_slice_mut().unwrap());
                wigner2m2.vector_d(theta, lpoly2m2_row.as_slice_mut().unwrap());
                wigner02.vector_d(theta, lpoly02_row.as_slice_mut().unwrap());
            });

        Mie {
            mie: sk_mie::Mie::new().with_cos_angles(cos_angles.as_array().to_vec()),
            lpoly_00: lpoly_00,
            lpoly_22: lpoly_22,
            lpoly_2m2: lpoly_2m2,
            lpoly_02: lpoly_02,
        }
    }

    pub fn integrate<'py>(
        &mut self,
        wavelength: f64,
        refractive_index_real: f64,
        refractive_index_imag: f64,
        size_param: PyReadonlyArray1<'py, f64>, // [size_param]
        pdf: PyReadonlyArray2<'py, f64>,        // [distribution, size_Param]
        size_weights: PyReadonlyArray1<'py, f64>, // [size_param]
        angle_weights: PyReadonlyArray1<'py, f64>, // [angle]
        mut xs_total: PyReadwriteArray1<'py, f64>,  // [distribution]
        mut xs_scattering: PyReadwriteArray1<'py, f64>, // [distribution]
        p11: PyReadwriteArray2<'py, f64>,       // [distribution, angle]
        p12: PyReadwriteArray2<'py, f64>,       // [distribution, angle]
        p33: PyReadwriteArray2<'py, f64>,       // [distribution, angle]
        p34: PyReadwriteArray2<'py, f64>,       // [distribution, angle]
        lm_a1: PyReadwriteArray2<'py, f64>,     // [distribution, legendre]
        lm_a2: PyReadwriteArray2<'py, f64>,     // [distribution, legendre]
        lm_a3: PyReadwriteArray2<'py, f64>,     // [distribution, legendre]
        lm_a4: PyReadwriteArray2<'py, f64>,     // [distribution, legendre]
        lm_b1: PyReadwriteArray2<'py, f64>,     // [distribution, legendre]
        lm_b2: PyReadwriteArray2<'py, f64>,     // [distribution, legendre]
    ) -> PyResult<()> {
        let complex_refractive = Complex64::new(refractive_index_real, refractive_index_imag);

        let mut xs_total = xs_total.as_array_mut();
        let mut xs_scattering = xs_scattering.as_array_mut();

        let mut Qext = Array1::zeros(xs_total.len());
        let mut Qsca = Array1::zeros(xs_total.len());

        Zip::from(size_param.as_array())
            .and(Qext.view_mut())
            .and(Qsca.view_mut())
            .for_each(|&size_param, Qe, Qs| {
                (*Qe, *Qs) = self.mie.calculate(size_param, complex_refractive);
            });

        for (i, pdf_row) in pdf.as_array().rows().into_iter().enumerate() {
            let pdf_row = pdf_row.as_slice().unwrap();

            Zip::from(Qext.view())
            .and(Qsca.view())
            .and(size_weights.as_array())
            .and(pdf_row)
            .for_each(|Qe, Qs, size_weight, pdf_val| {
                xs_total[i] += (Qe * size_weight * pdf_val);
                xs_scattering[i] += (Qs * size_weight * pdf_val);
            });

        }




        Ok(())
    }
}
