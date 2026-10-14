#![allow(non_snake_case)]

//! Mie scattering by homogeneous spheres, computed by `scatrs`.
//!
//! sasktran2 writes the refractive index as m = n - ik (k >= 0 for absorbing particles), the
//! convention of Wiscombe's MIEV0, while scatrs uses m = n + ik. The refractive index is
//! conjugated on the way in, and the amplitudes S1, S2 on the way out, so that every output
//! here keeps the sasktran2 convention.

use ndarray::{Array1, Array2, ArrayView1, ArrayViewMut1, Zip};
use num::complex::Complex64;
use rayon::prelude::*;
use scatrs::mie::{
    AmplitudeBuffers, AngularTable, MAX_MIE_ORDER, MIN_SIZE_PARAMETER, Mie, MieOptions,
    default_n_max,
};

use crate::prelude::*;

pub struct MieOutput {
    pub Qext: Array1<f64>,
    pub Qsca: Array1<f64>,
    pub S1: Array2<Complex64>,
    pub S2: Array2<Complex64>,
    pub size_param: Array1<f64>,
    pub refractive_index: Complex64,
    pub cos_angles: Array1<f64>,
}

impl MieOutput {
    pub fn new(
        size_param: ArrayView1<f64>,
        refractive_index: Complex64,
        cos_angles: ArrayView1<f64>,
    ) -> Self {
        MieOutput {
            Qext: Array1::zeros(size_param.len()),
            Qsca: Array1::zeros(size_param.len()),
            S1: Array2::zeros((size_param.len(), cos_angles.len())),
            S2: Array2::zeros((size_param.len(), cos_angles.len())),
            size_param: size_param.to_owned(),
            refractive_index,
            cos_angles: cos_angles.to_owned(),
        }
    }
}

/// sasktran2 has never restricted the inputs to the domain where scatrs has validated its
/// accuracy, so neither does this.
fn options() -> MieOptions {
    MieOptions::new().with_allow_unvalidated(true)
}

/// One sphere, writing Qext, Qsca and the amplitudes, all in the sasktran2 convention.
#[allow(clippy::too_many_arguments)]
fn calculate_sphere(
    mie: &mut Mie<f64>,
    table: &AngularTable<f64>,
    size_param: f64,
    refractive_index: Complex64,
    Qext: &mut f64,
    Qsca: &mut f64,
    mut S1: ArrayViewMut1<Complex64>,
    mut S2: ArrayViewMut1<Complex64>,
) -> Result<()> {
    S1.fill(Complex64::new(0.0, 0.0));
    S2.fill(Complex64::new(0.0, 0.0));
    // Below scatrs's smallest size parameter the efficiencies (~x^4) underflow; a sphere of
    // zero size does not scatter.
    if (0.0..MIN_SIZE_PARAMETER).contains(&size_param) {
        *Qext = 0.0;
        *Qsca = 0.0;
        return Ok(());
    }
    let sphere = mie
        .compute(size_param, refractive_index.conj(), &options())
        .map_err(|err| {
            anyhow!(
                "Mie calculation failed at size parameter {size_param} and refractive index \
                 {refractive_index} (sasktran2 convention m = n - ik): {err}"
            )
        })?;
    let efficiencies = sphere.efficiencies();
    *Qext = efficiencies.q_ext;
    *Qsca = efficiencies.q_sca;

    let (Some(s1), Some(s2)) = (S1.as_slice_mut(), S2.as_slice_mut()) else {
        anyhow::bail!("Mie amplitude outputs must be contiguous");
    };
    sphere.amplitudes_from_table(table, AmplitudeBuffers { s1, s2 })?;
    S1.mapv_inplace(|s| s.conj());
    S2.mapv_inplace(|s| s.conj());
    Ok(())
}

pub fn mie(
    size_param: ArrayView1<f64>,
    refractive_index: Complex64,
    cos_angles: ArrayView1<f64>,
) -> Result<MieOutput> {
    mie_with_threads(size_param, refractive_index, cos_angles, false)
}

/// Parallel work uses the caller's Rayon pool; the standalone Mie API stays serial.
pub(super) fn mie_with_threads(
    size_param: ArrayView1<f64>,
    refractive_index: Complex64,
    cos_angles: ArrayView1<f64>,
    parallel: bool,
) -> Result<MieOutput> {
    let mut output = MieOutput::new(size_param, refractive_index, cos_angles);
    // The angular functions are shared by every size parameter; the table needs the terms of
    // the largest sphere, and at least one.
    let num_orders = size_param
        .iter()
        .filter(|x| x.is_finite())
        .map(|&x| default_n_max(x))
        .max()
        .unwrap_or(1)
        .clamp(1, MAX_MIE_ORDER);
    let table = AngularTable::new(&cos_angles.to_vec(), num_orders)?;

    let rows = Zip::from(size_param)
        .and(&mut output.Qext)
        .and(&mut output.Qsca)
        .and(output.S1.outer_iter_mut())
        .and(output.S2.outer_iter_mut());
    let calculate = |mie: &mut Mie<f64>, (&x, Qext, Qsca, S1, S2)| {
        calculate_sphere(mie, &table, x, refractive_index, Qext, Qsca, S1, S2)
    };
    if parallel {
        rows.into_par_iter()
            .with_min_len(16)
            .try_for_each_init(Mie::new, calculate)?;
    } else {
        let mut mie = Mie::new();
        let mut result = Ok(());
        rows.for_each(|x, Qext, Qsca, S1, S2| {
            if result.is_ok() {
                result = calculate(&mut mie, (x, Qext, Qsca, S1, S2));
            }
        });
        result?;
    }

    Ok(output)
}

#[cfg(test)]
mod tests {
    use super::*;
    use ndarray::array;

    fn default_cos_angles() -> Array1<f64> {
        (0..7)
            .map(|i| (i as f64 * std::f64::consts::PI / 6.0).cos())
            .collect()
    }

    fn efficiencies(size_param: f64, refractive_index: Complex64) -> (f64, f64) {
        let empty = Array1::zeros(0);
        let output = mie(array![size_param].view(), refractive_index, empty.view()).unwrap();
        (output.Qext[0], output.Qsca[0])
    }

    /// MIEV0 prints 7 significant digits: allow one unit in the last, or `atol`.
    fn assert_amplitudes(output: &MieOutput, expected: &[[f64; 4]], atol: f64) {
        for (j, row) in expected.iter().enumerate() {
            let (s1, s2) = (output.S1[[0, j]], output.S2[[0, j]]);
            for (value, expected) in [s1.re, s1.im, s2.re, s2.im].into_iter().zip(row) {
                let tolerance = (1e-6 * expected.abs()).max(atol);
                assert!(
                    (value - expected).abs() < tolerance,
                    "angle {j}: S1 = {s1}, S2 = {s2}, expected {row:?}"
                );
            }
        }
    }

    #[test]
    fn test_miev0_efficiencies() {
        // Wiscombe's MIEV0 test cases, (size parameter, refractive index, Qsca).
        let cases = [
            (10.0, Complex64::new(1.5, 0.0), 2.8820, 1e-4),
            (100.0, Complex64::new(1.5, 0.0), 2.0944, 1e-4),
            (1000.0, Complex64::new(1.5, 0.0), 2.0139, 1e-4),
            (5000.0, Complex64::new(1.5, 0.0), 2.0086, 1e-4),
            (0.101, Complex64::new(0.75, 0.0), 0.000008, 1e-6),
            (10.0, Complex64::new(0.75, 0.0), 2.232265, 1e-6),
            (1000.0, Complex64::new(0.75, 0.0), 1.997908, 1e-6),
        ];
        for (size_param, refractive_index, expected, tolerance) in cases {
            let (Qext, Qsca) = efficiencies(size_param, refractive_index);
            assert!(
                (Qsca - expected).abs() < tolerance,
                "x = {size_param}: {Qsca}"
            );
            // Non-absorbing: extinction is scattering.
            assert!((Qext - Qsca).abs() < 1e-12 * Qext, "x = {size_param}");
        }
    }

    #[test]
    fn test_absorbing_uses_n_minus_ik() {
        let size_param = 2.0 * std::f64::consts::PI * 0.525 / 0.6328;

        let (Qext, Qsca) = efficiencies(size_param, Complex64::new(1.55, -0.1));
        assert!((Qext - 2.86165188243).abs() < 1e-7);
        assert!((Qsca - 1.66424911991).abs() < 1e-7);

        let (Qext, Qsca) = efficiencies(size_param, Complex64::new(1.55, 0.0));
        assert!((Qext - 3.10543).abs() < 1e-5);
        assert!((Qsca - 3.10543).abs() < 1e-5);
    }

    #[test]
    fn test_gain_medium_is_an_error() {
        let empty = Array1::zeros(0);
        let result = mie(array![1.0].view(), Complex64::new(1.5, 0.01), empty.view());
        assert!(result.is_err());
    }

    #[test]
    fn test_zero_size_does_not_scatter() {
        let cos_angles = default_cos_angles();
        let output = mie(
            array![0.0, 1.0].view(),
            Complex64::new(1.5, -0.01),
            cos_angles.view(),
        )
        .unwrap();
        assert_eq!(output.Qext[0], 0.0);
        assert_eq!(output.Qsca[0], 0.0);
        assert!(output.S1.row(0).iter().all(|s| s.norm() == 0.0));
        assert!(output.Qext[1] > 0.0);
    }

    #[test]
    fn test_parallel_matches_serial() {
        let cos_angles = default_cos_angles();
        let size_param = Array1::geomspace(1e-3, 50.0, 100).unwrap();
        let refractive_index = Complex64::new(1.5, -0.01);
        let serial = mie(size_param.view(), refractive_index, cos_angles.view()).unwrap();
        let pool = crate::util::create_pool(4).unwrap();
        let parallel = pool
            .install(|| {
                mie_with_threads(size_param.view(), refractive_index, cos_angles.view(), true)
            })
            .unwrap();
        assert_eq!(serial.Qext, parallel.Qext);
        assert_eq!(serial.Qsca, parallel.Qsca);
        assert_eq!(serial.S1, parallel.S1);
        assert_eq!(serial.S2, parallel.S2);
    }

    #[test]
    fn test_miev0_case_5() {
        let cos_angles = default_cos_angles();
        let output = mie(
            array![0.099].view(),
            Complex64::new(0.75, 0.0),
            cos_angles.view(),
        )
        .unwrap();

        assert!((output.Qsca[0] - 0.000007).abs() < 1e-5);
        assert_amplitudes(
            &output,
            &[
                [1.817558e-8, -1.654225e-4, 1.817558e-8, -1.654225e-4],
                [1.817558e-8, -1.653815e-4, 1.574051e-8, -1.432172e-4],
                [1.817558e-8, -1.652694e-4, 9.087788e-9, -8.261265e-5],
                [1.817558e-8, -1.651163e-4, 9.797186e-23, 2.938374e-8],
                [1.817558e-8, -1.649631e-4, -9.087788e-9, 8.250360e-5],
                [1.817558e-8, -1.648510e-4, -1.574051e-8, 1.427725e-4],
                [1.817558e-8, -1.648100e-4, -1.817558e-8, 1.648100e-4],
            ],
            1e-8,
        );
    }

    #[test]
    fn test_miev0_case_8() {
        let cos_angles = default_cos_angles();
        let output = mie(
            array![1000.0].view(),
            Complex64::new(0.75, 0.0),
            cos_angles.view(),
        )
        .unwrap();

        assert!((output.Qsca[0] - 1.997908).abs() < 1e-6);
        assert_amplitudes(
            &output,
            &[
                [4.994770e5, -1.336502e4, 4.994770e5, -1.336502e4],
                [-3.999296e2, -3.316361e2, -3.946018e2, -1.147791e2],
                [-5.209852e2, -5.776614e2, -1.970767e2, -6.937470e2],
                [-1.600887e2, 1.348013e2, -4.152365e1, 1.143000e2],
                [8.431720e1, -1.209493e2, -4.261732e1, 5.535055e1],
                [-7.556092e1, -8.134810e1, 4.218303e1, 9.100831e1],
                [1.705778e1, 4.842510e2, -1.705778e1, -4.842510e2],
            ],
            0.0,
        );
    }
}
