//! Mie optical properties of size distributions, integrated by `scatrs`.
//!
//! [`integrate_batch`] averages over many distributions of one material at once, sharing the
//! size nodes (`scatrs::ensemble::compute_batch`): the tables of a `MieDatabase`. The size
//! integral is adaptive, refined until the cross sections and the asymmetry parameter meet
//! `rtol`.
//!
//! The outputs are in the sasktran2 conventions: the refractive index is m = n - ik (scatrs
//! uses n + ik), and the Greek coefficients beta1 and beta2, and the phase-matrix element p12,
//! have the opposite sign to scatrs's.

use ndarray::{Array2, Array3, ArrayView1};
use num::complex::Complex64;
use scatrs::ensemble::{EnsembleOptions, GreekOrder, compute_batch};

use crate::prelude::*;

pub use scatrs::ensemble::SizeDistribution;

/// The averages over each distribution at each wavelength.
pub struct EnsembleTables {
    /// Cross sections per particle, in the square of the unit of the radii, [wavelength, distribution].
    pub xs_total: Array2<f64>,
    pub xs_scattering: Array2<f64>,
    pub xs_absorption: Array2<f64>,
    /// Normalized Greek coefficients a1, a2, a3, a4, b1, b2 (a1 = 1 at order 0),
    /// [wavelength, distribution, legendre].
    pub greek: [Array3<f64>; 6],
    /// Whether the size integration met `rtol`, [wavelength, distribution].
    pub converged: Array2<bool>,
    /// The normalized phase-matrix elements p11, p12, p33, p34 (half the integral of p11 over
    /// the cosine is 1) at the requested cosines, [wavelength, distribution, angle]; `None`
    /// if no cosines were requested.
    pub phase: Option<[Array3<f64>; 4]>,
}

/// Averages over each of `distributions` of spheres with refractive index
/// `refractive_index[i]` (sasktran2 convention, m = n - ik) at `wavelengths[i]`, in the unit
/// of the radii, with `num_coefficients` Greek coefficients and, if `cos_angles` is not empty,
/// the phase matrix at those cosines of the scattering angle. Uses a private pool of
/// `num_threads` (0 for Rayon's automatic count); the results do not depend on it.
#[allow(clippy::too_many_arguments)]
pub fn integrate_batch(
    distributions: &[SizeDistribution],
    wavelengths: ArrayView1<f64>,
    refractive_index: ArrayView1<Complex64>,
    num_coefficients: usize,
    cos_angles: ArrayView1<f64>,
    rtol: f64,
    max_panels: usize,
    num_threads: usize,
) -> Result<EnsembleTables> {
    anyhow::ensure!(
        num_coefficients >= 1,
        "at least one Greek coefficient is needed"
    );
    anyhow::ensure!(
        wavelengths.len() == refractive_index.len(),
        "{} refractive indices for {} wavelengths",
        refractive_index.len(),
        wavelengths.len()
    );
    if let Some(m) = refractive_index.iter().find(|m| m.im > 0.0) {
        anyhow::bail!(
            "refractive index {m} has a positive imaginary part, a gain medium in the sasktran2 \
             convention m = n - ik"
        );
    }
    let wavelengths = wavelengths.to_vec();
    let refractive_index: Vec<_> = refractive_index.iter().map(|m| m.conj()).collect();
    // sasktran2 has never restricted the sizes to scatrs's validated domain.
    let options = EnsembleOptions::new()
        .with_rtol(rtol)
        .with_max_panels(max_panels)
        .with_greek(GreekOrder::Order(num_coefficients - 1))
        .with_mu(cos_angles.to_vec())
        .with_allow_unvalidated(true);

    // Even one thread needs a pool: scatrs otherwise runs on the global pool.
    let pool = crate::util::create_pool(num_threads)?;
    let results =
        pool.install(|| compute_batch(distributions, &wavelengths, &refractive_index, &options))?;

    let shape = (wavelengths.len(), distributions.len());
    let mut tables = EnsembleTables {
        xs_total: Array2::zeros(shape),
        xs_scattering: Array2::zeros(shape),
        xs_absorption: Array2::zeros(shape),
        greek: std::array::from_fn(|_| Array3::zeros((shape.0, shape.1, num_coefficients))),
        converged: Array2::from_elem(shape, false),
        phase: (!cos_angles.is_empty())
            .then(|| std::array::from_fn(|_| Array3::zeros((shape.0, shape.1, cos_angles.len())))),
    };
    for (d, by_wavelength) in results.iter().enumerate() {
        for (w, result) in by_wavelength.iter().enumerate() {
            let values = &result.values;
            tables.xs_total[[w, d]] = values.c_ext;
            tables.xs_scattering[[w, d]] = values.c_sca;
            tables.xs_absorption[[w, d]] = values.c_abs;
            tables.converged[[w, d]] = result.converged;
            let greek = values
                .greek
                .as_ref()
                .ok_or_else(|| anyhow!("scatrs returned no Greek coefficients"))?;
            let coefficients = [
                (&greek.alpha1, 1.0),
                (&greek.alpha2, 1.0),
                (&greek.alpha3, 1.0),
                (&greek.alpha4, 1.0),
                (&greek.beta1, -1.0),
                (&greek.beta2, -1.0),
            ];
            for (table, (values, sign)) in tables.greek.iter_mut().zip(coefficients) {
                for (l, value) in values.iter().take(num_coefficients).enumerate() {
                    table[[w, d, l]] = sign * value;
                }
            }
            if let Some(tables) = &mut tables.phase {
                let phase = values
                    .phase
                    .as_ref()
                    .ok_or_else(|| anyhow!("scatrs returned no phase matrix"))?;
                let elements = [
                    (&phase.f11, 1.0),
                    (&phase.f12, -1.0),
                    (&phase.f33, 1.0),
                    (&phase.f34, 1.0),
                ];
                for (table, (values, sign)) in tables.iter_mut().zip(elements) {
                    for (j, value) in values.iter().enumerate() {
                        table[[w, d, j]] = sign * value;
                    }
                }
            }
        }
    }
    Ok(tables)
}

#[cfg(test)]
mod tests {
    use super::*;
    use ndarray::{Array1, array};

    fn lognormal(median_radius: f64, width: f64) -> SizeDistribution {
        SizeDistribution::LogNormal {
            median_radius,
            width,
        }
    }

    #[test]
    fn test_batch_layout_and_conventions() {
        let distributions = [lognormal(100.0, 1.5), lognormal(200.0, 1.6)];
        let wavelengths = array![400.0, 700.0, 1000.0];
        let m = Complex64::new(1.45, -0.01);
        let tables = integrate_batch(
            &distributions,
            wavelengths.view(),
            array![m, m, m].view(),
            16,
            array![0.5, -0.5].view(),
            1e-6,
            4096,
            1,
        )
        .unwrap();
        assert_eq!(tables.xs_total.dim(), (3, 2));
        assert_eq!(tables.greek[0].dim(), (3, 2, 16));
        assert!(tables.converged.iter().all(|&c| c));

        // Each entry is the single distribution computed on its own, within the tolerance,
        // with beta flipped and the refractive index conjugated.
        let single = scatrs::ensemble::compute(
            &distributions[1],
            &[700.0],
            &[m.conj()],
            &EnsembleOptions::new()
                .with_greek(GreekOrder::Order(15))
                .with_mu(vec![0.5, -0.5]),
        )
        .unwrap();
        let single = &single[0].values;
        let greek = single.greek.as_ref().unwrap();
        assert!((tables.xs_total[[1, 1]] / single.c_ext - 1.0).abs() < 1e-5);
        assert!((tables.xs_scattering[[1, 1]] / single.c_sca - 1.0).abs() < 1e-5);
        assert!((tables.xs_absorption[[1, 1]] / single.c_abs - 1.0).abs() < 1e-5);
        assert_eq!(tables.greek[0][[1, 1, 0]], 1.0);
        assert!((tables.greek[0][[1, 1, 1]] - greek.alpha1[1]).abs() < 1e-5);
        assert!((tables.greek[4][[1, 1, 2]] + greek.beta1[2]).abs() < 1e-5);
        assert!((tables.greek[5][[1, 1, 2]] + greek.beta2[2]).abs() < 1e-5);
        let phase = tables.phase.as_ref().unwrap();
        assert_eq!(phase[0].dim(), (3, 2, 2));
        let single_phase = single.phase.as_ref().unwrap();
        for j in 0..2 {
            let p11 = single_phase.f11[j];
            assert!((phase[0][[1, 1, j]] / p11 - 1.0).abs() < 1e-5);
            assert!((phase[1][[1, 1, j]] + single_phase.f12[j]).abs() < 1e-5 * p11);
            assert!((phase[2][[1, 1, j]] - single_phase.f33[j]).abs() < 1e-5 * p11);
            assert!((phase[3][[1, 1, j]] - single_phase.f34[j]).abs() < 1e-5 * p11);
        }
    }

    #[test]
    fn test_threads_do_not_change_results() {
        let distributions = [lognormal(150.0, 1.4)];
        let wavelengths = array![350.0, 750.0];
        let m = array![Complex64::new(1.43, -1e-7), Complex64::new(1.42, -1e-7)];
        let run = |threads| {
            integrate_batch(
                &distributions,
                wavelengths.view(),
                m.view(),
                8,
                Array1::zeros(0).view(),
                1e-4,
                4096,
                threads,
            )
            .unwrap()
        };
        let (serial, parallel) = (run(1), run(4));
        assert_eq!(serial.xs_total, parallel.xs_total);
        assert_eq!(serial.greek[0], parallel.greek[0]);
    }

    #[test]
    fn test_gain_medium_is_an_error() {
        let result = integrate_batch(
            &[lognormal(100.0, 1.5)],
            array![500.0].view(),
            array![Complex64::new(1.5, 0.01)].view(),
            4,
            Array1::zeros(0).view(),
            1e-4,
            4096,
            1,
        );
        assert!(result.is_err());
    }
}
