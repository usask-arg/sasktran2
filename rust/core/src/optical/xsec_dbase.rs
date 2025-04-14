use super::*;
use crate::interpolation::OutOfBoundsMode;
use crate::interpolation::linear::Grid;
use crate::interpolation::linear::Interp1Weights;
use crate::optical::read_fwf_xs::XsecListAtConditions;
use crate::util::*;
use anyhow::Result;
use ndarray::{Data, DataMut};
use ndarray::{Zip, prelude::*};

/// Trait for interpolating cross sections from a database.  This is used to interpolate
/// cross sections from a database to a given wavenumber and set of parameters.
pub trait XsecDatabaseInterp {
    fn xs_emplace<S1, S2, S3>(
        &self,
        wvnum: &ArrayBase<S1, Ix1>,
        params: &ArrayBase<S2, Ix1>,
        xs: &mut ArrayBase<S3, Ix1>,
        d_xs: Option<ArrayViewMut2<'_, f64>>,
    ) -> Result<(), String>
    where
        S1: Data<Elem = f64>,
        S2: Data<Elem = f64>,
        S3: DataMut<Elem = f64>;
}

/// A cross sectional database in the format that SASKTRAN requires it in.  This is
/// a cross sectional array of shape (len(param1), len(param2), ..., len(paramN), len(wvnum))
/// Where parameters are things like pressure, temperature, or foregin broadenerers, wvnum is always
/// The last dimension
pub struct SKXsecDatabase<D: Dimension> {
    xsec: Array<f64, D>,
    wvnum: Grid,
    params: Vec<Array1<f64>>,
}

impl<D: Dimension> SKXsecDatabase<D> {
    pub fn new(xsec: Array<f64, D>, wvnum: Grid, params: Vec<Array1<f64>>) -> Option<Self> {
        assert_eq!(params.len(), D::NDIM? - 1);

        Some(SKXsecDatabase {
            xsec,
            wvnum,
            params,
        })
    }
}

impl From<XsecListAtConditions> for SKXsecDatabase<Ix2> {
    fn from(db: XsecListAtConditions) -> Self {
        // We need to figure out how many dimensions we have
        let mut params: Vec<Array1<f64>> = vec![];
        let mut sidx: Vec<Vec<usize>> = vec![];
        for (_, vals) in &db.params {
            let unique_vals = unique_values(&Array1::from_vec(vals.clone()));

            if unique_vals.len() > 1 {
                params.push(Array1::from(unique_vals));
                sidx.push(argsort_f64(&vals));
            }
        }

        // Next we need to create a common wavenumber grid, concatenate all of the wavenumbers together
        // and then sort them
        let mut all_wvnum: Vec<f64> = vec![];
        for wvnum in &db.wvnum {
            all_wvnum.extend(wvnum.iter());
        }
        let unique_wvnum = unique_values(&Array1::from(all_wvnum));

        let mut xsec = Array2::zeros((params[0].len(), unique_wvnum.len()));

        for i in 0..params[0].len() {
            xsec.slice_mut(s![i, ..]).assign(&db.xsec[sidx[0][i]]);
        }

        SKXsecDatabase::<Ix2>::new(xsec, Grid::new(unique_wvnum), params).unwrap()
    }
}

impl XsecDatabaseInterp for SKXsecDatabase<Ix2> {
    fn xs_emplace<S1, S2, S3>(
        &self,
        wvnum: &ArrayBase<S1, Ix1>,
        params: &ArrayBase<S2, Ix1>,
        xs: &mut ArrayBase<S3, Ix1>,
        d_xs: Option<ArrayViewMut2<'_, f64>>,
    ) -> Result<(), String>
    where
        S1: Data<Elem = f64>,
        S2: Data<Elem = f64>,
        S3: DataMut<Elem = f64>,
    {
        let params_slice = &self.params[0];

        let weights = params_slice.interp1_weights(params[0], OutOfBoundsMode::Extend);

        for (i, weight, _) in weights.iter() {
            let internal_xs = self.xsec.slice(s![*i, ..]);

            Zip::indexed(wvnum).for_each(|j, wv| {
                let wvnum_weights = &self.wvnum.interp1_weights(*wv, OutOfBoundsMode::Zero);
                let local_xs = internal_xs[wvnum_weights[0].0] * wvnum_weights[0].1
                    + internal_xs[wvnum_weights[1].0] * wvnum_weights[1].1;

                xs[j] += local_xs * *weight;
            });
        }

        if let Some(mut d_xs_mut) = d_xs {
            for (i, _, d_weight) in weights.iter() {
                let internal_xs = self.xsec.slice(s![*i, ..]);
                Zip::indexed(wvnum).for_each(|j, wv| {
                    let wvnum_weights = &self.wvnum.interp1_weights(*wv, OutOfBoundsMode::Zero);
                    let local_xs = internal_xs[wvnum_weights[0].0] * wvnum_weights[0].1
                        + internal_xs[wvnum_weights[1].0] * wvnum_weights[1].1;

                    d_xs_mut[[0, j]] += local_xs * *d_weight;
                });
            }
        }

        // Have to iterate through all of the
        Ok(())
    }
}

impl XsecDatabaseInterp for SKXsecDatabase<Ix3> {
    fn xs_emplace<S1, S2, S3>(
        &self,
        wvnum: &ArrayBase<S1, Ix1>,
        params: &ArrayBase<S2, Ix1>,
        xs: &mut ArrayBase<S3, Ix1>,
        d_xs: Option<ArrayViewMut2<'_, f64>>,
    ) -> Result<(), String>
    where
        S1: Data<Elem = f64>,
        S2: Data<Elem = f64>,
        S3: DataMut<Elem = f64>,
    {
        let weights_0 = &self.params[0].interp1_weights(params[0], OutOfBoundsMode::Extend);
        let weights_1 = &self.params[1].interp1_weights(params[1], OutOfBoundsMode::Extend);

        for (i, weight_0, _) in weights_0.iter() {
            for (j, weight_1, _) in weights_1.iter() {
                let internal_xs = self.xsec.slice(s![*i, *j, ..]);
                Zip::indexed(wvnum).for_each(|k, wv| {
                    let wvnum_weights = &self.wvnum.interp1_weights(*wv, OutOfBoundsMode::Zero);
                    let local_xs = internal_xs[wvnum_weights[0].0] * wvnum_weights[0].1
                        + internal_xs[wvnum_weights[1].0] * wvnum_weights[1].1;

                    xs[k] += local_xs * *weight_0 * *weight_1;
                });
            }
        }

        if let Some(mut d_xs_mut) = d_xs {
            for (i, weight_0, d_weight_0) in weights_0.iter() {
                for (j, weight_1, d_weight_1) in weights_1.iter() {
                    let internal_xs = self.xsec.slice(s![*i, *j, ..]);
                    Zip::indexed(wvnum).for_each(|k, wv| {
                        let wvnum_weights = &self.wvnum.interp1_weights(*wv, OutOfBoundsMode::Zero);
                        let xs_interp = internal_xs[wvnum_weights[0].0] * wvnum_weights[0].1
                            + internal_xs[wvnum_weights[1].0] * wvnum_weights[1].1;

                        d_xs_mut[[0, k]] += xs_interp * *d_weight_0 * *weight_1;
                        d_xs_mut[[1, k]] += xs_interp * *weight_0 * *d_weight_1;
                    });
                }
            }
        }

        // Have to iterate through all of the
        Ok(())
    }
}

impl OpticalProperty for SKXsecDatabase<Ix2> {
    fn optical_quantities(
        &self,
        inputs: &dyn StorageInputs,
    ) -> Result<OpticalQuantities<OwnedRepr<f64>>> {
        let wvnum = inputs
            .wavenumbers_cminv()
            .ok_or_else(|| anyhow::anyhow!("Atmosphere does not have wavenumber_cminv"))?;

        let param = inputs
            .temperature_k()
            .ok_or_else(|| anyhow::anyhow!("Atmosphere does not have temperature_k"))?;

        // Reshape param to dim (1, param.len())
        let param = param.to_shape((1, param.len()))?;

        let mut xs = Array2::zeros((param.len(), wvnum.len()));

        Zip::from(xs.rows_mut())
            .and(param.columns())
            .par_for_each(|mut xs_row, param| {
                let _ = self.xs_emplace(&wvnum, &param, &mut xs_row, None);
            });

        let ssa = Array2::zeros(xs.dim());

        Ok(OpticalQuantities {
            cross_section: xs,
            ssa: ssa,
            a1: None,
            a2: None,
            a3: None,
            a4: None,
            b1: None,
            b2: None,
        })
    }
}

impl OpticalProperty for SKXsecDatabase<Ix3> {
    fn optical_quantities(
        &self,
        inputs: &dyn StorageInputs,
    ) -> Result<OpticalQuantities<OwnedRepr<f64>>> {
        let wvnum = inputs
            .wavenumbers_cminv()
            .ok_or_else(|| anyhow::anyhow!("Atmosphere does not have wavenumber_cminv"))?;

        let temperature = inputs
            .temperature_k()
            .ok_or_else(|| anyhow::anyhow!("Atmosphere does not have temperature_k"))?;

        let pressure = inputs
            .pressure_pa()
            .ok_or_else(|| anyhow::anyhow!("Atmosphere does not have pressure_pa"))?;

        let params = stack![Axis(0), pressure, temperature];

        let mut xs = Array2::zeros((params.shape()[1], wvnum.len()));

        Zip::from(xs.rows_mut())
            .and(params.columns())
            .par_for_each(|mut xs_row, param| {
                let _ = self.xs_emplace(&wvnum, &param, &mut xs_row, None);
            });

        let ssa = Array2::zeros(xs.dim());

        Ok(OpticalQuantities {
            cross_section: xs,
            ssa: ssa,
            a1: None,
            a2: None,
            a3: None,
            a4: None,
            b1: None,
            b2: None,
        })
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_construct_skxsec_database_ix2() {
        // Create a simple wavenumber grid
        let wavenumbers = Array1::linspace(300.0, 310.0, 11);
        let wvnum_grid = Grid::new(wavenumbers.clone());

        // Create a single parameter (temperature) with 3 values
        let temperatures = Array1::from_vec(vec![250.0, 273.0, 296.0]);
        let params = vec![temperatures.clone()];

        // Create cross-section data with shape (3, 11)
        // Each row represents cross-sections at a specific temperature
        // Columns are different wavelengths
        let mut xsec_data = Array2::zeros((3, 11));

        // Fill with some test data - let's make it depend on both temperature and wavelength
        // This makes a simple pattern where each temperature has a different baseline and slope
        for i in 0..3 {
            let temp = temperatures[i];
            for j in 0..11 {
                let wvnum = wavenumbers[j];
                // Create pattern: higher temps = higher baseline, and steeper slope with wavelength
                xsec_data[[i, j]] =
                    (temp - 200.0) * 1e-24 + (wvnum - 300.0) * 5e-26 * (i as f64 + 1.0);
            }
        }

        // Construct the database
        let db = SKXsecDatabase::new(xsec_data, wvnum_grid, params).unwrap();

        // Test interpolation at a few points
        let test_wvnum = Array1::from_vec(vec![302.5, 307.8]);
        let test_params = Array1::from_vec(vec![260.0]); // Between 250K and 273K

        let mut xs_result = Array1::zeros(test_wvnum.len());
        db.xs_emplace(&test_wvnum, &test_params, &mut xs_result, None)
            .unwrap();

        // Test interpolation with derivatives
        let mut xs_with_deriv = Array1::zeros(test_wvnum.len());
        let mut d_xs = Array2::zeros((1, test_wvnum.len()));
        db.xs_emplace(
            &test_wvnum,
            &test_params,
            &mut xs_with_deriv,
            Some(d_xs.view_mut()),
        )
        .unwrap();

        // Values should match regardless of whether we calculate derivatives
        assert!((xs_result[0] - xs_with_deriv[0]).abs() < 1e-10);
        assert!((xs_result[1] - xs_with_deriv[1]).abs() < 1e-10);
    }
}
