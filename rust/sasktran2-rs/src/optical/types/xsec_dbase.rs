use crate::atmosphere::traits::*;
use crate::interpolation::OutOfBoundsMode;
use crate::interpolation::grid1d::Grid1D;
use crate::interpolation::linear::Interp1Weights;
use crate::optical::storage::*;
use crate::optical::traits::*;
use ndarray::{Data, DataMut};
use ndarray::{Zip, prelude::*};
use std::collections::HashMap;

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

    fn d_xs_emplace<S1, S2, S3>(
        &self,
        wvnum: &ArrayBase<S1, Ix1>,
        params: &ArrayBase<S2, Ix1>,
        d_xs: &mut Vec<ArrayBase<S3, Ix1>>,
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
    wvnum: Grid1D,
    params: Vec<Array1<f64>>,
    param_names: Vec<String>,
}

impl<D: Dimension> SKXsecDatabase<D> {
    pub fn new(
        xsec: Array<f64, D>,
        wvnum: Grid1D,
        params: Vec<Array1<f64>>,
        param_names: Vec<String>,
    ) -> Option<Self> {
        assert_eq!(params.len(), D::NDIM? - 1);

        Some(SKXsecDatabase {
            xsec,
            wvnum,
            params,
            param_names,
        })
    }
}

impl XsecDatabaseInterp for SKXsecDatabase<Ix1> {
    fn xs_emplace<S1, S2, S3>(
        &self,
        wvnum: &ArrayBase<S1, Ix1>,
        _params: &ArrayBase<S2, Ix1>,
        xs: &mut ArrayBase<S3, Ix1>,
        _d_xs: Option<ArrayViewMut2<'_, f64>>,
    ) -> Result<(), String>
    where
        S1: Data<Elem = f64>,
        S2: Data<Elem = f64>,
        S3: DataMut<Elem = f64>,
    {
        // Only have to interpolate in wavenumber
        Zip::indexed(wvnum).for_each(|j, wv| {
            let wvnum_weights = &self.wvnum.interp1_weights(*wv, OutOfBoundsMode::Zero);
            let local_xs = self.xsec[wvnum_weights[0].0] * wvnum_weights[0].1
                + self.xsec[wvnum_weights[1].0] * wvnum_weights[1].1;

            xs[j] += local_xs;
        });

        Ok(())
    }

    fn d_xs_emplace<S1, S2, S3>(
        &self,
        _wvnum: &ArrayBase<S1, Ix1>,
        _params: &ArrayBase<S2, Ix1>,
        _d_xs: &mut Vec<ArrayBase<S3, Ix1>>,
    ) -> Result<(), String>
    where
        S1: Data<Elem = f64>,
        S2: Data<Elem = f64>,
        S3: DataMut<Elem = f64>,
    {
        // No derivatives
        Ok(())
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

    fn d_xs_emplace<S1, S2, S3>(
        &self,
        wvnum: &ArrayBase<S1, Ix1>,
        params: &ArrayBase<S2, Ix1>,
        d_xs: &mut Vec<ArrayBase<S3, Ix1>>,
    ) -> Result<(), String>
    where
        S1: Data<Elem = f64>,
        S2: Data<Elem = f64>,
        S3: DataMut<Elem = f64>,
    {
        let params_slice = &self.params[0];

        let weights = params_slice.interp1_weights(params[0], OutOfBoundsMode::Extend);

        let mut d_xs_0 = d_xs[0].view_mut();

        for (i, _, d_weight) in weights.iter() {
            let internal_xs = self.xsec.slice(s![*i, ..]);
            Zip::indexed(wvnum).for_each(|j, wv| {
                let wvnum_weights = &self.wvnum.interp1_weights(*wv, OutOfBoundsMode::Zero);
                let local_xs = internal_xs[wvnum_weights[0].0] * wvnum_weights[0].1
                    + internal_xs[wvnum_weights[1].0] * wvnum_weights[1].1;

                d_xs_0[[j]] += local_xs * *d_weight;
            });
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

    fn d_xs_emplace<S1, S2, S3>(
        &self,
        wvnum: &ArrayBase<S1, Ix1>,
        params: &ArrayBase<S2, Ix1>,
        d_xs: &mut Vec<ArrayBase<S3, Ix1>>,
    ) -> Result<(), String>
    where
        S1: Data<Elem = f64>,
        S2: Data<Elem = f64>,
        S3: DataMut<Elem = f64>,
    {
        let weights_0 = &self.params[0].interp1_weights(params[0], OutOfBoundsMode::Extend);
        let weights_1 = &self.params[1].interp1_weights(params[1], OutOfBoundsMode::Extend);

        let (d_xs_0, d_xs_1) = d_xs.split_at_mut(1);
        let d_xs_0 = &mut d_xs_0[0].view_mut();
        let d_xs_1 = &mut d_xs_1[0].view_mut();

        for (i, weight_0, d_weight_0) in weights_0.iter() {
            for (j, weight_1, d_weight_1) in weights_1.iter() {
                let internal_xs = self.xsec.slice(s![*i, *j, ..]);
                Zip::indexed(wvnum).for_each(|k, wv| {
                    let wvnum_weights = &self.wvnum.interp1_weights(*wv, OutOfBoundsMode::Zero);
                    let xs_interp = internal_xs[wvnum_weights[0].0] * wvnum_weights[0].1
                        + internal_xs[wvnum_weights[1].0] * wvnum_weights[1].1;

                    d_xs_0[[k]] += xs_interp * *d_weight_0 * *weight_1;
                    d_xs_1[[k]] += xs_interp * *weight_0 * *d_weight_1;
                });
            }
        }

        // Have to iterate through all of the
        Ok(())
    }
}

impl OpticalProperty for SKXsecDatabase<Ix1> {
    fn optical_quantities_emplace(
        &self,
        inputs: &dyn StorageInputs,
        aux_inputs: &dyn AuxOpticalInputs,
        optical_quantities: &mut OpticalQuantities,
    ) -> anyhow::Result<()> {
        let wavenumber_cminv = param_from_storage_or_aux(inputs, aux_inputs, "wavenumbers_cminv")?;

        // Just grab this to get the number of geometry points, we don't actually interpolate in this dimension
        let altitudes_m = param_from_storage_or_aux(inputs, aux_inputs, "altitude_m")?;

        let _ = optical_quantities.resize(altitudes_m.len(), wavenumber_cminv.len());

        let xs = &mut optical_quantities.cross_section;

        Zip::from(xs.rows_mut())
            .and(altitudes_m.view())
            .par_for_each(|mut row, param| {
                // Pass in altitude for no reason, but once again it's not used
                let params = Array1::from(vec![*param]);

                let _ = self.xs_emplace(&wavenumber_cminv, &params, &mut row, None);
            });

        Ok(())
    }

    fn optical_derivatives_emplace(
        &self,
        _inputs: &dyn StorageInputs,
        _aux_inputs: &dyn AuxOpticalInputs,
        _d_optical_quantities: &mut HashMap<String, OpticalQuantities>,
    ) -> anyhow::Result<()> {
        Ok(())
    }

    fn is_scatterer(&self) -> bool {
        false
    }
}

impl OpticalProperty for SKXsecDatabase<Ix2> {
    fn optical_quantities_emplace(
        &self,
        inputs: &dyn StorageInputs,
        aux_inputs: &dyn AuxOpticalInputs,
        optical_quantities: &mut OpticalQuantities,
    ) -> anyhow::Result<()> {
        let wavenumber_cminv = param_from_storage_or_aux(inputs, aux_inputs, "wavenumbers_cminv")?;
        let param = param_from_storage_or_aux(inputs, aux_inputs, &self.param_names[0])?;

        let _ = optical_quantities.resize(param.len(), wavenumber_cminv.len());

        let xs = &mut optical_quantities.cross_section;

        Zip::from(xs.rows_mut())
            .and(param.view())
            .par_for_each(|mut row, param| {
                let params = Array1::from(vec![*param]);

                let _ = self.xs_emplace(&wavenumber_cminv, &params, &mut row, None);
            });

        Ok(())
    }

    fn optical_derivatives_emplace(
        &self,
        inputs: &dyn StorageInputs,
        aux_inputs: &dyn AuxOpticalInputs,
        d_optical_quantities: &mut HashMap<String, OpticalQuantities>,
    ) -> anyhow::Result<()> {
        let wavenumber_cminv = param_from_storage_or_aux(inputs, aux_inputs, "wavenumbers_cminv")?;
        let param = param_from_storage_or_aux(inputs, aux_inputs, &self.param_names[0])?;

        if d_optical_quantities.contains_key(&self.param_names[0]) {
            d_optical_quantities
                .get_mut(&self.param_names[0])
                .unwrap()
                .resize(param.len(), wavenumber_cminv.len());
        } else {
            d_optical_quantities.insert(
                self.param_names[0].clone(),
                OpticalQuantities::new(param.len(), wavenumber_cminv.len(), false),
            );
        }

        let d_xs: &mut ArrayBase<ndarray::OwnedRepr<f64>, Dim<[usize; 2]>> =
            &mut d_optical_quantities
                .get_mut(&self.param_names[0])
                .unwrap()
                .cross_section;

        Zip::from(param.view())
            .and(d_xs.rows_mut())
            .par_for_each(|param, mut row| {
                let params = Array1::from(vec![*param]);

                let mut d_xs_scratch = vec![row.view_mut()];

                let _ = self.d_xs_emplace(&wavenumber_cminv, &params, &mut d_xs_scratch);
            });

        Ok(())
    }

    fn is_scatterer(&self) -> bool {
        false
    }
}

impl OpticalProperty for SKXsecDatabase<Ix3> {
    fn optical_quantities_emplace(
        &self,
        inputs: &dyn StorageInputs,
        aux_inputs: &dyn AuxOpticalInputs,
        optical_quantities: &mut OpticalQuantities,
    ) -> anyhow::Result<()> {
        let wavenumber_cminv = param_from_storage_or_aux(inputs, aux_inputs, "wavenumbers_cminv")?;
        let param_0 = param_from_storage_or_aux(inputs, aux_inputs, &self.param_names[0])?;
        let param_1 = param_from_storage_or_aux(inputs, aux_inputs, &self.param_names[1])?;

        let _ = optical_quantities.resize(param_0.len(), wavenumber_cminv.len());

        let xs = &mut optical_quantities.cross_section;

        Zip::from(xs.rows_mut())
            .and(param_0.view())
            .and(param_1.view())
            .par_for_each(|mut row, param_0, param_1| {
                let params = Array1::from(vec![*param_0, *param_1]);

                let _ = self.xs_emplace(&wavenumber_cminv, &params, &mut row, None);
            });

        Ok(())
    }

    fn optical_derivatives_emplace(
        &self,
        inputs: &dyn StorageInputs,
        aux_inputs: &dyn AuxOpticalInputs,
        d_optical_quantities: &mut HashMap<String, OpticalQuantities>,
    ) -> anyhow::Result<()> {
        let wavenumber_cminv = param_from_storage_or_aux(inputs, aux_inputs, "wavenumbers_cminv")?;
        let param_0 = param_from_storage_or_aux(inputs, aux_inputs, &self.param_names[0])?;
        let param_1 = param_from_storage_or_aux(inputs, aux_inputs, &self.param_names[1])?;

        for i in 0..2 {
            if d_optical_quantities.contains_key(&self.param_names[i]) {
                d_optical_quantities
                    .get_mut(&self.param_names[i])
                    .unwrap()
                    .resize(param_0.len(), wavenumber_cminv.len());
            } else {
                d_optical_quantities.insert(
                    self.param_names[i].clone(),
                    OpticalQuantities::new(param_0.len(), wavenumber_cminv.len(), false),
                );
            }
        }

        let map_ptr = d_optical_quantities as *mut HashMap<String, OpticalQuantities>;

        // TODO: This is fixed in rust 1.86 with get_disjoint_mut
        unsafe {
            let d_xs_0 = &mut (*map_ptr)
                .get_mut(&self.param_names[0])
                .unwrap()
                .cross_section;
            let d_xs_1 = &mut (*map_ptr)
                .get_mut(&self.param_names[1])
                .unwrap()
                .cross_section;

            Zip::from(d_xs_0.rows_mut())
                .and(d_xs_1.rows_mut())
                .and(param_0.view())
                .and(param_1.view())
                .par_for_each(|mut row_0, mut row_1, param_0, param_1| {
                    let params = Array1::from(vec![*param_0, *param_1]);

                    let mut d_xs_scratch = vec![row_0.view_mut(), row_1.view_mut()];

                    let _ = self.d_xs_emplace(&wavenumber_cminv, &params, &mut d_xs_scratch);
                });
        }

        Ok(())
    }

    fn is_scatterer(&self) -> bool {
        false
    }
}
