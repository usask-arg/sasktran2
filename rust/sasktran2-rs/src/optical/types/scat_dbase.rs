use crate::atmosphere::traits::*;
use crate::interpolation::{OutOfBoundsMode, grid1d::*, linear::*};
use crate::optical::storage::*;
use crate::optical::traits::*;
use crate::prelude::*;
use ndarray::*;

#[inline(always)]
fn assign_legendre(
    wvnum_weights: &[(usize, f64, f64); 2],
    leg_order: usize,
    num_legendre: usize,
    mut legendre_result: ArrayViewMut1<f64>,
    leg_db: ArrayView2<f64>,
    num_stokes: usize,
    weight: f64,
) {
    for l in 0..leg_order {
        let a1_index_db = l * 6;
        let a1_index_result = l * num_legendre;

        let local_a1 = leg_db[[wvnum_weights[0].0, a1_index_db]] * wvnum_weights[0].1
            + leg_db[[wvnum_weights[1].0, a1_index_db]] * wvnum_weights[1].1;

        legendre_result[[a1_index_result]] += local_a1 * weight;

        if num_stokes == 3 {
            let a2_index_db = a1_index_db + 1;
            let a2_index_result = a1_index_result + 1;

            let local_a2 = leg_db[[wvnum_weights[0].0, a2_index_db]] * wvnum_weights[0].1
                + leg_db[[wvnum_weights[1].0, a2_index_db]] * wvnum_weights[1].1;

            legendre_result[[a2_index_result]] += local_a2 * weight;

            let a3_index_db = a1_index_db + 2;
            let a3_index_result = a1_index_result + 2;

            let local_a3 = leg_db[[wvnum_weights[0].0, a3_index_db]] * wvnum_weights[0].1
                + leg_db[[wvnum_weights[1].0, a3_index_db]] * wvnum_weights[1].1;

            legendre_result[[a3_index_result]] += local_a3 * weight;

            // result is [a1, a2, a3, b1]
            // db is [a1, a2, a3, a4, b1, b2]
            let b1_index_db = a1_index_db + 4;
            let b1_index_result = a1_index_result + 3;

            let local_b1 = leg_db[[wvnum_weights[0].0, b1_index_db]] * wvnum_weights[0].1
                + leg_db[[wvnum_weights[1].0, b1_index_db]] * wvnum_weights[1].1;

            legendre_result[[b1_index_result]] += local_b1 * weight;
        }
    }
}

/// Trait for interpolating Scattering properties.
pub trait ScatteringDatabaseInterp {
    fn scat_prop_emplace<S1, S2, S3, S4, S5>(
        &self,
        wvnum: &ArrayBase<S1, Ix1>,
        params: &ArrayBase<S2, Ix1>,
        xs: &mut ArrayBase<S3, Ix1>,
        ssa: &mut ArrayBase<S4, Ix1>,
        legendre: &mut ArrayBase<S5, Ix2>,
        num_stokes: usize,
    ) -> Result<()>
    where
        S1: Data<Elem = f64>,
        S2: Data<Elem = f64>,
        S3: DataMut<Elem = f64>,
        S4: DataMut<Elem = f64>,
        S5: DataMut<Elem = f64>;

    fn d_scat_prop_emplace<S1, S2, S3, S4, S5>(
        &self,
        wvnum: &ArrayBase<S1, Ix1>,
        params: &ArrayBase<S2, Ix1>,
        d_xs: &mut Vec<ArrayBase<S3, Ix1>>,
        d_ssa: &mut Vec<ArrayBase<S4, Ix1>>,
        d_leg: &mut Vec<ArrayBase<S5, Ix2>>,
        num_stokes: usize,
    ) -> Result<()>
    where
        S1: Data<Elem = f64>,
        S2: Data<Elem = f64>,
        S3: DataMut<Elem = f64>,
        S4: DataMut<Elem = f64>,
        S5: DataMut<Elem = f64>;
}

pub struct ScatteringDatabase<D1: Dimension, D2: Dimension> {
    xsec: Array<f64, D1>,
    ssa: Array<f64, D1>,
    legendre: Array<f64, D2>,
    wvnum: Grid1D,
    params: Vec<Array1<f64>>,
    param_names: Vec<String>,
}

impl<D1: Dimension, D2: Dimension> ScatteringDatabase<D1, D2> {
    pub fn new(
        xsec: Array<f64, D1>,
        ssa: Array<f64, D1>,
        legendre: Array<f64, D2>,
        wvnum: Grid1D,
        params: Vec<Array1<f64>>,
        param_names: Vec<String>,
    ) -> Self {
        Self {
            xsec,
            ssa,
            legendre,
            wvnum,
            params,
            param_names,
        }
    }
}

impl ScatteringDatabaseInterp for ScatteringDatabase<Ix1, Ix2> {
    fn scat_prop_emplace<S1, S2, S3, S4, S5>(
        &self,
        wvnum: &ArrayBase<S1, Ix1>,
        _params: &ArrayBase<S2, Ix1>,
        xs: &mut ArrayBase<S3, Ix1>,
        ssa: &mut ArrayBase<S4, Ix1>,
        legendre: &mut ArrayBase<S5, Ix2>,
        num_stokes: usize,
    ) -> Result<()>
    where
        S1: Data<Elem = f64>,
        S2: Data<Elem = f64>,
        S3: DataMut<Elem = f64>,
        S4: DataMut<Elem = f64>,
        S5: DataMut<Elem = f64>,
    {
        let num_legendre = match num_stokes {
            1 => 1,
            3 => 4,
            4 => 6,
            _ => panic!("Invalid number of Stokes parameters"),
        };

        let leg_order = (self.legendre.dim().1 / 6).min(legendre.dim().1 / num_legendre);

        Zip::indexed(wvnum).for_each(|j, wv| {
            let wvnum_weights = &self.wvnum.interp1_weights(*wv, OutOfBoundsMode::Zero);

            let local_xs = self.xsec[wvnum_weights[0].0] * wvnum_weights[0].1
                + self.xsec[wvnum_weights[1].0] * wvnum_weights[1].1;

            let local_ssa = self.ssa[wvnum_weights[0].0] * wvnum_weights[0].1
                + self.ssa[wvnum_weights[1].0] * wvnum_weights[1].1;

            xs[j] += local_xs;
            ssa[j] += local_ssa;

            let legendre_result = legendre.index_axis_mut(Axis(0), j);
            assign_legendre(
                wvnum_weights,
                leg_order,
                num_legendre,
                legendre_result,
                self.legendre.view(),
                num_stokes,
                1.0,
            );
        });

        Ok(())
    }

    fn d_scat_prop_emplace<S1, S2, S3, S4, S5>(
        &self,
        _wvnum: &ArrayBase<S1, Ix1>,
        _params: &ArrayBase<S2, Ix1>,
        _d_xs: &mut Vec<ArrayBase<S3, Ix1>>,
        _d_ssa: &mut Vec<ArrayBase<S4, Ix1>>,
        _d_leg: &mut Vec<ArrayBase<S5, Ix2>>,
        _num_stokes: usize,
    ) -> Result<()>
    where
        S1: Data<Elem = f64>,
        S2: Data<Elem = f64>,
        S3: DataMut<Elem = f64>,
        S4: DataMut<Elem = f64>,
        S5: DataMut<Elem = f64>,
    {
        // No derivatives
        Ok(())
    }
}

impl OpticalProperty for ScatteringDatabase<Ix1, Ix2> {
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

        let num_stokes = inputs.num_stokes();
        let _ = optical_quantities.with_scatterer(inputs.num_singlescatter_moments(), num_stokes);

        let xs = &mut optical_quantities.cross_section;
        let ssa = &mut optical_quantities.ssa;
        let legendre = optical_quantities
            .legendre
            .as_mut()
            .ok_or_else(|| anyhow::anyhow!("Legendre coefficients not initialized"))?;

        Zip::from(xs.rows_mut())
            .and(ssa.rows_mut())
            .and(legendre.axis_iter_mut(Axis(0)))
            .and(altitudes_m.view())
            .par_for_each(|mut xs, mut ssa, mut legendre, param| {
                // Pass in altitude for no reason, but once again it's not used
                let params = Array1::from(vec![*param]);

                let _ = self.scat_prop_emplace(
                    &wavenumber_cminv,
                    &params,
                    &mut xs,
                    &mut ssa,
                    &mut legendre,
                    num_stokes,
                );
            });

        Ok(())
    }

    fn optical_derivatives_emplace(
        &self,
        _inputs: &dyn StorageInputs,
        _aux_inputs: &dyn AuxOpticalInputs,
        _d_optical_quantities: &mut std::collections::HashMap<String, OpticalQuantities>,
    ) -> anyhow::Result<()> {
        // No derivatives in this case

        Ok(())
    }
}

impl ScatteringDatabaseInterp for ScatteringDatabase<Ix2, Ix3> {
    fn scat_prop_emplace<S1, S2, S3, S4, S5>(
        &self,
        wvnum: &ArrayBase<S1, Ix1>,
        params: &ArrayBase<S2, Ix1>,
        xs: &mut ArrayBase<S3, Ix1>,
        ssa: &mut ArrayBase<S4, Ix1>,
        legendre: &mut ArrayBase<S5, Ix2>,
        num_stokes: usize,
    ) -> Result<()>
    where
        S1: Data<Elem = f64>,
        S2: Data<Elem = f64>,
        S3: DataMut<Elem = f64>,
        S4: DataMut<Elem = f64>,
        S5: DataMut<Elem = f64>,
    {
        let num_legendre = match num_stokes {
            1 => 1,
            3 => 4,
            4 => 6,
            _ => panic!("Invalid number of Stokes parameters"),
        };
        let leg_order = (self.legendre.dim().2 / 6).min(legendre.dim().1 / num_legendre);

        let weights_0 = &self.params[0].interp1_weights(params[0], OutOfBoundsMode::Extend);

        for (i0, weight0, _) in weights_0.iter() {
            let i0 = *i0;
            Zip::indexed(wvnum).for_each(|j, wv| {
                let wvnum_weights = &self.wvnum.interp1_weights(*wv, OutOfBoundsMode::Zero);

                let local_xs = self.xsec[[i0, wvnum_weights[0].0]] * wvnum_weights[0].1
                    + self.xsec[[i0, wvnum_weights[1].0]] * wvnum_weights[1].1;

                let local_ssa = self.ssa[[i0, wvnum_weights[0].0]] * wvnum_weights[0].1
                    + self.ssa[[i0, wvnum_weights[1].0]] * wvnum_weights[1].1;

                xs[j] += local_xs * (*weight0);
                ssa[j] += local_ssa * (*weight0);

                let legendre_result = legendre.index_axis_mut(Axis(0), j);

                let leg_db = self.legendre.index_axis(Axis(0), i0);

                assign_legendre(
                    wvnum_weights,
                    leg_order,
                    num_legendre,
                    legendre_result,
                    leg_db,
                    num_stokes,
                    *weight0,
                );
            });
        }

        Ok(())
    }

    fn d_scat_prop_emplace<S1, S2, S3, S4, S5>(
        &self,
        wvnum: &ArrayBase<S1, Ix1>,
        params: &ArrayBase<S2, Ix1>,
        d_xs: &mut Vec<ArrayBase<S3, Ix1>>,
        d_ssa: &mut Vec<ArrayBase<S4, Ix1>>,
        d_leg: &mut Vec<ArrayBase<S5, Ix2>>,
        num_stokes: usize,
    ) -> Result<()>
    where
        S1: Data<Elem = f64>,
        S2: Data<Elem = f64>,
        S3: DataMut<Elem = f64>,
        S4: DataMut<Elem = f64>,
        S5: DataMut<Elem = f64>,
    {
        // When we only have one derivative it's fairly straightforward
        let xs = &mut d_xs[0];
        let ssa = &mut d_ssa[0];
        let legendre = &mut d_leg[0];

        let num_legendre = match num_stokes {
            1 => 1,
            3 => 4,
            4 => 6,
            _ => panic!("Invalid number of Stokes parameters"),
        };
        let leg_order = (self.legendre.dim().2 / 6).min(legendre.dim().1 / num_legendre);

        let weights_0 = &self.params[0].interp1_weights(params[0], OutOfBoundsMode::Extend);

        for (i0, _, d_weight0) in weights_0.iter() {
            let i0 = *i0;
            Zip::indexed(wvnum).for_each(|j, wv| {
                let wvnum_weights = &self.wvnum.interp1_weights(*wv, OutOfBoundsMode::Zero);

                let local_xs = self.xsec[[i0, wvnum_weights[0].0]] * wvnum_weights[0].1
                    + self.xsec[[i0, wvnum_weights[1].0]] * wvnum_weights[1].1;

                let local_ssa = self.ssa[[i0, wvnum_weights[0].0]] * wvnum_weights[0].1
                    + self.ssa[[i0, wvnum_weights[1].0]] * wvnum_weights[1].1;

                xs[j] += local_xs * (*d_weight0);
                ssa[j] += local_ssa * (*d_weight0);

                let legendre_result = legendre.index_axis_mut(Axis(0), j);

                let leg_db = self.legendre.index_axis(Axis(0), i0);

                assign_legendre(
                    wvnum_weights,
                    leg_order,
                    num_legendre,
                    legendre_result,
                    leg_db,
                    num_stokes,
                    *d_weight0,
                );
            });
        }

        Ok(())
    }
}

impl OpticalProperty for ScatteringDatabase<Ix2, Ix3> {
    fn optical_quantities_emplace(
        &self,
        inputs: &dyn StorageInputs,
        aux_inputs: &dyn AuxOpticalInputs,
        optical_quantities: &mut OpticalQuantities,
    ) -> Result<()> {
        let wavenumber_cminv = param_from_storage_or_aux(inputs, aux_inputs, "wavenumbers_cminv")?;
        let param_0 = param_from_storage_or_aux(inputs, aux_inputs, &self.param_names[0])?;

        let _ = optical_quantities.resize(param_0.len(), wavenumber_cminv.len());
        let num_stokes = inputs.num_stokes();
        let _ = optical_quantities.with_scatterer(inputs.num_singlescatter_moments(), num_stokes);

        let xs = &mut optical_quantities.cross_section;
        let ssa = &mut optical_quantities.ssa;
        let legendre = optical_quantities
            .legendre
            .as_mut()
            .ok_or_else(|| anyhow::anyhow!("Legendre coefficients not initialized"))?;

        Zip::from(xs.rows_mut())
            .and(ssa.rows_mut())
            .and(legendre.axis_iter_mut(Axis(0)))
            .and(param_0.view())
            .par_for_each(|mut xs, mut ssa, mut legendre, param_0| {
                let params = Array1::from(vec![*param_0]);

                let _ = self.scat_prop_emplace(
                    &wavenumber_cminv,
                    &params,
                    &mut xs,
                    &mut ssa,
                    &mut legendre,
                    num_stokes,
                );
            });

        Ok(())
    }

    fn optical_derivatives_emplace(
        &self,
        inputs: &dyn StorageInputs,
        aux_inputs: &dyn AuxOpticalInputs,
        d_optical_quantities: &mut HashMap<String, OpticalQuantities>,
    ) -> Result<()> {
        let wavenumber_cminv = param_from_storage_or_aux(inputs, aux_inputs, "wavenumbers_cminv")?;
        let param_0 = param_from_storage_or_aux(inputs, aux_inputs, &self.param_names[0])?;

        if d_optical_quantities.contains_key(&self.param_names[0]) {
            d_optical_quantities
                .get_mut(&self.param_names[0])
                .unwrap()
                .resize(param_0.len(), wavenumber_cminv.len());
        } else {
            d_optical_quantities.insert(
                self.param_names[0].clone(),
                OpticalQuantities::new(param_0.len(), wavenumber_cminv.len(), false),
            );
        }

        let optical_quantities = d_optical_quantities.get_mut(&self.param_names[0]).unwrap();

        let _ = optical_quantities.resize(param_0.len(), wavenumber_cminv.len());
        let num_stokes = inputs.num_stokes();
        let _ = optical_quantities.with_scatterer(inputs.num_singlescatter_moments(), num_stokes);

        let xs = &mut optical_quantities.cross_section;
        let ssa = &mut optical_quantities.ssa;
        let legendre = optical_quantities
            .legendre
            .as_mut()
            .ok_or_else(|| anyhow::anyhow!("Legendre coefficients not initialized"))?;

        Zip::from(xs.rows_mut())
            .and(ssa.rows_mut())
            .and(legendre.axis_iter_mut(Axis(0)))
            .and(param_0.view())
            .par_for_each(|xs, ssa, legendre, param_0| {
                let params = Array1::from(vec![*param_0]);

                let _ = self.d_scat_prop_emplace(
                    &wavenumber_cminv,
                    &params,
                    &mut vec![xs],
                    &mut vec![ssa],
                    &mut vec![legendre],
                    num_stokes,
                );
            });

        Ok(())
    }
}

impl ScatteringDatabaseInterp for ScatteringDatabase<Ix3, Ix4> {
    fn scat_prop_emplace<S1, S2, S3, S4, S5>(
        &self,
        wvnum: &ArrayBase<S1, Ix1>,
        params: &ArrayBase<S2, Ix1>,
        xs: &mut ArrayBase<S3, Ix1>,
        ssa: &mut ArrayBase<S4, Ix1>,
        legendre: &mut ArrayBase<S5, Ix2>,
        num_stokes: usize,
    ) -> Result<()>
    where
        S1: Data<Elem = f64>,
        S2: Data<Elem = f64>,
        S3: DataMut<Elem = f64>,
        S4: DataMut<Elem = f64>,
        S5: DataMut<Elem = f64>,
    {
        let num_legendre = match num_stokes {
            1 => 1,
            3 => 4,
            4 => 6,
            _ => panic!("Invalid number of Stokes parameters"),
        };
        let leg_order = (self.legendre.dim().3 / 6).min(legendre.dim().1 / num_legendre);

        let weights_0 = &self.params[0].interp1_weights(params[0], OutOfBoundsMode::Extend);
        let weights_1 = &self.params[1].interp1_weights(params[1], OutOfBoundsMode::Extend);

        for (i0, weight0, _) in weights_0.iter() {
            let i0 = *i0;

            for (i1, weight1, _) in weights_1.iter() {
                let i1 = *i1;
                Zip::indexed(wvnum).for_each(|j, wv| {
                    let wvnum_weights = &self.wvnum.interp1_weights(*wv, OutOfBoundsMode::Zero);

                    let local_xs = self.xsec[[i0, i1, wvnum_weights[0].0]] * wvnum_weights[0].1
                        + self.xsec[[i0, i1, wvnum_weights[1].0]] * wvnum_weights[1].1;

                    let local_ssa = self.ssa[[i0, i1, wvnum_weights[0].0]] * wvnum_weights[0].1
                        + self.ssa[[i0, i1, wvnum_weights[1].0]] * wvnum_weights[1].1;

                    xs[j] += local_xs * (*weight0 * *weight1);
                    ssa[j] += local_ssa * (*weight0 * *weight1);

                    let legendre_result = legendre.index_axis_mut(Axis(0), j);

                    let leg_db = self.legendre.index_axis(Axis(0), i0);
                    let leg_db = leg_db.index_axis(Axis(0), i1);

                    assign_legendre(
                        wvnum_weights,
                        leg_order,
                        num_legendre,
                        legendre_result,
                        leg_db,
                        num_stokes,
                        *weight0 * *weight1,
                    );
                });
            }
        }

        Ok(())
    }

    fn d_scat_prop_emplace<S1, S2, S3, S4, S5>(
        &self,
        wvnum: &ArrayBase<S1, Ix1>,
        params: &ArrayBase<S2, Ix1>,
        d_xs: &mut Vec<ArrayBase<S3, Ix1>>,
        d_ssa: &mut Vec<ArrayBase<S4, Ix1>>,
        d_leg: &mut Vec<ArrayBase<S5, Ix2>>,
        num_stokes: usize,
    ) -> Result<()>
    where
        S1: Data<Elem = f64>,
        S2: Data<Elem = f64>,
        S3: DataMut<Elem = f64>,
        S4: DataMut<Elem = f64>,
        S5: DataMut<Elem = f64>,
    {
        let num_legendre = match num_stokes {
            1 => 1,
            3 => 4,
            4 => 6,
            _ => panic!("Invalid number of Stokes parameters"),
        };

        let weights_0 = &self.params[0].interp1_weights(params[0], OutOfBoundsMode::Extend);
        let weights_1 = &self.params[1].interp1_weights(params[1], OutOfBoundsMode::Extend);

        // Do two passes, one for the first param
        let xs = &mut d_xs[0];
        let ssa = &mut d_ssa[0];
        let legendre = &mut d_leg[0];
        let leg_order = (self.legendre.dim().3 / 6).min(legendre.dim().1 / num_legendre);

        for (i0, _, weight0) in weights_0.iter() {
            let i0 = *i0;

            for (i1, weight1, _) in weights_1.iter() {
                let i1 = *i1;
                Zip::indexed(wvnum).for_each(|j, wv| {
                    let wvnum_weights = &self.wvnum.interp1_weights(*wv, OutOfBoundsMode::Zero);

                    let local_xs = self.xsec[[i0, i1, wvnum_weights[0].0]] * wvnum_weights[0].1
                        + self.xsec[[i0, i1, wvnum_weights[1].0]] * wvnum_weights[1].1;

                    let local_ssa = self.ssa[[i0, i1, wvnum_weights[0].0]] * wvnum_weights[0].1
                        + self.ssa[[i0, i1, wvnum_weights[1].0]] * wvnum_weights[1].1;

                    xs[j] += local_xs * (*weight0 * *weight1);
                    ssa[j] += local_ssa * (*weight0 * *weight1);

                    let legendre_result = legendre.index_axis_mut(Axis(0), j);

                    let leg_db = self.legendre.index_axis(Axis(0), i0);
                    let leg_db = leg_db.index_axis(Axis(0), i1);

                    assign_legendre(
                        wvnum_weights,
                        leg_order,
                        num_legendre,
                        legendre_result,
                        leg_db,
                        num_stokes,
                        *weight0 * *weight1,
                    );
                });
            }
        }

        // And again for the second parameter
        let xs = &mut d_xs[1];
        let ssa = &mut d_ssa[1];
        let legendre = &mut d_leg[1];

        for (i0, weight0, _) in weights_0.iter() {
            let i0 = *i0;

            for (i1, _, weight1) in weights_1.iter() {
                let i1 = *i1;
                Zip::indexed(wvnum).for_each(|j, wv| {
                    let wvnum_weights = &self.wvnum.interp1_weights(*wv, OutOfBoundsMode::Zero);

                    let local_xs = self.xsec[[i0, i1, wvnum_weights[0].0]] * wvnum_weights[0].1
                        + self.xsec[[i0, i1, wvnum_weights[1].0]] * wvnum_weights[1].1;

                    let local_ssa = self.ssa[[i0, i1, wvnum_weights[0].0]] * wvnum_weights[0].1
                        + self.ssa[[i0, i1, wvnum_weights[1].0]] * wvnum_weights[1].1;

                    xs[j] += local_xs;
                    ssa[j] += local_ssa;

                    let legendre_result = legendre.index_axis_mut(Axis(0), j);

                    let leg_db = self.legendre.index_axis(Axis(0), i0);
                    let leg_db = leg_db.index_axis(Axis(0), i1);

                    assign_legendre(
                        wvnum_weights,
                        leg_order,
                        num_legendre,
                        legendre_result,
                        leg_db,
                        num_stokes,
                        *weight0 * *weight1,
                    );
                });
            }
        }

        Ok(())
    }
}

impl OpticalProperty for ScatteringDatabase<Ix3, Ix4> {
    fn optical_quantities_emplace(
        &self,
        inputs: &dyn StorageInputs,
        aux_inputs: &dyn AuxOpticalInputs,
        optical_quantities: &mut OpticalQuantities,
    ) -> Result<()> {
        let wavenumber_cminv = param_from_storage_or_aux(inputs, aux_inputs, "wavenumbers_cminv")?;
        let param_0 = param_from_storage_or_aux(inputs, aux_inputs, &self.param_names[0])?;
        let param_1 = param_from_storage_or_aux(inputs, aux_inputs, &self.param_names[1])?;

        let _ = optical_quantities.resize(param_0.len(), wavenumber_cminv.len());
        let num_stokes = inputs.num_stokes();
        let _ = optical_quantities.with_scatterer(inputs.num_singlescatter_moments(), num_stokes);

        let xs = &mut optical_quantities.cross_section;
        let ssa = &mut optical_quantities.ssa;
        let legendre = optical_quantities
            .legendre
            .as_mut()
            .ok_or_else(|| anyhow::anyhow!("Legendre coefficients not initialized"))?;

        Zip::from(xs.rows_mut())
            .and(ssa.rows_mut())
            .and(legendre.axis_iter_mut(Axis(0)))
            .and(param_0.view())
            .and(param_1.view())
            .par_for_each(|mut xs, mut ssa, mut legendre, param_0, param_1| {
                let params = Array1::from(vec![*param_0, *param_1]);

                let _ = self.scat_prop_emplace(
                    &wavenumber_cminv,
                    &params,
                    &mut xs,
                    &mut ssa,
                    &mut legendre,
                    num_stokes,
                );
            });

        Ok(())
    }

    fn optical_derivatives_emplace(
        &self,
        inputs: &dyn StorageInputs,
        aux_inputs: &dyn AuxOpticalInputs,
        d_optical_quantities: &mut HashMap<String, OpticalQuantities>,
    ) -> Result<()> {
        let wavenumber_cminv = param_from_storage_or_aux(inputs, aux_inputs, "wavenumbers_cminv")?;
        let param_0 = param_from_storage_or_aux(inputs, aux_inputs, &self.param_names[0])?;
        let param_1 = param_from_storage_or_aux(inputs, aux_inputs, &self.param_names[1])?;

        if d_optical_quantities.contains_key(&self.param_names[0]) {
            d_optical_quantities
                .get_mut(&self.param_names[0])
                .unwrap()
                .resize(param_0.len(), wavenumber_cminv.len());
        } else {
            d_optical_quantities.insert(
                self.param_names[0].clone(),
                OpticalQuantities::new(param_0.len(), wavenumber_cminv.len(), false),
            );
        }

        if d_optical_quantities.contains_key(&self.param_names[1]) {
            d_optical_quantities
                .get_mut(&self.param_names[1])
                .unwrap()
                .resize(param_0.len(), wavenumber_cminv.len());
        } else {
            d_optical_quantities.insert(
                self.param_names[1].clone(),
                OpticalQuantities::new(param_0.len(), wavenumber_cminv.len(), false),
            );
        }

        let optical_quantities = d_optical_quantities
            .get_disjoint_mut([self.param_names[0].as_str(), self.param_names[1].as_str()]);

        let [optical_quantities0, optical_quantities1] = optical_quantities;

        let optical_quantities0 = optical_quantities0.unwrap();
        let optical_quantities1 = optical_quantities1.unwrap();

        let _ = optical_quantities0.resize(param_0.len(), wavenumber_cminv.len());
        let num_stokes = inputs.num_stokes();
        let _ = optical_quantities0.with_scatterer(inputs.num_singlescatter_moments(), num_stokes);

        let xs0 = &mut optical_quantities0.cross_section;
        let ssa0 = &mut optical_quantities0.ssa;
        let legendre0 = optical_quantities0
            .legendre
            .as_mut()
            .ok_or_else(|| anyhow::anyhow!("Legendre coefficients not initialized"))?;

        let _ = optical_quantities1.resize(param_0.len(), wavenumber_cminv.len());
        let _ = optical_quantities1.with_scatterer(inputs.num_singlescatter_moments(), num_stokes);

        let xs1 = &mut optical_quantities1.cross_section;
        let ssa1 = &mut optical_quantities1.ssa;
        let legendre1 = optical_quantities1
            .legendre
            .as_mut()
            .ok_or_else(|| anyhow::anyhow!("Legendre coefficients not initialized"))?;

        Zip::indexed(param_0.view())
            .and(param_1.view())
            .for_each(|i, param_0, param_1| {
                let xs0 = xs0.index_axis_mut(Axis(0), i);
                let ssa0 = ssa0.index_axis_mut(Axis(0), i);
                let legendre0 = legendre0.index_axis_mut(Axis(0), i);

                let xs1 = xs1.index_axis_mut(Axis(0), i);
                let ssa1 = ssa1.index_axis_mut(Axis(0), i);
                let legendre1 = legendre1.index_axis_mut(Axis(0), i);

                let params = Array1::from(vec![*param_0, *param_1]);

                let _ = self.d_scat_prop_emplace(
                    &wavenumber_cminv,
                    &params,
                    &mut vec![xs0, xs1],
                    &mut vec![ssa0, ssa1],
                    &mut vec![legendre0, legendre1],
                    num_stokes,
                );
            });

        Ok(())
    }
}
