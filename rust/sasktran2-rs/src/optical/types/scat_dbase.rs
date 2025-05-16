use super::util::*;
use crate::atmosphere::traits::*;
use crate::interpolation::{OutOfBoundsMode, grid1d::*, linear::*};
use crate::optical::storage::*;
use crate::optical::traits::*;
use crate::prelude::*;
use ndarray::*;

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

            let mut legendre_result = legendre.index_axis_mut(Axis(0), j);

            for l in 0..leg_order {
                let a1_index_db = l * 6;
                let a1_index_result = l * num_legendre;

                let local_a1 = self.legendre[[wvnum_weights[0].0, a1_index_db]]
                    * wvnum_weights[0].1
                    + self.legendre[[wvnum_weights[1].0, a1_index_db]] * wvnum_weights[1].1;

                legendre_result[[a1_index_result]] += local_a1;
            }
        });

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
        todo!()
    }
}
