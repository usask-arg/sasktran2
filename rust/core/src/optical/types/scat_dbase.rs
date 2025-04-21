use crate::atmosphere::traits::*;
use crate::interpolation::{OutOfBoundsMode, grid1d::*, linear::*};
use crate::optical::storage::*;
use crate::optical::traits::*;
use ndarray::*;

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

impl OpticalProperty for ScatteringDatabase<Ix2, Ix3> {
    fn optical_quantities_emplace(
        &self,
        inputs: &dyn StorageInputs,
        aux_inputs: &dyn AuxOpticalInputs,
        optical_quantities: &mut OpticalQuantities,
    ) -> anyhow::Result<()> {
        let wavenumber_cminv = param_from_storage_or_aux(inputs, aux_inputs, "wavenumbers_cminv")?;

        let params: Vec<ArrayBase<CowRepr<'_, f64>, Dim<[usize; 1]>>> = self
            .param_names
            .iter()
            .map(|name| param_from_storage_or_aux(inputs, aux_inputs, name))
            .collect::<anyhow::Result<Vec<_>>>()?;

        let num_atmo_legendre = 16;

        let _ = optical_quantities.resize(params[0].len(), wavenumber_cminv.len());
        let _ = optical_quantities.with_scatterer(num_atmo_legendre, inputs.num_stokes());

        let (xs, ssa, legendre) = optical_quantities.mut_split();

        let legendre = legendre.unwrap();

        Zip::from(xs.rows_mut())
            .and(ssa.rows_mut())
            .and(legendre.axis_iter_mut(Axis(1)))
            .and(params[0].view())
            .par_for_each(|mut xs_row, mut ssa_row, mut leg_row, param| {
                let weights = &self.params[0].interp1_weights(*param, OutOfBoundsMode::Extend);

                for (i, weight, _) in weights.iter() {
                    let internal_xs = self.xsec.index_axis(Axis(0), *i);
                    let internal_ssa = self.ssa.index_axis(Axis(0), *i);
                    let internal_legendre = self.legendre.index_axis(Axis(1), *i);

                    Zip::indexed(&wavenumber_cminv).for_each(|j, wv| {
                        let wvnum_weights = self.wvnum.interp1_weights(*wv, OutOfBoundsMode::Zero);

                        let local_xs = internal_xs[wvnum_weights[0].0] * wvnum_weights[0].1
                            + internal_xs[wvnum_weights[1].0] * wvnum_weights[1].1;
                        let local_ssa = internal_ssa[wvnum_weights[0].0] * wvnum_weights[0].1
                            + internal_ssa[wvnum_weights[1].0] * wvnum_weights[1].1;

                        xs_row[[j]] += local_xs * *weight;
                        ssa_row[[j]] += local_ssa * *weight;

                        let internal_legendre = internal_legendre.index_axis(Axis(1), j);

                        Zip::indexed(&internal_legendre).for_each(|l, leg| {
                            let local_leg = leg * wvnum_weights[0].1 + leg * wvnum_weights[1].1;

                            leg_row[[l, j]] += local_leg * *weight;
                        })
                    });
                }
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
