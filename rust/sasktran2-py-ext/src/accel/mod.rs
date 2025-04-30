use crate::constituent::deriv_mapping::PyDerivMapping;
use ndarray::Zip;
use numpy::*;
use pyo3::prelude::*;
use sasktran2_rs::atmosphere::DerivMapping;

pub mod broadening;

#[pyfunction]
pub fn assign_absorber_derivatives<'py>(
    deriv_mapping: Bound<'py, PyAny>,
    d_extinction: PyReadonlyArray2<'py, f64>,
    d_ssa: PyReadonlyArray2<'py, f64>,
    atmo_ssa: PyReadonlyArray2<'py, f64>,
    atmo_extinction: PyReadonlyArray2<'py, f64>,
) -> PyResult<()> {
    let mut py_mapping = PyDerivMapping::new(deriv_mapping);
    let mut mapping = py_mapping.mut_view();

    Zip::from(mapping.d_extinction.columns_mut())
        .and(mapping.d_ssa.columns_mut())
        .and(d_extinction.as_array().columns())
        .and(d_ssa.as_array().columns())
        .and(atmo_ssa.as_array().columns())
        .and(atmo_extinction.as_array().columns())
        .par_for_each(
            |mut mapping_d_extinction_row,
             mut mapping_d_ssa_row,
             d_extinction_row,
             d_ssa_row,
             atmo_ssa_row,
             atmo_extinction_row| {
                mapping_d_extinction_row.assign(&d_extinction_row);
                mapping_d_ssa_row.assign(&d_ssa_row);

                Zip::from(mapping_d_ssa_row)
                    .and(atmo_ssa_row)
                    .and(d_extinction_row)
                    .and(atmo_extinction_row)
                    .for_each(
                        |mapping_d_ssa, atmo_ssa_val, d_extinction_val, atmo_extinction_val| {
                            *mapping_d_ssa = (*mapping_d_ssa - atmo_ssa_val) * d_extinction_val
                                / atmo_extinction_val;
                        },
                    );
            },
        );

    Ok(())
}
