use std::collections::HashMap;

use anyhow::{Result, anyhow};
use ndarray::Array2;
use numpy::PyReadonlyArray2;
use pyo3::types::PyDictMethods;
use pyo3::{prelude::*, types::PyDict};
use sasktran2_rs::atmosphere::traits::StorageInputs;
use sasktran2_rs::optical::storage::*;
use sasktran2_rs::optical::traits::*;

use crate::optical::line_absorber::PyLineAbsorber;
use crate::optical::xsec_dbase::*;

fn with_optical_downcast(
    bound_optical_property: &Bound<'_, PyAny>,
    f: impl FnOnce(&dyn OpticalProperty) -> Result<()>,
) -> Result<()> {
    let rust_obj = bound_optical_property
        .call_method0("_into_rust_object")
        .map_err(|_| anyhow!("Failed to call _into_rust_object"))?;

    macro_rules! try_downcast {
        ($($ty:ty),*) => {
            $(
                if let Ok(obj) = rust_obj.downcast::<$ty>() {
                    return f(obj.borrow().db());
                }
            )*
        };
    }

    try_downcast!(
        AbsorberDatabaseDim1,
        AbsorberDatabaseDim2,
        AbsorberDatabaseDim3,
        PyLineAbsorber
    );

    Err(anyhow!("Failed to downcast to any of the known types"))
}

pub struct PyOpticalProperty {
    py_optical_property: Py<PyAny>,
    py_atmosphere: Py<PyAny>,
}

impl PyOpticalProperty {
    pub fn new(py_optical_property: Py<PyAny>, py_atmosphere: Py<PyAny>) -> Self {
        Self {
            py_optical_property,
            py_atmosphere,
        }
    }
}

impl OpticalProperty for PyOpticalProperty {
    fn optical_quantities_emplace(
        &self,
        inputs: &dyn StorageInputs,
        aux_inputs: &dyn AuxOpticalInputs,
        optical_quantities: &mut OpticalQuantities,
    ) -> Result<()> {
        Python::with_gil(|py| {
            let bound_optical_property = self.py_optical_property.bind(py);

            if with_optical_downcast(bound_optical_property, |db| {
                db.optical_quantities_emplace(inputs, aux_inputs, optical_quantities)
            })
            .is_ok()
            {
                return Ok(());
            }

            let bound_atmosphere = self.py_atmosphere.bind(py);

            let aq = bound_optical_property
                .call_method1("atmosphere_quantities", (bound_atmosphere,))?;

            let cross_section: PyReadonlyArray2<f64> = aq.getattr("extinction")?.extract()?;
            let cross_section = cross_section.as_array().to_owned();

            let ssa = Array2::zeros(cross_section.dim());

            optical_quantities.cross_section = cross_section;
            optical_quantities.ssa = ssa;
            Ok(())
        })
    }

    fn optical_derivatives_emplace(
        &self,
        inputs: &dyn StorageInputs,
        aux_inputs: &dyn AuxOpticalInputs,
        d_optical_quantities: &mut HashMap<String, OpticalQuantities>,
    ) -> Result<()> {
        Python::with_gil(|py| {
            let bound_optical_property = self.py_optical_property.bind(py);

            if with_optical_downcast(bound_optical_property, |db| {
                db.optical_derivatives_emplace(inputs, aux_inputs, d_optical_quantities)
            })
            .is_ok()
            {
                return Ok(());
            }

            let bound_atmosphere = self.py_atmosphere.bind(py);

            let d_aq: Bound<'_, PyDict> = bound_optical_property
                .call_method1("optical_derivatives", (bound_atmosphere,))?
                .downcast_into()
                .unwrap();

            // Convert to PyDict
            d_aq.iter().for_each(|(key, value)| {
                let cross_section: PyReadonlyArray2<f64> =
                    value.getattr("d_extinction").unwrap().extract().unwrap();
                let cross_section = cross_section.as_array().to_owned();

                let ssa = Array2::zeros(cross_section.dim());

                let k = key.extract::<String>().unwrap();

                let fortran_ordering = !cross_section.is_standard_layout();

                d_optical_quantities.insert(
                    k,
                    OpticalQuantities {
                        cross_section,
                        ssa,
                        legendre: None,
                        fortran_ordering,
                    },
                );
            });

            Ok(())
        })
    }
}
