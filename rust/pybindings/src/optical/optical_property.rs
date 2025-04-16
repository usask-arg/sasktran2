use std::collections::HashMap;

use anyhow::Result;
use ndarray::Array2;
use numpy::PyReadonlyArray2;
use pyo3::types::PyDictMethods;
use pyo3::{prelude::*, types::PyDict};
use sk_core::{constituent::StorageInputs, optical::*};

use super::{AbsorberDatabaseDim2, AbsorberDatabaseDim3};

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
        optical_quantities: &mut OpticalQuantities,
    ) -> Result<()> {
        Python::with_gil(|py| {
            let bound_optical_property = self.py_optical_property.bind(py);
            let rust_optical = bound_optical_property
                .call_method0("_into_rust_object")
                .ok();

            if let Some(rust_optical) = rust_optical {
                if let Ok(absorber) = rust_optical.downcast::<AbsorberDatabaseDim2>() {
                    absorber
                        .borrow()
                        .db
                        .optical_quantities_emplace(inputs, optical_quantities)?;
                    return Ok(());
                }

                if let Ok(absorber) = rust_optical.downcast::<AbsorberDatabaseDim3>() {
                    absorber
                        .borrow()
                        .db
                        .optical_quantities_emplace(inputs, optical_quantities)?;
                    return Ok(());
                }
            }

            let bound_atmosphere = self.py_atmosphere.bind(py);

            let aq = bound_optical_property
                .call_method1("atmosphere_quantities", (bound_atmosphere,))?;

            let cross_section: PyReadonlyArray2<f64> = aq.getattr("extinction")?.extract()?;
            let cross_section = cross_section.as_array().to_owned();

            let ssa = Array2::zeros(cross_section.dim());

            optical_quantities.cross_section = cross_section.into();
            optical_quantities.ssa = ssa.into();
            Ok(())
        })
    }

    fn optical_derivatives_emplace(
        &self,
        inputs: &dyn StorageInputs,
        d_optical_quantities: &mut HashMap<String, OpticalQuantities>,
    ) -> Result<()> {
        Python::with_gil(|py| {
            let bound_optical_property = self.py_optical_property.bind(py);

            let rust_optical = bound_optical_property
                .call_method0("_into_rust_object")
                .ok();

            if let Some(rust_optical) = rust_optical {
                if let Ok(absorber) = rust_optical.downcast::<AbsorberDatabaseDim2>() {
                    absorber
                        .borrow()
                        .db
                        .optical_derivatives_emplace(inputs, d_optical_quantities)?;
                    return Ok(());
                }

                if let Ok(absorber) = rust_optical.downcast::<AbsorberDatabaseDim3>() {
                    absorber
                        .borrow()
                        .db
                        .optical_derivatives_emplace(inputs, d_optical_quantities)?;
                    return Ok(());
                }
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

                d_optical_quantities.insert(
                    k,
                    OpticalQuantities {
                        cross_section: cross_section.into(),
                        ssa: ssa.into(),
                        a1: None,
                        a2: None,
                        a3: None,
                        a4: None,
                        b1: None,
                        b2: None,
                    },
                );
            });

            Ok(())
        })
    }
}
