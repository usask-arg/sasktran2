use std::collections::HashMap;

use anyhow::Result;
use ndarray::Array2;
use numpy::PyReadonlyArray2;
use pyo3::types::PyDictMethods;
use pyo3::{prelude::*, types::PyDict};
use sk_core::{constituent::StorageInputs, optical::*};

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
    fn optical_quantities<'a>(
        &'a self,
        _inputs: &dyn StorageInputs,
    ) -> Result<OpticalQuantities<'a>> {
        Python::with_gil(|py| {
            let bound_optical_property = self.py_optical_property.bind(py);
            let bound_atmosphere = self.py_atmosphere.bind(py);

            let aq = bound_optical_property
                .call_method1("atmosphere_quantities", (bound_atmosphere,))?;

            let cross_section: PyReadonlyArray2<f64> = aq.getattr("extinction")?.extract()?;
            let cross_section = cross_section.as_array().to_owned();

            let ssa = Array2::zeros(cross_section.dim());

            Ok(OpticalQuantities {
                cross_section: cross_section.into(),
                ssa: ssa.into(),
                a1: None,
                a2: None,
                a3: None,
                a4: None,
                b1: None,
                b2: None,
            })
        })
    }

    fn optical_derivatives<'a>(
        &'a self,
        _inputs: &dyn StorageInputs,
    ) -> Result<std::collections::HashMap<String, OpticalQuantities<'a>>> {
        Python::with_gil(|py| {
            let bound_optical_property = self.py_optical_property.bind(py);
            let bound_atmosphere = self.py_atmosphere.bind(py);

            let d_aq: Bound<'_, PyDict> = bound_optical_property
                .call_method1("optical_derivatives", (bound_atmosphere,))?
                .downcast_into()
                .unwrap();

            let mut hash: HashMap<String, OpticalQuantities<'a>> = std::collections::HashMap::new();
            // Convert to PyDict
            d_aq.iter().for_each(|(key, value)| {
                let cross_section: PyReadonlyArray2<f64> =
                    value.getattr("d_extinction").unwrap().extract().unwrap();
                let cross_section = cross_section.as_array().to_owned();

                let ssa = Array2::zeros(cross_section.dim());

                let k = key.extract::<String>().unwrap();

                hash.insert(
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

            Ok(hash)
        })
    }
}
