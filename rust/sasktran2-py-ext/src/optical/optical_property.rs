use std::collections::HashMap;

use anyhow::{Result, anyhow};
use ndarray::{Array2, Array3};
use numpy::{IntoPyArray, PyReadonlyArray2, PyReadonlyArray3};
use pyo3::types::PyDictMethods;
use pyo3::{prelude::*, types::PyDict};
use sasktran2_rs::atmosphere::traits::StorageInputs;
use sasktran2_rs::optical::storage::*;
use sasktran2_rs::optical::traits::*;

use crate::optical::line_absorber::PyLineAbsorber;
use crate::optical::scat_dbase::{
    PyScatteringDatabaseDim1, PyScatteringDatabaseDim2, PyScatteringDatabaseDim3,
};
use crate::optical::xsec_absorber::PyXsecAbsorber;
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
                if let Ok(obj) = rust_obj.cast::<$ty>() {
                    return f(obj.borrow().db());
                }
            )*
        };
    }

    try_downcast!(
        AbsorberDatabaseDim1,
        AbsorberDatabaseDim2,
        AbsorberDatabaseDim3,
        PyScatteringDatabaseDim1,
        PyScatteringDatabaseDim2,
        PyScatteringDatabaseDim3,
        PyLineAbsorber,
        PyXsecAbsorber
    );

    Err(anyhow!("Failed to downcast to any of the known types"))
}

pub struct PyOpticalProperty {
    py_optical_property: Py<PyAny>,
    py_atmosphere: Py<PyAny>,
    aux_names: Vec<String>,
}

impl PyOpticalProperty {
    pub fn new(py_optical_property: Py<PyAny>, py_atmosphere: Py<PyAny>) -> Self {
        Self::new_with_aux_names(py_optical_property, py_atmosphere, Vec::new())
    }

    pub fn new_with_aux_names(
        py_optical_property: Py<PyAny>,
        py_atmosphere: Py<PyAny>,
        aux_names: Vec<String>,
    ) -> Self {
        Self {
            py_optical_property,
            py_atmosphere,
            aux_names,
        }
    }

    fn aux_kwargs<'py>(
        &self,
        py: Python<'py>,
        aux_inputs: &dyn AuxOpticalInputs,
    ) -> PyResult<Bound<'py, PyDict>> {
        let py_kwargs = PyDict::new(py);

        for name in &self.aux_names {
            if let Some(val) = aux_inputs.get_parameter(name) {
                py_kwargs.set_item(name, val.to_owned().into_pyarray(py))?;
            }
        }

        Ok(py_kwargs)
    }
}

fn extract_optional_array2(
    value: Option<Bound<'_, PyAny>>,
    dims: (usize, usize),
) -> Result<Array2<f64>> {
    if value.as_ref().is_none_or(|value| value.is_none()) {
        Ok(Array2::zeros(dims))
    } else {
        let array: PyReadonlyArray2<f64> = value
            .unwrap()
            .extract()
            .map_err(|e| anyhow!("Failed to extract array2 from optical quantities: {e}"))?;
        Ok(array.as_array().to_owned())
    }
}

fn extract_optional_legendre(value: Option<Bound<'_, PyAny>>) -> Result<Option<Array3<f64>>> {
    if value.as_ref().is_none_or(|value| value.is_none()) {
        Ok(None)
    } else {
        let array: PyReadonlyArray3<f64> = value.unwrap().extract().map_err(|e| {
            anyhow!("Failed to extract legendre array from optical quantities: {e}")
        })?;
        Ok(Some(array.as_array().permuted_axes([1, 2, 0]).to_owned()))
    }
}

impl OpticalProperty for PyOpticalProperty {
    fn optical_quantities_emplace(
        &self,
        inputs: &dyn StorageInputs,
        aux_inputs: &dyn AuxOpticalInputs,
        optical_quantities: &mut OpticalQuantities,
    ) -> Result<()> {
        Python::attach(|py| {
            let bound_optical_property = self.py_optical_property.bind(py);

            if with_optical_downcast(bound_optical_property, |db| {
                db.optical_quantities_emplace(inputs, aux_inputs, optical_quantities)
            })
            .is_ok()
            {
                return Ok(());
            }

            let bound_atmosphere = self.py_atmosphere.bind(py);
            let py_kwargs = self.aux_kwargs(py, aux_inputs)?;

            let aq = bound_optical_property.call_method(
                "atmosphere_quantities",
                (bound_atmosphere,),
                Some(&py_kwargs),
            )?;

            let cross_section: PyReadonlyArray2<f64> =
                aq.getattr("extinction")?.extract().map_err(|e| {
                    anyhow!("Failed to extract extinction array from atmosphere_quantities: {e}")
                })?;
            let cross_section = cross_section.as_array().to_owned();
            let dims = cross_section.dim();

            let ssa = extract_optional_array2(aq.getattr("ssa").ok(), dims)?;
            let legendre = extract_optional_legendre(aq.getattr("leg_coeff").ok())?;

            optical_quantities.cross_section = cross_section;
            optical_quantities.ssa = ssa;
            optical_quantities.legendre = legendre;
            Ok(())
        })
    }

    fn optical_derivatives_emplace(
        &self,
        inputs: &dyn StorageInputs,
        aux_inputs: &dyn AuxOpticalInputs,
        d_optical_quantities: &mut HashMap<String, OpticalQuantities>,
    ) -> Result<()> {
        Python::attach(|py| {
            let bound_optical_property = self.py_optical_property.bind(py);

            if with_optical_downcast(bound_optical_property, |db| {
                db.optical_derivatives_emplace(inputs, aux_inputs, d_optical_quantities)
            })
            .is_ok()
            {
                return Ok(());
            }

            let bound_atmosphere = self.py_atmosphere.bind(py);
            let py_kwargs = self.aux_kwargs(py, aux_inputs)?;

            let d_aq: Bound<'_, PyDict> = bound_optical_property
                .call_method("optical_derivatives", (bound_atmosphere,), Some(&py_kwargs))?
                .cast_into()
                .unwrap();

            // Convert to PyDict
            d_aq.iter().for_each(|(key, value)| {
                let cross_section: PyReadonlyArray2<f64> =
                    value.getattr("d_extinction").unwrap().extract().unwrap();
                let cross_section = cross_section.as_array().to_owned();
                let dims = cross_section.dim();

                let ssa = extract_optional_array2(value.getattr("d_ssa").ok(), dims).unwrap();
                let legendre =
                    extract_optional_legendre(value.getattr("d_leg_coeff").ok()).unwrap();

                let k = key.extract::<String>().unwrap();

                let fortran_ordering = !cross_section.is_standard_layout();

                d_optical_quantities.insert(
                    k,
                    OpticalQuantities {
                        cross_section,
                        ssa,
                        legendre,
                        fortran_ordering,
                    },
                );
            });

            Ok(())
        })
    }

    fn is_scatterer(&self) -> bool {
        // TODO: grab from python?
        false
    }
}
