#![allow(unused)]

use std::collections::HashMap;

use numpy::ndarray::{Array1, ArrayView1, ArrayView2, ArrayView3, ArrayViewMut2, ArrayViewMut3};
use numpy::{Element, PyUntypedArrayMethods};
use numpy::{PyReadonlyArray1, PyReadwriteArray2, PyReadwriteArray3};
use pyo3::prelude::*;
use pyo3::types::PyDict;
use sk_core::constituent::{AtmosphereStorageOutputImmutView, AtmosphereStorageOutputView};

use sk_core::constituent::{DerivMapping, DerivMappingGenerator, StorageInputs, StorageOutputs};

use super::deriv_mapping::PyDerivMapping;

fn get_optional_array1<'py, T: Element>(
    obj: &Bound<'py, PyAny>,
    attr_name: &str,
) -> PyResult<Option<PyReadonlyArray1<'py, T>>> {
    match obj.getattr(attr_name) {
        Ok(value) => {
            if value.is_none() {
                Ok(None)
            } else {
                let array: PyReadonlyArray1<T> = value.extract()?;
                Ok(Some(array))
            }
        }
        Err(err) => {
            if err.is_instance_of::<pyo3::exceptions::PyAttributeError>(obj.py()) {
                Ok(None)
            } else {
                Err(err)
            }
        }
    }
}

pub struct AtmosphereStorageOutputs<'py> {
    num_stokes: usize,
    pub py_total_extinction: PyReadwriteArray2<'py, f64>,
    pub py_ssa: PyReadwriteArray2<'py, f64>,
    pub py_legendre: PyReadwriteArray3<'py, f64>,
}

impl StorageOutputs for AtmosphereStorageOutputs<'_> {
    fn mut_view(&mut self) -> AtmosphereStorageOutputView {
        AtmosphereStorageOutputView {
            total_extinction: self.py_total_extinction.as_array_mut(),
            ssa: self.py_ssa.as_array_mut(),
            legendre: self.py_legendre.as_array_mut(),
        }
    }

    fn view<'a>(&'a self) -> AtmosphereStorageOutputImmutView<'a> {
        AtmosphereStorageOutputImmutView {
            total_extinction: self.py_total_extinction.as_array(),
            ssa: self.py_ssa.as_array(),
            legendre: self.py_legendre.as_array(),
        }
    }
}

pub struct AtmosphereStorageInputs<'py> {
    pub num_stokes: usize,
    pub py_pressure_pa: Option<PyReadonlyArray1<'py, f64>>,
    pub py_temperature_k: Option<PyReadonlyArray1<'py, f64>>,
    pub py_wavelength_nm: Option<PyReadonlyArray1<'py, f64>>,
    pub py_wavenumber_cminv: Option<PyReadonlyArray1<'py, f64>>,
    pub py_equation_of_state: Bound<'py, PyAny>,
}

impl<'py> StorageInputs for AtmosphereStorageInputs<'py> {
    fn num_stokes(&self) -> usize {
        self.num_stokes
    }

    fn pressure_pa(&self) -> Option<ArrayView1<f64>> {
        if let Some(array) = &self.py_pressure_pa {
            Some(array.as_array())
        } else {
            None
        }
    }

    fn temperature_k(&self) -> Option<ArrayView1<f64>> {
        if let Some(array) = &self.py_temperature_k {
            Some(array.as_array())
        } else {
            None
        }
    }

    fn wavelengths_nm(&self) -> Option<ArrayView1<f64>> {
        if let Some(array) = &self.py_wavelength_nm {
            Some(array.as_array())
        } else {
            None
        }
    }

    fn wavenumbers_cminv(&self) -> Option<ArrayView1<f64>> {
        if let Some(array) = &self.py_wavenumber_cminv {
            Some(array.as_array())
        } else {
            None
        }
    }

    fn air_numberdensity_dict(&self) -> HashMap<String, Array1<f64>> {
        let state_eqn = &self.py_equation_of_state;

        let air_numden: Bound<'py, PyDict> = state_eqn
            .getattr("air_numberdensity")
            .unwrap()
            .downcast_into()
            .unwrap();

        let mut dict = HashMap::new();
        for (key, value) in air_numden.iter() {
            let array: PyReadonlyArray1<f64> = value.extract().unwrap();
            dict.insert(key.to_string(), array.as_array().to_owned());
        }
        dict
    }
}

pub struct PyDerivativeGenerator<'py> {
    storage: Bound<'py, PyAny>,
}

pub struct AtmosphereStorage<'py> {
    pub num_stokes: usize,
    pub deriv_generator: PyDerivativeGenerator<'py>,
    pub inputs: AtmosphereStorageInputs<'py>,
    pub outputs: AtmosphereStorageOutputs<'py>,
}

impl<'py> AtmosphereStorage<'py> {
    pub fn new(atmo: &Bound<'py, PyAny>) -> Self {
        let pressure_pa_array = get_optional_array1::<f64>(atmo, "pressure_pa").unwrap();
        let temperature_k_array = get_optional_array1::<f64>(atmo, "temperature_k").unwrap();
        let wavelengths_nm_array = get_optional_array1::<f64>(atmo, "wavelengths_nm").unwrap();
        let wavenumber_cminv_array = get_optional_array1::<f64>(atmo, "wavenumbers_cminv").unwrap();

        let storage = atmo.getattr("storage").unwrap();
        let total_extinction_obj = storage.getattr("total_extinction").unwrap();
        let total_extinction: PyReadwriteArray2<f64> = total_extinction_obj.extract().unwrap();

        let ssa_obj = storage.getattr("ssa").unwrap();
        let ssa: PyReadwriteArray2<f64> = ssa_obj.extract().unwrap();

        let legendre_obj = storage.getattr("leg_coeff").unwrap();
        let legendre: PyReadwriteArray3<f64> = legendre_obj.extract().unwrap();

        let state_eqn_obj = atmo.getattr("state_equation").unwrap();

        let num_stokes_obj = atmo.getattr("nstokes").unwrap();
        let num_stokes: usize = num_stokes_obj.extract().unwrap();

        AtmosphereStorage {
            num_stokes,
            deriv_generator: PyDerivativeGenerator { storage },
            inputs: AtmosphereStorageInputs {
                num_stokes,
                py_pressure_pa: pressure_pa_array,
                py_temperature_k: temperature_k_array,
                py_wavelength_nm: wavelengths_nm_array,
                py_wavenumber_cminv: wavenumber_cminv_array,
                py_equation_of_state: state_eqn_obj,
            },
            outputs: AtmosphereStorageOutputs {
                num_stokes,
                py_total_extinction: total_extinction,
                py_ssa: ssa,
                py_legendre: legendre,
            },
        }
    }

    pub fn num_geometry(&self) -> usize {
        self.outputs.py_total_extinction.shape()[0]
    }
    pub fn num_wavelengths(&self) -> usize {
        self.outputs.py_total_extinction.shape()[1]
    }
}

impl<'py> DerivMappingGenerator<'py> for PyDerivativeGenerator<'py> {
    fn get_derivative_mapping(&self, name: &str) -> impl DerivMapping<'py> {
        let args = (name,);

        let deriv = self
            .storage
            .call_method1("get_derivative_mapping", args)
            .unwrap();

        PyDerivMapping::new(deriv)
    }
}
