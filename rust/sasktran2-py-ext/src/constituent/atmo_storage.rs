#![allow(unused)]

use std::collections::HashMap;

use numpy::ndarray::{Array1, ArrayView1, ArrayView2, ArrayView3, ArrayViewMut2, ArrayViewMut3};
use numpy::{Element, PyUntypedArrayMethods};
use numpy::{PyReadonlyArray1, PyReadwriteArray2, PyReadwriteArray3};
use pyo3::prelude::*;
use pyo3::types::PyDict;
use sasktran2_rs::atmosphere::AtmosphereStorageAccess;
use sasktran2_rs::atmosphere::traits::*;
use sasktran2_rs::bindings::config;
use sasktran2_rs::constituent::traits::*;

use crate::config::PyConfig;

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
    pub py_emission_source: PyReadwriteArray2<'py, f64>,
}

impl StorageOutputs for AtmosphereStorageOutputs<'_> {
    fn mut_view(&mut self) -> AtmosphereStorageOutputView<'_> {
        AtmosphereStorageOutputView {
            total_extinction: self.py_total_extinction.as_array_mut(),
            ssa: self.py_ssa.as_array_mut(),
            legendre: self.py_legendre.as_array_mut(),
            emission_source: self.py_emission_source.as_array_mut(),
        }
    }

    fn view<'a>(&'a self) -> AtmosphereStorageOutputImmutView<'a> {
        AtmosphereStorageOutputImmutView {
            total_extinction: self.py_total_extinction.as_array(),
            ssa: self.py_ssa.as_array(),
            legendre: self.py_legendre.as_array(),
            emission_source: self.py_emission_source.as_array(),
        }
    }
}

pub struct AtmosphereStorageInputs<'py> {
    pub num_stokes: usize,
    pub num_legendre: usize,
    pub calculate_pressure_derivative: bool,
    pub calculate_temperature_derivative: bool,
    pub calculate_specific_humidity_derivative: bool,
    pub py_altitude_m: PyReadonlyArray1<'py, f64>,
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

    fn pressure_pa(&self) -> Option<ArrayView1<'_, f64>> {
        self.py_pressure_pa.as_ref().map(|array| array.as_array())
    }

    fn temperature_k(&self) -> Option<ArrayView1<'_, f64>> {
        self.py_temperature_k.as_ref().map(|array| array.as_array())
    }

    fn wavelengths_nm(&self) -> Option<ArrayView1<'_, f64>> {
        self.py_wavelength_nm.as_ref().map(|array| array.as_array())
    }

    fn wavenumbers_cminv(&self) -> Option<ArrayView1<'_, f64>> {
        self.py_wavenumber_cminv
            .as_ref()
            .map(|array| array.as_array())
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

    fn dry_air_numberdensity_dict(&self) -> HashMap<String, Array1<f64>> {
        let state_eqn = &self.py_equation_of_state;

        let air_numden: Bound<'py, PyDict> = state_eqn
            .getattr("dry_air_numberdensity")
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

    fn altitude_m(&self) -> ArrayView1<'_, f64> {
        self.py_altitude_m.as_array()
    }

    fn calculate_pressure_derivative(&self) -> bool {
        self.calculate_pressure_derivative
    }

    fn calculate_temperature_derivative(&self) -> bool {
        self.calculate_temperature_derivative
    }

    fn calculate_specific_humidity_derivative(&self) -> bool {
        self.calculate_specific_humidity_derivative
    }

    fn num_singlescatter_moments(&self) -> usize {
        self.num_legendre
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
    pub fn new(atmo: &Bound<'py, PyAny>) -> PyResult<Self> {
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

        let num_legendre = legendre.shape()[0];

        let emission_source_obj = storage.getattr("emission_source").unwrap();
        let emission_source: PyReadwriteArray2<f64> = emission_source_obj.extract().unwrap();

        let state_eqn_obj = atmo.getattr("state_equation").unwrap();

        let num_stokes_obj = atmo.getattr("nstokes").unwrap();
        let num_stokes: usize = num_stokes_obj.extract().unwrap();

        let geometry_obj = atmo.getattr("model_geometry").unwrap();
        let altitudes_m = geometry_obj.call_method0("altitudes").unwrap();
        let altitudes_m: PyReadonlyArray1<f64> = altitudes_m.extract().unwrap();

        let calculate_pressure_derivative = atmo
            .getattr("calculate_pressure_derivative")
            .unwrap()
            .extract()
            .unwrap_or(false);

        let calculate_temperature_derivative = atmo
            .getattr("calculate_temperature_derivative")
            .unwrap()
            .extract()
            .unwrap_or(false);

        let calculate_specific_humidity_derivative = atmo
            .getattr("calculate_specific_humidity_derivative")
            .unwrap()
            .extract()
            .unwrap_or(false);

        Ok(AtmosphereStorage {
            num_stokes,
            deriv_generator: PyDerivativeGenerator { storage },
            inputs: AtmosphereStorageInputs {
                num_stokes,
                num_legendre,
                calculate_pressure_derivative,
                calculate_temperature_derivative,
                calculate_specific_humidity_derivative,
                py_altitude_m: altitudes_m,
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
                py_emission_source: emission_source,
            },
        })
    }

    pub fn num_geometry(&self) -> usize {
        self.outputs.py_total_extinction.shape()[0]
    }
    pub fn num_wavelengths(&self) -> usize {
        self.outputs.py_total_extinction.shape()[1]
    }
}

impl AtmosphereStorageAccess for AtmosphereStorage<'_> {
    fn split_inputs_outputs(&mut self) -> (&impl StorageInputs, &mut impl StorageOutputs) {
        let inputs = &self.inputs;
        let outputs = &mut self.outputs;

        (inputs, outputs)
    }

    fn split_inputs_outputs_deriv<'a>(
        &'a self,
    ) -> (
        &'a impl StorageInputs,
        &'a impl StorageOutputs,
        &'a impl DerivMappingGenerator<'a>,
    ) {
        (&self.inputs, &self.outputs, &self.deriv_generator)
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
