use ndarray::Array2;
use numpy::*;
use pyo3::prelude::*;
use sasktran2_rs::atmosphere::{DerivMapping, DerivMappingView};

pub struct PyDerivMapping<'py> {
    py_mapping: Bound<'py, PyAny>,
    pub d_extinction: PyReadwriteArray2<'py, f64>,
    pub d_ssa: PyReadwriteArray2<'py, f64>,
    pub scat_factor: Option<PyReadwriteArray2<'py, f64>>,
    pub d_legendre: Option<PyReadwriteArray3<'py, f64>>,
    pub d_emission: PyReadwriteArray2<'py, f64>,
}

impl<'py> PyDerivMapping<'py> {
    pub fn new(py_mapping: Bound<'py, PyAny>) -> Self {
        let d_extinction: PyReadwriteArray2<'py, f64> = py_mapping
            .getattr("d_extinction")
            .unwrap()
            .extract()
            .unwrap();

        let d_ssa: PyReadwriteArray2<'py, f64> =
            py_mapping.getattr("d_ssa").unwrap().extract().unwrap();

        let d_emission: PyReadwriteArray2<'py, f64> =
            py_mapping.getattr("d_emission").unwrap().extract().unwrap();

        PyDerivMapping {
            py_mapping,
            d_extinction,
            d_ssa,
            scat_factor: None,
            d_legendre: None,
            d_emission,
        }
    }
}

impl<'py> DerivMapping<'_> for PyDerivMapping<'py> {
    fn with_scatterer(mut self) -> Self {
        self.scat_factor = Some(
            self.py_mapping
                .getattr("scat_factor")
                .unwrap()
                .extract()
                .unwrap(),
        );
        self.d_legendre = Some(
            self.py_mapping
                .getattr("d_leg_coeff")
                .unwrap()
                .extract()
                .unwrap(),
        );
        self
    }

    fn mut_view(&mut self) -> DerivMappingView<'_> {
        DerivMappingView {
            d_extinction: self.d_extinction.as_array_mut(),
            d_ssa: self.d_ssa.as_array_mut(),
            d_legendre: self.d_legendre.as_mut().map(|arr| arr.as_array_mut()),
            scat_factor: self.scat_factor.as_mut().map(|arr| arr.as_array_mut()),
            d_emission: self.d_emission.as_array_mut(),
        }
    }

    fn set_interpolator(&mut self, interpolator: &Array2<f64>) {
        let pyarray = interpolator.clone().into_pyarray(self.py_mapping.py());

        self.py_mapping
            .setattr("interpolator", pyarray)
            .expect("Failed to set interpolator");
    }

    fn set_interp_dim(&mut self, interp_dim: &str) {
        self.py_mapping
            .setattr("interp_dim", interp_dim)
            .expect("Failed to set interp_dim");
    }

    fn set_assign_name(&mut self, assign_name: &str) {
        self.py_mapping
            .setattr("assign_name", assign_name)
            .expect("Failed to set assign_name");
    }
}
