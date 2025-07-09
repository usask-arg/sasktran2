use numpy::{PyArray1, PyArrayMethods, PyReadonlyArray1};
use pyo3::prelude::*;
use sasktran2_rs::bindings::geodetic;

use crate::prelude::IntoPyResult;

fn pyarray_to_slice(pyarray: PyReadonlyArray1<f64>) -> PyResult<[f64; 3]> {
    let mut result = [0.0; 3];

    // Ensure the array has exactly 3 elements
    if pyarray.len()? != 3 {
        return Err(PyErr::new::<pyo3::exceptions::PyValueError, _>(
            "Array must have exactly 3 elements",
        ));
    }

    // Copy into result
    let slice = pyarray.as_slice()?;
    result.copy_from_slice(slice);

    Ok(result)
}

fn slice_to_pyarray<'py>(py: Python<'py>, slice: &[f64; 3]) -> PyResult<Bound<'py, PyArray1<f64>>> {
    let result = PyArray1::zeros(py, (3,), false);

    unsafe {
        result.as_slice_mut()?.copy_from_slice(slice);
    }

    Ok(result)
}

#[pyclass(unsendable)]
pub struct PyGeodetic {
    pub output: geodetic::Geodetic,
}

#[pymethods]
impl PyGeodetic {
    #[new]
    fn new(equatorial_radius: f64, flattening_factor: f64) -> PyResult<Self> {
        let output =
            geodetic::Geodetic::new(equatorial_radius, flattening_factor).into_pyresult()?;
        Ok(Self { output })
    }

    #[getter]
    fn get_altitude(&self) -> PyResult<f64> {
        self.output.altitude().into_pyresult()
    }

    #[getter]
    fn get_latitude(&self) -> PyResult<f64> {
        self.output.latitude().into_pyresult()
    }

    #[getter]
    fn get_longitude(&self) -> PyResult<f64> {
        self.output.longitude().into_pyresult()
    }

    #[getter]
    fn get_location<'py>(&self, py: Python<'py>) -> PyResult<Bound<'py, PyArray1<f64>>> {
        let location = self.output.location().into_pyresult()?;
        slice_to_pyarray(py, &location)
    }

    #[getter]
    fn get_local_south<'py>(&self, py: Python<'py>) -> PyResult<Bound<'py, PyArray1<f64>>> {
        let local_south = self.output.local_south().into_pyresult()?;
        slice_to_pyarray(py, &local_south)
    }

    #[getter]
    fn get_local_up<'py>(&self, py: Python<'py>) -> PyResult<Bound<'py, PyArray1<f64>>> {
        let local_up = self.output.local_up().into_pyresult()?;
        slice_to_pyarray(py, &local_up)
    }

    #[getter]
    fn get_local_west<'py>(&self, py: Python<'py>) -> PyResult<Bound<'py, PyArray1<f64>>> {
        let local_west = self.output.local_west().into_pyresult()?;
        slice_to_pyarray(py, &local_west)
    }

    #[allow(clippy::type_complexity)]
    fn altitude_intercepts<'py>(
        &self,
        py: Python<'py>,
        altitude: f64,
        observer: PyReadonlyArray1<f64>,
        look_vector: PyReadonlyArray1<f64>,
    ) -> PyResult<(Bound<'py, PyArray1<f64>>, Bound<'py, PyArray1<f64>>)> {
        let observer = pyarray_to_slice(observer)?;
        let look_vector = pyarray_to_slice(look_vector)?;

        let (altitude_intercepts, look_vector_intercepts) = self
            .output
            .altitude_intercepts(altitude, observer, look_vector)
            .into_pyresult()?;

        let altitude_intercepts = slice_to_pyarray(py, &altitude_intercepts)?;
        let look_vector_intercepts = slice_to_pyarray(py, &look_vector_intercepts)?;

        Ok((altitude_intercepts, look_vector_intercepts))
    }

    #[allow(clippy::wrong_self_convention)]
    fn from_lat_lon_alt(&mut self, latitude: f64, longitude: f64, altitude: f64) -> PyResult<()> {
        self.output
            .from_lat_lon_alt(latitude, longitude, altitude)
            .map_err(|e| PyErr::new::<pyo3::exceptions::PyRuntimeError, _>(format!("{e}")))
    }

    #[allow(clippy::wrong_self_convention)]
    fn from_tangent_altitude(
        &mut self,
        tangent_altitude: f64,
        observer: PyReadonlyArray1<f64>,
        look_vector: PyReadonlyArray1<f64>,
    ) -> PyResult<[f64; 3]> {
        let observer = pyarray_to_slice(observer)?;
        let look_vector = pyarray_to_slice(look_vector)?;

        let result = self
            .output
            .from_tangent_altitude(tangent_altitude, observer, look_vector)
            .map_err(|e| PyErr::new::<pyo3::exceptions::PyRuntimeError, _>(format!("{e}")))?;

        Ok(result)
    }

    #[allow(clippy::wrong_self_convention)]
    fn from_tangent_point(
        &mut self,
        observer: PyReadonlyArray1<f64>,
        look_vector: PyReadonlyArray1<f64>,
    ) -> PyResult<()> {
        let observer = pyarray_to_slice(observer)?;
        let look_vector = pyarray_to_slice(look_vector)?;

        self.output
            .from_tangent_point(observer, look_vector)
            .map_err(|e| PyErr::new::<pyo3::exceptions::PyRuntimeError, _>(format!("{e}")))
    }

    #[allow(clippy::wrong_self_convention)]
    fn from_xyz(&mut self, location: PyReadonlyArray1<f64>) -> PyResult<()> {
        let location = pyarray_to_slice(location)?;
        let x = location[0];
        let y = location[1];
        let z = location[2];
        self.output
            .from_xyz(x, y, z)
            .map_err(|e| PyErr::new::<pyo3::exceptions::PyRuntimeError, _>(format!("{e}")))
    }

    fn is_valid(&self) -> PyResult<bool> {
        self.output
            .is_valid()
            .map_err(|e| PyErr::new::<pyo3::exceptions::PyRuntimeError, _>(format!("{e}")))
    }

    fn osculating_spheroid(&mut self) -> PyResult<(f64, [f64; 3])> {
        self.output
            .osculating_spheroid()
            .map_err(|e| PyErr::new::<pyo3::exceptions::PyRuntimeError, _>(format!("{e}")))
    }
}
