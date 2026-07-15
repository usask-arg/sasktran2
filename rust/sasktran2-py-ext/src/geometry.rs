use crate::prelude::*;
use numpy::PyArray1;
use numpy::{PyReadonlyArray1, ToPyArray};
use pyo3::prelude::*;
use sasktran2_rs::bindings::geometry;

#[pyclass(eq, eq_int)]
#[derive(PartialEq, Clone)]
pub enum GeometryType {
    Spherical,
    PlaneParallel,
    PseudoSpherical,
    Ellipsoidal,
}

#[pyclass(eq, eq_int)]
#[derive(PartialEq, Clone)]
#[allow(clippy::enum_variant_names)]
pub enum InterpolationMethod {
    LinearInterpolation,
    ShellInterpolation,
    LowerInterpolation,
}

#[pyclass(unsendable)]
pub struct PyGeometry1D {
    pub geometry: geometry::Geometry1D,
}

#[pyclass(unsendable)]
pub struct PyGeometry2D {
    pub geometry: geometry::Geometry2D,
}

#[pymethods]
impl PyGeometry1D {
    #[new]
    fn new(
        cos_sza: f64,
        solar_azimuth: f64,
        earth_radius_m: f64,
        altitude_grid_m: PyReadonlyArray1<f64>,
        interpolation_method: PyRef<'_, InterpolationMethod>,
        geometry_type: PyRef<'_, GeometryType>,
    ) -> PyResult<Self> {
        let interpolation_method = match *interpolation_method {
            InterpolationMethod::LinearInterpolation => geometry::InterpolationMethod::Linear,
            InterpolationMethod::ShellInterpolation => geometry::InterpolationMethod::Shell,
            InterpolationMethod::LowerInterpolation => geometry::InterpolationMethod::Lower,
        };

        let geometry_type = match *geometry_type {
            GeometryType::Spherical => geometry::GeometryType::Spherical,
            GeometryType::PlaneParallel => geometry::GeometryType::PlaneParallel,
            GeometryType::PseudoSpherical => geometry::GeometryType::PseudoSpherical,
            GeometryType::Ellipsoidal => geometry::GeometryType::Ellipsoidal,
        };

        let altitude_grid_m = altitude_grid_m.as_slice()?;
        let geometry = geometry::Geometry1D::new(
            cos_sza,
            solar_azimuth,
            earth_radius_m,
            altitude_grid_m.to_vec(),
            interpolation_method,
            geometry_type,
        );

        Ok(PyGeometry1D { geometry })
    }

    fn altitudes<'py>(&self, py: Python<'py>) -> PyResult<Bound<'py, PyArray1<f64>>> {
        let altitudes = self.geometry.altitudes_m().into_pyresult()?;

        Ok(altitudes.to_pyarray(py).to_owned())
    }

    #[getter]
    fn get_refractive_index<'py>(this: Bound<'py, Self>) -> PyResult<Bound<'py, PyArray1<f64>>> {
        let binding = &this.borrow().geometry;
        let array = binding.refractive_index_mut().into_pyresult()?;

        unsafe { Ok(PyArray1::borrow_from_array(&array, this.into_any())) }
    }

    #[setter]
    fn set_refractive_index(&mut self, refractive_index: PyReadonlyArray1<f64>) -> PyResult<()> {
        let binding = &mut self.geometry;
        let refractive_index = refractive_index.as_array();

        let mut view = binding.refractive_index_mut().into_pyresult()?;
        view.assign(&refractive_index);

        Ok(())
    }
}

#[pymethods]
impl PyGeometry2D {
    #[new]
    fn new(
        cos_sza: f64,
        solar_azimuth: f64,
        earth_radius_m: f64,
        altitude_grid_m: PyReadonlyArray1<f64>,
        horizontal_angle_grid_radians: PyReadonlyArray1<f64>,
        interpolation_method: PyRef<'_, InterpolationMethod>,
    ) -> PyResult<Self> {
        let interpolation_method = match *interpolation_method {
            InterpolationMethod::LinearInterpolation => geometry::InterpolationMethod::Linear,
            InterpolationMethod::ShellInterpolation => geometry::InterpolationMethod::Shell,
            InterpolationMethod::LowerInterpolation => geometry::InterpolationMethod::Lower,
        };
        let altitude_grid_m = altitude_grid_m.as_slice()?;
        let horizontal_angle_grid_radians = horizontal_angle_grid_radians.as_slice()?;
        let geometry = geometry::Geometry2D::new(
            cos_sza,
            solar_azimuth,
            earth_radius_m,
            altitude_grid_m.to_vec(),
            horizontal_angle_grid_radians.to_vec(),
            interpolation_method,
        )
        .into_pyresult()?;

        Ok(Self { geometry })
    }

    fn altitudes<'py>(&self, py: Python<'py>) -> PyResult<Bound<'py, PyArray1<f64>>> {
        let altitudes = self.geometry.altitudes_m().into_pyresult()?;
        Ok(altitudes.to_pyarray(py).to_owned())
    }

    fn horizontal_angles<'py>(&self, py: Python<'py>) -> PyResult<Bound<'py, PyArray1<f64>>> {
        let horizontal_angles = self.geometry.horizontal_angles().into_pyresult()?;
        Ok(horizontal_angles.to_pyarray(py).to_owned())
    }

    fn location_shape(&self) -> PyResult<(usize, usize)> {
        self.geometry.location_shape().into_pyresult()
    }

    fn location_index(&self, altitude_index: usize, horizontal_index: usize) -> PyResult<usize> {
        self.geometry
            .location_index(altitude_index, horizontal_index)
            .into_pyresult()
    }
}
