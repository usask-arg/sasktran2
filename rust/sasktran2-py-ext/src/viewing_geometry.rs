use crate::prelude::IntoPyResult;
use pyo3::exceptions::PyValueError;
use pyo3::prelude::*;
use sasktran2_rs::bindings::viewing_geometry;

#[pyclass]
pub struct PyGroundViewingSolar {
    cos_sza: f64,
    relative_azimuth: f64,
    observer_altitude_m: f64,
    cos_viewing_zenith: f64,
}

#[pymethods]
impl PyGroundViewingSolar {
    #[new]
    fn new(
        cos_sza: f64,
        relative_azimuth: f64,
        cos_viewing_zenith: f64,
        observer_altitude_m: f64,
    ) -> Self {
        Self {
            cos_sza,
            relative_azimuth,
            observer_altitude_m,
            cos_viewing_zenith,
        }
    }
}

#[pyclass]
pub struct PyTangentAltitudeSolar {
    tangent_altitude_m: f64,
    relative_azimuth: f64,
    observer_altitude_m: f64,
    cos_sza: f64,
}

#[pyclass]
pub struct PyTangentAltitude {
    tangent_altitude_m: f64,
    observer_altitude_m: f64,
    horizontal_angle_radians: f64,
    viewing_azimuth_radians: f64,
}

#[pymethods]
impl PyTangentAltitude {
    #[new]
    fn new(
        tangent_altitude_m: f64,
        observer_altitude_m: f64,
        horizontal_angle_radians: f64,
        viewing_azimuth_radians: f64,
    ) -> PyResult<Self> {
        if !tangent_altitude_m.is_finite() || tangent_altitude_m < 0.0 {
            return Err(PyValueError::new_err(
                "tangent_altitude_m must be finite and non-negative",
            ));
        }
        if !observer_altitude_m.is_finite() || observer_altitude_m < tangent_altitude_m {
            return Err(PyValueError::new_err(
                "observer_altitude_m must be finite and greater than or equal to tangent_altitude_m",
            ));
        }
        if !horizontal_angle_radians.is_finite() || !viewing_azimuth_radians.is_finite() {
            return Err(PyValueError::new_err("viewing angles must be finite"));
        }
        Ok(Self {
            tangent_altitude_m,
            observer_altitude_m,
            horizontal_angle_radians,
            viewing_azimuth_radians,
        })
    }
}

#[pymethods]
impl PyTangentAltitudeSolar {
    #[new]
    fn new(
        tangent_altitude_m: f64,
        relative_azimuth: f64,
        observer_altitude_m: f64,
        cos_sza: f64,
    ) -> Self {
        Self {
            tangent_altitude_m,
            relative_azimuth,
            observer_altitude_m,
            cos_sza,
        }
    }
}

#[pyclass]
pub struct PySolarAnglesObserverLocation {
    cos_sza: f64,
    relative_azimuth: f64,
    cos_viewing_zenith: f64,
    observer_altitude_m: f64,
}

#[pymethods]
impl PySolarAnglesObserverLocation {
    #[new]
    fn new(
        cos_sza: f64,
        relative_azimuth: f64,
        cos_viewing_zenith: f64,
        observer_altitude_m: f64,
    ) -> Self {
        Self {
            cos_sza,
            relative_azimuth,
            cos_viewing_zenith,
            observer_altitude_m,
        }
    }
}

#[pyclass]
pub struct PyFluxObserverSolar {
    cos_sza: f64,
    observer_altitude_m: f64,
}

#[pymethods]
impl PyFluxObserverSolar {
    #[new]
    fn new(cos_sza: f64, observer_altitude_m: f64) -> Self {
        Self {
            cos_sza,
            observer_altitude_m,
        }
    }
}

#[pyclass(unsendable)]
pub struct PyViewingGeometry {
    pub viewing_geometry: viewing_geometry::ViewingGeometry,
}

#[pymethods]
impl PyViewingGeometry {
    #[new]
    fn new() -> Self {
        Self {
            viewing_geometry: viewing_geometry::ViewingGeometry::new(),
        }
    }

    fn add_ray<'py>(&mut self, ray: Bound<'py, PyAny>) -> PyResult<()> {
        if let Ok(ray) = ray.cast::<PyTangentAltitude>() {
            let ray = ray.borrow();
            self.viewing_geometry
                .add_tangent_altitude(
                    ray.tangent_altitude_m,
                    ray.observer_altitude_m,
                    ray.horizontal_angle_radians,
                    ray.viewing_azimuth_radians,
                )
                .into_pyresult()?;
        }

        if let Ok(ray) = ray.cast::<PyGroundViewingSolar>() {
            let ray = ray.borrow();
            self.viewing_geometry.add_ground_viewing_solar(
                ray.cos_sza,
                ray.relative_azimuth,
                ray.observer_altitude_m,
                ray.cos_viewing_zenith,
            );
        }

        if let Ok(ray) = ray.cast::<PyTangentAltitudeSolar>() {
            let ray = ray.borrow();
            self.viewing_geometry.add_tangent_altitude_solar(
                ray.tangent_altitude_m,
                ray.relative_azimuth,
                ray.observer_altitude_m,
                ray.cos_sza,
            );
        }

        if let Ok(ray) = ray.cast::<PySolarAnglesObserverLocation>() {
            let ray = ray.borrow();
            self.viewing_geometry.add_solar_angles_observer_location(
                ray.cos_sza,
                ray.relative_azimuth,
                ray.cos_viewing_zenith,
                ray.observer_altitude_m,
            );
        }
        Ok(())
    }

    fn add_flux_observer<'py>(&mut self, observer: Bound<'py, PyAny>) {
        if let Ok(observer) = observer.cast::<PyFluxObserverSolar>() {
            let observer = observer.borrow();
            self.viewing_geometry
                .add_flux_observer_solar(observer.cos_sza, observer.observer_altitude_m);
        }
    }
}
