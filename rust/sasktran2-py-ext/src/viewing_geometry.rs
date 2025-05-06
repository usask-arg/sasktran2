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

    fn add_ray<'py>(&mut self, ray: Bound<'py, PyAny>) {
        if let Ok(ray) = ray.downcast::<PyGroundViewingSolar>() {
            let ray = ray.borrow();
            self.viewing_geometry.add_ground_viewing_solar(
                ray.cos_sza,
                ray.relative_azimuth,
                ray.observer_altitude_m,
                ray.cos_viewing_zenith,
            );
        }

        if let Ok(ray) = ray.downcast::<PyTangentAltitudeSolar>() {
            let ray = ray.borrow();
            self.viewing_geometry.add_tangent_altitude_solar(
                ray.tangent_altitude_m,
                ray.relative_azimuth,
                ray.observer_altitude_m,
                ray.cos_sza,
            );
        }

        if let Ok(ray) = ray.downcast::<PySolarAnglesObserverLocation>() {
            let ray = ray.borrow();
            self.viewing_geometry.add_solar_angles_observer_location(
                ray.cos_sza,
                ray.relative_azimuth,
                ray.cos_viewing_zenith,
                ray.observer_altitude_m,
            );
        }
    }
}
