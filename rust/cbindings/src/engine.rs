use crate::atmosphere::Atmosphere;
use crate::config::Config;
use crate::ffi;
use crate::geometry::Geometry1D;
use crate::output::Output;
use crate::viewing_geometry::ViewingGeometry;
use crate::prelude::*;

pub struct Engine<'a> {
    pub engine: *mut ffi::Engine,
    pub config: &'a Config,
    pub geometry: &'a Geometry1D,
    pub viewing_geometry: &'a ViewingGeometry,
}
impl<'a> Engine<'a> {
    pub fn new(
        config: &'a Config,
        geometry: &'a Geometry1D,
        viewing_geometry: &'a ViewingGeometry,
    ) -> Self {
        Engine {
            engine: unsafe {
                ffi::sk_engine_create(
                    config.config,
                    geometry.geometry,
                    viewing_geometry.viewing_geometry,
                )
            },
            config: config,
            geometry: geometry,
            viewing_geometry: viewing_geometry,
        }
    }

    pub fn calculate_radiance(&self, atmosphere: &Atmosphere) -> Result<Output> {
        let num_stokes = self.config.num_stokes()?;
        let num_los = self.viewing_geometry.num_rays()?;
        let num_wavel = atmosphere.num_wavel();

        let output = Output::new(num_stokes, num_los, num_wavel);
        unsafe {
            ffi::sk_engine_calculate_radiance(self.engine, atmosphere.atmosphere, output.output);
        }
        Ok(output)
    }
}
impl<'a> Drop for Engine<'a> {
    fn drop(&mut self) {
        unsafe {
            ffi::sk_engine_destroy(self.engine);
        }
    }
}

#[cfg(test)]
mod tests {
    use crate::geometry::GeometryType;
    use crate::geometry::InterpolationMethod;

    use super::*;

    #[test]
    fn test_engine() {
        let config = Config::new();
        let cos_sza = 0.5;
        let saa = 0.0;
        let earth_radius = 6371000.0;
        let grid_values = vec![1.0, 2.0, 3.0];
        let interp_method = InterpolationMethod::Linear;
        let geotype = GeometryType::Spherical;

        let geometry = Geometry1D::new(
            cos_sza,
            saa,
            earth_radius,
            grid_values,
            interp_method,
            geotype,
        );
        let viewing_geometry = ViewingGeometry::new();

        let engine = Engine::new(&config, &geometry, &viewing_geometry);

        assert!(!engine.engine.is_null());
    }
}
