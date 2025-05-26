use super::prelude::*;
use crate::prelude::*;
use sasktran2_sys::ffi;

#[repr(i32)]
#[derive(Debug, Clone, Copy, PartialEq)]
pub enum MultipleScatterSource {
    DiscreteOrdinates = 0,
    SuccessiveOrders = 1,
    TwoStream = 2,
    None = 3,
}

#[repr(i32)]
#[derive(Debug, Clone, Copy, PartialEq)]
pub enum SingleScatterSource {
    Exact = 0,
    SolarTable = 1,
    DiscreteOrdinates = 2,
    None = 3,
}

#[repr(i32)]
#[derive(Debug, Clone, Copy, PartialEq)]
pub enum OccultationSource {
    None = 1,
    Standard = 0,
}

#[repr(i32)]
#[derive(Debug, Clone, Copy, PartialEq)]
pub enum EmissionSource {
    None = 1,
    Standard = 0,
}

#[repr(i32)]
#[derive(Debug, Clone, Copy, PartialEq)]
pub enum StokesBasis {
    Standard = 0,
    Solar = 1,
    Observer = 2,
}

#[repr(i32)]
#[derive(Debug, Clone, Copy, PartialEq)]
pub enum ThreadingModel {
    Wavelength = 0,
    Source = 1,
}

#[derive(Debug, Clone, Copy, PartialEq)]
pub enum ThreadingLib {
    Rayon,
    OpenMP,
}

#[repr(i32)]
#[derive(Debug, Clone, Copy, PartialEq)]
pub enum InputValidationMode {
    Strict = 0,
    Standard = 1,
    Disabled = 2,
}

/// A wrapper around the c++ Config object, implemented on the c++ side, see
/// cpp/include/sasktran2/config.h
pub struct Config {
    // c++ object
    pub config: *mut ffi::Config,
    // rust specific config settings
    threading_lib: ThreadingLib,
}

impl Default for Config {
    fn default() -> Self {
        Config::new()
    }
}

impl Config {
    pub fn new() -> Self {
        Config {
            config: unsafe { ffi::sk_config_create() },
            threading_lib: ThreadingLib::OpenMP,
        }
    }

    pub fn num_threads(&self) -> Result<usize> {
        let mut num_threads = 0i32;
        let error_code = unsafe { ffi::sk_config_get_num_threads(self.config, &mut num_threads) };

        if error_code != 0 {
            Err(anyhow!(
                "Error getting number of threads: error code {}",
                error_code
            ))
        } else {
            Ok(num_threads as usize)
        }
    }

    pub fn with_num_threads(&mut self, num_threads: usize) -> Result<&mut Self> {
        let error_code = unsafe { ffi::sk_config_set_num_threads(self.config, num_threads as i32) };

        if num_threads < 1 {
            return Err(anyhow!("Number of threads must be greater than 0"));
        }

        threading::set_num_threads(num_threads)?;

        if error_code != 0 {
            Err(anyhow!(
                "Error setting number of threads: error code {}",
                error_code
            ))
        } else {
            Ok(self)
        }
    }

    pub fn threading_model(&self) -> Result<ThreadingModel> {
        let mut threading_model = 0i32;
        let error_code =
            unsafe { ffi::sk_config_get_threading_model(self.config, &mut threading_model) };

        if error_code != 0 {
            Err(anyhow!(
                "Error getting threading model: error code {}",
                error_code
            ))
        } else {
            Ok(unsafe { std::mem::transmute::<i32, ThreadingModel>(threading_model) })
        }
    }

    pub fn with_threading_model(&mut self, threading_model: ThreadingModel) -> Result<&mut Self> {
        let error_code =
            unsafe { ffi::sk_config_set_threading_model(self.config, threading_model as i32) };

        if error_code != 0 {
            Err(anyhow!(
                "Error setting threading model: error code {}",
                error_code
            ))
        } else {
            Ok(self)
        }
    }

    pub fn input_validation_mode(&self) -> Result<InputValidationMode> {
        let mut input_validation_mode = 0i32;
        let error_code = unsafe {
            ffi::sk_config_get_input_validation_mode(self.config, &mut input_validation_mode)
        };

        if error_code != 0 {
            Err(anyhow!(
                "Error getting input validation mode: error code {}",
                error_code
            ))
        } else {
            Ok(unsafe { std::mem::transmute::<i32, InputValidationMode>(input_validation_mode) })
        }
    }

    pub fn with_input_validation_mode(
        &mut self,
        input_validation_mode: InputValidationMode,
    ) -> Result<&mut Self> {
        let error_code = unsafe {
            ffi::sk_config_set_input_validation_mode(self.config, input_validation_mode as i32)
        };

        if error_code != 0 {
            Err(anyhow!(
                "Error setting input validation mode: error code {}",
                error_code
            ))
        } else {
            Ok(self)
        }
    }

    pub fn num_stokes(&self) -> Result<usize> {
        let mut num_stokes = 0i32;
        let error_code = unsafe { ffi::sk_config_get_num_stokes(self.config, &mut num_stokes) };

        if error_code != 0 {
            Err(anyhow!(
                "Error getting number of stokes: error code {}",
                error_code
            ))
        } else {
            Ok(num_stokes as usize)
        }
    }

    pub fn with_num_stokes(&mut self, num_stokes: usize) -> Result<&mut Self> {
        let error_code = unsafe { ffi::sk_config_set_num_stokes(self.config, num_stokes as i32) };

        if error_code != 0 {
            Err(anyhow!(
                "Error setting number of stokes: error code {}",
                error_code
            ))
        } else {
            Ok(self)
        }
    }

    pub fn multiple_scatter_source(&self) -> Result<MultipleScatterSource> {
        let mut multiple_scatter_source = 0i32;
        let error_code = unsafe {
            ffi::sk_config_get_multiple_scatter_source(self.config, &mut multiple_scatter_source)
        };

        if error_code != 0 {
            Err(anyhow!(
                "Error getting multiple scatter source: error code {}",
                error_code
            ))
        } else {
            Ok(unsafe {
                std::mem::transmute::<i32, MultipleScatterSource>(multiple_scatter_source)
            })
        }
    }

    pub fn with_multiple_scatter_source(
        &mut self,
        source: MultipleScatterSource,
    ) -> Result<&mut Self> {
        let error_code =
            unsafe { ffi::sk_config_set_multiple_scatter_source(self.config, source as i32) };

        if error_code != 0 {
            Err(anyhow!(
                "Error setting multiple scatter source: error code {}",
                error_code
            ))
        } else {
            Ok(self)
        }
    }

    pub fn single_scatter_source(&self) -> Result<SingleScatterSource> {
        let mut single_scatter_source = 0i32;
        let error_code = unsafe {
            ffi::sk_config_get_single_scatter_source(self.config, &mut single_scatter_source)
        };

        if error_code != 0 {
            Err(anyhow!(
                "Error getting single scatter source: error code {}",
                error_code
            ))
        } else {
            Ok(unsafe { std::mem::transmute::<i32, SingleScatterSource>(single_scatter_source) })
        }
    }

    pub fn with_single_scatter_source(&mut self, source: SingleScatterSource) -> Result<&mut Self> {
        let error_code =
            unsafe { ffi::sk_config_set_single_scatter_source(self.config, source as i32) };

        if error_code != 0 {
            Err(anyhow!(
                "Error setting single scatter source: error code {}",
                error_code
            ))
        } else {
            Ok(self)
        }
    }

    pub fn occultation_source(&self) -> Result<OccultationSource> {
        let mut occultation_source = 0i32;
        let error_code =
            unsafe { ffi::sk_config_get_occultation_source(self.config, &mut occultation_source) };

        if error_code != 0 {
            Err(anyhow!(
                "Error getting occultation source: error code {}",
                error_code
            ))
        } else {
            Ok(unsafe { std::mem::transmute::<i32, OccultationSource>(occultation_source) })
        }
    }

    pub fn with_occultation_source(&mut self, source: OccultationSource) -> Result<&mut Self> {
        let error_code =
            unsafe { ffi::sk_config_set_occultation_source(self.config, source as i32) };

        if error_code != 0 {
            Err(anyhow!(
                "Error setting occultation source: error code {}",
                error_code
            ))
        } else {
            Ok(self)
        }
    }

    pub fn emission_source(&self) -> Result<EmissionSource> {
        let mut emission_source = 0i32;
        let error_code =
            unsafe { ffi::sk_config_get_emission_source(self.config, &mut emission_source) };

        if error_code != 0 {
            Err(anyhow!(
                "Error getting emission source: error code {}",
                error_code
            ))
        } else {
            Ok(unsafe { std::mem::transmute::<i32, EmissionSource>(emission_source) })
        }
    }

    pub fn with_emission_source(&mut self, source: EmissionSource) -> Result<&mut Self> {
        let error_code = unsafe { ffi::sk_config_set_emission_source(self.config, source as i32) };

        if error_code != 0 {
            Err(anyhow!(
                "Error setting emission source: error code {}",
                error_code
            ))
        } else {
            Ok(self)
        }
    }
    pub fn stokes_basis(&self) -> Result<StokesBasis> {
        let mut stokes_basis = 0i32;
        let error_code = unsafe { ffi::sk_config_get_stokes_basis(self.config, &mut stokes_basis) };

        if error_code != 0 {
            Err(anyhow!(
                "Error getting stokes basis: error code {}",
                error_code
            ))
        } else {
            Ok(unsafe { std::mem::transmute::<i32, StokesBasis>(stokes_basis) })
        }
    }

    pub fn with_stokes_basis(&mut self, basis: StokesBasis) -> Result<&mut Self> {
        let error_code = unsafe { ffi::sk_config_set_stokes_basis(self.config, basis as i32) };

        if error_code != 0 {
            Err(anyhow!(
                "Error setting stokes basis: error code {}",
                error_code
            ))
        } else {
            Ok(self)
        }
    }

    pub fn delta_m_scaling(&self) -> Result<bool> {
        let mut delta_m_scaling = 0i32;
        let error_code =
            unsafe { ffi::sk_config_get_apply_delta_scaling(self.config, &mut delta_m_scaling) };

        if error_code != 0 {
            Err(anyhow!(
                "Error getting delta m scaling: error code {}",
                error_code
            ))
        } else {
            Ok(delta_m_scaling != 0)
        }
    }

    pub fn with_delta_m_scaling(&mut self, delta_m_scaling: bool) -> Result<&mut Self> {
        let error_code = unsafe {
            ffi::sk_config_set_apply_delta_scaling(self.config, if delta_m_scaling { 1 } else { 0 })
        };

        if error_code != 0 {
            Err(anyhow!(
                "Error setting delta m scaling: error code {}",
                error_code
            ))
        } else {
            Ok(self)
        }
    }

    pub fn los_refraction(&self) -> Result<bool> {
        let mut los_refraction = 0i32;
        let error_code =
            unsafe { ffi::sk_config_get_los_refraction(self.config, &mut los_refraction) };

        if error_code != 0 {
            Err(anyhow!(
                "Error getting los refraction: error code {}",
                error_code
            ))
        } else {
            Ok(los_refraction != 0)
        }
    }

    pub fn with_los_refraction(&mut self, los_refraction: bool) -> Result<&mut Self> {
        let error_code = unsafe {
            ffi::sk_config_set_los_refraction(self.config, if los_refraction { 1 } else { 0 })
        };

        if error_code != 0 {
            Err(anyhow!(
                "Error setting los refraction: error code {}",
                error_code
            ))
        } else {
            Ok(self)
        }
    }

    pub fn output_los_optical_depth(&self) -> Result<bool> {
        let mut output_los_optical_depth = 0i32;
        let error_code = unsafe {
            ffi::sk_config_get_output_los_optical_depth(self.config, &mut output_los_optical_depth)
        };

        if error_code != 0 {
            Err(anyhow!(
                "Error getting output los optical depth: error code {}",
                error_code
            ))
        } else {
            Ok(output_los_optical_depth != 0)
        }
    }

    pub fn with_output_los_optical_depth(
        &mut self,
        output_los_optical_depth: bool,
    ) -> Result<&mut Self> {
        let error_code = unsafe {
            ffi::sk_config_set_output_los_optical_depth(
                self.config,
                if output_los_optical_depth { 1 } else { 0 },
            )
        };

        if error_code != 0 {
            Err(anyhow!(
                "Error setting output los optical depth: error code {}",
                error_code
            ))
        } else {
            Ok(self)
        }
    }

    pub fn solar_refraction(&self) -> Result<bool> {
        let mut solar_refraction = 0i32;
        let error_code =
            unsafe { ffi::sk_config_get_solar_refraction(self.config, &mut solar_refraction) };

        if error_code != 0 {
            Err(anyhow!(
                "Error getting solar refraction: error code {}",
                error_code
            ))
        } else {
            Ok(solar_refraction != 0)
        }
    }

    pub fn with_solar_refraction(&mut self, solar_refraction: bool) -> Result<&mut Self> {
        let error_code = unsafe {
            ffi::sk_config_set_solar_refraction(self.config, if solar_refraction { 1 } else { 0 })
        };

        if error_code != 0 {
            Err(anyhow!(
                "Error setting solar refraction: error code {}",
                error_code
            ))
        } else {
            Ok(self)
        }
    }

    pub fn multiple_scatter_refraction(&self) -> Result<bool> {
        let mut multiple_scatter_refraction = 0i32;
        let error_code = unsafe {
            ffi::sk_config_get_multiple_scatter_refraction(
                self.config,
                &mut multiple_scatter_refraction,
            )
        };

        if error_code != 0 {
            Err(anyhow!(
                "Error getting multiple scatter refraction: error code {}",
                error_code
            ))
        } else {
            Ok(multiple_scatter_refraction != 0)
        }
    }

    pub fn with_multiple_scatter_refraction(
        &mut self,
        multiple_scatter_refraction: bool,
    ) -> Result<&mut Self> {
        let error_code = unsafe {
            ffi::sk_config_set_multiple_scatter_refraction(
                self.config,
                if multiple_scatter_refraction { 1 } else { 0 },
            )
        };

        if error_code != 0 {
            Err(anyhow!(
                "Error setting multiple scatter refraction: error code {}",
                error_code
            ))
        } else {
            Ok(self)
        }
    }

    pub fn num_sza(&self) -> Result<usize> {
        let mut num_sza = 0i32;
        let error_code = unsafe { ffi::sk_config_get_num_do_sza(self.config, &mut num_sza) };

        if error_code != 0 {
            Err(anyhow!(
                "Error getting number of sza: error code {}",
                error_code
            ))
        } else {
            Ok(num_sza as usize)
        }
    }

    pub fn with_num_sza(&mut self, num_sza: usize) -> Result<&mut Self> {
        let error_code = unsafe { ffi::sk_config_set_num_do_sza(self.config, num_sza as i32) };

        if error_code != 0 {
            Err(anyhow!(
                "Error setting number of sza: error code {}",
                error_code
            ))
        } else {
            Ok(self)
        }
    }

    pub fn num_successive_orders_iterations(&self) -> Result<usize> {
        let mut num_iterations = 0i32;
        let error_code = unsafe {
            ffi::sk_config_get_num_hr_spherical_iterations(self.config, &mut num_iterations)
        };

        if error_code != 0 {
            Err(anyhow!(
                "Error getting number of successive orders iterations: error code {}",
                error_code
            ))
        } else {
            Ok(num_iterations as usize)
        }
    }

    pub fn with_num_successive_orders_iterations(
        &mut self,
        num_iterations: usize,
    ) -> Result<&mut Self> {
        let error_code = unsafe {
            ffi::sk_config_set_num_hr_spherical_iterations(self.config, num_iterations as i32)
        };

        if error_code != 0 {
            Err(anyhow!(
                "Error setting number of successive orders iterations: error code {}",
                error_code
            ))
        } else {
            Ok(self)
        }
    }

    pub fn init_successive_orders_with_discrete_ordinates(&self) -> Result<bool> {
        let mut init = 0i32;
        let error_code =
            unsafe { ffi::sk_config_get_initialize_hr_with_do(self.config, &mut init) };

        if error_code != 0 {
            Err(anyhow!(
                "Error getting init successive orders with discrete ordinates: error code {}",
                error_code
            ))
        } else {
            Ok(init != 0)
        }
    }

    pub fn with_init_successive_orders_with_discrete_ordinates(
        &mut self,
        init: bool,
    ) -> Result<&mut Self> {
        let error_code = unsafe {
            ffi::sk_config_set_initialize_hr_with_do(self.config, if init { 1 } else { 0 })
        };

        if error_code != 0 {
            Err(anyhow!(
                "Error setting init successive orders with discrete ordinates: error code {}",
                error_code
            ))
        } else {
            Ok(self)
        }
    }

    pub fn num_streams(&self) -> Result<usize> {
        let mut num_streams = 0i32;
        let error_code = unsafe { ffi::sk_config_get_num_streams(self.config, &mut num_streams) };

        if error_code != 0 {
            Err(anyhow!(
                "Error getting number of streams: error code {}",
                error_code
            ))
        } else {
            Ok(num_streams as usize)
        }
    }

    pub fn with_num_streams(&mut self, num_streams: usize) -> Result<&mut Self> {
        let error_code = unsafe { ffi::sk_config_set_num_streams(self.config, num_streams as i32) };

        if error_code != 0 {
            Err(anyhow!(
                "Error setting number of streams: error code {}",
                error_code
            ))
        } else {
            Ok(self)
        }
    }

    pub fn num_forced_azimuth(&self) -> Result<usize> {
        let mut num_azimuth = 0i32;
        let error_code =
            unsafe { ffi::sk_config_get_num_do_forced_azimuth(self.config, &mut num_azimuth) };

        if error_code != 0 {
            Err(anyhow!(
                "Error getting number of forced azimuth: error code {}",
                error_code
            ))
        } else {
            Ok(num_azimuth as usize)
        }
    }

    pub fn with_num_forced_azimuth(&mut self, num_azimuth: usize) -> Result<&mut Self> {
        let error_code =
            unsafe { ffi::sk_config_set_num_do_forced_azimuth(self.config, num_azimuth as i32) };

        if error_code != 0 {
            Err(anyhow!(
                "Error setting number of forced azimuth: error code {}",
                error_code
            ))
        } else {
            Ok(self)
        }
    }

    pub fn do_backprop(&self) -> Result<bool> {
        let mut do_backprop = 0i32;
        let error_code = unsafe { ffi::sk_config_get_do_backprop(self.config, &mut do_backprop) };

        if error_code != 0 {
            Err(anyhow!(
                "Error getting do backprop: error code {}",
                error_code
            ))
        } else {
            Ok(do_backprop != 0)
        }
    }

    pub fn with_do_backprop(&mut self, do_backprop: bool) -> Result<&mut Self> {
        let error_code =
            unsafe { ffi::sk_config_set_do_backprop(self.config, if do_backprop { 1 } else { 0 }) };

        if error_code != 0 {
            Err(anyhow!(
                "Error setting do backprop: error code {}",
                error_code
            ))
        } else {
            Ok(self)
        }
    }

    pub fn num_successive_orders_points(&self) -> Result<usize> {
        let mut num_points = 0i32;
        let error_code =
            unsafe { ffi::sk_config_get_num_hr_full_incoming_points(self.config, &mut num_points) };

        if error_code != 0 {
            Err(anyhow!(
                "Error getting number of successive orders points: error code {}",
                error_code
            ))
        } else {
            Ok(num_points as usize)
        }
    }

    pub fn with_num_successive_orders_points(&mut self, num_points: usize) -> Result<&mut Self> {
        let error_code = unsafe {
            ffi::sk_config_set_num_hr_full_incoming_points(self.config, num_points as i32)
        };

        if error_code != 0 {
            Err(anyhow!(
                "Error setting number of successive orders points: error code {}",
                error_code
            ))
        } else {
            Ok(self)
        }
    }

    pub fn num_singlescatter_moments(&self) -> Result<usize> {
        let mut num_moments = 0i32;
        let error_code =
            unsafe { ffi::sk_config_get_num_singlescatter_moments(self.config, &mut num_moments) };

        if error_code != 0 {
            Err(anyhow!(
                "Error getting number of single scatter moments: error code {}",
                error_code
            ))
        } else {
            Ok(num_moments as usize)
        }
    }

    pub fn with_num_singlescatter_moments(&mut self, num_moments: usize) -> Result<&mut Self> {
        let error_code = unsafe {
            ffi::sk_config_set_num_singlescatter_moments(self.config, num_moments as i32)
        };

        if error_code != 0 {
            Err(anyhow!(
                "Error setting number of single scatter moments: error code {}",
                error_code
            ))
        } else {
            Ok(self)
        }
    }

    pub fn num_successive_orders_incoming(&self) -> Result<usize> {
        let mut num_incoming = 0i32;
        let error_code =
            unsafe { ffi::sk_config_get_num_hr_incoming(self.config, &mut num_incoming) };

        if error_code != 0 {
            Err(anyhow!(
                "Error getting number of successive orders incoming: error code {}",
                error_code
            ))
        } else {
            Ok(num_incoming as usize)
        }
    }

    pub fn with_num_successive_orders_incoming(
        &mut self,
        num_incoming: usize,
    ) -> Result<&mut Self> {
        let error_code =
            unsafe { ffi::sk_config_set_num_hr_incoming(self.config, num_incoming as i32) };

        if error_code != 0 {
            Err(anyhow!(
                "Error setting number of successive orders incoming: error code {}",
                error_code
            ))
        } else {
            Ok(self)
        }
    }

    pub fn num_successive_orders_outgoing(&self) -> Result<usize> {
        let mut num_outgoing = 0i32;
        let error_code =
            unsafe { ffi::sk_config_get_num_hr_outgoing(self.config, &mut num_outgoing) };

        if error_code != 0 {
            Err(anyhow!(
                "Error getting number of successive orders outgoing: error code {}",
                error_code
            ))
        } else {
            Ok(num_outgoing as usize)
        }
    }

    pub fn with_num_successive_orders_outgoing(
        &mut self,
        num_outgoing: usize,
    ) -> Result<&mut Self> {
        let error_code =
            unsafe { ffi::sk_config_set_num_hr_outgoing(self.config, num_outgoing as i32) };

        if error_code != 0 {
            Err(anyhow!(
                "Error setting number of successive orders outgoing: error code {}",
                error_code
            ))
        } else {
            Ok(self)
        }
    }

    pub fn with_threading_lib(&mut self, threading_lib: ThreadingLib) -> Result<&mut Self> {
        self.threading_lib = threading_lib;
        Ok(self)
    }

    pub fn threading_lib(&self) -> ThreadingLib {
        self.threading_lib
    }
}

impl Drop for Config {
    fn drop(&mut self) {
        unsafe {
            ffi::sk_config_destroy(self.config);
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_config() {
        let _config = Config::new();
    }

    #[test]
    fn test_config_defaults() {
        let config = Config::new();

        // Test default values
        assert_eq!(config.num_threads().unwrap(), 1);
        assert_eq!(
            config.threading_model().unwrap(),
            ThreadingModel::Wavelength
        );
        assert_eq!(
            config.input_validation_mode().unwrap(),
            InputValidationMode::Strict
        );
        assert_eq!(config.num_stokes().unwrap(), 1);
        assert_eq!(
            config.multiple_scatter_source().unwrap(),
            MultipleScatterSource::None
        );
        assert_eq!(
            config.single_scatter_source().unwrap(),
            SingleScatterSource::Exact
        );
        assert_eq!(config.stokes_basis().unwrap(), StokesBasis::Standard);
    }

    #[test]
    fn test_config_setters() {
        let mut config = Config::new();

        // Test setters and getters
        config.with_num_threads(4).unwrap();
        assert_eq!(config.num_threads().unwrap(), 4);

        config.with_threading_model(ThreadingModel::Source).unwrap();
        assert_eq!(config.threading_model().unwrap(), ThreadingModel::Source);

        config
            .with_input_validation_mode(InputValidationMode::Strict)
            .unwrap();
        assert_eq!(
            config.input_validation_mode().unwrap(),
            InputValidationMode::Strict
        );

        config.with_num_stokes(3).unwrap();
        assert_eq!(config.num_stokes().unwrap(), 3);

        config.with_threading_lib(ThreadingLib::Rayon).unwrap();
        assert_eq!(config.threading_lib(), ThreadingLib::Rayon);

        config
            .with_multiple_scatter_source(MultipleScatterSource::SuccessiveOrders)
            .unwrap();
        assert_eq!(
            config.multiple_scatter_source().unwrap(),
            MultipleScatterSource::SuccessiveOrders
        );

        config
            .with_single_scatter_source(SingleScatterSource::SolarTable)
            .unwrap();
        assert_eq!(
            config.single_scatter_source().unwrap(),
            SingleScatterSource::SolarTable
        );

        config.with_stokes_basis(StokesBasis::Observer).unwrap();
        assert_eq!(config.stokes_basis().unwrap(), StokesBasis::Observer);
    }

    #[test]
    fn test_boolean_options() {
        let mut config = Config::new();

        config.with_delta_m_scaling(true).unwrap();
        assert!(config.delta_m_scaling().unwrap());

        config.with_delta_m_scaling(false).unwrap();
        assert!(!config.delta_m_scaling().unwrap());

        config.with_los_refraction(true).unwrap();
        assert!(config.los_refraction().unwrap());

        config.with_output_los_optical_depth(true).unwrap();
        assert!(config.output_los_optical_depth().unwrap());

        config.with_solar_refraction(true).unwrap();
        assert!(config.solar_refraction().unwrap());

        config.with_multiple_scatter_refraction(true).unwrap();
        assert!(config.multiple_scatter_refraction().unwrap());

        config.with_do_backprop(true).unwrap();
        assert!(config.do_backprop().unwrap());
    }

    #[test]
    fn test_numerical_parameters() {
        let mut config = Config::new();

        config.with_num_sza(10).unwrap();
        assert_eq!(config.num_sza().unwrap(), 10);

        config.with_num_streams(8).unwrap();
        assert_eq!(config.num_streams().unwrap(), 8);

        config.with_num_forced_azimuth(6).unwrap();
        assert_eq!(config.num_forced_azimuth().unwrap(), 6);

        config.with_num_successive_orders_iterations(5).unwrap();
        assert_eq!(config.num_successive_orders_iterations().unwrap(), 5);

        config.with_num_successive_orders_points(100).unwrap();
        assert_eq!(config.num_successive_orders_points().unwrap(), 100);

        config.with_num_singlescatter_moments(16).unwrap();
        assert_eq!(config.num_singlescatter_moments().unwrap(), 16);

        config.with_num_successive_orders_incoming(20).unwrap();
        assert_eq!(config.num_successive_orders_incoming().unwrap(), 20);

        config.with_num_successive_orders_outgoing(30).unwrap();
        assert_eq!(config.num_successive_orders_outgoing().unwrap(), 30);
    }

    #[test]
    fn test_init_successive_orders() {
        let mut config = Config::new();

        config
            .with_init_successive_orders_with_discrete_ordinates(true)
            .unwrap();
        assert!(
            config
                .init_successive_orders_with_discrete_ordinates()
                .unwrap()
        );

        config
            .with_init_successive_orders_with_discrete_ordinates(false)
            .unwrap();
        assert!(
            !config
                .init_successive_orders_with_discrete_ordinates()
                .unwrap()
        );
    }
}
