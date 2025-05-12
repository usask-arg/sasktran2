use crate::prelude::*;
use pyo3::prelude::*;
use sasktran2_rs::bindings::config;

#[pyclass(eq, eq_int)]
#[derive(PartialEq, Clone)]
pub enum MultipleScatterSource {
    DiscreteOrdinates,
    SuccessiveOrders,
    TwoStream,
    NoSource,
}

#[pyclass(eq, eq_int)]
#[derive(PartialEq, Clone)]
pub enum SingleScatterSource {
    NoSource,
    Exact,
    Table,
    DiscreteOrdinates,
}

#[pyclass(eq, eq_int)]
#[derive(PartialEq, Clone)]
pub enum OccultationSource {
    NoSource,
    Standard,
}

#[pyclass(eq, eq_int)]
#[derive(PartialEq, Clone)]
pub enum EmissionSource {
    NoSource,
    Standard,
}

#[pyclass(eq, eq_int)]
#[derive(PartialEq, Clone)]
pub enum StokesBasis {
    Standard,
    Solar,
    Observer,
}

#[pyclass(eq, eq_int)]
#[derive(PartialEq, Clone)]
pub enum ThreadingModel {
    Wavelength,
    Source,
}

#[pyclass(eq, eq_int)]
#[derive(PartialEq, Clone)]
pub enum InputValidationMode {
    Strict,
    Standard,
    Disabled,
}

#[pyclass(eq, eq_int)]
#[derive(PartialEq, Clone)]
pub enum ThreadingLib {
    Rayon,
    OpenMP,
}

#[pyclass(unsendable)]
pub struct PyConfig {
    pub config: config::Config,
}

#[pymethods]
impl PyConfig {
    #[new]
    fn new() -> Self {
        let config = config::Config::new();
        Self { config }
    }

    #[getter]
    fn get_multiple_scatter_source(&self) -> MultipleScatterSource {
        match self.config.multiple_scatter_source().unwrap() {
            config::MultipleScatterSource::DiscreteOrdinates => {
                MultipleScatterSource::DiscreteOrdinates
            }
            config::MultipleScatterSource::SuccessiveOrders => {
                MultipleScatterSource::SuccessiveOrders
            }
            config::MultipleScatterSource::TwoStream => MultipleScatterSource::TwoStream,
            config::MultipleScatterSource::None => MultipleScatterSource::NoSource,
        }
    }

    #[setter]
    fn set_multiple_scatter_source(
        &mut self,
        source: PyRef<'_, MultipleScatterSource>,
    ) -> PyResult<()> {
        let source = match *source {
            MultipleScatterSource::DiscreteOrdinates => {
                config::MultipleScatterSource::DiscreteOrdinates
            }
            MultipleScatterSource::SuccessiveOrders => {
                config::MultipleScatterSource::SuccessiveOrders
            }
            MultipleScatterSource::TwoStream => config::MultipleScatterSource::TwoStream,
            MultipleScatterSource::NoSource => config::MultipleScatterSource::None,
        };
        self.config
            .with_multiple_scatter_source(source)
            .into_pyresult()?;

        Ok(())
    }

    #[getter]
    fn get_num_threads(&self) -> PyResult<usize> {
        self.config.num_threads().into_pyresult()
    }

    #[setter]
    fn set_num_threads(&mut self, num_threads: usize) -> PyResult<()> {
        self.config.with_num_threads(num_threads).into_pyresult()?;

        Ok(())
    }

    #[getter]
    fn threading_model(&self) -> PyResult<ThreadingModel> {
        let model = self.config.threading_model().into_pyresult()?;
        match model {
            config::ThreadingModel::Wavelength => Ok(ThreadingModel::Wavelength),
            config::ThreadingModel::Source => Ok(ThreadingModel::Source),
        }
    }

    #[setter]
    fn set_threading_model(&mut self, model: PyRef<'_, ThreadingModel>) -> PyResult<()> {
        let model = match *model {
            ThreadingModel::Wavelength => config::ThreadingModel::Wavelength,
            ThreadingModel::Source => config::ThreadingModel::Source,
        };
        self.config.with_threading_model(model).into_pyresult()?;

        Ok(())
    }

    #[getter]
    fn get_input_validation_mode(&self) -> PyResult<InputValidationMode> {
        let mode = self.config.input_validation_mode().into_pyresult()?;
        match mode {
            config::InputValidationMode::Strict => Ok(InputValidationMode::Strict),
            config::InputValidationMode::Standard => Ok(InputValidationMode::Standard),
            config::InputValidationMode::Disabled => Ok(InputValidationMode::Disabled),
        }
    }

    #[setter]
    fn set_input_validation_mode(&mut self, mode: PyRef<'_, InputValidationMode>) -> PyResult<()> {
        let mode = match *mode {
            InputValidationMode::Strict => config::InputValidationMode::Strict,
            InputValidationMode::Standard => config::InputValidationMode::Standard,
            InputValidationMode::Disabled => config::InputValidationMode::Disabled,
        };
        self.config
            .with_input_validation_mode(mode)
            .into_pyresult()?;

        Ok(())
    }

    #[getter]
    fn get_num_stokes(&self) -> PyResult<usize> {
        self.config.num_stokes().into_pyresult()
    }

    #[setter]
    fn set_num_stokes(&mut self, num_stokes: usize) -> PyResult<()> {
        self.config.with_num_stokes(num_stokes).into_pyresult()?;

        Ok(())
    }

    #[getter]
    fn get_single_scatter_source(&self) -> PyResult<SingleScatterSource> {
        let source = self.config.single_scatter_source().into_pyresult()?;

        match source {
            config::SingleScatterSource::Exact => Ok(SingleScatterSource::Exact),
            config::SingleScatterSource::SolarTable => Ok(SingleScatterSource::Table),
            config::SingleScatterSource::DiscreteOrdinates => {
                Ok(SingleScatterSource::DiscreteOrdinates)
            }
            config::SingleScatterSource::None => Ok(SingleScatterSource::NoSource),
        }
    }

    #[setter]
    fn set_single_scatter_source(
        &mut self,
        source: PyRef<'_, SingleScatterSource>,
    ) -> PyResult<()> {
        let source = match *source {
            SingleScatterSource::Exact => config::SingleScatterSource::Exact,
            SingleScatterSource::Table => config::SingleScatterSource::SolarTable,
            SingleScatterSource::DiscreteOrdinates => {
                config::SingleScatterSource::DiscreteOrdinates
            }
            SingleScatterSource::NoSource => config::SingleScatterSource::None,
        };
        self.config
            .with_single_scatter_source(source)
            .into_pyresult()?;

        Ok(())
    }

    #[getter]
    fn get_occultation_source(&self) -> PyResult<OccultationSource> {
        let source = self.config.occultation_source().into_pyresult()?;

        match source {
            config::OccultationSource::Standard => Ok(OccultationSource::Standard),
            config::OccultationSource::None => Ok(OccultationSource::NoSource),
        }
    }

    #[setter]
    fn set_occultation_source(&mut self, source: PyRef<'_, OccultationSource>) -> PyResult<()> {
        let source = match *source {
            OccultationSource::Standard => config::OccultationSource::Standard,
            OccultationSource::NoSource => config::OccultationSource::None,
        };
        self.config
            .with_occultation_source(source)
            .into_pyresult()?;

        Ok(())
    }

    #[getter]
    fn get_emission_source(&self) -> PyResult<EmissionSource> {
        let source = self.config.emission_source().into_pyresult()?;

        match source {
            config::EmissionSource::Standard => Ok(EmissionSource::Standard),
            config::EmissionSource::None => Ok(EmissionSource::NoSource),
        }
    }

    #[setter]
    fn set_emission_source(&mut self, source: PyRef<'_, EmissionSource>) -> PyResult<()> {
        let source = match *source {
            EmissionSource::Standard => config::EmissionSource::Standard,
            EmissionSource::NoSource => config::EmissionSource::None,
        };
        self.config.with_emission_source(source).into_pyresult()?;

        Ok(())
    }

    #[getter]
    fn get_stokes_basis(&self) -> PyResult<StokesBasis> {
        let basis = self.config.stokes_basis().into_pyresult()?;

        match basis {
            config::StokesBasis::Standard => Ok(StokesBasis::Standard),
            config::StokesBasis::Solar => Ok(StokesBasis::Solar),
            config::StokesBasis::Observer => Ok(StokesBasis::Observer),
        }
    }

    #[setter]
    fn set_stokes_basis(&mut self, basis: PyRef<'_, StokesBasis>) -> PyResult<()> {
        let basis = match *basis {
            StokesBasis::Standard => config::StokesBasis::Standard,
            StokesBasis::Solar => config::StokesBasis::Solar,
            StokesBasis::Observer => config::StokesBasis::Observer,
        };
        self.config.with_stokes_basis(basis).into_pyresult()?;

        Ok(())
    }

    #[getter]
    fn delta_m_scaling(&self) -> PyResult<bool> {
        self.config.delta_m_scaling().into_pyresult()
    }

    #[setter]
    fn set_delta_m_scaling(&mut self, delta_m_scaling: bool) -> PyResult<()> {
        self.config
            .with_delta_m_scaling(delta_m_scaling)
            .into_pyresult()?;

        Ok(())
    }

    #[getter]
    fn los_refraction(&self) -> PyResult<bool> {
        self.config.los_refraction().into_pyresult()
    }

    #[setter]
    fn set_los_refraction(&mut self, los_refraction: bool) -> PyResult<()> {
        self.config
            .with_los_refraction(los_refraction)
            .into_pyresult()?;

        Ok(())
    }

    #[getter]
    fn get_output_los_optical_depth(&self) -> PyResult<bool> {
        self.config.output_los_optical_depth().into_pyresult()
    }

    #[setter]
    fn set_output_los_optical_depth(&mut self, output_los_optical_depth: bool) -> PyResult<()> {
        self.config
            .with_output_los_optical_depth(output_los_optical_depth)
            .into_pyresult()?;

        Ok(())
    }

    #[getter]
    fn get_solar_refraction(&self) -> PyResult<bool> {
        self.config.solar_refraction().into_pyresult()
    }

    #[setter]
    fn set_solar_refraction(&mut self, solar_refraction: bool) -> PyResult<()> {
        self.config
            .with_solar_refraction(solar_refraction)
            .into_pyresult()?;

        Ok(())
    }

    #[getter]
    fn get_multiple_scatter_refraction(&self) -> PyResult<bool> {
        self.config.multiple_scatter_refraction().into_pyresult()
    }

    #[setter]
    fn set_multiple_scatter_refraction(
        &mut self,
        multiple_scatter_refraction: bool,
    ) -> PyResult<()> {
        self.config
            .with_multiple_scatter_refraction(multiple_scatter_refraction)
            .into_pyresult()?;

        Ok(())
    }

    #[getter]
    fn get_num_sza(&self) -> PyResult<usize> {
        self.config.num_sza().into_pyresult()
    }

    #[setter]
    fn set_num_sza(&mut self, num_sza: usize) -> PyResult<()> {
        self.config.with_num_sza(num_sza).into_pyresult()?;

        Ok(())
    }

    #[getter]
    fn get_num_successive_orders_iterations(&self) -> PyResult<usize> {
        self.config
            .num_successive_orders_iterations()
            .into_pyresult()
    }

    #[setter]
    fn set_num_successive_orders_iterations(&mut self, num_iterations: usize) -> PyResult<()> {
        self.config
            .with_num_successive_orders_iterations(num_iterations)
            .into_pyresult()?;

        Ok(())
    }

    #[getter]
    fn init_successive_orders_with_discrete_ordinates(&self) -> PyResult<bool> {
        self.config
            .init_successive_orders_with_discrete_ordinates()
            .into_pyresult()
    }

    #[setter]
    fn set_init_successive_orders_with_discrete_ordinates(&mut self, init: bool) -> PyResult<()> {
        self.config
            .with_init_successive_orders_with_discrete_ordinates(init)
            .into_pyresult()?;

        Ok(())
    }

    #[getter]
    fn get_num_streams(&self) -> PyResult<usize> {
        self.config.num_streams().into_pyresult()
    }

    #[setter]
    fn set_num_streams(&mut self, num_streams: usize) -> PyResult<()> {
        self.config.with_num_streams(num_streams).into_pyresult()?;

        Ok(())
    }

    #[getter]
    fn get_num_forced_azimuth(&self) -> PyResult<usize> {
        self.config.num_forced_azimuth().into_pyresult()
    }

    #[setter]
    fn set_num_forced_azimuth(&mut self, num_forced_azimuth: usize) -> PyResult<()> {
        self.config
            .with_num_forced_azimuth(num_forced_azimuth)
            .into_pyresult()?;

        Ok(())
    }

    #[getter]
    fn get_do_backprop(&self) -> PyResult<bool> {
        self.config.do_backprop().into_pyresult()
    }

    #[setter]
    fn set_do_backprop(&mut self, do_backprop: bool) -> PyResult<()> {
        self.config.with_do_backprop(do_backprop).into_pyresult()?;

        Ok(())
    }

    #[getter]
    fn get_num_successive_orders_points(&self) -> PyResult<usize> {
        self.config.num_successive_orders_points().into_pyresult()
    }

    #[setter]
    fn set_num_successive_orders_points(&mut self, num_points: usize) -> PyResult<()> {
        self.config
            .with_num_successive_orders_points(num_points)
            .into_pyresult()?;

        Ok(())
    }

    #[getter]
    fn get_num_singlescatter_moments(&self) -> PyResult<usize> {
        self.config.num_singlescatter_moments().into_pyresult()
    }

    #[setter]
    fn set_num_singlescatter_moments(&mut self, num_moments: usize) -> PyResult<()> {
        self.config
            .with_num_singlescatter_moments(num_moments)
            .into_pyresult()?;

        Ok(())
    }

    #[getter]
    fn get_num_successive_orders_incoming(&self) -> PyResult<usize> {
        self.config.num_successive_orders_incoming().into_pyresult()
    }

    #[setter]
    fn set_num_successive_orders_incoming(&mut self, num_incoming: usize) -> PyResult<()> {
        self.config
            .with_num_successive_orders_incoming(num_incoming)
            .into_pyresult()?;

        Ok(())
    }

    #[getter]
    fn get_num_successive_orders_outgoing(&self) -> PyResult<usize> {
        self.config.num_successive_orders_outgoing().into_pyresult()
    }

    #[setter]
    fn set_num_successive_orders_outgoing(&mut self, num_outgoing: usize) -> PyResult<()> {
        self.config
            .with_num_successive_orders_outgoing(num_outgoing)
            .into_pyresult()?;

        Ok(())
    }

    #[getter]
    fn get_threading_lib(&self) -> PyResult<ThreadingLib> {
        let lib = self.config.threading_lib();
        match lib {
            config::ThreadingLib::Rayon => Ok(ThreadingLib::Rayon),
            config::ThreadingLib::OpenMP => Ok(ThreadingLib::OpenMP),
        }
    }

    #[setter]
    fn set_threading_lib(&mut self, lib: PyRef<'_, ThreadingLib>) -> PyResult<()> {
        let lib = match *lib {
            ThreadingLib::Rayon => config::ThreadingLib::Rayon,
            ThreadingLib::OpenMP => config::ThreadingLib::OpenMP,
        };
        self.config.with_threading_lib(lib).into_pyresult()?;

        Ok(())
    }
}
