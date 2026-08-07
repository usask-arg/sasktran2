#include <sasktran2/config.h>
#include <sasktran2/validation/validation.h>
#include <cmath>

namespace sasktran2 {
    Config::Config()
        : m_nthreads(1), m_wavelength_batch_size(1), m_nstokes(1),
          m_ndostreams(16), m_enable_wfs(true), m_apply_delta_scaling(false),
          m_los_refraction(false), m_ms_refraction(false),
          m_solar_refraction(false),
          m_wf_precision(WeightingFunctionPrecision::full),
          m_nsinglescatter_moments(16), m_ndosza(1),
          m_ndosphericaliterations(0), m_hr_nincoming(110), m_hr_noutgoing(110),
          m_hr_nspherical_iterations(50), m_hr_num_incoming_points(-1),
          m_successive_orders_max_iterations(50),
          m_successive_orders_relative_tolerance(1.0e-6),
          m_successive_orders_absolute_tolerance(1.0e-12),
          m_successive_orders_anderson_depth(3),
          m_successive_orders_damping(1.0),
          m_successive_orders_initialization(
              SuccessiveOrdersInitialization::none),
          m_do_forced_azimuth(-1), m_do_backprop(false),
          m_singlescatter_phasemode(SingleScatterPhaseMode::from_legendre),
          m_threading_model(ThreadingModel::wavelength),
#ifdef SKTRAN_RUST_SUPPORT
          m_two_stream_backend(TwoStreamBackend::rust),
#else
          m_two_stream_backend(TwoStreamBackend::cpp),
#endif
          m_output_los_optical_depth(false),
          m_input_validation_mode(InputValidationMode::strict),
          m_log_level(spdlog::level::warn) {
        set_multiple_scatter_source(MultipleScatterSource::none);
        set_single_scatter_source(SingleScatterSource::exact);
        set_occultation_source(OccultationSource::none);
        set_emission_source(EmissionSource::none);
        set_stokes_basis(StokesBasis::standard);

        // Set initial log level
        spdlog::set_level(m_log_level);

        m_flux_types.push_back(FluxType::upwelling);
        m_flux_types.push_back(FluxType::downwelling);
    }

    void Config::validate_config() const {
#ifndef SKTRAN_RUST_SUPPORT
        if (m_two_stream_backend == TwoStreamBackend::rust &&
            (m_multiple_scatter_source == MultipleScatterSource::twostream ||
             m_emission_source == EmissionSource::twostream)) {
            spdlog::critical(
                "two_stream_backend=rust requires a build with Rust support");
            sasktran2::validation::throw_configuration_error();
        }
        if (m_multiple_scatter_source ==
            MultipleScatterSource::successive_orders_rust) {
            spdlog::critical(
                "multiple_scatter_source=successive_orders_rust requires a "
                "build with Rust support");
            sasktran2::validation::throw_configuration_error();
        }
#endif
        if (!m_successive_orders_altitude_grid_m.empty() &&
            m_multiple_scatter_source !=
                MultipleScatterSource::successive_orders_rust) {
            spdlog::critical(
                "successive_orders_altitude_grid_m is only supported by "
                "multiple_scatter_source=successive_orders_rust");
            sasktran2::validation::throw_configuration_error();
        }
        if (input_validation_mode() == InputValidationMode::disabled) {
            return;
        }
        // Check that the number of stokes is valid
        if (m_nstokes != 1 && m_nstokes != 3) {
            spdlog::critical("Invalid number of stokes: {}, must be 1 or 3",
                             m_nstokes);

            sasktran2::validation::throw_configuration_error();
        }

        // Check that the number of threads is valid
        if (m_nthreads < 1) {
            spdlog::critical(
                "Invalid number of threads: {}, must be at least 1",
                m_nthreads);

            sasktran2::validation::throw_configuration_error();
        }

        if (m_wavelength_batch_size < 1) {
            spdlog::critical("Invalid wavelength batch size: {}, must be at "
                             "least 1",
                             m_wavelength_batch_size);
            sasktran2::validation::throw_configuration_error();
        }

        // Check that the number of streams is valid
        if (m_ndostreams < 2) {
            spdlog::critical(
                "Invalid number of streams: {}, must be at least 2",
                m_ndostreams);

            sasktran2::validation::throw_configuration_error();
        }

#ifndef SKTRAN_RUST_SUPPORT
        if (m_ndostreams > 40) {
            spdlog::critical(
                "Invalid number of streams: {}, must be less than 40. Compile "
                "with rust support enabled to use more streams.",
                m_ndostreams);

            sasktran2::validation::throw_configuration_error();
        }
#endif

        if (m_ndostreams % 2 != 0) {
            spdlog::critical("Invalid number of streams: {}, must be even",
                             m_ndostreams);

            sasktran2::validation::throw_configuration_error();
        }

        if (m_hr_num_incoming_points != -1 && m_hr_num_incoming_points < 2) {
            spdlog::critical(
                "num_successive_order_points must be -1 or at least 2");
            sasktran2::validation::throw_configuration_error();
        }

        if (m_multiple_scatter_source ==
            MultipleScatterSource::successive_orders_rust) {
            if (m_successive_orders_max_iterations < 1) {
                spdlog::critical(
                    "successive_orders_max_iterations must be at least 1");
                sasktran2::validation::throw_configuration_error();
            }
            if (!std::isfinite(m_successive_orders_relative_tolerance) ||
                m_successive_orders_relative_tolerance < 0 ||
                !std::isfinite(m_successive_orders_absolute_tolerance) ||
                m_successive_orders_absolute_tolerance < 0) {
                spdlog::critical(
                    "successive-orders tolerances must be finite and "
                    "non-negative");
                sasktran2::validation::throw_configuration_error();
            }
            if (m_successive_orders_anderson_depth < 0) {
                spdlog::critical(
                    "successive_orders_anderson_depth must be non-negative");
                sasktran2::validation::throw_configuration_error();
            }
            if (!std::isfinite(m_successive_orders_damping) ||
                m_successive_orders_damping <= 0 ||
                m_successive_orders_damping > 1) {
                spdlog::critical(
                    "successive_orders_damping must be in the interval (0, "
                    "1]");
                sasktran2::validation::throw_configuration_error();
            }
            if (m_successive_orders_initialization ==
                SuccessiveOrdersInitialization::discrete_ordinates) {
                spdlog::critical(
                    "successive_orders_initialization=discrete_ordinates is "
                    "not supported by successive_orders_rust");
                sasktran2::validation::throw_configuration_error();
            }
        } else if (m_successive_orders_initialization ==
                   SuccessiveOrdersInitialization::twostream) {
            spdlog::critical(
                "successive_orders_initialization=twostream requires "
                "multiple_scatter_source=successive_orders_rust");
            sasktran2::validation::throw_configuration_error();
        }

        for (std::size_t altitude_index = 0;
             altitude_index < m_successive_orders_altitude_grid_m.size();
             ++altitude_index) {
            if (!std::isfinite(
                    m_successive_orders_altitude_grid_m[altitude_index]) ||
                (altitude_index > 0 &&
                 m_successive_orders_altitude_grid_m[altitude_index] <=
                     m_successive_orders_altitude_grid_m[altitude_index - 1])) {
                spdlog::critical(
                    "successive_orders_altitude_grid_m must contain finite, "
                    "strictly increasing altitudes");
                sasktran2::validation::throw_configuration_error();
            }
        }

        // Check that the number of single scatter moments is valid
        if (m_nsinglescatter_moments < 1) {
            spdlog::critical("Invalid number of single scatter moments: {}, "
                             "must be at least 1",
                             m_nsinglescatter_moments);

            sasktran2::validation::throw_configuration_error();
        }

        if (m_nsinglescatter_moments < m_ndostreams &&
            m_multiple_scatter_source != MultipleScatterSource::none &&
            (m_single_scatter_source !=
             SingleScatterSource::discrete_ordinates)) {
            spdlog::critical("Invalid number of single scatter moments: {}, "
                             "must be at least the number of streams, {}",
                             m_nsinglescatter_moments, m_ndostreams);

            sasktran2::validation::throw_configuration_error();
        }

        // Very num sza
        if (m_ndosza < 1) {
            spdlog::critical("Invalid number of dosza: {}, must be at least 1",
                             m_ndosza);

            sasktran2::validation::throw_configuration_error();
        }

        // Special twostream options
        if (m_multiple_scatter_source == MultipleScatterSource::twostream) {
            if (m_ndostreams != 2) {
                spdlog::critical("Invalid number of streams: {}, must be 2 for "
                                 "twostream multiple scatter source",
                                 m_ndostreams);
                sasktran2::validation::throw_configuration_error();
            }
        }

        // Validation for discrete_ordinates emission source
        if (m_emission_source == EmissionSource::discrete_ordinates) {
            if (m_single_scatter_source !=
                SingleScatterSource::discrete_ordinates) {
                spdlog::critical("emission_source=discrete_ordinates requires "
                                 "single_scatter_source=discrete_ordinates");
                sasktran2::validation::throw_configuration_error();
            }
            if (m_multiple_scatter_source !=
                MultipleScatterSource::discrete_ordinates) {
                spdlog::critical("emission_source=discrete_ordinates requires "
                                 "multiple_scatter_source=discrete_ordinates");
                sasktran2::validation::throw_configuration_error();
            }
        }
    }

    void Config::validate_config_geometry(int geotype) const {
        if (input_validation_mode() == InputValidationMode::disabled) {
            return;
        }

        // emission_source=discrete_ordinates only works with plane parallel or
        // spherical geometry (pseudospherical is treated as plane parallel for
        // DO purposes)
        if (m_emission_source == EmissionSource::discrete_ordinates) {
            // geotype values: 0=planeparallel, 1=pseudospherical, 2=spherical,
            // 3=ellipsoidal
            if (geotype != 0 && geotype != 1 && geotype != 2) {
                spdlog::critical(
                    "emission_source=discrete_ordinates requires plane "
                    "parallel, pseudospherical, or spherical geometry");
                sasktran2::validation::throw_configuration_error();
            }
        }
    }

} // namespace sasktran2
