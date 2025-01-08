#include <sasktran2/config.h>
#include <sasktran2/validation/validation.h>

namespace sasktran2 {
    Config::Config()
        : m_nstokes(1), m_nthreads(1), m_ndostreams(16), m_enable_wfs(true),
          m_apply_delta_scaling(false), m_los_refraction(false),
          m_ms_refraction(false), m_solar_refraction(false),
          m_wf_precision(WeightingFunctionPrecision::full),
          m_nsinglescatter_moments(16), m_ndosza(1),
          m_ndosphericaliterations(0), m_hr_nincoming(110), m_hr_noutgoing(110),
          m_hr_nspherical_iterations(50), m_hr_num_incoming_points(-1),
          m_do_forced_azimuth(-1), m_do_backprop(false),
          m_singlescatter_phasemode(SingleScatterPhaseMode::from_legendre),
          m_threading_model(ThreadingModel::wavelength),
          m_initialize_hr_with_do_solution(false),
          m_input_validation_mode(InputValidationMode::strict) {
        set_multiple_scatter_source(MultipleScatterSource::none);
        set_single_scatter_source(SingleScatterSource::exact);
        set_occultation_source(OccultationSource::none);
        set_emission_source(EmissionSource::none);
        set_stokes_basis(StokesBasis::standard);
    }

    void Config::validate_config() const {
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

        // Check that the number of streams is valid
        if (m_ndostreams < 2) {
            spdlog::critical(
                "Invalid number of streams: {}, must be at least 2",
                m_ndostreams);

            sasktran2::validation::throw_configuration_error();
        }

        if (m_ndostreams > 40) {
            spdlog::critical(
                "Invalid number of streams: {}, must be less than 40",
                m_ndostreams);

            sasktran2::validation::throw_configuration_error();
        }

        if (m_ndostreams % 2 != 0) {
            spdlog::critical("Invalid number of streams: {}, must be even",
                             m_ndostreams);

            sasktran2::validation::throw_configuration_error();
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
    }

} // namespace sasktran2
