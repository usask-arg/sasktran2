#include "sasktran2/atmosphere/grid_storage.h"
#include "sasktran2/do_source.h"
#include "sasktran2/geometry.h"
#include "sasktran2/raytracing.h"
#include "sasktran2/source_interface.h"
#include "sktran_disco/twostream/meta.h"
#ifdef SKTRAN_RUST_SUPPORT
#include "sktran_disco/twostream/rust_source.h"
#endif
#include <memory>
#include <sasktran2.h>
#include <sasktran2/validation/validation.h>
#ifdef SKTRAN_OPENMP_SUPPORT
#include <omp.h>
#endif

template <int NSTOKES> void Sasktran2<NSTOKES>::initialize() {
    m_config.validate_config();

    if (m_geometry_1d != nullptr) {
        m_geometry_1d->validate();
    } else {
        m_geometry_2d->validate();

        if ((m_config.single_scatter_source() !=
                 sasktran2::Config::SingleScatterSource::none &&
             m_config.single_scatter_source() !=
                 sasktran2::Config::SingleScatterSource::exact) ||
            m_config.multiple_scatter_source() !=
                sasktran2::Config::MultipleScatterSource::none ||
            (m_config.emission_source() !=
                 sasktran2::Config::EmissionSource::none &&
             m_config.emission_source() !=
                 sasktran2::Config::EmissionSource::standard &&
             m_config.emission_source() !=
                 sasktran2::Config::EmissionSource::volume_emission_rate)) {
            throw std::invalid_argument(
                "Geometry2D currently supports exact single scattering, "
                "occultation, standard emission, and volume emission rate "
                "sources, with multiple scattering disabled");
        }
        if (!m_viewing_geometry.flux_observers().empty()) {
            throw std::invalid_argument(
                "Geometry2D does not yet support flux observers");
        }
        if (m_config.los_refraction()) {
            throw std::invalid_argument(
                "Geometry2D engine integration does not yet accept per-ray "
                "refractive-index profiles");
        }
    }

    m_config.validate_config_geometry(
        m_geometry->coordinates().geometry_type());
    construct_raytracer();
    construct_integrator();
    construct_source_terms();
    calculate_geometry();
}

template <int NSTOKES> void Sasktran2<NSTOKES>::construct_raytracer() {
    if (m_geometry_2d != nullptr) {
#ifdef SKTRAN_RUST_SUPPORT
        m_raytracer_2d =
            std::make_unique<sasktran2::raytracing::RustRayTracer2D>(
                *m_geometry_2d);
        return;
#else
        throw std::invalid_argument(
            "Geometry2D engine integration requires Rust support");
#endif
    }

    if (m_geometry->coordinates().geometry_type() ==
        sasktran2::geometrytype::spherical) {
#if defined(SKTRAN_USE_RUST_RAYTRACER) && defined(SKTRAN_RUST_SUPPORT)
        m_raytracer = std::make_unique<sasktran2::raytracing::RustRayTracer>(
            *m_geometry_1d);
#else
        m_raytracer =
            std::make_unique<sasktran2::raytracing::SphericalShellRayTracer>(
                *m_geometry_1d);
#endif
    } else if (m_geometry->coordinates().geometry_type() ==
                   sasktran2::geometrytype::planeparallel ||
               m_geometry->coordinates().geometry_type() ==
                   sasktran2::geometrytype::pseudospherical) {
        m_raytracer =
            std::make_unique<sasktran2::raytracing::PlaneParallelRayTracer>(
                *m_geometry_1d);
    } else {
        spdlog::error("Requested geometry type is not yet supported.");
    }
}

template <int NSTOKES> void Sasktran2<NSTOKES>::construct_integrator() {
    m_source_integrator =
        std::make_unique<sasktran2::SourceIntegrator<NSTOKES>>(true);
}

template <int NSTOKES> void Sasktran2<NSTOKES>::construct_source_terms() {

    if (m_geometry_2d != nullptr) {
        if (m_config.single_scatter_source() ==
            sasktran2::Config::SingleScatterSource::exact) {
#ifdef SKTRAN_RUST_SUPPORT
            m_source_terms.emplace_back(
                std::make_unique<
                    sasktran2::solartransmission::SingleScatterSource<
                        sasktran2::solartransmission::SolarTransmissionExact,
                        NSTOKES>>(*m_geometry_2d, *m_raytracer_2d));
            m_los_source_terms.push_back(m_source_terms.back().get());
#else
            throw std::invalid_argument(
                "Geometry2D exact single scattering requires Rust support");
#endif
        }
        if (m_config.occultation_source() ==
            sasktran2::Config::OccultationSource::standard) {
            m_source_terms.emplace_back(
                std::make_unique<sasktran2::solartransmission::
                                     OccultationSource<NSTOKES>>());
            m_los_source_terms.push_back(m_source_terms.back().get());
        }
        if (m_config.emission_source() ==
            sasktran2::Config::EmissionSource::standard) {
            m_source_terms.emplace_back(
                std::make_unique<sasktran2::emission::EmissionSource<
                    NSTOKES, sasktran2::Config::EmissionSource::standard>>());
            m_los_source_terms.push_back(m_source_terms.back().get());
        }
        if (m_config.emission_source() ==
            sasktran2::Config::EmissionSource::volume_emission_rate) {
            m_source_terms.emplace_back(
                std::make_unique<sasktran2::emission::EmissionSource<
                    NSTOKES, sasktran2::Config::EmissionSource::
                                 volume_emission_rate>>());
            m_los_source_terms.push_back(m_source_terms.back().get());
        }
        for (auto& source : m_source_terms) {
            source->initialize_config(m_config);
        }
        return;
    }

    if (m_config.single_scatter_source() ==
        sasktran2::Config::SingleScatterSource::exact) {
        m_source_terms.emplace_back(
            std::make_unique<sasktran2::solartransmission::SingleScatterSource<
                sasktran2::solartransmission::SolarTransmissionExact, NSTOKES>>(
                *m_geometry_1d, *m_raytracer));

        m_los_source_terms.push_back(
            m_source_terms[m_source_terms.size() - 1].get());
    }

    if (m_config.single_scatter_source() ==
        sasktran2::Config::SingleScatterSource::solartable) {
        m_source_terms.emplace_back(
            std::make_unique<sasktran2::solartransmission::SingleScatterSource<
                sasktran2::solartransmission::SolarTransmissionTable, NSTOKES>>(
                *m_geometry_1d, *m_raytracer));

        m_los_source_terms.push_back(
            m_source_terms[m_source_terms.size() - 1].get());
    }

    if (m_config.occultation_source() ==
        sasktran2::Config::OccultationSource::standard) {
        m_source_terms.emplace_back(
            std::make_unique<
                sasktran2::solartransmission::OccultationSource<NSTOKES>>());

        m_los_source_terms.push_back(
            m_source_terms[m_source_terms.size() - 1].get());
    }

    if (m_config.emission_source() ==
        sasktran2::Config::EmissionSource::standard) {
        m_source_terms.emplace_back(
            std::make_unique<sasktran2::emission::EmissionSource<
                NSTOKES, sasktran2::Config::EmissionSource::standard>>());

        m_los_source_terms.push_back(
            m_source_terms[m_source_terms.size() - 1].get());
    }

    if (m_config.emission_source() ==
        sasktran2::Config::EmissionSource::volume_emission_rate) {
        m_source_terms.emplace_back(
            std::make_unique<sasktran2::emission::EmissionSource<
                NSTOKES,
                sasktran2::Config::EmissionSource::volume_emission_rate>>());

        m_los_source_terms.push_back(
            m_source_terms[m_source_terms.size() - 1].get());
    }

    if (m_config.emission_source() ==
        sasktran2::Config::EmissionSource::twostream) {
        if constexpr (NSTOKES == 1) {
#ifdef SKTRAN_RUST_SUPPORT
            if (m_config.two_stream_backend() ==
                sasktran2::Config::TwoStreamBackend::rust) {
                m_source_terms.emplace_back(
                    std::make_unique<RustTwoStreamSourceAdapter<
                        sasktran2::twostream::SourceType::ONLY_THERMAL>>(
                        *m_geometry_1d));
            } else
#endif
            {
                m_source_terms.emplace_back(
                    std::make_unique<TwoStreamSource<
                        1, sasktran2::twostream::SourceType::ONLY_THERMAL>>(
                        *m_geometry_1d));
            }

            m_los_source_terms.push_back(
                m_source_terms[m_source_terms.size() - 1].get());
        } else {
            spdlog::error(
                "TwoStreamSource is only implemented for NSTOKES = 1");
        }
    }

    if (m_config.multiple_scatter_source() ==
        sasktran2::Config::MultipleScatterSource::discrete_ordinates) {

#ifdef SASKTRAN_DISCO_FULL_COMPILE
        if constexpr (NSTOKES == 1) {
            if (m_config.num_do_streams() == 2) {
                if (m_geometry->coordinates().geometry_type() ==
                    sasktran2::geometrytype::spherical) {
                    m_source_terms.emplace_back(
                        std::make_unique<
                            sasktran2::DOSourceInterpolatedPostProcessing<
                                NSTOKES, 2>>(*m_geometry_1d, *m_raytracer));
                } else if (m_config.emission_source() ==
                           sasktran2::Config::EmissionSource::
                               discrete_ordinates) {
                    // The hand-specialized CNSTR=2 plane-parallel
                    // postprocessor only implements the solar particular
                    // solution. Use the generic implementation when thermal
                    // emission is enabled so the atmospheric thermal source is
                    // integrated as well.
                    m_source_terms.emplace_back(
                        std::make_unique<
                            sasktran2::DOSourcePlaneParallelPostProcessing<
                                NSTOKES, -1>>(*m_geometry_1d));
                } else {
                    m_source_terms.emplace_back(
                        std::make_unique<
                            sasktran2::DOSourcePlaneParallelPostProcessing<
                                NSTOKES, 2>>(*m_geometry_1d));
                }
            } else if (m_config.num_do_streams() == 4) {
                if (m_geometry->coordinates().geometry_type() ==
                    sasktran2::geometrytype::spherical) {
                    m_source_terms.emplace_back(
                        std::make_unique<
                            sasktran2::DOSourceInterpolatedPostProcessing<
                                NSTOKES, 4>>(*m_geometry_1d, *m_raytracer));
                } else {
                    m_source_terms.emplace_back(
                        std::make_unique<
                            sasktran2::DOSourcePlaneParallelPostProcessing<
                                NSTOKES, 4>>(*m_geometry_1d));
                }
            } else {
                if (m_geometry->coordinates().geometry_type() ==
                    sasktran2::geometrytype::spherical) {
                    m_source_terms.emplace_back(
                        std::make_unique<
                            sasktran2::DOSourceInterpolatedPostProcessing<
                                NSTOKES, -1>>(*m_geometry_1d, *m_raytracer));
                } else {
                    m_source_terms.emplace_back(
                        std::make_unique<
                            sasktran2::DOSourcePlaneParallelPostProcessing<
                                NSTOKES, -1>>(*m_geometry_1d));
                }
            }
        } else {
            if (m_geometry->coordinates().geometry_type() ==
                sasktran2::geometrytype::spherical) {
                m_source_terms.emplace_back(
                    std::make_unique<
                        sasktran2::DOSourceInterpolatedPostProcessing<NSTOKES,
                                                                      -1>>(
                        *m_geometry_1d, *m_raytracer));
            } else {
                m_source_terms.emplace_back(
                    std::make_unique<
                        sasktran2::DOSourcePlaneParallelPostProcessing<NSTOKES,
                                                                       -1>>(
                        *m_geometry_1d));
            }
        }
#else
        if (m_geometry->coordinates().geometry_type() ==
            sasktran2::geometrytype::spherical) {
            m_source_terms.emplace_back(
                std::make_unique<
                    sasktran2::DOSourceInterpolatedPostProcessing<NSTOKES, -1>>(
                    *m_geometry_1d, *m_raytracer));
        } else {
            m_source_terms.emplace_back(
                std::make_unique<sasktran2::DOSourcePlaneParallelPostProcessing<
                    NSTOKES, -1>>(*m_geometry_1d));
        }

#endif

        m_los_source_terms.push_back(
            m_source_terms[m_source_terms.size() - 1].get());
    } else if (m_config.multiple_scatter_source() ==
               sasktran2::Config::MultipleScatterSource::hr) {
        m_source_terms.emplace_back(
            std::make_unique<sasktran2::hr::DiffuseTable<NSTOKES>>(
                *m_raytracer, *m_geometry_1d));
        m_los_source_terms.push_back(
            m_source_terms[m_source_terms.size() - 1].get());
    } else if (m_config.multiple_scatter_source() ==
               sasktran2::Config::MultipleScatterSource::twostream) {
        if constexpr (NSTOKES == 1) {
#ifdef SKTRAN_RUST_SUPPORT
            if (m_config.two_stream_backend() ==
                sasktran2::Config::TwoStreamBackend::rust) {
                m_source_terms.emplace_back(
                    std::make_unique<RustTwoStreamSourceAdapter<
                        sasktran2::twostream::SourceType::ONLY_SOLAR>>(
                        *m_geometry_1d));
            } else
#endif
            {
                m_source_terms.emplace_back(
                    std::make_unique<TwoStreamSource<
                        NSTOKES, sasktran2::twostream::SourceType::ONLY_SOLAR>>(
                        *m_geometry_1d));
            }
            m_los_source_terms.push_back(
                m_source_terms[m_source_terms.size() - 1].get());
        } else {
            spdlog::error(
                "TwoStreamSource is only implemented for NSTOKES = 1");
        }
    }

    for (auto& source : m_source_terms) {
        source->initialize_config(m_config);
    }
}

template <int NSTOKES> void Sasktran2<NSTOKES>::calculate_geometry() {
    FrameMarkStart("Geometry");
    ZoneScopedN("calculate_geometry");

    m_internal_viewing_geometry.traced_rays.clear();
    m_internal_viewing_geometry.flux_observers.clear();

    if (m_geometry_2d != nullptr) {
        m_internal_viewing_geometry.traced_rays.resize(
            m_viewing_geometry.observer_rays().size());

        for (int i = 0; i < m_viewing_geometry.observer_rays().size(); ++i) {
            const auto& viewing_ray = m_viewing_geometry.observer_rays()[i];
            auto ray = viewing_ray->construct_ray(m_geometry->coordinates());
#ifdef SKTRAN_RUST_SUPPORT
            m_raytracer_2d->trace_ray(
                ray, m_internal_viewing_geometry.traced_rays[i]);
#endif
        }

        m_source_integrator->initialize_geometry(
            m_internal_viewing_geometry.traced_rays, *m_geometry_2d);
        for (auto& source : m_source_terms) {
            source->initialize_geometry(m_internal_viewing_geometry);
        }

        FrameMarkEnd("Geometry");
        return;
    }

    // Trace every ray that we are given
    m_internal_viewing_geometry.traced_rays.resize(
        m_viewing_geometry.observer_rays().size());

    for (int i = 0; i < m_viewing_geometry.observer_rays().size(); ++i) {
        ZoneScopedN("Trace LOS Ray");
        const auto& viewing_ray = m_viewing_geometry.observer_rays()[i];
        sasktran2::viewinggeometry::ViewingRay ray =
            viewing_ray->construct_ray(m_geometry->coordinates());

#ifdef SASKTRAN_DEBUG_ASSERTS
        if (!ray.look_away.allFinite() || !ray.observer.position.allFinite()) {
            spdlog::error("Error constructing LOS ray: {}", i);

            // spdlog::error("ray look: {}", ray.look_away);
            // spdlog::error("obs: {}", ray.observer.position);

            spdlog::error("LOS Type: {}", typeid(viewing_ray).name());
        }
#endif

        m_raytracer->trace_ray(ray, m_internal_viewing_geometry.traced_rays[i],
                               m_config.los_refraction());
    }

    // Construct the flux observer points
    m_internal_viewing_geometry.flux_observers.resize(
        m_viewing_geometry.flux_observers().size());
    for (int i = 0; i < m_viewing_geometry.flux_observers().size(); ++i) {
        m_internal_viewing_geometry.flux_observers[i] =
            m_viewing_geometry.flux_observers()[i]->construct_flux_observer(
                m_geometry->coordinates());
    }

    // Initialize the integrator
    m_source_integrator->initialize_geometry(
        m_internal_viewing_geometry.traced_rays, *m_geometry_1d);

    for (auto& source : m_source_terms) {
        ZoneScopedN("Source Term Geometry Init");
        source->initialize_geometry(m_internal_viewing_geometry);
    }

    FrameMarkEnd("Geometry");
}

template <int NSTOKES>
void Sasktran2<NSTOKES>::validate_input_atmosphere(
    const sasktran2::atmosphere::Atmosphere<NSTOKES>& atmosphere) const {
    if (m_config.input_validation_mode() ==
        sasktran2::Config::InputValidationMode::disabled) {
        return;
    }

    // Check that we have the required legendre information for the number of
    // NSTOKES requested

    // Check that the atmosphere parameters have the correct dimensions
    if (atmosphere.storage().total_extinction.rows() != m_geometry->size()) {
        spdlog::error(
            "Atmosphere total extinction does not have the correct dimensions");
        throw std::runtime_error(
            "Invalid input. Check log for more information");
    }

    if (atmosphere.storage().ssa.rows() != m_geometry->size()) {
        spdlog::error("Atmosphere single scatter albedo does not have the "
                      "correct dimensions");
        throw std::runtime_error(
            "Invalid input. Check log for more information");
    }

    if (atmosphere.storage().total_extinction.cols() !=
        atmosphere.num_wavel()) {
        spdlog::error(
            "Atmosphere total extinction does not have the correct dimensions");
        throw std::runtime_error(
            "Invalid input. Check log for more information");
    }

    if (atmosphere.storage().ssa.cols() != atmosphere.num_wavel()) {
        spdlog::error("Atmosphere single scatter albedo does not have the "
                      "correct dimensions");
        throw std::runtime_error(
            "Invalid input. Check log for more information");
    }

    // Verify that all extinction values are finite and greater than 0
    sasktran2::validation::verify_finite(atmosphere.storage().total_extinction,
                                         "Atmosphere total extinction");
    sasktran2::validation::verify_greater_than(
        atmosphere.storage().total_extinction, "Atmosphere total extinction",
        0.0);

    // Verify that the SSA values are finite and betwen 0 and 1
    sasktran2::validation::verify_finite(atmosphere.storage().ssa,
                                         "Atmosphere single scatter albedo");
    sasktran2::validation::verify_greater_than(
        atmosphere.storage().ssa, "Atmosphere single scatter albedo", 0.0);
    sasktran2::validation::verify_less_than(
        atmosphere.storage().ssa, "Atmosphere single scatter albedo", 1.0);

    // Check that the atmosphere geometry matches the global geometry
    if (atmosphere.num_wavel() != atmosphere.surface().brdf_args().cols()) {
        spdlog::error("Atmosphere albedo does not have the correct dimensions");
    }
}

template <int NSTOKES>
void Sasktran2<NSTOKES>::calculate_radiance(
    const sasktran2::atmosphere::Atmosphere<NSTOKES>& atmosphere,
    sasktran2::Output<NSTOKES>& output, bool only_initialize) const {

#ifdef SKTRAN_OPENMP_SUPPORT
    omp_set_num_threads(m_config.num_threads());
    Eigen::setNbThreads(m_config.num_source_threads());
#endif

    validate_input_atmosphere(atmosphere);

    const_cast<sasktran2::atmosphere::AtmosphereGridStorageFull<NSTOKES>&>(
        atmosphere.storage())
        .determine_maximum_order();

    // Use this method for observer geometries and make a different method for
    // interior fluxes?

    const int wavelength_batch_size =
        effective_wavelength_batch_size(atmosphere.num_wavel());
    m_source_integrator->initialize_thread_storage(m_config.num_threads(),
                                                   wavelength_batch_size);

    // Initialize each source term with the atmosphere
    for (auto& source : m_source_terms) {
        source->set_wavelength_block_capacity(wavelength_batch_size);
        source->initialize_atmosphere(atmosphere);
    }

    m_source_integrator->initialize_atmosphere(atmosphere);
    m_source_integrator->initialize_derivative_sparsity(m_los_source_terms);

    auto& radiance =
        const_cast<std::vector<sasktran2::WavelengthBlockDual<NSTOKES>>&>(
            m_thread_radiance);
    radiance.resize(m_config.num_threads());
    for (auto& thread_radiance : radiance) {
        thread_radiance.resize(wavelength_batch_size, atmosphere.num_deriv(),
                               false);
    }
    auto& flux = const_cast<std::vector<
        sasktran2::Dual<double, sasktran2::dualstorage::dense, 1>>&>(
        m_thread_flux);
    flux.resize(m_config.num_threads());
    for (auto& thread_flux : flux) {
        thread_flux.resize(1, atmosphere.num_deriv(), false);
    }

    output.set_wavelength_block_capacity(wavelength_batch_size);
    output.initialize(m_config, *m_geometry, m_internal_viewing_geometry,
                      atmosphere);

    // Optical depth depends only on the initialized atmosphere and geometry.
    // Calculate it before the early return used by the Rust wavelength
    // threading path so threaded and non-threaded outputs are identical.
    if (m_config.output_los_optical_depth()) {
        m_source_integrator->integrate_optical_depth(
            output.los_optical_depth());
    }

    if (only_initialize) {
        return;
    }

    const int num_blocks =
        (atmosphere.num_wavel() + wavelength_batch_size - 1) /
        wavelength_batch_size;
#pragma omp parallel for num_threads(m_config.num_wavelength_threads())
    for (int block_index = 0; block_index < num_blocks; ++block_index) {
#ifdef SKTRAN_OPENMP_SUPPORT
        const int thread_idx = omp_get_thread_num();
#else
        const int thread_idx = 0;
#endif
        const int start = block_index * wavelength_batch_size;
        const sasktran2::WavelengthBlock<> block{
            start,
            std::min(wavelength_batch_size, atmosphere.num_wavel() - start)};
        calculate_radiance_block_thread(output, block, thread_idx);
    }
}

template <int NSTOKES>
int Sasktran2<NSTOKES>::effective_wavelength_batch_size(
    int num_wavelengths) const {
    const int requested_block_size = std::max(
        1, std::min(m_config.wavelength_batch_size(), num_wavelengths));
    int block_size = requested_block_size;
    if (!m_internal_viewing_geometry.flux_observers.empty()) {
        block_size = 1;
    }
    for (const auto& source : m_source_terms) {
        const int source_block_size = source->maximum_wavelength_block_size();
        if (source_block_size < 1) {
            throw std::logic_error(
                "Source reported an invalid maximum wavelength block size");
        }
        block_size = std::min(block_size, source_block_size);
    }
    if (block_size < requested_block_size) {
        spdlog::debug(
            "Reducing wavelength batch size from {} to {} for active sources "
            "and observers",
            requested_block_size, block_size);
    }
    return block_size;
}

template <int NSTOKES>
void Sasktran2<NSTOKES>::calculate_radiance_block_thread(
    sasktran2::Output<NSTOKES>& output,
    const sasktran2::WavelengthBlock<>& batch, int thread_idx) const {
    if (thread_idx < 0 || thread_idx >= m_thread_radiance.size() ||
        batch.count < 1 ||
        batch.count > m_thread_radiance[thread_idx].block_capacity()) {
        throw std::invalid_argument(
            "Invalid wavelength block for engine thread storage");
    }

    FrameMarkStart("WavelengthBlock");
    for (auto& source : m_source_terms) {
        ZoneScopedN("Batch Source Calculation");
        source->calculate(batch, thread_idx);
    }

    auto& radiance_storage =
        const_cast<std::vector<sasktran2::WavelengthBlockDual<NSTOKES>>&>(
            m_thread_radiance);
    auto& flux_storage = const_cast<std::vector<
        sasktran2::Dual<double, sasktran2::dualstorage::dense, 1>>&>(
        m_thread_flux);
#pragma omp parallel for num_threads(m_config.num_source_threads())            \
    schedule(dynamic)
    for (int ray_index = 0; ray_index < m_internal_viewing_geometry.num_rays();
         ++ray_index) {
#ifdef SKTRAN_OPENMP_SUPPORT
        const int ray_threadidx = omp_get_thread_num() + thread_idx;
#else
        const int ray_threadidx = thread_idx;
#endif
        auto& radiance = radiance_storage[ray_threadidx];
        radiance.set_zero(batch.count);
        m_source_integrator->integrate(radiance, m_los_source_terms, batch,
                                       ray_index, thread_idx, ray_threadidx);

        for (const auto* source : m_los_source_terms) {
            source->start_of_ray_source(batch, ray_index, thread_idx,
                                        ray_threadidx, radiance);
        }

        output.assign(batch, radiance, ray_index, ray_threadidx);
    }

#pragma omp parallel for num_threads(m_config.num_source_threads())            \
    schedule(dynamic)
    for (int i = 0; i < m_internal_viewing_geometry.flux_observers.size();
         ++i) {
#ifdef SKTRAN_OPENMP_SUPPORT
        int ray_threadidx = omp_get_thread_num() + thread_idx;
#else
        int ray_threadidx = thread_idx;
#endif
        auto& flux = flux_storage[ray_threadidx];
        for (int lane = 0; lane < batch.count; ++lane) {
            const int wavelength = batch.wavelength(lane);
            for (int flux_type_idx = 0;
                 flux_type_idx < m_config.get_flux_types().size();
                 ++flux_type_idx) {
                auto flux_type = m_config.get_flux_types()[flux_type_idx];
                flux.value.setZero();
                flux.deriv.setZero();

                for (const SourceTermInterface<NSTOKES>* source :
                     m_los_source_terms) {
                    source->flux(wavelength, i, thread_idx, ray_threadidx, flux,
                                 flux_type);
                }

                output.assign_flux(flux, i, wavelength, ray_threadidx,
                                   flux_type_idx);
            }
        }
    }
    FrameMarkEnd("WavelengthBlock");
}

template class Sasktran2<1>;
template class Sasktran2<3>;
