#include "sasktran2/atmosphere/grid_storage.h"
#include "sasktran2/do_source.h"
#include "sasktran2/geometry.h"
#include "sasktran2/raytracing.h"
#include "sasktran2/source_interface.h"
#include <memory>
#include <sasktran2.h>
#include <sasktran2/validation/validation.h>
#ifdef SKTRAN_OPENMP_SUPPORT
#include <omp.h>
#endif

template <int NSTOKES> void Sasktran2<NSTOKES>::construct_raytracer() {
    if (m_geometry->coordinates().geometry_type() ==
        sasktran2::geometrytype::spherical) {
        m_raytracer =
            std::make_unique<sasktran2::raytracing::SphericalShellRayTracer>(
                *m_geometry);
    } else if (m_geometry->coordinates().geometry_type() ==
                   sasktran2::geometrytype::planeparallel ||
               m_geometry->coordinates().geometry_type() ==
                   sasktran2::geometrytype::pseudospherical) {
        m_raytracer =
            std::make_unique<sasktran2::raytracing::PlaneParallelRayTracer>(
                *m_geometry);
    } else {
        spdlog::error("Requested geometry type is not yet supported.");
    }
}

template <int NSTOKES> void Sasktran2<NSTOKES>::construct_integrator() {
    m_source_integrator =
        std::make_unique<sasktran2::SourceIntegrator<NSTOKES>>(true);
}

template <int NSTOKES> void Sasktran2<NSTOKES>::construct_source_terms() {

    if (m_config.single_scatter_source() ==
        sasktran2::Config::SingleScatterSource::exact) {
        m_source_terms.emplace_back(
            std::make_unique<sasktran2::solartransmission::SingleScatterSource<
                sasktran2::solartransmission::SolarTransmissionExact, NSTOKES>>(
                *m_geometry, *m_raytracer));

        m_los_source_terms.push_back(
            m_source_terms[m_source_terms.size() - 1].get());
    }

    if (m_config.single_scatter_source() ==
        sasktran2::Config::SingleScatterSource::solartable) {
        m_source_terms.emplace_back(
            std::make_unique<sasktran2::solartransmission::SingleScatterSource<
                sasktran2::solartransmission::SolarTransmissionTable, NSTOKES>>(
                *m_geometry, *m_raytracer));

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
            std::make_unique<sasktran2::emission::EmissionSource<NSTOKES>>());

        m_los_source_terms.push_back(
            m_source_terms[m_source_terms.size() - 1].get());
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
                                NSTOKES, 2>>(*m_geometry, *m_raytracer));
                } else {
                    m_source_terms.emplace_back(
                        std::make_unique<
                            sasktran2::DOSourcePlaneParallelPostProcessing<
                                NSTOKES, 2>>(*m_geometry));
                }
            } else if (m_config.num_do_streams() == 4) {
                if (m_geometry->coordinates().geometry_type() ==
                    sasktran2::geometrytype::spherical) {
                    m_source_terms.emplace_back(
                        std::make_unique<
                            sasktran2::DOSourceInterpolatedPostProcessing<
                                NSTOKES, 4>>(*m_geometry, *m_raytracer));
                } else {
                    m_source_terms.emplace_back(
                        std::make_unique<
                            sasktran2::DOSourcePlaneParallelPostProcessing<
                                NSTOKES, 4>>(*m_geometry));
                }
            } else {
                if (m_geometry->coordinates().geometry_type() ==
                    sasktran2::geometrytype::spherical) {
                    m_source_terms.emplace_back(
                        std::make_unique<
                            sasktran2::DOSourceInterpolatedPostProcessing<
                                NSTOKES, -1>>(*m_geometry, *m_raytracer));
                } else {
                    m_source_terms.emplace_back(
                        std::make_unique<
                            sasktran2::DOSourcePlaneParallelPostProcessing<
                                NSTOKES, -1>>(*m_geometry));
                }
            }
        } else {
            if (m_geometry->coordinates().geometry_type() ==
                sasktran2::geometrytype::spherical) {
                m_source_terms.emplace_back(
                    std::make_unique<
                        sasktran2::DOSourceInterpolatedPostProcessing<NSTOKES,
                                                                      -1>>(
                        *m_geometry, *m_raytracer));
            } else {
                m_source_terms.emplace_back(
                    std::make_unique<
                        sasktran2::DOSourcePlaneParallelPostProcessing<NSTOKES,
                                                                       -1>>(
                        *m_geometry));
            }
        }
#else
        if (m_geometry->coordinates().geometry_type() ==
            sasktran2::geometrytype::spherical) {
            m_source_terms.emplace_back(
                std::make_unique<
                    sasktran2::DOSourceInterpolatedPostProcessing<NSTOKES, -1>>(
                    *m_geometry, *m_raytracer));
        } else {
            m_source_terms.emplace_back(
                std::make_unique<sasktran2::DOSourcePlaneParallelPostProcessing<
                    NSTOKES, -1>>(*m_geometry));
        }

#endif

        m_los_source_terms.push_back(
            m_source_terms[m_source_terms.size() - 1].get());
    } else if (m_config.multiple_scatter_source() ==
               sasktran2::Config::MultipleScatterSource::hr) {
        m_source_terms.emplace_back(
            std::make_unique<sasktran2::hr::DiffuseTable<NSTOKES>>(
                *m_raytracer, *m_geometry));
        m_los_source_terms.push_back(
            m_source_terms[m_source_terms.size() - 1].get());
    } else if (m_config.multiple_scatter_source() ==
               sasktran2::Config::MultipleScatterSource::twostream) {
        if constexpr (NSTOKES == 1) {
            m_source_terms.emplace_back(
                std::make_unique<TwoStreamSource<NSTOKES>>(*m_geometry));
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
    // Trace every ray that we are given
    m_traced_rays.clear();
    m_traced_rays.resize(m_viewing_geometry.observer_rays().size());

    for (int i = 0; i < m_viewing_geometry.observer_rays().size(); ++i) {
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

        m_raytracer->trace_ray(ray, m_traced_rays[i],
                               m_config.los_refraction());
    }

    // Initialize the integrator
    m_source_integrator->initialize_geometry(m_traced_rays, *m_geometry);

    for (auto& source : m_source_terms) {
        source->initialize_geometry(m_traced_rays);
    }
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
    sasktran2::Output<NSTOKES>& output) const {
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

    // Initialize each source term with the atmosphere
    for (auto& source : m_source_terms) {
        source->initialize_atmosphere(atmosphere);
    }

    m_source_integrator->initialize_atmosphere(atmosphere);

    // Allocate memory, should be moved to thread storage?
    std::vector<sasktran2::Dual<double, sasktran2::dualstorage::dense, NSTOKES>>
        radiance(m_config.num_threads(),
                 {NSTOKES, atmosphere.num_deriv(), true});

    output.initialize(m_config, *m_geometry, m_traced_rays, atmosphere);

#pragma omp parallel for num_threads(m_config.num_wavelength_threads())
    for (int w = 0; w < atmosphere.num_wavel(); ++w) {
#ifdef SKTRAN_OPENMP_SUPPORT
        int thread_idx = omp_get_thread_num();
#else
        int thread_idx = 0;
#endif

        // Trigger source term generation for this wavelength
        for (auto& source : m_source_terms) {
            source->calculate(w, thread_idx);
        }

#pragma omp parallel for num_threads(m_config.num_source_threads())            \
    schedule(dynamic)
        for (int i = 0; i < m_traced_rays.size(); ++i) {
#ifdef SKTRAN_OPENMP_SUPPORT
            int ray_threadidx = omp_get_thread_num() + thread_idx;
#else
            int ray_threadidx = 0;
#endif

            // Set the radiance thread storage to 0
            radiance[ray_threadidx].value.setZero();
            radiance[ray_threadidx].deriv.setZero();

            // Integrate all of the sources for the ray
            m_source_integrator->integrate(radiance[ray_threadidx],
                                           m_los_source_terms, w, i, thread_idx,
                                           ray_threadidx);

            // Add on any start of ray sources
            for (const SourceTermInterface<NSTOKES>* source :
                 m_los_source_terms) {
                source->start_of_ray_source(w, i, thread_idx, ray_threadidx,
                                            radiance[ray_threadidx]);
            }

            // And assign it to the output
            output.assign(radiance[ray_threadidx], i, w, ray_threadidx);
        }

        // TODO: Is this where we should generate fluxes or other quantities
        // that aren't through the integrator?
    }
}

template class Sasktran2<1>;
template class Sasktran2<3>;
