#include "sasktran2/atmosphere/grid_storage.h"
#include "sasktran2/do_source.h"
#include "sasktran2/geometry.h"
#include "sasktran2/raytracing.h"
#include "sasktran2/runtime_backend_tuner.h"
#include "sasktran2/source_interface.h"
#include "../successive_orders/factory.h"
#include "sktran_disco/twostream/cpp_source.h"
#include "sktran_disco/twostream/meta.h"
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
                 sasktran2::Config::SingleScatterSource::exact &&
             m_config.single_scatter_source() !=
                 sasktran2::Config::SingleScatterSource::solartable) ||
            (m_config.multiple_scatter_source() !=
                 sasktran2::Config::MultipleScatterSource::none &&
             m_config.multiple_scatter_source() !=
                 sasktran2::Config::MultipleScatterSource::successive_orders) ||
            (m_config.emission_source() !=
                 sasktran2::Config::EmissionSource::none &&
             m_config.emission_source() !=
                 sasktran2::Config::EmissionSource::standard &&
             m_config.emission_source() !=
                 sasktran2::Config::EmissionSource::volume_emission_rate)) {
            throw std::invalid_argument(
                "Geometry2D currently supports exact and table single "
                "scattering, "
                "successive-orders multiple scattering, occultation, "
                "standard emission, and volume emission rate sources");
        }
        if (!m_viewing_geometry.flux_observers().empty()) {
            throw std::invalid_argument(
                "Geometry2D does not yet support flux observers");
        }
        if (m_config.multiple_scatter_source() ==
                sasktran2::Config::MultipleScatterSource::successive_orders &&
            m_config.multiple_scatter_refraction()) {
            throw std::invalid_argument(
                "Geometry2D successive orders does not support diffuse-ray "
                "refraction");
        }
    }

    m_config.validate_config_geometry(
        m_geometry->coordinates().geometry_type());
    if (m_geometry_1d != nullptr) {
        sasktran2::detail::RuntimeBackendTuner::resolve(
            m_config,
            static_cast<int>(m_geometry_1d->altitude_grid().grid().size()) - 1);
    }
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
#ifdef SKTRAN_RUST_SUPPORT
        std::shared_ptr<sasktran2::solartransmission::SolarTransmissionTable2D>
            shared_solar_table;
        const bool needs_solar_table =
            m_config.single_scatter_source() ==
                sasktran2::Config::SingleScatterSource::solartable ||
            (m_config.single_scatter_source() ==
                 sasktran2::Config::SingleScatterSource::exact &&
             m_config.solar_refraction()) ||
            m_config.multiple_scatter_source() ==
                sasktran2::Config::MultipleScatterSource::successive_orders;
        if (needs_solar_table) {
            shared_solar_table = std::make_shared<
                sasktran2::solartransmission::SolarTransmissionTable2D>(
                *m_geometry_2d, *m_raytracer_2d);
        }
#endif
        if (m_config.single_scatter_source() ==
            sasktran2::Config::SingleScatterSource::exact) {
#ifdef SKTRAN_RUST_SUPPORT
            m_source_terms.emplace_back(
                std::make_unique<
                    sasktran2::solartransmission::SingleScatterSource<
                        sasktran2::solartransmission::SolarTransmissionExact,
                        NSTOKES>>(*m_geometry_2d, *m_raytracer_2d,
                                  m_config.solar_refraction()
                                      ? shared_solar_table
                                      : nullptr));
            m_los_source_terms.push_back(m_source_terms.back().get());
#else
            throw std::invalid_argument(
                "Geometry2D exact single scattering requires Rust support");
#endif
        }
        if (m_config.single_scatter_source() ==
            sasktran2::Config::SingleScatterSource::solartable) {
#ifdef SKTRAN_RUST_SUPPORT
            m_source_terms.emplace_back(
                std::make_unique<
                    sasktran2::solartransmission::SingleScatterSource<
                        sasktran2::solartransmission::SolarTransmissionTable2D,
                        NSTOKES>>(*m_geometry_2d, *m_raytracer_2d,
                                  shared_solar_table));
            m_los_source_terms.push_back(m_source_terms.back().get());
#else
            throw std::invalid_argument(
                "Geometry2D table single scattering requires Rust support");
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
        if (m_config.multiple_scatter_source() ==
            sasktran2::Config::MultipleScatterSource::successive_orders) {
#ifdef SKTRAN_RUST_SUPPORT
            m_source_terms.emplace_back(
                sasktran2::successive_orders::make_successive_orders_source<
                    NSTOKES>(*m_raytracer_2d, *m_geometry_2d,
                             shared_solar_table));
            m_los_source_terms.push_back(m_source_terms.back().get());
#else
            throw std::invalid_argument(
                "Geometry2D successive orders requires Rust support");
#endif
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
            m_source_terms.emplace_back(
                std::make_unique<CppTwoStreamSourceAdapter<
                    sasktran2::twostream::SourceType::ONLY_THERMAL>>(
                    *m_geometry_1d));

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
               sasktran2::Config::MultipleScatterSource::
                   successive_orders_legacy) {
        m_source_terms.emplace_back(
            std::make_unique<sasktran2::hr::DiffuseTable<NSTOKES>>(
                *m_raytracer, *m_geometry_1d));
        m_los_source_terms.push_back(
            m_source_terms[m_source_terms.size() - 1].get());
    } else if (m_config.multiple_scatter_source() ==
               sasktran2::Config::MultipleScatterSource::successive_orders) {
        m_source_terms.emplace_back(
            sasktran2::successive_orders::make_successive_orders_source<
                NSTOKES>(*m_raytracer, *m_geometry_1d));
        m_los_source_terms.push_back(m_source_terms.back().get());
    } else if (m_config.multiple_scatter_source() ==
               sasktran2::Config::MultipleScatterSource::twostream) {
        if constexpr (NSTOKES == 1) {
            m_source_terms.emplace_back(
                std::make_unique<CppTwoStreamSourceAdapter<
                    sasktran2::twostream::SourceType::ONLY_SOLAR>>(
                    *m_geometry_1d));
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

template <int NSTOKES>
void Sasktran2<NSTOKES>::calculate_geometry(bool refresh_los_sources) {
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
            const bool refract_ray =
                !m_refractive_profiles_2d.empty() &&
                !m_refractive_profiles_2d[i].array().isApproxToConstant(1.0,
                                                                        0.0);
            if (refract_ray) {
                m_raytracer_2d->trace_ray(
                    ray, m_refractive_profiles_2d[i],
                    m_internal_viewing_geometry.traced_rays[i]);
            } else {
                m_raytracer_2d->trace_ray(
                    ray, m_internal_viewing_geometry.traced_rays[i]);
            }
#endif
        }

        m_source_integrator->initialize_geometry(
            m_internal_viewing_geometry.traced_rays, *m_geometry_2d);
        for (auto& source : m_source_terms) {
            if (refresh_los_sources) {
                source->refresh_los_geometry(m_internal_viewing_geometry);
            } else {
                source->initialize_geometry(m_internal_viewing_geometry);
            }
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
void Sasktran2<NSTOKES>::set_2d_refractive_profiles(
    const Eigen::Ref<const Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic,
                                         Eigen::RowMajor>>& profiles) {
    if (m_geometry_2d == nullptr) {
        throw std::invalid_argument(
            "Per-ray refractive profiles are only valid for Geometry2D");
    }
    const Eigen::Index expected_rays =
        static_cast<Eigen::Index>(m_viewing_geometry.observer_rays().size());
    if (profiles.rows() != expected_rays ||
        profiles.cols() != m_geometry_2d->num_altitudes()) {
        throw std::invalid_argument(
            "Geometry2D refractive profiles must have shape (num_los, "
            "num_altitudes)");
    }
    if (!profiles.allFinite() || (profiles.array() <= 0.0).any()) {
        throw std::invalid_argument(
            "Geometry2D refractive profiles must be finite and positive");
    }

    m_refractive_profiles_2d.resize(static_cast<std::size_t>(profiles.rows()));
    for (Eigen::Index ray = 0; ray < profiles.rows(); ++ray) {
        m_refractive_profiles_2d[static_cast<std::size_t>(ray)] =
            profiles.row(ray).transpose();
    }
    calculate_geometry(true);
}

template <int NSTOKES>
void Sasktran2<NSTOKES>::assign_2d_surface_interpolation_weights(
    Eigen::Ref<RowMajorMatrix> weights) const {
    if (m_geometry_2d == nullptr) {
        throw std::logic_error(
            "Surface interpolation weights require a Geometry2D engine");
    }
    if (weights.rows() != m_internal_viewing_geometry.num_rays() ||
        weights.cols() != m_geometry_2d->num_horizontal_locations()) {
        throw std::invalid_argument(
            "Surface interpolation weight buffer has incorrect dimensions");
    }
    weights.setZero();
    std::vector<std::pair<int, double>> stencil;
    for (int ray_index = 0; ray_index < m_internal_viewing_geometry.num_rays();
         ++ray_index) {
        const auto& ray = m_internal_viewing_geometry.traced_rays[ray_index];
        if (!ray.ground_is_hit || ray.layers.empty()) {
            continue;
        }
        m_geometry_2d->assign_horizontal_interpolation_weights(
            ray.layers.front().exit, stencil);
        for (const auto& [horizontal_index, weight] : stencil) {
            weights(ray_index, horizontal_index) += weight;
        }
    }
}

template <int NSTOKES>
void Sasktran2<NSTOKES>::assign_2d_horizontal_edge_usage(
    Eigen::Ref<RowMajorIntMatrix> usage) const {
    if (m_geometry_2d == nullptr) {
        throw std::logic_error(
            "Horizontal edge usage requires a Geometry2D engine");
    }
    if (usage.rows() != m_internal_viewing_geometry.num_rays() ||
        usage.cols() != 2) {
        throw std::invalid_argument(
            "Horizontal edge usage buffer has incorrect dimensions");
    }

    usage.setZero();
    const auto& horizontal_grid = m_geometry_2d->horizontal_angle_grid();
    const double lower = horizontal_grid[0];
    const double upper = horizontal_grid[horizontal_grid.size() - 1];
    constexpr double angular_tolerance_rad = 1e-8;

    for (int ray_index = 0; ray_index < m_internal_viewing_geometry.num_rays();
         ++ray_index) {
        const auto& ray = m_internal_viewing_geometry.traced_rays[ray_index];
        for (const auto& layer : ray.layers) {
            for (const auto* location : {&layer.entrance, &layer.exit}) {
                const double angle =
                    m_geometry_2d->horizontal_angle_at(*location);
                if (angle < lower - angular_tolerance_rad) {
                    usage(ray_index, 0) = 1;
                }
                if (angle > upper + angular_tolerance_rad) {
                    usage(ray_index, 1) = 1;
                }
            }
        }
    }
}

template <int NSTOKES>
void Sasktran2<NSTOKES>::validate_input_atmosphere(
    const sasktran2::atmosphere::Atmosphere<NSTOKES>& atmosphere) const {
    if (m_geometry_2d != nullptr && m_config.los_refraction() &&
        m_refractive_profiles_2d.empty()) {
        throw std::runtime_error(
            "Geometry2D line-of-sight refraction requires one refractive-"
            "index profile per viewing ray");
    }
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
void Sasktran2<NSTOKES>::initialize_jvp(
    const sasktran2::atmosphere::Atmosphere<NSTOKES>& atmosphere,
    sasktran2::OutputJVP<NSTOKES>& output) const {
#ifdef SKTRAN_OPENMP_SUPPORT
    omp_set_num_threads(m_config.num_threads());
    Eigen::setNbThreads(m_config.num_source_threads());
#endif

    if (!supports_linearization(sasktran2::LinearizationMode::JVP)) {
        throw std::logic_error(
            "The configured sources do not implement a native JVP");
    }
    validate_input_atmosphere(atmosphere);
    if (atmosphere.num_deriv() == 0) {
        throw std::invalid_argument(
            "JVP calculation requires an atmosphere constructed with "
            "calculate_derivatives=true");
    }
    const_cast<sasktran2::atmosphere::AtmosphereGridStorageFull<NSTOKES>&>(
        atmosphere.storage())
        .determine_maximum_order();

    // The native source JVP hooks currently operate on one wavelength at a
    // time. Feeding them from batched source caches adds dispatch/cache
    // overhead without vectorizing the product itself, so retain the scalar
    // source path until a native block-JVP hook is available.
    constexpr int wavelength_batch_size = 1;
    m_source_integrator->initialize_thread_storage(m_config.num_threads(), 1);
    for (auto& source : m_source_terms) {
        source->set_wavelength_block_capacity(wavelength_batch_size);
        source->initialize_atmosphere_native(atmosphere);
    }
    m_source_integrator->initialize_atmosphere(atmosphere);
    output.set_wavelength_block_capacity(wavelength_batch_size);
    output.initialize(m_config, *m_geometry, m_internal_viewing_geometry,
                      atmosphere);

    m_native_jvp_tangents.resize(m_config.num_threads());
    for (auto& tangent : m_native_jvp_tangents) {
        tangent.resize(atmosphere.num_deriv());
    }
}

template <int NSTOKES>
void Sasktran2<NSTOKES>::calculate_jvp_wavelength_thread(
    sasktran2::OutputJVP<NSTOKES>& output, int wavelength,
    int wavelength_threadidx) const {
    if (wavelength < 0 || wavelength >= output.num_wavel() ||
        wavelength_threadidx < 0 ||
        wavelength_threadidx >= m_config.num_wavelength_threads() ||
        wavelength_threadidx >=
            static_cast<int>(m_native_jvp_tangents.size())) {
        throw std::invalid_argument(
            "Invalid wavelength or worker for native JVP");
    }
    const sasktran2::WavelengthBlock<> block{wavelength, 1};
    auto& native_tangent = m_native_jvp_tangents[wavelength_threadidx];
    output.native_tangent(wavelength, native_tangent);
    // Source-wide directional state must be complete before LOS workers
    // consume it concurrently.
    for (auto& source : m_source_terms) {
        source->calculate_jvp(block, wavelength_threadidx, native_tangent);
    }
#pragma omp parallel for num_threads(m_config.num_source_threads())            \
    schedule(dynamic)
    for (int ray = 0; ray < m_internal_viewing_geometry.num_rays(); ++ray) {
#ifdef SKTRAN_OPENMP_SUPPORT
        const int ray_threadidx = omp_get_thread_num() + wavelength_threadidx;
#else
        const int ray_threadidx = wavelength_threadidx;
#endif
        sasktran2::RadianceJVP<NSTOKES> radiance;
        m_source_integrator->integrate_jvp(
            radiance, m_los_source_terms, wavelength, ray, wavelength_threadidx,
            ray_threadidx, native_tangent);
        output.assign_native(ray, wavelength, radiance.value, radiance.jvp);
    }
}

template <int NSTOKES>
void Sasktran2<NSTOKES>::calculate_jvp(
    const sasktran2::atmosphere::Atmosphere<NSTOKES>& atmosphere,
    sasktran2::OutputJVP<NSTOKES>& output) const {
    initialize_jvp(atmosphere, output);
    const int num_wavelength_threads = m_config.num_wavelength_threads();
#pragma omp parallel for num_threads(num_wavelength_threads)
    for (int wavelength = 0; wavelength < atmosphere.num_wavel();
         ++wavelength) {
#ifdef SKTRAN_OPENMP_SUPPORT
        const int wavelength_threadidx = omp_get_thread_num();
#else
        const int wavelength_threadidx = 0;
#endif
        calculate_jvp_wavelength_thread(output, wavelength,
                                        wavelength_threadidx);
    }
}

template <int NSTOKES>
void Sasktran2<NSTOKES>::initialize_vjp(
    const sasktran2::atmosphere::Atmosphere<NSTOKES>& atmosphere,
    sasktran2::OutputVJP<NSTOKES>& output) const {
#ifdef SKTRAN_OPENMP_SUPPORT
    omp_set_num_threads(m_config.num_threads());
    Eigen::setNbThreads(m_config.num_source_threads());
#endif

    if (!supports_linearization(sasktran2::LinearizationMode::VJP)) {
        throw std::logic_error(
            "The configured sources do not implement a native VJP");
    }
    validate_input_atmosphere(atmosphere);
    if (atmosphere.num_deriv() == 0) {
        throw std::invalid_argument(
            "VJP calculation requires an atmosphere constructed with "
            "calculate_derivatives=true");
    }
    const_cast<sasktran2::atmosphere::AtmosphereGridStorageFull<NSTOKES>&>(
        atmosphere.storage())
        .determine_maximum_order();

    const int wavelength_batch_size =
        effective_wavelength_batch_size(atmosphere.num_wavel());
    m_source_integrator->initialize_thread_storage(m_config.num_threads(),
                                                   wavelength_batch_size);
    m_source_integrator->initialize_vjp_storage(
        m_config.num_threads(), static_cast<int>(m_los_source_terms.size()));
    for (auto& source : m_source_terms) {
        source->set_wavelength_block_capacity(wavelength_batch_size);
        source->initialize_atmosphere_native(atmosphere);
    }
    m_source_integrator->initialize_atmosphere(atmosphere);
    output.set_wavelength_block_capacity(wavelength_batch_size);
    output.initialize(m_config, *m_geometry, m_internal_viewing_geometry,
                      atmosphere);

    m_native_vjp_gradients.resize(m_config.num_threads());
    m_native_vjp_radiance.resize(m_config.num_threads());
    m_native_vjp_cotangent.resize(m_config.num_threads());
    for (int thread = 0; thread < m_config.num_threads(); ++thread) {
        m_native_vjp_gradients[thread].resize(atmosphere.num_deriv(),
                                              wavelength_batch_size);
        m_native_vjp_radiance[thread].resize(NSTOKES, wavelength_batch_size);
        m_native_vjp_cotangent[thread].resize(NSTOKES, wavelength_batch_size);
    }
}

template <int NSTOKES>
void Sasktran2<NSTOKES>::calculate_vjp_block_thread(
    sasktran2::OutputVJP<NSTOKES>& output,
    const sasktran2::WavelengthBlock<>& block, int wavelength_threadidx) const {
    if (block.start < 0 || block.count < 1 ||
        block.end() > output.num_wavel() || wavelength_threadidx < 0 ||
        wavelength_threadidx >= m_config.num_wavelength_threads() ||
        wavelength_threadidx >=
            static_cast<int>(m_native_vjp_gradients.size()) ||
        wavelength_threadidx >=
            static_cast<int>(m_native_vjp_radiance.size()) ||
        wavelength_threadidx >=
            static_cast<int>(m_native_vjp_cotangent.size()) ||
        block.count > m_native_vjp_gradients[wavelength_threadidx].cols()) {
        throw std::invalid_argument(
            "Invalid wavelength block or worker for native VJP");
    }
    for (auto& source : m_source_terms) {
        source->calculate_vjp(block, wavelength_threadidx);
    }
    const int num_source_threads = m_config.num_source_threads();
    if (num_source_threads == 1) {
        m_native_vjp_gradients[wavelength_threadidx]
            .leftCols(block.count)
            .setZero();
    } else {
        for (int thread = 0; thread < num_source_threads; ++thread) {
            m_native_vjp_gradients[thread].leftCols(block.count).setZero();
        }
    }
#pragma omp parallel for num_threads(num_source_threads) schedule(dynamic)
    for (int ray = 0; ray < m_internal_viewing_geometry.num_rays(); ++ray) {
#ifdef SKTRAN_OPENMP_SUPPORT
        const int ray_threadidx = omp_get_thread_num() + wavelength_threadidx;
#else
        const int ray_threadidx = wavelength_threadidx;
#endif
        auto& radiance = m_native_vjp_radiance[ray_threadidx];
        auto& cotangent = m_native_vjp_cotangent[ray_threadidx];
        for (int lane = 0; lane < block.count; ++lane) {
            cotangent.col(lane) =
                output.native_cotangent(ray, block.wavelength(lane));
        }
        m_source_integrator->integrate_vjp_block(
            radiance, m_los_source_terms, block, ray, wavelength_threadidx,
            ray_threadidx, cotangent, m_native_vjp_gradients[ray_threadidx]);
        for (int lane = 0; lane < block.count; ++lane) {
            output.assign_native_value(ray, block.wavelength(lane),
                                       radiance.col(lane));
        }
    }
    auto& reduced_gradient = m_native_vjp_gradients[wavelength_threadidx];
    if (num_source_threads > 1) {
        for (int thread = 1; thread < num_source_threads; ++thread) {
            reduced_gradient.leftCols(block.count) +=
                m_native_vjp_gradients[thread].leftCols(block.count);
        }
    }
    // Globally coupled sources complete one adjoint after all LOS cotangents
    // have reached their thread-local accumulators.
    auto reduced_gradient_block = reduced_gradient.leftCols(block.count);
    for (const auto& source : m_source_terms) {
        source->finalize_vjp(block, wavelength_threadidx,
                             reduced_gradient_block);
    }
    output.accumulate_native_gradient(block, wavelength_threadidx,
                                      reduced_gradient);
}

template <int NSTOKES>
void Sasktran2<NSTOKES>::calculate_vjp(
    const sasktran2::atmosphere::Atmosphere<NSTOKES>& atmosphere,
    sasktran2::OutputVJP<NSTOKES>& output) const {
    initialize_vjp(atmosphere, output);
    const int wavelength_batch_size =
        effective_wavelength_batch_size(atmosphere.num_wavel());
    const int num_wavelength_threads = m_config.num_wavelength_threads();
    const int num_blocks =
        (atmosphere.num_wavel() + wavelength_batch_size - 1) /
        wavelength_batch_size;
#pragma omp parallel for num_threads(num_wavelength_threads)
    for (int block_index = 0; block_index < num_blocks; ++block_index) {
#ifdef SKTRAN_OPENMP_SUPPORT
        const int wavelength_threadidx = omp_get_thread_num();
#else
        const int wavelength_threadidx = 0;
#endif
        const int start = block_index * wavelength_batch_size;
        const sasktran2::WavelengthBlock<> block{
            start,
            std::min(wavelength_batch_size, atmosphere.num_wavel() - start)};
        calculate_vjp_block_thread(output, block, wavelength_threadidx);
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
