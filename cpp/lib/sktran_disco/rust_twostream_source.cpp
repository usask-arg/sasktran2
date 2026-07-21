#include "sktran_disco/twostream/rust_source.h"

#include <algorithm>
#include <cmath>
#include <cstdint>

#ifdef SKTRAN_RUST_SUPPORT

namespace {
    template <typename T>
    ::rust::Slice<const T> as_slice(const std::vector<T>& values) {
        return {values.data(), values.size()};
    }

    ::rust::Slice<const double> as_slice(const double* values,
                                         std::size_t size) {
        return {values, size};
    }

    std::pair<double, double>
    normalized_source_weights(const sasktran2::raytracing::TracedLayer& layer) {
        const double sum =
            layer.od_quad_start_fraction + layer.od_quad_end_fraction;
        if (std::isfinite(layer.od_quad_start_fraction) &&
            std::isfinite(layer.od_quad_end_fraction) && sum > 0) {
            return {layer.od_quad_start_fraction / sum,
                    layer.od_quad_end_fraction / sum};
        }
        return {0.5, 0.5};
    }
} // namespace

template <sasktran2::twostream::SourceType SOURCE_TYPE>
RustTwoStreamSourceAdapter<SOURCE_TYPE>::RustTwoStreamSourceAdapter(
    const sasktran2::Geometry1D& geometry)
    : m_geometry(geometry) {
    m_spec.configure(2, geometry.size() - 1);
}

template <sasktran2::twostream::SourceType SOURCE_TYPE>
void RustTwoStreamSourceAdapter<SOURCE_TYPE>::initialize_config(
    const sasktran2::Config& config) {
    m_config = &config;
}

template <sasktran2::twostream::SourceType SOURCE_TYPE>
void RustTwoStreamSourceAdapter<SOURCE_TYPE>::initialize_geometry(
    const sasktran2::viewinggeometry::InternalViewingGeometry&
        internal_viewing) {
    ZoneScopedN("Rust Twostream Initialize Geometry");
    const int nlyr = m_geometry.size() - 1;
    m_los_rays = &internal_viewing.traced_rays;

    const bool spherical = m_geometry.coordinates().geometry_type() ==
                           sasktran2::geometrytype::spherical;
    std::vector<double> sza_grid;
    const double reference_cos_sza =
        m_geometry.coordinates().cos_sza_at_reference();
    if constexpr (sasktran2::twostream::has_solar<SOURCE_TYPE>()) {
        const int requested_sza = m_config->num_do_sza();
        if (spherical && requested_sza > 1) {
            const auto bounds =
                sasktran2::raytracing::min_max_cos_sza_of_all_rays(*m_los_rays);
            if (bounds.first <= bounds.second &&
                bounds.second - bounds.first > 1e-14) {
                sza_grid.resize(requested_sza);
                for (int index = 0; index < requested_sza; ++index) {
                    sza_grid[index] =
                        bounds.first + (bounds.second - bounds.first) * index /
                                           (requested_sza - 1);
                }
            }
        }
    }
    if (sza_grid.empty()) {
        sza_grid.push_back(reference_cos_sza);
    }

    m_pconfigs.clear();
    m_geometry_layers.clear();
    m_pconfigs.reserve(sza_grid.size());
    m_geometry_layers.reserve(sza_grid.size());
    for (const double cos_sza : sza_grid) {
        auto pconfig =
            std::make_unique<sasktran_disco::PersistentConfiguration<1>>();
        pconfig->configure(m_spec, *m_config, cos_sza, nlyr,
                           internal_viewing.traced_rays);
        auto layers = std::make_unique<sasktran_disco::GeometryLayerArray<1>>(
            *pconfig, m_geometry);
        m_pconfigs.push_back(std::move(pconfig));
        m_geometry_layers.push_back(std::move(layers));
    }

    std::vector<double> layer_thickness(nlyr);
    std::vector<double> chapman_factors(sza_grid.size() * nlyr * nlyr);
    for (int layer = 0; layer < nlyr; ++layer) {
        layer_thickness[layer] = m_geometry_layers[0]->layer_ceiling()(layer) -
                                 m_geometry_layers[0]->layer_floor()(layer);
    }
    for (std::size_t sza = 0; sza < sza_grid.size(); ++sza) {
        for (int boundary = 0; boundary < nlyr; ++boundary) {
            for (int layer = 0; layer < nlyr; ++layer) {
                chapman_factors[sza * nlyr * nlyr + boundary * nlyr + layer] =
                    m_geometry_layers[sza]->chapman_factors()(boundary, layer);
            }
        }
    }

    constexpr int source_mode =
        sasktran2::twostream::has_solar<SOURCE_TYPE>() ? 0 : 1;
    m_rust_source.emplace(sasktran2::rust::twostream::new_rust_twostream_source(
        as_slice(layer_thickness),
        as_slice(chapman_factors.data(), nlyr * nlyr), sza_grid[0], source_mode,
        static_cast<std::size_t>(m_config->num_threads())));

    if (spherical) {
        const auto& altitude_grid = m_geometry.altitude_grid().grid();
        const std::size_t nlevel = altitude_grid.size();
        const double earth_radius = m_geometry.coordinates().earth_radius();
        std::vector<std::size_t> ray_offsets;
        std::vector<std::uint8_t> ground_hit;
        std::vector<double> ground_cos_sza;
        std::vector<std::size_t> segment_layers;
        std::vector<double> segment_fractions;
        std::vector<double> segment_cosines;
        std::vector<double> segment_relative_azimuths;
        std::vector<double> segment_cos_sza;
        std::vector<std::size_t> od_offsets;
        std::vector<std::size_t> od_indices;
        std::vector<double> od_weights;
        ray_offsets.reserve(m_los_rays->size() + 1);
        ground_hit.reserve(m_los_rays->size());
        ground_cos_sza.reserve(m_los_rays->size());
        ray_offsets.push_back(0);
        od_offsets.push_back(0);

        for (const auto& ray : *m_los_rays) {
            ground_hit.push_back(ray.ground_is_hit ? 1 : 0);
            ground_cos_sza.push_back(ray.layers.empty()
                                         ? reference_cos_sza
                                         : ray.layers.front().cos_sza_exit);
            for (std::size_t layer_index = 0; layer_index < ray.layers.size();
                 ++layer_index) {
                const auto& layer = ray.layers[layer_index];
                const auto source_weights = normalized_source_weights(layer);
                const Eigen::Vector3d position =
                    source_weights.first * layer.entrance.position +
                    source_weights.second * layer.exit.position;
                const double altitude = position.norm() - earth_radius;
                const double* upper =
                    std::upper_bound(altitude_grid.data(),
                                     altitude_grid.data() + nlevel, altitude);
                std::size_t bottom_layer =
                    upper == altitude_grid.data()
                        ? 0
                        : static_cast<std::size_t>(upper -
                                                   altitude_grid.data() - 1);
                bottom_layer = std::min(bottom_layer, nlevel - 2);
                const double bottom = altitude_grid(bottom_layer);
                const double top = altitude_grid(bottom_layer + 1);
                const double fraction_from_top =
                    std::clamp((top - altitude) / (top - bottom), 0.0, 1.0);
                const Eigen::Vector3d local_up = position.normalized();
                const double cosine = std::clamp(
                    -layer.average_look_away.dot(local_up), -1.0, 1.0);
                const double sin_azimuth =
                    source_weights.first * std::sin(layer.saz_entrance) +
                    source_weights.second * std::sin(layer.saz_exit);
                const double cos_azimuth =
                    source_weights.first * std::cos(layer.saz_entrance) +
                    source_weights.second * std::cos(layer.saz_exit);

                segment_layers.push_back(static_cast<std::size_t>(nlyr - 1) -
                                         bottom_layer);
                segment_fractions.push_back(fraction_from_top);
                segment_cosines.push_back(cosine);
                segment_relative_azimuths.push_back(
                    std::atan2(sin_azimuth, cos_azimuth));
                segment_cos_sza.push_back(
                    std::clamp(source_weights.first * layer.cos_sza_entrance +
                                   source_weights.second * layer.cos_sza_exit,
                               -1.0, 1.0));

                const auto stencil = ray.optical_depth_weights(layer_index);
                for (std::size_t entry = 0; entry < stencil.size(); ++entry) {
                    const auto [index, weight] = stencil[entry];
                    if (index < 0 ||
                        static_cast<std::size_t>(index) >= nlevel) {
                        throw sasktran_disco::InternalRuntimeError(
                            "Traced ray optical-depth stencil is outside the "
                            "one-dimensional atmosphere grid");
                    }
                    if (weight != 0) {
                        od_indices.push_back(nlevel - 1 - index);
                        od_weights.push_back(weight);
                    }
                }
                od_offsets.push_back(od_indices.size());
            }
            ray_offsets.push_back(segment_layers.size());
        }

        sasktran2::rust::twostream::set_spherical_geometry(
            **m_rust_source, as_slice(sza_grid), as_slice(chapman_factors),
            as_slice(ray_offsets), as_slice(ground_hit),
            as_slice(ground_cos_sza), as_slice(segment_layers),
            as_slice(segment_fractions), as_slice(segment_cosines),
            as_slice(segment_relative_azimuths), as_slice(segment_cos_sza),
            as_slice(od_offsets), as_slice(od_indices), as_slice(od_weights));
        return;
    }

    std::vector<double> viewing_cosines;
    std::vector<double> relative_azimuths;
    viewing_cosines.reserve(m_los_rays->size());
    relative_azimuths.reserve(m_los_rays->size());
    for (const auto& ray : *m_los_rays) {
        const double viewing_cosine = -ray.observer_and_look.look_away.z();
        if (viewing_cosine <= 0) {
            throw sasktran_disco::InternalRuntimeError(
                "Rust two-stream currently only supports upwelling plane "
                "parallel radiances");
        }
        viewing_cosines.push_back(viewing_cosine);
        relative_azimuths.push_back(-ray.observer_and_look.relative_azimuth);
    }
    sasktran2::rust::twostream::set_views(**m_rust_source,
                                          as_slice(viewing_cosines),
                                          as_slice(relative_azimuths));
}

template <sasktran2::twostream::SourceType SOURCE_TYPE>
void RustTwoStreamSourceAdapter<SOURCE_TYPE>::initialize_atmosphere(
    const sasktran2::atmosphere::Atmosphere<1>& atmosphere) {
    ZoneScopedN("Rust Twostream Atmosphere Batch");
    m_atmosphere = &atmosphere;
    if (m_los_rays->empty()) {
        return;
    }
    const auto& storage = atmosphere.storage();
    const std::size_t nwavel = atmosphere.num_wavel();
    const std::size_t nlevel = storage.total_extinction.rows();
    const std::size_t num_legendre = storage.leg_coeff.dimension(0);

    std::vector<double> surface_albedo(nwavel);
    for (std::size_t wave = 0; wave < nwavel; ++wave) {
        surface_albedo[wave] =
            atmosphere.surface().brdf(wave, 0, 0, 0)(0, 0) * EIGEN_PI;
    }

    sasktran2::rust::twostream::solve(
        **m_rust_source, nwavel, nlevel,
        as_slice(storage.total_extinction.data(), nlevel * nwavel),
        as_slice(storage.ssa.data(), nlevel * nwavel),
        as_slice(storage.leg_coeff.data(), num_legendre * nlevel * nwavel),
        num_legendre, as_slice(storage.f.data(), nlevel * nwavel),
        as_slice(storage.emission_source.data(), nlevel * nwavel),
        as_slice(surface_albedo),
        as_slice(atmosphere.surface().emission().data(), nwavel),
        as_slice(storage.solar_irradiance.data(), nwavel),
        atmosphere.num_deriv() > 0);
}

template <sasktran2::twostream::SourceType SOURCE_TYPE>
void RustTwoStreamSourceAdapter<SOURCE_TYPE>::start_of_ray_source(
    const sasktran2::WavelengthBlock<>& block, int losidx, int wavel_threadidx,
    int threadidx, sasktran2::WavelengthBlockDual<1>& source) const {
    ZoneScopedN("Rust Twostream Cached Source");
    const std::size_t nwavel = m_atmosphere->num_wavel();
    const std::size_t nlevel = m_atmosphere->storage().total_extinction.rows();
    const std::size_t wavelength_start = block.start;
    const std::size_t surface_index = losidx * nwavel + wavelength_start;
    const auto radiance = sasktran2::rust::twostream::radiance(**m_rust_source);
    source.value.row(0).head(block.count) +=
        Eigen::Map<const Eigen::RowVectorXd>(radiance.data() + surface_index,
                                             block.count);

    if (m_atmosphere->num_deriv() == 0) {
        return;
    }

    const auto extinction =
        sasktran2::rust::twostream::extinction_jacobian(**m_rust_source);
    const auto ssa = sasktran2::rust::twostream::ssa_jacobian(**m_rust_source);
    const auto b1 = sasktran2::rust::twostream::b1_jacobian(**m_rust_source);
    const auto emission =
        sasktran2::rust::twostream::emission_jacobian(**m_rust_source);
    const std::size_t view_level_offset = losidx * nlevel * nwavel;
    for (std::size_t rust_level = 0; rust_level < nlevel; ++rust_level) {
        const std::size_t cpp_level = nlevel - 1 - rust_level;
        const std::size_t index =
            view_level_offset + rust_level * nwavel + wavelength_start;
        source.deriv.row(cpp_level).head(block.count) +=
            Eigen::Map<const Eigen::RowVectorXd>(extinction.data() + index,
                                                 block.count);
        source.deriv.row(m_atmosphere->ssa_deriv_start_index() + cpp_level)
            .head(block.count) += Eigen::Map<const Eigen::RowVectorXd>(
            ssa.data() + index, block.count);
        for (int group = 0; group < m_atmosphere->num_scattering_deriv_groups();
             ++group) {
            const int derivative_index =
                m_atmosphere->scat_deriv_start_index() + group * nlevel +
                cpp_level;
            for (int lane = 0; lane < block.count; ++lane) {
                const int wavelength = block.wavelength(lane);
                double phase_derivative = m_atmosphere->storage().d_leg_coeff(
                    1, cpp_level, wavelength, group);
                if (m_atmosphere->storage().applied_f_order > 0) {
                    const double delta_m =
                        m_atmosphere->storage().f(cpp_level, wavelength);
                    // Reverse b1* - 3 f / (1 - f), matching the standard DO
                    // layer derivative for the two-stream l = 1 moment.
                    phase_derivative -= 3.0 *
                                        m_atmosphere->storage().d_f(
                                            cpp_level, wavelength, group) /
                                        ((1.0 - delta_m) * (1.0 - delta_m));
                }
                source.deriv(derivative_index, lane) +=
                    phase_derivative * b1[index + lane];
            }
        }
        if constexpr (sasktran2::twostream::has_thermal<SOURCE_TYPE>()) {
            if (m_atmosphere->include_emission_derivatives()) {
                source.deriv
                    .row(m_atmosphere->emission_deriv_start_index() + cpp_level)
                    .head(block.count) += Eigen::Map<const Eigen::RowVectorXd>(
                    emission.data() + index, block.count);
            }
        }
    }

    const auto albedo_jacobian =
        sasktran2::rust::twostream::surface_albedo_jacobian(**m_rust_source);
    for (int deriv = 0; deriv < m_atmosphere->surface().num_deriv(); ++deriv) {
        source.deriv.row(m_atmosphere->surface_deriv_start_index() + deriv)
            .head(block.count) += Eigen::Map<const Eigen::RowVectorXd>(
            albedo_jacobian.data() + surface_index, block.count);
    }
    if constexpr (sasktran2::twostream::has_thermal<SOURCE_TYPE>()) {
        if (m_atmosphere->include_emission_derivatives()) {
            const auto surface_emission =
                sasktran2::rust::twostream::surface_emission_jacobian(
                    **m_rust_source);
            source.deriv.row(m_atmosphere->surface_emission_deriv_start_index())
                .head(block.count) += Eigen::Map<const Eigen::RowVectorXd>(
                surface_emission.data() + surface_index, block.count);
        }
    }
}

template class RustTwoStreamSourceAdapter<
    sasktran2::twostream::SourceType::ONLY_SOLAR>;
template class RustTwoStreamSourceAdapter<
    sasktran2::twostream::SourceType::ONLY_THERMAL>;

#endif
