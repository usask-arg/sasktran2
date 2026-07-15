#include "sktran_disco/twostream/rust_source.h"

#ifdef SKTRAN_RUST_SUPPORT

namespace {
    ::rust::Slice<const double> as_slice(const std::vector<double>& values) {
        return {values.data(), values.size()};
    }

    ::rust::Slice<const double> as_slice(const double* values,
                                         std::size_t size) {
        return {values, size};
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
    m_pconfig.configure(m_spec, *m_config,
                        m_geometry.coordinates().cos_sza_at_reference(), nlyr,
                        internal_viewing.traced_rays);
    m_los_rays = &internal_viewing.traced_rays;
    m_geometry_layers = std::make_unique<sasktran_disco::GeometryLayerArray<1>>(
        m_pconfig, m_geometry);

    std::vector<double> layer_thickness(nlyr);
    std::vector<double> chapman_factors(nlyr * nlyr);
    for (int layer = 0; layer < nlyr; ++layer) {
        layer_thickness[layer] = m_geometry_layers->layer_ceiling()(layer) -
                                 m_geometry_layers->layer_floor()(layer);
        for (int boundary = 0; boundary < nlyr; ++boundary) {
            chapman_factors[boundary * nlyr + layer] =
                m_geometry_layers->chapman_factors()(boundary, layer);
        }
    }

    constexpr int source_mode =
        sasktran2::twostream::has_solar<SOURCE_TYPE>() ? 0 : 1;
    m_rust_source.emplace(sasktran2::rust::twostream::new_rust_twostream_source(
        as_slice(layer_thickness), as_slice(chapman_factors),
        m_geometry.coordinates().cos_sza_at_reference(), source_mode,
        static_cast<std::size_t>(m_config->num_threads())));

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
        num_legendre, as_slice(storage.emission_source.data(), nlevel * nwavel),
        as_slice(surface_albedo),
        as_slice(atmosphere.surface().emission().data(), nwavel),
        as_slice(storage.solar_irradiance.data(), nwavel),
        atmosphere.num_deriv() > 0);
}

template <sasktran2::twostream::SourceType SOURCE_TYPE>
void RustTwoStreamSourceAdapter<SOURCE_TYPE>::start_of_ray_source(
    int wavelidx, int losidx, int wavel_threadidx, int threadidx,
    sasktran2::Dual<double, sasktran2::dualstorage::dense, 1>& source) const {
    ZoneScopedN("Rust Twostream Cached Source");
    const std::size_t nwavel = m_atmosphere->num_wavel();
    const std::size_t nlevel = m_atmosphere->storage().total_extinction.rows();
    const std::size_t surface_index = losidx * nwavel + wavelidx;
    source.value(0) +=
        sasktran2::rust::twostream::radiance(**m_rust_source)[surface_index];

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
            view_level_offset + rust_level * nwavel + wavelidx;
        source.deriv(cpp_level) += extinction[index];
        source.deriv(m_atmosphere->ssa_deriv_start_index() + cpp_level) +=
            ssa[index];
        for (int group = 0; group < m_atmosphere->num_scattering_deriv_groups();
             ++group) {
            source.deriv(m_atmosphere->scat_deriv_start_index() +
                         group * nlevel + cpp_level) +=
                m_atmosphere->storage().d_leg_coeff(1, cpp_level, wavelidx,
                                                    group) *
                b1[index];
        }
        if constexpr (sasktran2::twostream::has_thermal<SOURCE_TYPE>()) {
            if (m_atmosphere->include_emission_derivatives()) {
                source.deriv(m_atmosphere->emission_deriv_start_index() +
                             cpp_level) += emission[index];
            }
        }
    }

    const double albedo = sasktran2::rust::twostream::surface_albedo_jacobian(
        **m_rust_source)[surface_index];
    for (int deriv = 0; deriv < m_atmosphere->surface().num_deriv(); ++deriv) {
        source.deriv(m_atmosphere->surface_deriv_start_index() + deriv) +=
            albedo;
    }
    if constexpr (sasktran2::twostream::has_thermal<SOURCE_TYPE>()) {
        if (m_atmosphere->include_emission_derivatives()) {
            source.deriv(m_atmosphere->surface_emission_deriv_start_index()) +=
                sasktran2::rust::twostream::surface_emission_jacobian(
                    **m_rust_source)[surface_index];
        }
    }
}

template class RustTwoStreamSourceAdapter<
    sasktran2::twostream::SourceType::ONLY_SOLAR>;
template class RustTwoStreamSourceAdapter<
    sasktran2::twostream::SourceType::ONLY_THERMAL>;

#endif
