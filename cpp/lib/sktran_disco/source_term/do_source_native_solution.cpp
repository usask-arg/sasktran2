#include "sasktran2/do_source.h"
#include "sktran_disco/sktran_do.h"
#include <memory>

namespace sasktran2 {
    template <int NSTOKES, int CNSTR>
    DOSourceNativeSolution<NSTOKES, CNSTR>::
        DOSourceNativeSolution(
            const sasktran2::Geometry1D& geometry,
            const sasktran2::raytracing::RayTracerBase& raytracer)
        : DOSource<NSTOKES, CNSTR>(geometry, raytracer) {}

    template <int NSTOKES, int CNSTR>
    void DOSourceNativeSolution<NSTOKES, CNSTR>::
        accumulate_solved_azimuth(
            sasktran_disco::OpticalLayerArray<NSTOKES, CNSTR>& optical_layer,
            DOSourceThreadStorage<NSTOKES, CNSTR>& thread_storage, int szaidx,
            sasktran_disco::AEOrder m, int threadidx) {
        m_radiance_storage->accumulate_sources(optical_layer, m, thread_storage,
                                              szaidx, threadidx);
    }

    template <int NSTOKES, int CNSTR>
    void DOSourceNativeSolution<NSTOKES, CNSTR>::calculate(
        int wavelidx, int threadidx) {
        DOSource<NSTOKES, CNSTR>::calculate(wavelidx, threadidx);
    }

    template <int NSTOKES, int CNSTR>
    void
    DOSourceNativeSolution<NSTOKES, CNSTR>::initialize_geometry(
        const std::vector<sasktran2::raytracing::TracedRay>& los_rays) {
        DOSource<NSTOKES, CNSTR>::initialize_geometry(los_rays);

        m_radiance_storage =
            std::make_unique<sasktran2::DORadianceStorage<NSTOKES, CNSTR>>(
                *(this->m_thread_storage)[0]
                     .sza_calculators[0]
                     .geometry_layers.get(),
                *(this->m_thread_storage)[0]
                     .sza_calculators[0]
                     .persistent_config.get(),
                *this->m_sza_grid, *this->m_config, this->m_geometry);
    }

    template <int NSTOKES, int CNSTR>
    void
    DOSourceNativeSolution<NSTOKES, CNSTR>::initialize_atmosphere(
        const sasktran2::atmosphere::Atmosphere<NSTOKES>& atmosphere) {
        m_atmosphere = &atmosphere;

        DOSource<NSTOKES, CNSTR>::initialize_atmosphere(atmosphere);

        m_radiance_storage->initialize_atmosphere(atmosphere);
    }

    template <int NSTOKES, int CNSTR>
    void DOSourceNativeSolution<NSTOKES, CNSTR>::initialize_config(
        const sasktran2::Config& config) {
        DOSource<NSTOKES, CNSTR>::initialize_config(config);
    }

    SASKTRAN_DISCO_INSTANTIATE_TEMPLATE(DOSourceNativeSolution);
}