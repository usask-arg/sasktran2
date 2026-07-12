#include "sasktran2/viewinggeometry_internal.h"
#include <sasktran2/solartransmission.h>
#include <sasktran2/dual.h>
#include <sasktran2/config.h>
#include <sasktran2/atmosphere/atmosphere.h>

namespace sasktran2::solartransmission {
    template <int NSTOKES>
    void OccultationSource<NSTOKES>::initialize_config(
        const sasktran2::Config& config) {}

    template <int NSTOKES>
    void OccultationSource<NSTOKES>::initialize_geometry(
        const sasktran2::viewinggeometry::InternalViewingGeometry&
            internal_viewing) {
        m_ground_is_hit.resize(internal_viewing.traced_rays.size());
        std::transform(internal_viewing.traced_rays.begin(),
                       internal_viewing.traced_rays.end(),
                       m_ground_is_hit.begin(),
                       [](const auto& ray) { return ray.ground_is_hit; });
    }

    template <int NSTOKES>
    void OccultationSource<NSTOKES>::initialize_geometry(
        const std::vector<sasktran2::raytracing::TracedRay2D>& traced_rays) {
        m_ground_is_hit.resize(traced_rays.size());
        std::transform(traced_rays.begin(), traced_rays.end(),
                       m_ground_is_hit.begin(),
                       [](const auto& ray) { return ray.ground_is_hit; });
    }

    template <int NSTOKES>
    void OccultationSource<NSTOKES>::integrated_source(
        int wavelidx, int losidx, int layeridx, int wavel_threadidx,
        int threadidx, const sasktran2::raytracing::SphericalLayer& layer,
        const sasktran2::SparseODDualView& shell_od,
        sasktran2::Dual<double, sasktran2::dualstorage::dense, NSTOKES>& source,
        typename SourceTermInterface<NSTOKES>::IntegrationDirection direction)
        const {}

    template <int NSTOKES>
    void OccultationSource<NSTOKES>::calculate(int wavelidx, int threadidx) {}

    template <int NSTOKES>
    void OccultationSource<NSTOKES>::initialize_atmosphere(
        const sasktran2::atmosphere::Atmosphere<NSTOKES>& atmosphere) {}

    template <int NSTOKES>
    void OccultationSource<NSTOKES>::end_of_ray_source(
        int wavelidx, int losidx, int wavel_threadidx, int threadidx,
        sasktran2::Dual<double, sasktran2::dualstorage::dense, NSTOKES>& source)
        const {
        if (m_ground_is_hit.at(losidx)) {
            return;
        }
        if constexpr (NSTOKES == 1) {
            source.value.array() += 1;
        } else {
            source.value(0) += 1;
        }
    }

    template class OccultationSource<1>;
    template class OccultationSource<3>;

} // namespace sasktran2::solartransmission
