#pragma once

#ifdef SKTRAN_RUST_SUPPORT

#include "sasktran2-core/src/twostream/cxx.rs.h"
#include "sasktran2/config.h"
#include "sasktran2/geometry.h"
#include "sasktran2/source_interface.h"
#include "sasktran2/viewinggeometry_internal.h"
#include "sktran_disco/sktran_do.h"
#include "sktran_disco/twostream/meta.h"
#include <optional>

template <sasktran2::twostream::SourceType SOURCE_TYPE>
class RustTwoStreamSourceAdapter final : public SourceTermInterface<1> {
  private:
    const sasktran2::Geometry1D& m_geometry;
    const sasktran2::Config* m_config = nullptr;
    const sasktran2::atmosphere::Atmosphere<1>* m_atmosphere = nullptr;
    sasktran_disco::SKTRAN_DO_UserSpec m_spec;
    std::vector<std::unique_ptr<sasktran_disco::PersistentConfiguration<1>>>
        m_pconfigs;
    std::vector<std::unique_ptr<sasktran_disco::GeometryLayerArray<1>>>
        m_geometry_layers;
    const std::vector<sasktran2::raytracing::TracedRay>* m_los_rays = nullptr;
    std::optional<::rust::Box<sasktran2::rust::twostream::RustTwoStreamSource>>
        m_rust_source;

  public:
    explicit RustTwoStreamSourceAdapter(const sasktran2::Geometry1D& geometry);

    bool requires_integration() const override { return false; }
    bool has_interior_source() const override { return false; }
    void initialize_config(const sasktran2::Config& config) override;
    void initialize_geometry(
        const sasktran2::viewinggeometry::InternalViewingGeometry&
            internal_viewing) override;
    void initialize_atmosphere(
        const sasktran2::atmosphere::Atmosphere<1>& atmosphere) override;
    void calculate(int wavelidx, int threadidx) override {}

    void integrated_source(
        int wavelidx, int losidx, int layeridx, int wavel_threadidx,
        int threadidx, const sasktran2::raytracing::TracedLayer& layer,
        const sasktran2::raytracing::GridWeightStencilView& entrance_weights,
        const sasktran2::raytracing::GridWeightStencilView& exit_weights,
        const sasktran2::SparseODDualView& shell_od,
        sasktran2::Dual<double, sasktran2::dualstorage::dense, 1>& source,
        IntegrationDirection direction =
            IntegrationDirection::none) const override {}

    void end_of_ray_source(
        int wavelidx, int losidx, int wavel_threadidx, int threadidx,
        sasktran2::Dual<double, sasktran2::dualstorage::dense, 1>& source)
        const override {}

    void start_of_ray_source(
        int wavelidx, int losidx, int wavel_threadidx, int threadidx,
        sasktran2::Dual<double, sasktran2::dualstorage::dense, 1>& source)
        const override;
};

#endif
