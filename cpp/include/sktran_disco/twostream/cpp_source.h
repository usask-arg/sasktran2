#pragma once

#include "sasktran2/config.h"
#include "sasktran2/geometry.h"
#include "sasktran2/source_interface.h"
#include "sasktran2/viewinggeometry_internal.h"
#include "sktran_disco/twostream/meta.h"

#include <limits>
#include <memory>

/**
 * Wavelength-batched C++ two-stream source.
 *
 * The numerical implementation lives out of line so including the source in
 * the engine does not instantiate the Eigen packet kernel in every translation
 * unit.  Independent wavelengths are stored in fixed-width Eigen arrays; Eigen
 * therefore selects the native SIMD instruction set at compile time.
 */
template <sasktran2::twostream::SourceType SOURCE_TYPE>
class CppTwoStreamSourceAdapter final : public SourceTermInterface<1> {
  private:
    struct Impl;
    std::unique_ptr<Impl> m_impl;

  public:
    explicit CppTwoStreamSourceAdapter(const sasktran2::Geometry1D& geometry);
    ~CppTwoStreamSourceAdapter() override;

    CppTwoStreamSourceAdapter(const CppTwoStreamSourceAdapter&) = delete;
    CppTwoStreamSourceAdapter&
    operator=(const CppTwoStreamSourceAdapter&) = delete;

    bool requires_integration() const override { return false; }
    bool has_interior_source() const override { return false; }

    void initialize_config(const sasktran2::Config& config) override;
    void initialize_geometry(
        const sasktran2::viewinggeometry::InternalViewingGeometry&
            internal_viewing) override;
    void initialize_atmosphere(
        const sasktran2::atmosphere::Atmosphere<1>& atmosphere) override;
    void set_wavelength_block_capacity(int block_capacity) override;
    int maximum_wavelength_block_size() const override {
        return std::numeric_limits<int>::max();
    }
    void calculate(const sasktran2::WavelengthBlock<>& block,
                   int threadidx) override;

    void integrated_source(
        const sasktran2::WavelengthBlock<>&, int, int, int, int,
        const sasktran2::raytracing::TracedLayer&,
        const sasktran2::raytracing::GridWeightStencilView&,
        const sasktran2::raytracing::GridWeightStencilView&,
        const sasktran2::WavelengthBlockODView&,
        sasktran2::WavelengthBlockDual<1>&,
        IntegrationDirection = IntegrationDirection::none) const override {}

    void end_of_ray_source(const sasktran2::WavelengthBlock<>&, int, int, int,
                           sasktran2::WavelengthBlockDual<1>&) const override {}

    void start_of_ray_source(
        const sasktran2::WavelengthBlock<>& block, int losidx,
        int wavel_threadidx, int threadidx,
        sasktran2::WavelengthBlockDual<1>& source) const override;
};
