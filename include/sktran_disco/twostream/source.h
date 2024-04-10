#pragma once

#include "sasktran2/config.h"
#include "sasktran2/geometry.h"
#include "sasktran2/math/scattering.h"
#include "sasktran2/source_interface.h"
#include "sktran_disco/sktran_do_geometrylayerarray.h"
#include "sktran_disco/sktran_do_pconfig.h"
#include "sktran_disco/sktran_do_specs.h"
#include "storage.h"
#include "solutions.h"
#include <memory>

template <int NSTOKES>
class TwoStreamSource : public SourceTermInterface<NSTOKES> {
  private:
    std::vector<sasktran2::twostream::Solution> m_solutions;
    std::vector<sasktran2::twostream::Input> m_inputs;
    mutable std::vector<sasktran2::twostream::Sources> m_sources;

    const sasktran2::Geometry1D& m_geometry;
    const sasktran2::Config* m_config;
    sasktran_disco::PersistentConfiguration<1> m_pconfig;
    sasktran_disco::SKTRAN_DO_UserSpec m_spec;

    std::unique_ptr<sasktran_disco::GeometryLayerArray<1>> m_geometry_layers;

    const std::vector<sasktran2::raytracing::TracedRay>* m_los_rays;

  public:
    TwoStreamSource(const sasktran2::Geometry1D& geometry)
        : m_geometry(geometry) {
        m_spec.configure(2, geometry.size() - 1);
    };

    virtual void initialize_config(const sasktran2::Config& config) {
        m_solutions.resize(config.num_threads());
        m_inputs.resize(config.num_threads());
        m_sources.resize(config.num_threads());

        m_config = &config;
    };

    virtual void initialize_geometry(
        const std::vector<sasktran2::raytracing::TracedRay>& los_rays) {
        m_pconfig.configure(m_spec, *m_config,
                            m_geometry.coordinates().cos_sza_at_reference(),
                            m_geometry.size() - 1, los_rays);
        m_los_rays = &los_rays;

        m_geometry_layers =
            std::make_unique<sasktran_disco::GeometryLayerArray<1>>(m_pconfig,
                                                                    m_geometry);

        for (auto& input : m_inputs) {
            input.geometry_layers = m_geometry_layers.get();
            input.init(m_geometry.size() - 1);

            input.csz = m_geometry.coordinates().cos_sza_at_reference();
            input.mu = 0.5;

            input.init(m_geometry.size() - 1);
        }

        for (auto& solution : m_solutions) {
            solution.init(m_geometry.size() - 1);
        }

        for (auto& source : m_sources) {
            source.init(m_geometry.size() - 1);
        }
    };

    virtual void initialize_atmosphere(
        const sasktran2::atmosphere::Atmosphere<NSTOKES>& atmosphere) {
        for (auto& input : m_inputs) {
            input.atmosphere = &atmosphere;
        }
    };

    virtual void calculate(int wavelidx, int threadidx) {
        auto& input = m_inputs[threadidx];
        auto& solution = m_solutions[threadidx];

        input.calculate(wavelidx);

        sasktran2::twostream::solve(input, solution);
    };

    virtual void integrated_source(
        int wavelidx, int losidx, int layeridx, int wavel_threadidx,
        int threadidx, const sasktran2::raytracing::SphericalLayer& layer,
        const sasktran2::SparseODDualView& shell_od,
        sasktran2::Dual<double, sasktran2::dualstorage::dense, NSTOKES>& source)
        const {};

    virtual void end_of_ray_source(
        int wavelidx, int losidx, int wavel_threadidx, int threadidx,
        sasktran2::Dual<double, sasktran2::dualstorage::dense, NSTOKES>& source)
        const {

        };

    /** Calculates the radiance at the start of the ray, i.e., the source term
     * has done the equivalent of the integration along the ray.  This is useful
     * if the source term has it's own way of performing integration that is
     * different than the standard method used in the model.  It can also be
     * used if the source term calculates quantities that can be used to get the
     * radiance directly instead of doing integration.
     *
     *  Typically source terms will only either implement start_of_ray_source,
     * or integrated_source + end_of_ray_source and not both
     *
     * @param wavelidx
     * @param losidx
     * @param wavel_threadidx
     * @param threadidx
     * @param source
     */
    virtual void start_of_ray_source(
        int wavelidx, int losidx, int wavel_threadidx, int threadidx,
        sasktran2::Dual<double, sasktran2::dualstorage::dense, NSTOKES>& source)
        const override {
        auto& sources = m_sources[threadidx];
        const auto& ray = (*m_los_rays)[losidx];
        const auto& solution = m_solutions[threadidx];
        const auto& input = m_inputs[threadidx];

        double viewing_zenith = -ray.observer_and_look.look_away.z();

        if (viewing_zenith < 0) {
            throw sasktran_disco::InternalRuntimeError(
                "Error, currently only calculation of upwelling radiances "
                "is supported in plane parallel mode");
        }

        double azimuth = -ray.observer_and_look.relative_azimuth;

        sources.final_weight_factors.setConstant(1.0);
        for (int i = 0; i < input.nlyr - 1; ++i) {
            sources.final_weight_factors(Eigen::seq(i + 1, Eigen::last)) *=
                sources.beamtrans.value(i);
        }

        sasktran2::twostream::post_process(input, viewing_zenith, azimuth,
                                           solution, sources);

        double integrated_source =
            (solution.particular[0].G_plus_bottom.value(Eigen::last) +
             solution.bvp_coeffs[0].rhs(Eigen::last - 1, 0) *
                 solution.homog[0].X_plus.value(Eigen::last) *
                 solution.homog[0].omega.value(Eigen::last) +
             solution.bvp_coeffs[0].rhs(Eigen::last, 0) *
                 solution.homog[0].X_minus.value(Eigen::last)) *
            2 * input.mu * input.albedo *
            sources.final_weight_factors(Eigen::last) *
            sources.beamtrans.value(Eigen::last);

        // TODO: Add albedo related backprop terms

        integrated_source +=
            sources.source.value.dot(sources.final_weight_factors);

        source.value(0) += integrated_source;
    };
};
