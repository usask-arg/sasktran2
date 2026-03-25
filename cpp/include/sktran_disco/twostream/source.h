#pragma once

#include "meta.h"
#include "sasktran2/config.h"
#include "sasktran2/geometry.h"
#include "sasktran2/math/scattering.h"
#include "sasktran2/source_interface.h"
#include "sasktran2/viewinggeometry_internal.h"
#include "sktran_disco/sktran_do_geometrylayerarray.h"
#include "sktran_disco/sktran_do_pconfig.h"
#include "sktran_disco/sktran_do_specs.h"
#include "sktran_disco/twostream/backprop.h"
#include "storage.h"
#include "solutions.h"
#include <memory>

template <int NSTOKES, sasktran2::twostream::SourceType SOURCE_TYPE>
class TwoStreamSource : public SourceTermInterface<NSTOKES> {
  private:
    mutable std::vector<sasktran2::twostream::Solution<SOURCE_TYPE>>
        m_solutions;
    std::vector<sasktran2::twostream::Input<SOURCE_TYPE>> m_inputs;
    mutable std::vector<sasktran2::twostream::Sources<SOURCE_TYPE>> m_sources;
    mutable std::vector<std::array<
        Eigen::MatrixXd, sasktran2::twostream::num_azimuth<SOURCE_TYPE>()>>
        m_bvp_backprop_storage;

    const sasktran2::Geometry1D& m_geometry;
    const sasktran2::Config* m_config;
    sasktran_disco::PersistentConfiguration<1> m_pconfig;
    sasktran_disco::SKTRAN_DO_UserSpec m_spec;

    std::unique_ptr<sasktran_disco::GeometryLayerArray<1>> m_geometry_layers;

    const std::vector<sasktran2::raytracing::TracedRay>* m_los_rays;

    std::vector<Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic,
                              Eigen::RowMajorBit>>
        m_los_attenuation_factors;
    mutable std::vector<Eigen::VectorXd> m_internal_gradients;

  public:
    TwoStreamSource(const sasktran2::Geometry1D& geometry)
        : m_geometry(geometry) {
        m_spec.configure(2, geometry.size() - 1);
    };

    virtual bool requires_integration() const override { return false; }

    virtual void initialize_config(const sasktran2::Config& config) {
        m_solutions.resize(config.num_threads());
        m_inputs.resize(config.num_threads());
        m_sources.resize(config.num_threads());
        m_bvp_backprop_storage.resize(config.num_threads());
        m_internal_gradients.resize(config.num_threads());

        m_config = &config;
    };

    virtual void initialize_geometry(
        const sasktran2::viewinggeometry::InternalViewingGeometry&
            internal_viewing) override {
        ZoneScopedN("Twostream Initialize Geometry");
        m_pconfig.configure(
            m_spec, *m_config, m_geometry.coordinates().cos_sza_at_reference(),
            m_geometry.size() - 1, internal_viewing.traced_rays);
        m_los_rays = &internal_viewing.traced_rays;

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

        for (auto& backprop : m_bvp_backprop_storage) {
            for (auto& matrix : backprop) {
                matrix.resize(2 * (m_geometry.size() - 1), 1);
            }
        }
        m_los_attenuation_factors.resize(m_los_rays->size());

        for (int i = 0; i < (*m_los_rays).size(); ++i) {
            const auto& ray = (*m_los_rays)[i];
            Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic,
                          Eigen::RowMajorBit>& atten_matrix =
                m_los_attenuation_factors[i];
            double viewing_secant = 1 / (-ray.observer_and_look.look_away.z());

            atten_matrix.resize(m_geometry.size() - 1, m_geometry.size() - 1);
            atten_matrix.setZero();

            // layer 0 = TOA, has no attenuation
            for (int j = 1; j < m_geometry.size() - 1; ++j) {
                atten_matrix(j, Eigen::seq(0, j - 1))
                    .setConstant(viewing_secant);
            }
        }

        for (auto& internal_gradient : m_internal_gradients) {
            int nlyr = m_geometry.size() - 1;
            if constexpr (sasktran2::twostream::has_thermal<SOURCE_TYPE>()) {
                // +2 for albedo, thermal surface
                internal_gradient.resize(5 * nlyr + 2);
            } else {
                // +1 for albedo
                internal_gradient.resize(3 * nlyr + 1);
            }
        }
    };

    virtual void initialize_atmosphere(
        const sasktran2::atmosphere::Atmosphere<NSTOKES>& atmosphere) {
        for (auto& input : m_inputs) {
            input.atmosphere = &atmosphere;
        }
    };

    virtual void calculate(int wavelidx, int threadidx) {
        ZoneScopedN("Twostream Calculate");
        auto& input = m_inputs[threadidx];
        auto& solution = m_solutions[threadidx];

        {
            ZoneScopedN("Twostream Input Calcultaion");
            input.calculate(wavelidx);
        }

        {
            ZoneScopedN("Twostream Solution");
            sasktran2::twostream::solve(input, solution);
        }
    };

    virtual void integrated_source(
        int wavelidx, int losidx, int layeridx, int wavel_threadidx,
        int threadidx, const sasktran2::raytracing::SphericalLayer& layer,
        const sasktran2::SparseODDualView& shell_od,
        sasktran2::Dual<double, sasktran2::dualstorage::dense, NSTOKES>& source,
        typename SourceTermInterface<NSTOKES>::IntegrationDirection direction =
            SourceTermInterface<NSTOKES>::IntegrationDirection::none) const {};

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
        ZoneScopedN("Twostream Start of Ray Source");
        auto& sources = m_sources[threadidx];
        const auto& ray = (*m_los_rays)[losidx];
        auto& solution = m_solutions[threadidx];
        const auto& input = m_inputs[threadidx];
        auto& internal_gradient = m_internal_gradients[threadidx];

        internal_gradient.setZero();

        double viewing_zenith = -ray.observer_and_look.look_away.z();

        if (viewing_zenith < 0) {
            throw sasktran_disco::InternalRuntimeError(
                "Error, currently only calculation of upwelling radiances "
                "is supported in plane parallel mode");
        }

        double azimuth = -ray.observer_and_look.relative_azimuth;

        sources.final_weight_factors.noalias() =
            (-m_los_attenuation_factors[losidx] * input.od);
        sources.final_weight_factors =
            sources.final_weight_factors.array().exp();

        sasktran2::twostream::post_process(input, viewing_zenith, azimuth,
                                           solution, sources);

        double integrated_source =
            (solution.particular[0].G_plus_bottom.value(
                 Eigen::placeholders::last) +
             solution.bvp_coeffs[0].rhs(Eigen::placeholders::last - 1, 0) *
                 solution.homog[0].X_plus.value(Eigen::placeholders::last) *
                 solution.homog[0].omega.value(Eigen::placeholders::last) +
             solution.bvp_coeffs[0].rhs(Eigen::placeholders::last, 0) *
                 solution.homog[0].X_minus.value(Eigen::placeholders::last)) *
            2 * input.mu * input.albedo *
            sources.final_weight_factors(Eigen::placeholders::last) *
            sources.beamtrans.value(Eigen::placeholders::last);

        if constexpr (sasktran2::twostream::has_thermal<SOURCE_TYPE>()) {
            integrated_source +=
                input.thermal_surf *
                sources.final_weight_factors(Eigen::placeholders::last) *
                sources.beamtrans.value(Eigen::placeholders::last);
        }

        double ground_source = integrated_source;

        integrated_source +=
            sources.source.value.dot(sources.final_weight_factors);

        source.value(0) += integrated_source;

        if (input.atmosphere->num_deriv() > 0) {
            double ground_weight =
                2 * input.mu * input.albedo *
                sources.final_weight_factors(Eigen::placeholders::last) *
                sources.beamtrans.value(Eigen::placeholders::last);

            sasktran2::twostream::backprop::GradientMap<SOURCE_TYPE> grad(
                *input.atmosphere, internal_gradient.data());

            // assign d_albedo
            grad.d_albedo =
                (solution.particular[0].G_plus_bottom.value(
                     Eigen::placeholders::last) +
                 solution.bvp_coeffs[0].rhs(Eigen::placeholders::last - 1, 0) *
                     solution.homog[0].X_plus.value(Eigen::placeholders::last) *
                     solution.homog[0].omega.value(Eigen::placeholders::last) +
                 solution.bvp_coeffs[0].rhs(Eigen::placeholders::last, 0) *
                     solution.homog[0].X_minus.value(
                         Eigen::placeholders::last)) *
                2 * input.mu *
                sources.final_weight_factors(Eigen::placeholders::last) *
                sources.beamtrans.value(Eigen::placeholders::last);

            if constexpr (sasktran2::twostream::has_thermal<SOURCE_TYPE>()) {
                grad.d_thermal_surf =
                    sources.final_weight_factors(Eigen::placeholders::last) *
                    sources.beamtrans.value(Eigen::placeholders::last);
            }

            sasktran2::twostream::backprop::full(
                input, solution, sources, sources.final_weight_factors,
                m_bvp_backprop_storage[threadidx], grad, ground_weight);

            {
                ZoneScopedN("Twostream Backprop Attenuation");
                // Backprop the attenuation factors
                for (int i = 0; i < input.nlyr; ++i) {
                    sasktran2::twostream::backprop::od(
                        input,
                        (-sources.final_weight_factors(i) *
                         sources.source.value(i) *
                         m_los_attenuation_factors[losidx].row(i)),
                        grad);
                }
            }

            // Backprop ground source LOS attenuation
            sasktran2::twostream::backprop::od(
                input,
                (-ground_source * Eigen::RowVectorXd::Ones(input.nlyr) /
                 viewing_zenith),
                grad);

            // Map the derivatives to the atmosphere
            sasktran2::twostream::backprop::map_to_atmosphere(input, grad,
                                                              source);
        }
    };
};
