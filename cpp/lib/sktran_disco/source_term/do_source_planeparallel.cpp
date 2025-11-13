#include "sasktran2/config.h"
#include "sasktran2/do_source.h"
#include "sasktran2/geometry.h"
#include "sasktran2/raytracing.h"
#include "sasktran2/viewinggeometry_internal.h"
#include "sktran_disco/sktran_do.h"
#include "sktran_disco/sktran_do_opticallayer.h"
#include "sktran_disco/sktran_do_polarization_types.h"
#include "sktran_disco/sktran_do_postprocessing.h"
#include "sktran_disco/sktran_do_types.h"

namespace sasktran2 {
    template <int NSTOKES, int CNSTR>
    DOSourcePlaneParallelPostProcessing<NSTOKES, CNSTR>::
        DOSourcePlaneParallelPostProcessing(
            const sasktran2::Geometry1D& geometry)
        : m_geometry(geometry) {}

    template <int NSTOKES, int CNSTR>
    void DOSourcePlaneParallelPostProcessing<NSTOKES, CNSTR>::calculate(
        int wavelidx, int threadidx) {
        bool include_single_scatter =
            m_config->single_scatter_source() ==
            sasktran2::Config::SingleScatterSource::discrete_ordinates;

        auto& solver = m_thread_storage[threadidx].sza_calculators[0];

        sasktran_disco::OpticalLayerArray<NSTOKES, CNSTR> optical_layer(
            *solver.persistent_config, wavelidx, m_do_los,
            *solver.geometry_layers, *m_atmosphere, *m_config);

        sasktran_disco::RTESolver<NSTOKES, CNSTR> rte(*solver.persistent_config,
                                                      optical_layer);

        auto& input_derivatives = optical_layer.inputDerivatives();

        // Have to set storage to 0 for each ray
        for (int j = 0; j < m_do_los.size(); ++j) {
            m_radiances[threadidx][j].value.setZero();

            if (m_radiances[threadidx][j].deriv.size() == 0) {
                m_radiances[threadidx][j].deriv.resize(
                    NSTOKES, m_atmosphere->num_deriv());
            }

            m_radiances[threadidx][j].deriv.setZero();
        }

        for (int j = 0; j < m_internal_viewing->flux_observers.size(); ++j) {
            m_flux[threadidx][j].resize(m_config->get_flux_types().size(),
                                        m_atmosphere->num_deriv(), true);
        }

        int num_azi = m_config->num_do_streams();

        if (m_config->num_do_forced_azimuth() > 0) {
            num_azi = m_config->num_do_forced_azimuth();
        }

        for (int m = 0; m < num_azi; ++m) {
            rte.solve(m);

            if (m == 0) {
                compute_flux(threadidx, optical_layer);
            }

            ZoneScopedN("DO PostProcessing");
            for (int j = 0; j < m_do_los.size(); ++j) {
                double observeraltitude = m_do_los[j].observeraltitude;
                double observer_opticaldepth = 0.0;

                if (observeraltitude >= 0) {
                    // We might be inside the atmosphere, so we have to
                    // adjust the LOS optical depth
                    observer_opticaldepth =
                        optical_layer.opticalDepthAt(observeraltitude);
                }

                // Start with assigning the radiance reflected from the
                // ground
                auto& ground_radiance =
                    optical_layer.reflectedIntensity(m, m_do_los[j]);

                m_component[threadidx].value = ground_radiance.value;
                m_component[threadidx].deriv = ground_radiance.deriv;

                // Calculate radiance recursively upwards through atmosphere
                for (auto layer = optical_layer.template iteratorAcross<
                                  sasktran_disco::Propagating::UP>();
                     layer.isValid(); ++layer) {
                    if (layer.entryOpticalDepth() < observer_opticaldepth) {
                        // Layer doesn't contribute
                        continue;
                    }
                    double layerfraction = 1;
                    if (layer.exitOpticalDepth() < observer_opticaldepth) {
                        layerfraction = (layer.entryOpticalDepth() -
                                         observer_opticaldepth) /
                                        (layer.entryOpticalDepth() -
                                         layer.exitOpticalDepth());
                    }

                    // Attenuation through the layer is exp(-opticaldepth *
                    // layerfraction)
                    auto& dual_optical_depth = layer.layer().dual_thickness();

                    m_component[threadidx].value *=
                        exp(-1.0 * dual_optical_depth.value * layerfraction /
                            m_do_los[j].coszenith);

                    input_derivatives.reverse_trace(j).multiply_by_constant(
                        exp(-1.0 * dual_optical_depth.value * layerfraction /
                            m_do_los[j].coszenith));

                    m_component[threadidx].deriv *=
                        exp(-1.0 * dual_optical_depth.value * layerfraction /
                            m_do_los[j].coszenith);
                    auto seq =
                        Eigen::seq(dual_optical_depth.layer_start,
                                   dual_optical_depth.layer_start +
                                       dual_optical_depth.deriv.size() - 1);
                    if constexpr (NSTOKES == 1) {
                        m_component[threadidx].deriv(seq) +=
                            -layerfraction / m_do_los[j].coszenith *
                            m_component[threadidx].value *
                            dual_optical_depth.deriv;
                    } else {
                        for (int k = 0; k < NSTOKES; ++k) {
                            m_component[threadidx].deriv(seq, k) +=
                                -layerfraction / m_do_los[j].coszenith *
                                m_component[threadidx].value(k) *
                                dual_optical_depth.deriv;
                        }
                    }

                    // Calculate and add in new source
                    layer.ptr()->integrate_source(
                        m, m_do_los[j].coszenith, observer_opticaldepth,
                        m_lp_coszen[j][m], m_integral[threadidx],
                        input_derivatives.reverse_trace(j), nullptr,
                        include_single_scatter);

                    m_component[threadidx].value += m_integral[threadidx].value;
                    m_component[threadidx].deriv += m_integral[threadidx].deriv;
                }

                if (m_config->do_backprop() &&
                    input_derivatives.numDerivative() > 0) {
                    rte.backprop(m, input_derivatives.reverse_trace(j),
                                 m_component[threadidx]);
                }

                if constexpr (NSTOKES == 1) {
                    double azimuthal_factor = cos(m * m_do_los[j].azimuth);

                    m_radiances[threadidx][j].value(0) +=
                        m_component[threadidx].value * azimuthal_factor;

                    for (int k = 0; k < m_component[threadidx].deriv.size();
                         ++k) {
                        for (int l = 0;
                             l < input_derivatives.layerDerivatives()[k]
                                     .group_and_triangle_fraction.size();
                             ++l) {
                            const std::pair<sasktran_disco::uint, double>&
                                group_fraction =
                                    input_derivatives.layerDerivatives()[k]
                                        .group_and_triangle_fraction[l];
                            const auto& extinction =
                                input_derivatives.layerDerivatives()[k]
                                    .extinctions[l];

                            m_radiances[threadidx][j].deriv(
                                group_fraction.first) +=
                                group_fraction.second *
                                m_component[threadidx].deriv(k) * extinction *
                                azimuthal_factor;
                        }
                    }

                } else {
                    for (int s = 0; s < NSTOKES; ++s) {
                        double azimuthal_factor;
                        if (s < 2) {
                            azimuthal_factor = cos(m * m_do_los[j].azimuth);
                        } else {
                            azimuthal_factor = sin(m * m_do_los[j].azimuth);
                        }
                        m_radiances[threadidx][j].value(s) +=
                            m_component[threadidx].value(s) * azimuthal_factor;

                        for (int k = 0; k < m_component[threadidx].deriv.rows();
                             ++k) {
                            for (int l = 0;
                                 l < input_derivatives.layerDerivatives()[k]
                                         .group_and_triangle_fraction.size();
                                 ++l) {
                                const std::pair<sasktran_disco::uint, double>&
                                    group_fraction =
                                        input_derivatives.layerDerivatives()[k]
                                            .group_and_triangle_fraction[l];
                                const auto& extinction =
                                    input_derivatives.layerDerivatives()[k]
                                        .extinctions[l];

                                m_radiances[threadidx][j].deriv(
                                    s, group_fraction.first) +=
                                    group_fraction.second *
                                    m_component[threadidx].deriv(k, s) *
                                    extinction * azimuthal_factor;
                            }
                        }
                    }
                }
            }
        }
    }

    template <int NSTOKES, int CNSTR>
    void DOSourcePlaneParallelPostProcessing<NSTOKES, CNSTR>::compute_flux(
        int threadidx,
        sasktran_disco::OpticalLayerArray<NSTOKES, CNSTR>& optical_layer) {
        // TODO: I think we can do this more efficiently by collapsing the
        // homogeneous solution vectors first?

        using MatrixView = Eigen::Map<Eigen::MatrixXd>;
        using ConstMatrixView = Eigen::Map<const Eigen::MatrixXd>;
        using ConstTensorView =
            Eigen::TensorMap<const Eigen::Tensor<double, 3>>; // For d_homog

        const int nstr = m_config->num_do_streams();

        const auto& input_derivatives = optical_layer.inputDerivatives();
        const int num_total_derivatives = input_derivatives.numDerivative();

        Eigen::VectorXd temp_deriv(num_total_derivatives);

        // Go through all of our flux types
        for (int flux_type_idx = 0;
             flux_type_idx < m_config->get_flux_types().size();
             ++flux_type_idx) {

            sasktran2::Config::FluxType flux_type =
                m_config->get_flux_types()[flux_type_idx];

            // And flux observers
            for (int flux_obs_idx = 0;
                 flux_obs_idx < m_internal_viewing->flux_observers.size();
                 ++flux_obs_idx) {
                auto& flux = m_flux[threadidx][flux_obs_idx];
                flux.value(flux_type_idx) = 0.0;

                temp_deriv.setZero();
                const auto& flux_observer =
                    m_internal_viewing->flux_observers[flux_obs_idx];

                // Then figure out what layer we are in
                double altitude = flux_observer.observer.radius() -
                                  m_geometry.coordinates().earth_radius();
                const sasktran_disco::OpticalLayer<NSTOKES, CNSTR>& layer =
                    *optical_layer.layerAtAltitude(altitude);
                int p = layer.index();

                const int num_layer_derivatives =
                    input_derivatives.numDerivativeLayer(p);
                const int layer_d_start = input_derivatives.layerStartIndex(p);

                double layer_fraction =
                    (layer.altitude(sasktran_disco::Location::CEILING) -
                     altitude) /
                    (layer.altitude(sasktran_disco::Location::CEILING) -
                     layer.altitude(sasktran_disco::Location::FLOOR));
                double x = layer_fraction * layer.dual_thickness().value;

                // For the layer, we need the homogeneous solutions and the BVP
                // coefficients
                const auto& solution = layer.solution(0);

                auto& Dm = layer.postprocessing_cache().Dm[0];
                auto& Dp = layer.postprocessing_cache().Dp[0];

                auto& Hp = layer.postprocessing_cache().hp[0];
                auto& Hm = layer.postprocessing_cache().hm[0];

                const auto& h_minus = solution.value.dual_homog_minus();
                const auto& h_plus = solution.value.dual_homog_plus();

                const auto& dual_Aplus = solution.value.dual_green_A_plus();
                const auto& dual_Aminus = solution.value.dual_green_A_minus();

                ConstMatrixView homog_plus_matrix(h_plus.value.data(),
                                                  nstr / 2 * NSTOKES,
                                                  nstr / 2 * NSTOKES);
                ConstMatrixView homog_minus_matrix(h_minus.value.data(),
                                                   nstr / 2 * NSTOKES,
                                                   nstr / 2 * NSTOKES);

                ConstTensorView d_homog_plus(
                    h_plus.deriv.data(), num_layer_derivatives,
                    nstr / 2 * NSTOKES, nstr / 2 * NSTOKES);

                ConstTensorView d_homog_minus(
                    h_minus.deriv.data(), num_layer_derivatives,
                    nstr / 2 * NSTOKES, nstr / 2 * NSTOKES);

                const auto& eigval = solution.value.dual_eigval();
                const auto& transmission = layer.ceiling_beam_transmittanc();
                const auto& average_secant = layer.dual_average_secant();
                const auto& dual_thickness = layer.dual_thickness();

                const auto& dual_L = solution.boundary.L_coeffs;
                const auto& dual_M = solution.boundary.M_coeffs;

                // j = SOLUTION INDEX
                for (int j = 0; j < nstr / 2 * NSTOKES; ++j) {
                    sasktran_disco::postprocessing::d_minus_sampled(
                        dual_thickness, eigval, j, layer_fraction, transmission,
                        average_secant, layer_d_start, Dm);
                    sasktran_disco::postprocessing::d_plus_sampled(
                        dual_thickness, eigval, j, layer_fraction, transmission,
                        average_secant, layer_d_start, Dp);

                    sasktran_disco::postprocessing::h_plus_sampled(
                        dual_thickness, eigval, j, layer_fraction, Hp);
                    sasktran_disco::postprocessing::h_minus_sampled(
                        dual_thickness, eigval, j, layer_fraction, Hm);

                    // l = QUADRATURE ANGLE INDEX
                    for (int l = 0; l < nstr / 2; ++l) {
                        double quadrature_weight =
                            (*m_thread_storage[threadidx]
                                  .sza_calculators[0]
                                  .persistent_config->quadrature_weights())[l] *
                            (*m_thread_storage[threadidx]
                                  .sza_calculators[0]
                                  .persistent_config
                                  ->quadrature_cos_angle())[l] *
                            EIGEN_PI * 2.0;

                        double L = dual_L.value(j);
                        double M = dual_M.value(j);

                        int h_lidx = l * NSTOKES;

                        if (flux_type ==
                            sasktran2::Config::FluxType::downwelling) {
                            double homog_plus_factor =
                                quadrature_weight *
                                (L * Hp.value + dual_Aplus.value(j) * Dm.value);

                            double homog_minus_factor =
                                quadrature_weight *
                                (M * Hm.value +
                                 dual_Aminus.value(j) * Dp.value);

                            m_flux[threadidx][flux_obs_idx].value(
                                flux_type_idx) +=
                                homog_minus_factor *
                                    homog_minus_matrix(h_lidx, j) +
                                homog_plus_factor *
                                    homog_plus_matrix(h_lidx, j);

                            // Derivatives
                            // L, M, Dm, Dp have total derivatives
                            for (int d = 0; d < num_total_derivatives; ++d) {
                                temp_deriv(d) +=
                                    quadrature_weight *
                                        homog_plus_matrix(h_lidx, j) *
                                        (dual_L.deriv(d, j) * Hp.value +
                                         Dm.deriv(d) * dual_Aplus.value(j)) +
                                    quadrature_weight *
                                        homog_minus_matrix(h_lidx, j) *
                                        (dual_M.deriv(d, j) * Hm.value +
                                         Dp.deriv(d) * dual_Aminus.value(j));
                            }
                            // And then Hm, Hp, Aplus, Aminus, and homog
                            // matrices have layer derivatives
                            for (int d = 0; d < num_layer_derivatives; ++d) {
                                temp_deriv(d + layer_d_start) +=
                                    quadrature_weight *
                                        homog_plus_matrix(h_lidx, j) *
                                        (Hp.deriv(d) * L +
                                         dual_Aplus.deriv(d, j) * Dm.value) +
                                    quadrature_weight *
                                        homog_minus_matrix(h_lidx, j) *
                                        (Hm.deriv(d) * M +
                                         dual_Aminus.deriv(d, j) * Dp.value) +
                                    homog_plus_factor *
                                        d_homog_plus(d, h_lidx, j) +
                                    homog_minus_factor *
                                        d_homog_minus(d, h_lidx, j);
                            }
                        }

                        if (flux_type ==
                            sasktran2::Config::FluxType::upwelling) {
                            double homog_minus_factor =
                                quadrature_weight *
                                (L * Hp.value + dual_Aplus.value(j) * Dm.value);

                            double homog_plus_factor =
                                quadrature_weight *
                                (M * Hm.value +
                                 dual_Aminus.value(j) * Dp.value);

                            m_flux[threadidx][flux_obs_idx].value(
                                flux_type_idx) +=
                                homog_minus_factor *
                                    homog_minus_matrix(h_lidx, j) +
                                homog_plus_factor *
                                    homog_plus_matrix(h_lidx, j);

                            // Derivatives
                            // L, M, Dm, Dp have total derivatives
                            for (int d = 0; d < num_total_derivatives; ++d) {
                                temp_deriv(d) +=
                                    quadrature_weight *
                                        homog_minus_matrix(h_lidx, j) *
                                        (dual_L.deriv(d, j) * Hp.value +
                                         Dm.deriv(d) * dual_Aplus.value(j)) +
                                    quadrature_weight *
                                        homog_plus_matrix(h_lidx, j) *
                                        (dual_M.deriv(d, j) * Hm.value +
                                         Dp.deriv(d) * dual_Aminus.value(j));
                            }
                            // And then Hm, Hp, Aplus, Aminus, and homog
                            // matrices have layer derivatives
                            for (int d = 0; d < num_layer_derivatives; ++d) {
                                temp_deriv(d + layer_d_start) +=
                                    quadrature_weight *
                                        homog_minus_matrix(h_lidx, j) *
                                        (Hp.deriv(d) * L +
                                         dual_Aplus.deriv(d, j) * Dm.value) +
                                    quadrature_weight *
                                        homog_plus_matrix(h_lidx, j) *
                                        (Hm.deriv(d) * M +
                                         dual_Aminus.deriv(d, j) * Dp.value) +
                                    homog_plus_factor *
                                        d_homog_plus(d, h_lidx, j) +
                                    homog_minus_factor *
                                        d_homog_minus(d, h_lidx, j);
                            }
                        }
                    }
                }

                // Now add in the downwelling solar flux
                if (m_config->single_scatter_source() ==
                    sasktran2::Config::SingleScatterSource::
                        discrete_ordinates) {
                    if (flux_type == sasktran2::Config::FluxType::downwelling) {
                        // The solar flux is just the transmission * cos_sza
                        double csz =
                            m_geometry.coordinates().cos_sza_at_reference();
                        m_flux[threadidx][flux_obs_idx].value(flux_type_idx) +=
                            transmission.value * csz;

                        // Derivatives
                        for (int d = 0; d < num_total_derivatives; ++d) {
                            temp_deriv(d) += transmission.deriv(d) * csz;
                        }
                    }
                }

                if (num_total_derivatives > 0) {
                    flux.deriv(flux_type_idx, Eigen::all).setZero();
                }
                for (int k = 0; k < num_total_derivatives; ++k) {
                    for (int l = 0; l < input_derivatives.layerDerivatives()[k]
                                            .group_and_triangle_fraction.size();
                         ++l) {
                        const std::pair<sasktran_disco::uint, double>&
                            group_fraction =
                                input_derivatives.layerDerivatives()[k]
                                    .group_and_triangle_fraction[l];
                        const auto& extinction =
                            input_derivatives.layerDerivatives()[k]
                                .extinctions[l];

                        flux.deriv(flux_type_idx, group_fraction.first) +=
                            group_fraction.second * temp_deriv(k) * extinction;
                    }
                }
            }
        }
    }

    template <int NSTOKES, int CNSTR>
    void
    DOSourcePlaneParallelPostProcessing<NSTOKES, CNSTR>::initialize_geometry(
        const sasktran2::viewinggeometry::InternalViewingGeometry&
            internal_viewing) {
        m_internal_viewing = &internal_viewing;

        m_do_los.resize(internal_viewing.traced_rays.size());
        m_lp_coszen.resize(internal_viewing.traced_rays.size());

        for (int i = 0; i < m_do_los.size(); ++i) {
            auto& do_los = m_do_los[i];
            const auto& ray = internal_viewing.traced_rays[i];

            do_los.coszenith = -ray.observer_and_look.look_away.z();

            if (do_los.coszenith < 0) {
                throw sasktran_disco::InternalRuntimeError(
                    "Error, currently only calculation of upwelling radiances "
                    "is supported in plane parallel mode");
            }

            do_los.azimuth = -ray.observer_and_look.relative_azimuth;

            do_los.cos_scattering_angle =
                m_geometry.coordinates().sun_unit().dot(
                    ray.observer_and_look.look_away);

            do_los.observeraltitude =
                ray.observer_and_look.observer.position.z() -
                m_geometry.coordinates().earth_radius();

            do_los.unsorted_index = i;

            for (int m = 0; m < m_nstr; ++m) {
                m_lp_coszen[i].push_back(
                    sasktran_disco::VectorDim1<
                        sasktran_disco::LegendrePhaseContainer<NSTOKES>>(
                        m_nstr));
                for (int l = 0; l < m_nstr; ++l) {
                    m_lp_coszen[i][m][l].fill(m, l, m_do_los[i].coszenith);
                }
            }
        }

        for (int thidx = 0; thidx < m_thread_storage.size(); ++thidx) {
            auto& sza_calculators = m_thread_storage[thidx].sza_calculators;
            for (int i = 0; i < sza_calculators.size(); ++i) {
                auto& sza_calculator = sza_calculators[i];
                double cos_sza =
                    m_geometry.coordinates().cos_sza_at_reference();
                sza_calculator.persistent_config->configure(
                    sza_calculator.userspec, *m_config, cos_sza,
                    (int)m_geometry.altitude_grid().grid().size() - 1,
                    internal_viewing.traced_rays);

                sza_calculator.geometry_layers = std::make_unique<
                    sasktran_disco::GeometryLayerArray<NSTOKES, CNSTR>>(
                    *sza_calculator.persistent_config, m_geometry);
            }
            // These quantities aren't needed for each SZA
            m_thread_storage[thidx].postprocessing_cache.resize(
                m_geometry.altitude_grid().grid().size() - 1);
        }
    }

    template <int NSTOKES, int CNSTR>
    void
    DOSourcePlaneParallelPostProcessing<NSTOKES, CNSTR>::initialize_atmosphere(
        const sasktran2::atmosphere::Atmosphere<NSTOKES>& atmosphere) {
        m_atmosphere = &atmosphere;

        m_component.resize(m_thread_storage.size(), m_atmosphere->num_deriv());
        m_radiances.resize(
            m_thread_storage.size(),
            std::vector<sasktran2::Dual<double, sasktran2::dualstorage::dense,
                                        NSTOKES>>(
                m_do_los.size(),
                sasktran2::Dual<double, sasktran2::dualstorage::dense, NSTOKES>(
                    NSTOKES, m_atmosphere->num_deriv())));

        m_flux.resize(
            m_thread_storage.size(),
            std::vector<
                sasktran2::Dual<double, sasktran2::dualstorage::dense, -1>>(
                m_internal_viewing->flux_observers.size(),
                sasktran2::Dual<double, sasktran2::dualstorage::dense, -1>(
                    m_config->get_flux_types().size(),
                    m_atmosphere->num_deriv())));

        m_integral.resize(m_thread_storage.size(), m_atmosphere->num_deriv());
    }

    template <int NSTOKES, int CNSTR>
    void DOSourcePlaneParallelPostProcessing<NSTOKES, CNSTR>::initialize_config(
        const sasktran2::Config& config) {
        m_config = &config;

        // Create the thread storage
        m_thread_storage.resize(config.num_threads());

        m_nstr = config.num_do_streams();

        // Always only one SZA for plane parallel
        int num_sza = 1;

        for (int i = 0; i < m_thread_storage.size(); ++i) {
            auto& sza_calculators = m_thread_storage[i].sza_calculators;
            sza_calculators.resize(num_sza);
            for (auto& calculator : sza_calculators) {
                calculator.persistent_config = std::make_unique<
                    sasktran_disco::PersistentConfiguration<NSTOKES, CNSTR>>();
            }
        }
    }

    template <int NSTOKES, int CNSTR>
    void
    DOSourcePlaneParallelPostProcessing<NSTOKES, CNSTR>::start_of_ray_source(
        int wavelidx, int losidx, int wavel_threadidx, int threadidx,
        sasktran2::Dual<double, sasktran2::dualstorage::dense, NSTOKES>& source)
        const {
        source.value += m_radiances[threadidx][losidx].value;
        source.deriv += m_radiances[threadidx][losidx].deriv;
    }

    template <int NSTOKES, int CNSTR>
    void DOSourcePlaneParallelPostProcessing<NSTOKES, CNSTR>::flux(
        int wavelidx, int fluxidx, int wavel_threadidx, int threadidx,
        sasktran2::Dual<double, sasktran2::dualstorage::dense, 1>& flux,
        sasktran2::Config::FluxType flux_type) const {
        flux.value.setConstant(
            m_flux[threadidx][fluxidx].value(static_cast<int>(flux_type)));

        flux.deriv =
            m_flux[threadidx][fluxidx].deriv.row(static_cast<int>(flux_type));
    }

    SASKTRAN_DISCO_INSTANTIATE_TEMPLATE(DOSourcePlaneParallelPostProcessing);
} // namespace sasktran2
