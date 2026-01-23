#include "sktran_disco/sktran_do.h"
#include "sktran_disco/sktran_do_layerarray.h"
#include "sktran_disco/sktran_do_types.h"

template <int NSTOKES, int CNSTR>
void sasktran_disco::OpticalLayerArray<NSTOKES, CNSTR>::
    computeReflectedIntensities(AEOrder m,
                                const sasktran_disco::LineOfSight& los) {
    // Here we evaluate the upwelling term at the surface.  Note that the sign
    // convention switches compared to the standard sign convention.  +<->- on W
    // and Z

    auto& trace = m_input_derivatives.reverse_trace(los.unsorted_index);

    m_reflection_computed[m][los.unsorted_index] = true;
    if (m >= m_surface.sk2_surface().max_azimuthal_order()) {
        return;
    }
    const uint N = this->M_NSTR / 2;

    const auto& layer = bottom();
    const auto& solution = layer.solution(m);
    const auto& input_deriv = m_input_derivatives;
    int layerStart = (int)input_deriv.layerStartIndex(layer.index());
    int numLayerDeriv = (int)input_deriv.numDerivativeLayer(layer.index());

    int l_offset = layer.index() * this->M_NSTR * NSTOKES;
    int m_offset = l_offset + this->M_NSTR * NSTOKES / 2;

    Radiance<NSTOKES> diffuse_contrib((int)input_deriv.numDerivative()),
        direct_contrib((int)input_deriv.numDerivative());
    diffuse_contrib.setzero();
    diffuse_contrib.setzero();

    direct_contrib.setzero();
    direct_contrib.setzero();

    const auto& rho = m_surface.storage().brdf.los_stream;

    // Construct the Dual quantity for the layer transmittance
    Dual<double> beam_transmittance =
        layer.dual_beamTransmittance(Location::FLOOR, input_deriv);
    Dual<double> stream_transmittance;

    for (StreamIndex i = 0; i < N * NSTOKES; ++i) {
        LayerDual<double> dual_rho(
            (uint)input_deriv.numDerivativeLayer(layer.index()), layer.index(),
            (uint)input_deriv.layerStartIndex(layer.index()));
        // TODO: something is definitely wrong here, but it might be fine for
        // scalar surface reflection...
        int s1 = i % NSTOKES;
        if (i % NSTOKES != 0) {
            dual_rho.value = 0.0;
            dual_rho.deriv.setZero();
        } else {
            dual_rho.value = rho(los.unsorted_index, i / NSTOKES);

            for (uint j = 0; j < input_deriv.numDerivativeLayer(layer.index());
                 ++j) {
                const auto& d_rho =
                    m_surface.storage()
                        .d_brdf[input_deriv
                                    .layerDerivatives()[input_deriv
                                                            .layerStartIndex(
                                                                layer.index()) +
                                                        j]
                                    .surface_deriv_index]
                        .los_stream;
                dual_rho.deriv(j) =
                    input_deriv
                        .layerDerivatives()
                            [input_deriv.layerStartIndex(layer.index()) + j]
                        .d_albedo *
                    d_rho(los.unsorted_index, i / NSTOKES);
            }
        }

        Radiance<NSTOKES> stream_contrib((int)input_deriv.numDerivative());
        if constexpr (NSTOKES == 1) {
            stream_contrib.value = solution.value.dual_Gplus_bottom().value(i);
        } else {
            stream_contrib.value(s1) =
                solution.value.dual_Gplus_bottom().value(i);
        }

        stream_contrib.deriv(Eigen::all, s1) =
            solution.value.dual_Gplus_bottom().deriv(Eigen::all, i);

        double factor = (1.0 + kronDelta(m, 0)) * (*this->M_MU)[i / NSTOKES] *
                        (*this->M_WT)[i / NSTOKES];

        // Positive homogeneous solutions
        for (uint j = 0; j < N * NSTOKES; ++j) {
            stream_transmittance = layer.dual_streamTransmittance(
                Location::INSIDE, m, j, input_deriv);
            uint homogIndex = j * N * NSTOKES + i;

            if constexpr (NSTOKES == 1) {
                stream_contrib.value +=
                    solution.boundary.L_coeffs.value(j) *
                    solution.value.dual_homog_plus().value(homogIndex) *
                    stream_transmittance.value;
            } else {
                stream_contrib.value(s1) +=
                    solution.boundary.L_coeffs.value(j) *
                    solution.value.dual_homog_plus().value(homogIndex) *
                    stream_transmittance.value;
            }

            if (this->M_BACKPROP_BVP) {
                if (s1 == 0) {
                    trace.bvp_coeff_weights()(l_offset + j, s1) +=
                        solution.value.dual_homog_plus().value(homogIndex) *
                        stream_transmittance.value * factor * dual_rho.value;
                }
            } else {
                // LCoeffs have full derivatives
                for (uint k = 0; k < input_deriv.numDerivative(); ++k) {
                    stream_contrib.deriv(k, s1) +=
                        solution.boundary.L_coeffs.deriv(k, j) *
                        solution.value.dual_homog_plus().value(homogIndex) *
                        stream_transmittance.value;
                }
            }

            // Homog only have layer derivs
            for (int k = 0; k < numLayerDeriv; ++k) {
                stream_contrib.deriv(k + layerStart, s1) +=
                    solution.boundary.L_coeffs.value(j) *
                    solution.value.dual_homog_plus().deriv(k, homogIndex) *
                    stream_transmittance.value;
                stream_contrib.deriv(k + layerStart, s1) +=
                    solution.boundary.L_coeffs.value(j) *
                    solution.value.dual_homog_plus().value(homogIndex) *
                    stream_transmittance.deriv(k + layerStart);
            }

            if constexpr (NSTOKES == 1) {
                stream_contrib.value +=
                    solution.boundary.M_coeffs.value(j) *
                    solution.value.dual_homog_minus().value(homogIndex);
            } else {
                stream_contrib.value(s1) +=
                    solution.boundary.M_coeffs.value(j) *
                    solution.value.dual_homog_minus().value(homogIndex);
            }

            if (this->M_BACKPROP_BVP) {
                if (s1 == 0) {
                    trace.bvp_coeff_weights()(m_offset + j, s1) +=
                        solution.value.dual_homog_minus().value(homogIndex) *
                        factor * dual_rho.value;
                }
            } else {
                // MCoeffs have full derivatives
                for (uint k = 0; k < input_deriv.numDerivative(); ++k) {
                    stream_contrib.deriv(k, s1) +=
                        solution.boundary.M_coeffs.deriv(k, j) *
                        solution.value.dual_homog_minus().value(homogIndex);
                }
            }
            // Homog only have layer derivs
            for (int k = 0; k < numLayerDeriv; ++k) {
                stream_contrib.deriv(k + layerStart, s1) +=
                    solution.boundary.M_coeffs.value(j) *
                    solution.value.dual_homog_minus().deriv(k, homogIndex);
            }
        }
        // Add stream i contribution

        if constexpr (NSTOKES == 1) {
            diffuse_contrib.value +=
                factor * stream_contrib.value * dual_rho.value;
        } else {
            diffuse_contrib.value(s1) +=
                factor * stream_contrib.value(s1) * dual_rho.value;
        }

        // stream_contrib has full derivatives
        for (uint k = 0; k < input_deriv.numDerivative(); ++k) {
            diffuse_contrib.deriv(k, s1) +=
                factor * stream_contrib.deriv(k, s1) * dual_rho.value;
        }
        // rho only have layer derivs
        for (int k = 0; k < numLayerDeriv; ++k) {
            if constexpr (NSTOKES == 1) {
                diffuse_contrib.deriv(k + layerStart) +=
                    factor * stream_contrib.value * dual_rho.deriv(k);
            } else {
                diffuse_contrib.deriv(k + layerStart, s1) +=
                    factor * stream_contrib.value(s1) * dual_rho.deriv(k);
            }
        }
    }
    LayerDual<double> dual_albedo_sun(
        (uint)input_deriv.numDerivativeLayer(layer.index()), layer.index(),
        (uint)input_deriv.layerStartIndex(layer.index()));

    dual_albedo_sun.value =
        m_surface.storage().brdf.los_solar(los.unsorted_index, 0);

    for (uint j = 0; j < input_deriv.numDerivativeLayer(layer.index()); ++j) {
        double d_sun =
            m_surface.storage()
                .d_brdf[input_deriv
                            .layerDerivatives()
                                [input_deriv.layerStartIndex(layer.index()) + j]
                            .surface_deriv_index]
                .los_solar(los.unsorted_index, 0);

        dual_albedo_sun.deriv(j) =
            input_deriv
                .layerDerivatives()[input_deriv.layerStartIndex(layer.index()) +
                                    j]
                .d_albedo *
            d_sun;
    }

    if (m_include_direct_bounce) {
        if constexpr (NSTOKES == 1) {
            direct_contrib.value = this->M_CSZ / PI * beam_transmittance.value *
                                   dual_albedo_sun.value;
        } else {
            direct_contrib.value(0) = this->M_CSZ / PI *
                                      beam_transmittance.value *
                                      dual_albedo_sun.value;
        }

        // beam has full deriv
        for (uint k = 0; k < input_deriv.numDerivative(); ++k) {
            direct_contrib.deriv(k, 0) += this->M_CSZ / PI *
                                          beam_transmittance.deriv(k) *
                                          dual_albedo_sun.value;
        }
        // albedo has layer deriv
        for (int k = 0; k < numLayerDeriv; ++k) {
            direct_contrib.deriv(k + layerStart, 0) +=
                this->M_CSZ / PI * beam_transmittance.value *
                dual_albedo_sun.deriv(k);
        }

        if (m_config.ss_only()) {
            m_ground_reflection[m][los.unsorted_index].value =
                direct_contrib.value;
            m_ground_reflection[m][los.unsorted_index].deriv =
                direct_contrib.deriv;
        } else {
            m_ground_reflection[m][los.unsorted_index].value =
                direct_contrib.value + diffuse_contrib.value;
            m_ground_reflection[m][los.unsorted_index].deriv =
                direct_contrib.deriv + diffuse_contrib.deriv;
        }

        // Thermal source
        if (m == 0) {
            double thermal_emission =
                m_surface.sk2_surface().emission()[m_wavel_index];

            if constexpr (NSTOKES == 1) {
                m_ground_reflection[m][los.unsorted_index].value +=
                    thermal_emission;
            } else {
                m_ground_reflection[m][los.unsorted_index].value(0) +=
                    thermal_emission;
            }

            // Derivatives wrt to surface emission parameters
            for (int k = 0; k < numLayerDeriv; ++k) {
                double d_thermal_emission =
                    1.0 * input_deriv
                              .layerDerivatives()[input_deriv.layerStartIndex(
                                                      layer.index()) +
                                                  k]
                              .d_surface_emission;

                m_ground_reflection[m][los.unsorted_index].deriv(
                    k + layerStart) += d_thermal_emission;
            }
        }

    } else {
        m_ground_reflection[m][los.unsorted_index].value =
            diffuse_contrib.value;
        m_ground_reflection[m][los.unsorted_index].deriv =
            diffuse_contrib.deriv;
    }
}

template <int NSTOKES, int CNSTR>
sasktran_disco::OpticalLayerArray<NSTOKES, CNSTR>::OpticalLayerArray(
    const PersistentConfiguration<NSTOKES, CNSTR>& config, int wavelidx,
    const std::vector<LineOfSight>& los,
    const GeometryLayerArray<NSTOKES, CNSTR>& geometry_layers,
    const sasktran2::atmosphere::Atmosphere<NSTOKES>& atmosphere,
    const sasktran2::Config& sk_config)
    : OpticalLayerArrayROP<NSTOKES>(config), m_config(config),
      m_include_direct_bounce(
          sk_config.single_scatter_source() ==
          sasktran2::Config::SingleScatterSource::discrete_ordinates),
      m_num_los(los.size()),
      m_chapman_factors(geometry_layers.chapman_factors()),
      m_optical_interpolator(geometry_layers.interpolating_matrix()),
      m_input_derivatives(config.pool().thread_data().input_derivatives()),
      m_transmission(config.pool().thread_data().transmission()),
      m_surface(config.pool().thread_data().surface_storage(),
                atmosphere.surface(), wavelidx) {

    bool include_thermal_emission =
        sk_config.emission_source() ==
        sasktran2::Config::EmissionSource::discrete_ordinates;

    m_wavel_index = wavelidx;
    m_direct_toa = atmosphere.storage().solar_irradiance(wavelidx);
    // Allocations
    m_layers.reserve(this->M_NLYR);
    m_surface.storage().resize(this->M_NSTR, m_num_los,
                               atmosphere.surface().num_deriv());
    m_surface.storage().csz = this->M_CSZ;
    m_surface.storage().mu = this->M_MU;

    m_surface.storage().los = &los;

    // Accumulation quantities for the layers
    double ceiling_depth = 0;
    double floor_depth = 0;

    const auto& d = atmosphere.storage().leg_coeff.dimensions();
    Eigen::Map<const Eigen::MatrixXd> phase(
        &atmosphere.storage().leg_coeff(0, 0, wavelidx), d[0], d[1]);

    for (LayerIndex p = 0; p < this->M_NLYR; ++p) {
        double ceil_h = geometry_layers.layer_ceiling()(p);
        double floor_h = geometry_layers.layer_floor()(p);
        double layer_dh = ceil_h - floor_h;

        double od = 0.0;
        double ssa = 0.0;
        double f = 0.0;

        // thermal emission: S(x) = b0 * exp(-b1 * x)
        // where x = optical depth from top of layer
        // b0 = emission at top (ceiling), b1 = ln(B_top/B_bot) / layer_od
        // Note: q index runs from bottom (q=0) to top of atmosphere
        // so higher q = higher altitude = layer ceiling
        int min_q = -1; // bottom of layer (lower altitude, lower q)
        int max_q = -1; // top of layer (higher altitude, higher q)
        for (int q = 0; q < m_optical_interpolator.cols(); ++q) {
            if (m_optical_interpolator(p, q) > 0) {
                if (min_q == -1 || q < min_q) {
                    min_q = q;
                }
                if (max_q == -1 || q > max_q) {
                    max_q = q;
                }
            }
        }

        // Get emission at top (ceiling) and bottom (floor) of layer
        double b0_top = 0.0;
        double b0_bot = 0.0;
        if (max_q >= 0) {
            b0_top =
                atmosphere.storage().emission_source(max_q, wavelidx); // top
        }
        if (min_q >= 0 && min_q != max_q) {
            b0_bot =
                atmosphere.storage().emission_source(min_q, wavelidx); // bottom
        } else {
            b0_bot = b0_top; // Single point layer
        }

        // Copy the legendre coefficients to a new vector
        std::unique_ptr<
            VectorDim1<sasktran_disco::LegendreCoefficient<NSTOKES>>>
            lephasef(
                new VectorDim1<sasktran_disco::LegendreCoefficient<NSTOKES>>);

        lephasef->resize(this->M_NSTR);

        for (int q = 0; q < m_optical_interpolator.cols(); ++q) {
            if (m_optical_interpolator(p, q) > 0) {
                double weight = m_optical_interpolator(p, q);

                double kext =
                    atmosphere.storage().total_extinction(q, wavelidx);
                double kscat = atmosphere.storage().ssa(q, wavelidx) * kext;

                // Store extinction in od during weighting
                od += kext * weight;

                // Store scattering extinction in ssa
                ssa += kscat * weight;

                // Numerator of f
                f = atmosphere.storage().f(q, wavelidx);

                for (uint k = 0; k < this->M_NSTR; ++k) {
                    auto& temp = (*lephasef)[k];

                    if constexpr (NSTOKES == 1) {
                        double test = phase(k, q);
                        temp.a1 += weight * kscat *
                                   (phase(k, q) - (2 * k + 1) * f / (1 - f));
                    } else if constexpr (NSTOKES == 3) {
                        auto stokes_seq = Eigen::seq(k * 4, (k + 1) * 4 - 1);

                        Eigen::Vector<double, 4> p_bot = phase(stokes_seq, q);

                        // a's have the f subtraction and b's don't
                        temp.a1 += weight * kscat *
                                   (p_bot(0) - f * (2 * k + 1) / (1 - f));
                        temp.a2 += weight * kscat *
                                   (p_bot(1) - f * (2 * k + 1) / (1 - f));
                        temp.a3 += weight * kscat *
                                   (p_bot(2) - f * (2 * k + 1) / (1 - f));

                        temp.b1 += -weight * kscat * p_bot(3);
                    } else {
                    }
                }
            }
        }

        if (ssa > 0) {
            // Divide the phase elements by the scattering extinction
            for (uint k = 0; k < this->M_NSTR; ++k) {
                auto& temp = (*lephasef)[k];

                if constexpr (NSTOKES == 1) {
                    temp.a1 /= ssa;
                } else if constexpr (NSTOKES == 3) {
                    temp.a1 /= ssa;
                    temp.a2 /= ssa;
                    temp.a3 /= ssa;
                    temp.b1 /= ssa;
                } else {
                }
            }
        } else {
            // If there is no scattering, we should set a1 to be 1
            auto& temp = (*lephasef)[0];
            temp.a1 = 0;
        }

        // Convert ssa from scattering extinction to ssa
        ssa /= od;

        // Then convert od to optical depth
        od *= layer_dh;

        floor_depth += od;

        double total_ext = od / (ceil_h - floor_h);
        double scat_ext = total_ext * ssa;
        scat_ext = std::max(scat_ext,
                            total_ext * this->m_userspec->getSSAEqual1Dither());

        // Calculate thermal emission coefficients
        // S(x) = b0 * exp(-b1 * x) where x is optical depth from layer top
        // At x=0: S(0) = b0 = b0_top (emission at ceiling)
        // At x=od: S(od) = b0 * exp(-b1 * od) = b0_bot
        // Therefore: b1 = ln(b0_top / b0_bot) / od
        double b0 = b0_top;
        double b1 = 0.0;
        if (od > 1e-10 && b0_top > 1e-30 && b0_bot > 1e-30 &&
            std::abs(b0_top - b0_bot) > 1e-15 * std::max(b0_top, b0_bot)) {
            b1 = std::log(b0_top / b0_bot) / od;
        }

        m_layers.push_back(std::unique_ptr<OpticalLayer<NSTOKES, CNSTR>>(
            new OpticalLayer<NSTOKES, CNSTR>(
                config, p, scat_ext, total_ext, std::move(lephasef),
                ceiling_depth, floor_depth, ceil_h, floor_h,
                m_input_derivatives, include_thermal_emission, b0, b1)));
        ceiling_depth = floor_depth;
    }

    m_input_derivatives.set_num_los(m_num_los, this->M_NSTR, this->M_NLYR);
    if (atmosphere.num_deriv() > 0 &&
        sk_config.wf_precision() ==
            sasktran2::Config::WeightingFunctionPrecision::full) {
        int numderiv = atmosphere.num_deriv();
        int num_atmo_grid = (int)atmosphere.storage().total_extinction.rows();
        int num_scattering_groups = atmosphere.num_scattering_deriv_groups();

        if (!m_input_derivatives.is_geometry_configured()) {
            const Eigen::MatrixXd& atmosphere_mapping = m_optical_interpolator;

            // Now avg_k in layers = atmosphere_mapping @ atmosphere_extinction
            for (LayerIndex p = 0; p < this->M_NLYR; ++p) {
                auto& ptrb_layer = layer(p);

                // Go through the atmosphere mapping and add non-zero elements
                // for SSA/ext
                for (int scatidx = 0; scatidx < num_scattering_groups;
                     ++scatidx) {
                    LayerInputDerivative<NSTOKES>& deriv_scat =
                        m_input_derivatives.addDerivative(this->M_NSTR, p);

                    for (int i = 0; i < num_atmo_grid; ++i) {
                        if (atmosphere_mapping(p, i) > 0) {
                            // How d atmosphere legendre influences layer
                            // legendre
                            deriv_scat.group_and_triangle_fraction.emplace_back(
                                atmosphere.scat_deriv_start_index() +
                                    num_atmo_grid * scatidx + i,
                                atmosphere_mapping(p, i));
                            deriv_scat.extinctions.emplace_back(1);
                        }
                    }
                }

                LayerInputDerivative<NSTOKES>& deriv_ext =
                    m_input_derivatives.addDerivative(this->M_NSTR, p);

                deriv_ext.d_optical_depth = 1;
                // Go through the atmosphere mapping and add non-zero elements
                // for SSA/ext
                for (int i = 0; i < num_atmo_grid; ++i) {
                    if (atmosphere_mapping(p, i) > 0) {

                        // How d atmosphere k influences layer od
                        deriv_ext.group_and_triangle_fraction.emplace_back(
                            i, atmosphere_mapping(p, i) *
                                   (ptrb_layer.altitude(Location::CEILING) -
                                    ptrb_layer.altitude(Location::FLOOR)));
                        deriv_ext.extinctions.emplace_back(1);
                    }
                }

                // if we are doing thermal emission each layer needs a
                // derivative for b0 and b1
                if (atmosphere.include_emission_derivatives()) {
                    // Find the min and max q indices for this layer
                    // q=0 is bottom of atmosphere, higher q = higher altitude
                    int min_q_deriv = -1;
                    int max_q_deriv = -1;
                    for (int i = 0; i < num_atmo_grid; ++i) {
                        if (atmosphere_mapping(p, i) > 0) {
                            if (min_q_deriv == -1 || i < min_q_deriv) {
                                min_q_deriv = i;
                            }
                            if (max_q_deriv == -1 || i > max_q_deriv) {
                                max_q_deriv = i;
                            }
                        }
                    }

                    // b0 = emission at top of layer (max_q)
                    // d_b0/d_emission_top = 1
                    LayerInputDerivative<NSTOKES>& deriv_b0 =
                        m_input_derivatives.addDerivative(this->M_NSTR, p);
                    deriv_b0.d_thermal_b0 = 1;
                    if (max_q_deriv >= 0) {
                        // Only link to the top (ceiling) emission source
                        deriv_b0.group_and_triangle_fraction.emplace_back(
                            atmosphere.emission_deriv_start_index() +
                                max_q_deriv,
                            1.0);
                        deriv_b0.extinctions.emplace_back(1);
                    }

                    // b1 = ln(b0_top / b0_bot) / od
                    // d_b1/d_emission_top = 1 / (b0_top * od)
                    // d_b1/d_emission_bot = -1 / (b0_bot * od)
                    // d_b1/d_od = -b1 / od = -ln(b0_top/b0_bot) / od^2
                    // Note: We store the partial derivative factors in
                    // extinctions and apply them during the derivative
                    // calculation phase
                    LayerInputDerivative<NSTOKES>& deriv_b1 =
                        m_input_derivatives.addDerivative(this->M_NSTR, p);
                    deriv_b1.d_thermal_b1 = 1;

                    // Link to top emission source (for d_b1/d_emission_top)
                    if (max_q_deriv >= 0) {
                        deriv_b1.group_and_triangle_fraction.emplace_back(
                            atmosphere.emission_deriv_start_index() +
                                max_q_deriv,
                            1.0); // Weight will be computed in reassignment
                                  // phase
                        deriv_b1.extinctions.emplace_back(1);
                    }
                    // Link to bottom emission source (for d_b1/d_emission_bot)
                    if (min_q_deriv >= 0 && min_q_deriv != max_q_deriv) {
                        deriv_b1.group_and_triangle_fraction.emplace_back(
                            atmosphere.emission_deriv_start_index() +
                                min_q_deriv,
                            1.0); // Weight will be computed in reassignment
                                  // phase
                        deriv_b1.extinctions.emplace_back(1);
                    }
                    // Link to optical depth (for d_b1/d_od)
                    // We need to add contributions from all atmosphere grid
                    // points that affect the layer optical depth
                    for (int i = 0; i < num_atmo_grid; ++i) {
                        if (atmosphere_mapping(p, i) > 0) {
                            deriv_b1.group_and_triangle_fraction.emplace_back(
                                i, atmosphere_mapping(p, i) *
                                       (ptrb_layer.altitude(Location::CEILING) -
                                        ptrb_layer.altitude(Location::FLOOR)));
                            deriv_b1.extinctions.emplace_back(1);
                        }
                    }
                }

                LayerInputDerivative<NSTOKES>& deriv_ssa =
                    m_input_derivatives.addDerivative(this->M_NSTR, p);
                deriv_ssa.d_SSA = 1;
                // Go through the atmosphere mapping and add non-zero elements
                // for SSA/ext
                for (int i = 0; i < num_atmo_grid; ++i) {
                    if (atmosphere_mapping(p, i) > 0) {
                        // How d atmosphere ssa influences layer ssa
                        deriv_ssa.group_and_triangle_fraction.emplace_back(
                            i + num_atmo_grid, atmosphere_mapping(p, i));
                        deriv_ssa.extinctions.emplace_back(1);

                        // How d atmosphere k influences layer ssa
                        deriv_ssa.group_and_triangle_fraction.emplace_back(
                            i, atmosphere_mapping(p, i));
                        deriv_ssa.extinctions.emplace_back(1);
                    }
                }
            }
            // Add one final derivative for the surface
            for (int k = 0; k < m_surface.sk2_surface().num_deriv(); ++k) {
                LayerInputDerivative<NSTOKES>& deriv_albedo =
                    m_input_derivatives.addDerivative(this->M_NSTR,
                                                      this->M_NLYR - 1);
                deriv_albedo.d_albedo = 1;
                deriv_albedo.surface_deriv_index = k;
                deriv_albedo.group_and_triangle_fraction.emplace_back(
                    atmosphere.surface_deriv_start_index() + k, 1);
                deriv_albedo.extinctions.emplace_back(1);
            }

            // And if we are including emission, one more derivative for the
            // emission surface
            if (atmosphere.include_emission_derivatives()) {
                LayerInputDerivative<NSTOKES>& deriv_albedo =
                    m_input_derivatives.addDerivative(this->M_NSTR,
                                                      this->M_NLYR - 1);
                deriv_albedo.d_surface_emission = 1;
                deriv_albedo.surface_deriv_index = 0;
                deriv_albedo.group_and_triangle_fraction.emplace_back(
                    atmosphere.surface_emission_deriv_start_index(), 1);
                deriv_albedo.extinctions.emplace_back(1);
            }

            m_input_derivatives.set_geometry_configured();
            m_input_derivatives.sort(this->M_NLYR);

            for (Dual<double>& trans : m_transmission) {
                trans.deriv.resize(m_input_derivatives.numDerivative());
            }
        }

        // Go through the derivatives and reassign a few things
        for (int k = 0; k < m_input_derivatives.numDerivative(); ++k) {
            LayerInputDerivative<NSTOKES>& deriv =
                m_input_derivatives.layerDerivatives()[k];
            const OpticalLayer<NSTOKES, CNSTR>* layer =
                m_layers[deriv.layer_index].get();

            // Have to check what kind of derivative this is
            if (deriv.d_optical_depth > 0) {
                // OD derivative, don't have to do anything here
            } else if (deriv.d_thermal_b0 > 0) {
                // b0 derivative - b0 = emission at top
                // d_b0/d_emission_top = 1, already set correctly
            } else if (deriv.d_thermal_b1 > 0) {
                // b1 = ln(b0_top / b0_bot) / od
                // d_b1/d_emission_top = 1 / (b0_top * od)
                // d_b1/d_emission_bot = -1 / (b0_bot * od)
                // d_b1/d_od = -b1 / od

                double b0_val = layer->dual_thermal_b0().value;
                double b1_val = layer->dual_thermal_b1().value;
                double od_val =
                    layer->opticalDepth(sasktran_disco::Location::INSIDE);

                // The group_and_triangle_fraction entries are:
                // [0]: top emission source (if exists)
                // [1]: bottom emission source (if different from top)
                // [2+]: optical depth contributions
                int emission_deriv_start =
                    atmosphere.emission_deriv_start_index();
                int idx = 0;
                for (int l = 0; l < deriv.group_and_triangle_fraction.size();
                     ++l) {
                    auto& group = deriv.group_and_triangle_fraction[l];
                    auto& extinction = deriv.extinctions[l];

                    if (group.first >= emission_deriv_start &&
                        group.first < emission_deriv_start + num_atmo_grid) {
                        // This is an emission source derivative
                        int atmo_idx = group.first - emission_deriv_start;
                        double emission_val =
                            atmosphere.storage().emission_source(atmo_idx,
                                                                 wavelidx);
                        // Check if this is the top or bottom
                        // Top has b0_val = emission_val (approximately)
                        if (emission_val > 0 && od_val > 1e-10) {
                            if (idx == 0) {
                                // First emission source is top
                                // d_b1/d_emission_top = 1/(b0_top * od)
                                extinction = 1.0 / (emission_val * od_val);
                            } else {
                                // Second emission source is bottom
                                // d_b1/d_emission_bot = -1/(b0_bot * od)
                                extinction = -1.0 / (emission_val * od_val);
                            }
                        } else {
                            extinction = 0.0;
                        }
                        idx++;
                    } else {
                        // This is an optical depth contribution
                        // d_b1/d_od = -b1 / od (already scaled by atmosphere
                        // mapping in group.second)
                        if (od_val > 1e-10) {
                            extinction = -b1_val / od_val;
                        } else {
                            extinction = 0.0;
                        }
                    }
                }
            } else if (deriv.d_SSA > 0) {

                // SSA derivative

                for (int l = 0; l < deriv.group_and_triangle_fraction.size();
                     ++l) {
                    auto& group = deriv.group_and_triangle_fraction[l];
                    auto& extinction = deriv.extinctions[l];

                    if (group.first >= num_atmo_grid) {
                        // Layer SSA contribution to dI/d atmosphere SSA, these
                        // are equal to the relative fraction of extinction
                        // contributions

                        int atmo_index = group.first - num_atmo_grid;

                        extinction = atmosphere.storage().total_extinction(
                                         atmo_index, wavelidx) /
                                     layer->totalExt();
                    } else {
                        // This is layer SSA contribution to dI/d atmosphere k,
                        int atmo_index = group.first;

                        extinction =
                            (atmosphere.storage().ssa(atmo_index, wavelidx) -
                             layer->ssa()) /
                            layer->totalExt();
                    }
                }

            } else if (deriv.d_albedo == 0 && deriv.d_surface_emission == 0) {
                // Scattering derivative
                for (int l = 0; l < deriv.group_and_triangle_fraction.size();
                     ++l) {
                    auto& group = deriv.group_and_triangle_fraction[l];
                    auto& extinction = deriv.extinctions[l];
                    int atmo_index = group.first % num_atmo_grid;

                    if (group.first >= 2 * num_atmo_grid) {
                        // Layer legendre contribution to dI / d atmosphere
                        // legendre
                        int group_index = (group.first - 2 * num_atmo_grid) /
                                          int(num_atmo_grid);

                        double f = atmosphere.storage().f(atmo_index, wavelidx);
                        extinction =
                            atmosphere.storage().ssa(atmo_index, wavelidx) *
                            atmosphere.storage().total_extinction(atmo_index,
                                                                  wavelidx) /
                            layer->scatExt();

                        if (atmosphere.storage().d_leg_coeff.size() > 0) {
                            for (int l = 0; l < (int)this->M_NSTR; ++l) {
                                if constexpr (NSTOKES == 1) {
                                    deriv.d_legendre_coeff[l].a1 =
                                        atmosphere.storage().d_leg_coeff(
                                            l, atmo_index, wavelidx,
                                            group_index) +
                                        (atmosphere.storage().leg_coeff(
                                             l, atmo_index, wavelidx) -
                                         (2 * l + 1) * f / (1 - f) -
                                         layer->legendre_coeff()[l].a1);

                                    if (atmosphere.storage().applied_f_order >
                                        0) {
                                        deriv.d_legendre_coeff[l].a1 +=
                                            -(2 * l + 1) / (1 - f) / (1 - f) *
                                            atmosphere.storage().d_f(
                                                atmo_index, wavelidx,
                                                group_index);
                                    }
                                }

                                if constexpr (NSTOKES == 3) {
                                    deriv.d_legendre_coeff[l].a1 =
                                        atmosphere.storage().d_leg_coeff(
                                            l * 4, atmo_index, wavelidx,
                                            group_index) +
                                        (atmosphere.storage().leg_coeff(
                                             l * 4, atmo_index, wavelidx) -
                                         (2 * l + 1) * f / (1 - f) -
                                         layer->legendre_coeff()[l].a1);

                                    deriv.d_legendre_coeff[l].a2 =
                                        atmosphere.storage().d_leg_coeff(
                                            l * 4 + 1, atmo_index, wavelidx,
                                            group_index) +
                                        (atmosphere.storage().leg_coeff(
                                             l * 4 + 1, atmo_index, wavelidx) -
                                         (2 * l + 1) * f / (1 - f) -
                                         layer->legendre_coeff()[l].a2);

                                    deriv.d_legendre_coeff[l].a3 =
                                        atmosphere.storage().d_leg_coeff(
                                            l * 4 + 2, atmo_index, wavelidx,
                                            group_index) +
                                        (atmosphere.storage().leg_coeff(
                                             l * 4 + 2, atmo_index, wavelidx) -
                                         (2 * l + 1) * f / (1 - f) -
                                         layer->legendre_coeff()[l].a3);

                                    if (atmosphere.storage().applied_f_order >
                                        0) {
                                        deriv.d_legendre_coeff[l].a1 +=
                                            -(2 * l + 1) / (1 - f) / (1 - f) *
                                            atmosphere.storage().d_f(
                                                atmo_index, wavelidx,
                                                group_index);

                                        deriv.d_legendre_coeff[l].a2 +=
                                            -(2 * l + 1) / (1 - f) / (1 - f) *
                                            atmosphere.storage().d_f(
                                                atmo_index, wavelidx,
                                                group_index);

                                        deriv.d_legendre_coeff[l].a3 +=
                                            -(2 * l + 1) / (1 - f) / (1 - f) *
                                            atmosphere.storage().d_f(
                                                atmo_index, wavelidx,
                                                group_index);
                                    }

                                    deriv.d_legendre_coeff[l].b1 =
                                        (-1) * atmosphere.storage().d_leg_coeff(
                                                   l * 4 + 3, atmo_index,
                                                   wavelidx, group_index) +
                                        ((-1) * atmosphere.storage().leg_coeff(
                                                    l * 4 + 3, atmo_index,
                                                    wavelidx) -
                                         layer->legendre_coeff()[l].b1);
                                }
                            }
                        }

                    } else if (group.first >= num_atmo_grid) {
                    } else {
                    }
                }
            }
        }
    }

    // Post configure the layers, PS beam transmittances and derivative
    // calculations
    for (auto& layer : m_layers) {
        // Start by telling each layer to calculate the derivative of thickness
        layer->configureDerivative();
    }

    configureTransmission();

    // resize ground reflection vector
    m_ground_reflection.resize(
        this->M_NSTR,
        std::vector<Radiance<NSTOKES>>(
            static_cast<uint>(los.size()),
            Radiance<NSTOKES>((int)m_input_derivatives.numDerivative())));
    m_reflection_computed.resize(this->M_NSTR,
                                 std::vector<bool>(los.size(), false));
}

template <int NSTOKES, int CNSTR>
void sasktran_disco::OpticalLayerArray<NSTOKES,
                                       CNSTR>::configureTransmission() {
    // This is a weird function where it is called so many times that we
    // basically have to have two separate paths for when derivates are on/off

    bool compute_deriv = m_input_derivatives.numDerivative() > 0;

    // Start at TOA and work down
    // Start by storing OD
    m_transmission[0].value = 0;
    m_transmission[0].deriv.setZero();

    int index = 1;
    for (auto& layer : m_layers) {
        m_transmission[index].value = 0.0;
        m_transmission[index].deriv.setZero();

        if (m_chapman_factors(layer->index(), 0) < 0) {
            // Transmission is 0,
            m_transmission[index].value = -1.0;
            layer->dual_average_secant().value = 1.0 / this->M_CSZ;
            if (compute_deriv) {
                layer->dual_average_secant().deriv =
                    Eigen::VectorXd::Zero(m_input_derivatives.numDerivative());
            }
            ++index;
            continue;
        }

        for (LayerIndex p = 0; p < m_layers.size(); ++p) {
            const auto& dual_thickness = m_layers[p]->dual_thickness();

            m_transmission[index].value +=
                m_chapman_factors(layer->index(), p) * dual_thickness.value;

            if (compute_deriv) {
                const auto seq =
                    Eigen::seq(dual_thickness.layer_start,
                               dual_thickness.layer_start +
                                   dual_thickness.deriv.size() - 1);
                m_transmission[index].deriv(seq) +=
                    m_chapman_factors(layer->index(), p) * dual_thickness.deriv;
            }
        }
        layer->dual_average_secant().value =
            (m_transmission[index].value - m_transmission[index - 1].value) /
            layer->dual_thickness().value;

        if (compute_deriv) {
            layer->dual_average_secant().deriv =
                (m_transmission[index].deriv -
                 m_transmission[index - 1].deriv) /
                layer->dual_thickness().value;

            const auto seq =
                Eigen::seq(layer->dual_thickness().layer_start,
                           layer->dual_thickness().layer_start +
                               layer->dual_thickness().deriv.size() - 1);
            layer->dual_average_secant().deriv(seq) -=
                layer->dual_thickness().deriv / layer->dual_thickness().value *
                layer->dual_average_secant().value;
        }
        ++index;
    }

    // Now convert OD to transmission
    m_transmission[0].value = directIntensityTOA();
    index = 1;
    for (auto& layer : m_layers) {
        if (m_transmission[index].value < 0) {
            // Transmission is 0
            m_transmission[index].value = 0.0;
            layer->dual_average_secant().value = 1e10;
        } else {
            m_transmission[index].value =
                std::exp(-m_transmission[index].value) * directIntensityTOA();
        }

        if (compute_deriv) {
            m_transmission[index].deriv =
                -m_transmission[index].deriv * m_transmission[index].value;
        }

        layer->set_transmittances(m_transmission[index - 1],
                                  m_transmission[index]);

        ++index;
    }
}

SASKTRAN_DISCO_INSTANTIATE_TEMPLATE(sasktran_disco::OpticalLayerArray);

template class sasktran_disco::OpticalLayerArrayIterator<
    sasktran_disco::Propagating::UP, 1>;
template class sasktran_disco::OpticalLayerArrayIterator<
    sasktran_disco::Propagating::DOWN, 1>;

template class sasktran_disco::OpticalLayerArrayIterator<
    sasktran_disco::Propagating::UP, 3>;
template class sasktran_disco::OpticalLayerArrayIterator<
    sasktran_disco::Propagating::DOWN, 3>;

template class sasktran_disco::OpticalLayerArrayIterator<
    sasktran_disco::Propagating::UP, 4>;
template class sasktran_disco::OpticalLayerArrayIterator<
    sasktran_disco::Propagating::DOWN, 4>;

template class sasktran_disco::OpticalLayerArrayIterator<
    sasktran_disco::Propagating::UP, 1, 2>;
template class sasktran_disco::OpticalLayerArrayIterator<
    sasktran_disco::Propagating::DOWN, 1, 2>;

template class sasktran_disco::OpticalLayerArrayIterator<
    sasktran_disco::Propagating::UP, 3, 2>;
template class sasktran_disco::OpticalLayerArrayIterator<
    sasktran_disco::Propagating::DOWN, 3, 2>;

template class sasktran_disco::OpticalLayerArrayIterator<
    sasktran_disco::Propagating::UP, 4, 2>;
template class sasktran_disco::OpticalLayerArrayIterator<
    sasktran_disco::Propagating::DOWN, 4, 2>;

template class sasktran_disco::OpticalLayerArrayIterator<
    sasktran_disco::Propagating::UP, 1, 4>;
template class sasktran_disco::OpticalLayerArrayIterator<
    sasktran_disco::Propagating::DOWN, 1, 4>;

template class sasktran_disco::OpticalLayerArrayIterator<
    sasktran_disco::Propagating::UP, 3, 4>;
template class sasktran_disco::OpticalLayerArrayIterator<
    sasktran_disco::Propagating::DOWN, 3, 4>;

template class sasktran_disco::OpticalLayerArrayIterator<
    sasktran_disco::Propagating::UP, 4, 4>;
template class sasktran_disco::OpticalLayerArrayIterator<
    sasktran_disco::Propagating::DOWN, 4, 4>;

template class sasktran_disco::OpticalLayerArrayIterator<
    sasktran_disco::Propagating::UP, 1, 16>;
template class sasktran_disco::OpticalLayerArrayIterator<
    sasktran_disco::Propagating::DOWN, 1, 16>;

template class sasktran_disco::OpticalLayerArrayIterator<
    sasktran_disco::Propagating::UP, 3, 16>;
template class sasktran_disco::OpticalLayerArrayIterator<
    sasktran_disco::Propagating::DOWN, 3, 16>;

template class sasktran_disco::OpticalLayerArrayIterator<
    sasktran_disco::Propagating::UP, 4, 16>;
template class sasktran_disco::OpticalLayerArrayIterator<
    sasktran_disco::Propagating::DOWN, 4, 16>;
