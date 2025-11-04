
    template <int NSTOKES, int CNSTR>
    void DORadianceStorage<NSTOKES, CNSTR>::accumulate_sources(
        sasktran_disco::OpticalLayerArray<NSTOKES, CNSTR>& optical_layer,
        sasktran_disco::AEOrder m,
        sasktran2::DOSourceThreadStorage<NSTOKES, CNSTR>& thread_storage,
        int szaidx, int thread_idx) {
        using MatrixView = Eigen::Map<Eigen::MatrixXd>;
        using ConstMatrixView = Eigen::Map<const Eigen::MatrixXd>;
        using ConstTensorView = Eigen::TensorMap<const Eigen::Tensor<double, 3>>; // For d_homog

        auto& storage = m_storage[thread_idx];
        if (m == 0) {
            m_storage[thread_idx].radiances.value.setZero();
            m_storage[thread_idx].radiances.deriv.setZero();
        }

        const int nstr = m_config.num_do_streams();
        const int num_azi = nstr;

        const auto& input_derivatives = optical_layer.inputDerivatives();
        const int num_total_derivatives = input_derivatives.numDerivative();

        // TODO: Move to cache
        Eigen::Matrix<double, -1, NSTOKES> temp_deriv_downw(
            input_derivatives.numDerivative(), NSTOKES);

        Eigen::Matrix<double, -1, NSTOKES> temp_deriv_upw(
            input_derivatives.numDerivative(), NSTOKES);

        // lidx = layer index
        for (int lidx = 0; lidx < m_altitude_grid->grid().size(); ++lidx) {
            // Pull out the layer and where we are inside the layer
            auto& layer = optical_layer.layer(
                (int)m_altitude_grid->grid().size() - lidx - 1);
            sasktran_disco::LayerIndex p = layer.index();
            double altitude = m_altitude_grid->grid()(lidx);

            const int num_layer_derivatives = input_derivatives.numDerivativeLayer(p);
            const int layer_d_start = input_derivatives.layerStartIndex(p);

            double layer_fraction =
                (layer.altitude(sasktran_disco::Location::CEILING) - altitude) /
                (layer.altitude(sasktran_disco::Location::CEILING) -
                 layer.altitude(sasktran_disco::Location::FLOOR));
            double x = layer_fraction * layer.dual_thickness().value;

            // For the layer, we need the homogeneous solutions and the BVP
            // coefficients
            const auto& solution = layer.solution(m);

            const auto& h_minus = solution.value.dual_homog_minus();
            const auto& h_plus = solution.value.dual_homog_plus();

            const auto& g_plus_bottom = solution.value.dual_Gplus_bottom();
            const auto& g_plus_top = solution.value.dual_Gplus_top();
            const auto& g_minus_bottom = solution.value.dual_Gminus_bottom();
            const auto& g_minus_top = solution.value.dual_Gminus_top();

            ConstMatrixView homog_plus_matrix(
                h_plus.value.data(), nstr / 2 * NSTOKES, nstr / 2 * NSTOKES);
            ConstMatrixView homog_minus_matrix(
                h_minus.value.data(), nstr / 2 * NSTOKES, nstr / 2 * NSTOKES);

            ConstTensorView d_homog_plus(
                h_plus.deriv.data(), num_layer_derivatives, nstr / 2 * NSTOKES, nstr / 2 * NSTOKES
            );

            ConstTensorView d_homog_minus(
                h_minus.deriv.data(), num_layer_derivatives, nstr / 2 * NSTOKES, nstr / 2 * NSTOKES
            );

            const auto& eigval = solution.value.dual_eigval();

            const auto& dual_L = solution.boundary.L_coeffs;
            const auto& dual_M = solution.boundary.M_coeffs;


            // l = QUADRATURE ANGLE INDEX
            for (int l = 0; l < nstr / 2; ++l) {
                temp_deriv_downw.setZero();
                temp_deriv_upw.setZero();

                auto [upwelling_index, downwelling_index] = do_to_grid_zenith_index(l);

                int source_index_upw = m * NSTOKES +
                                    upwelling_index * num_azi * NSTOKES +
                                    szaidx * num_azi * NSTOKES * nstr +
                                    lidx * num_azi * NSTOKES * nstr *
                                        m_sza_grid.grid().size();

                int source_index_downw = m * NSTOKES +
                                    downwelling_index * num_azi * NSTOKES +
                                    szaidx * num_azi * NSTOKES * nstr +
                                    lidx * num_azi * NSTOKES * nstr *
                                        m_sza_grid.grid().size();

                // Loop over the homogeneous solutions

                // j = SOLUTION INDEX
                for (int j = 0; j < nstr / 2 * NSTOKES; ++j) {
                    double L = dual_L.value(j);
                    double M = dual_M.value(j);

                    double exp_f = std::exp(-layer.dual_thickness().value * eigval.value(j));
                    // Notes:
                    // d_exp_f / d_thickness = -exp_f * eigval.value(j)
                    // d_exp_f / d_eigval = -exp_f * layer.dual_thickness().value

                    // Offset to the downwelling direction depends on the quadrature angle...
                    // l=0 corresponds to the smallest angle, so the offset is the smallest

                    // s = STOKES INDEX
                    for (int s = 0; s < NSTOKES; ++s) {
                        int h_lidx = l * NSTOKES + s;

                        double neg =
                            sasktran_disco::stokes_negation_factor<NSTOKES>(s);

                        storage.radiances.value[source_index_upw + s] +=
                            0.5 * (1 + exp_f) * 
                            (L * homog_plus_matrix(h_lidx, j) + M * homog_minus_matrix(h_lidx, j));
                        
                        storage.radiances
                            .value[source_index_downw + s] +=
                            0.5 * (1 + exp_f) *
                            (L * homog_minus_matrix(h_lidx, j) + M * homog_plus_matrix(h_lidx, j)) *
                            neg;

                        // Derivatives
                        // homog_plus, homog_minus, thickness, and eigval all have local derivatives
                        for (int d = 0; d < num_layer_derivatives; ++d) {
                            temp_deriv_upw(d + layer_d_start, s) += 0.5 * (1 + exp_f) *
                                (L * d_homog_plus(d, h_lidx, j) +
                                 M * d_homog_minus(d, h_lidx, j));

                            temp_deriv_upw(d + layer_d_start, s) += -0.5 * exp_f *
                                (L * homog_plus_matrix(h_lidx, j) +
                                 M * homog_minus_matrix(h_lidx, j)) *
                                (eigval.deriv(d, j) * layer.dual_thickness().value + layer.dual_thickness().deriv(d) * eigval.value(j));

                            temp_deriv_downw(d + layer_d_start, s) += 0.5 * (1 + exp_f) *
                                (L * d_homog_minus(d, h_lidx, j) +
                                 M * d_homog_plus(d, h_lidx, j)) * neg;

                            temp_deriv_downw(d + layer_d_start, s) += -0.5 * exp_f *
                                (L * homog_minus_matrix(h_lidx, j) +
                                 M * homog_plus_matrix(h_lidx, j)) *
                                (eigval.deriv(d, j) * layer.dual_thickness().value + layer.dual_thickness().deriv(d) * eigval.value(j)) * neg;
                        }
                        // L and M have cross derivatives
                        for(int d = 0; d < num_total_derivatives; ++d) {
                            temp_deriv_upw(d, s) += 0.5 * (1 + exp_f) *
                                (dual_L.deriv(d, j) * homog_plus_matrix(h_lidx, j) +
                                 dual_M.deriv(d, j) * homog_minus_matrix(h_lidx, j));

                            temp_deriv_downw(d, s) += 0.5 * (1 + exp_f) *
                                (dual_L.deriv(d, j) * homog_minus_matrix(h_lidx, j) +
                                 dual_M.deriv(d, j) * homog_plus_matrix(h_lidx, j)) * neg;
                        }

                    }
                }
                // s = STOKES INDEX
                for (int s = 0; s < NSTOKES; ++s) {
                    storage.radiances.value[source_index_upw + s] += 
                        0.5 * (g_minus_bottom.value(l * NSTOKES + s) + 
                               g_minus_top.value(l * NSTOKES + s));
                    storage.radiances.value[source_index_downw + s] +=
                        0.5 * (g_plus_bottom.value(l * NSTOKES + s)
                                 + g_plus_top.value(l * NSTOKES + s));

                    // derivatives, straightforward
                    for(int d = 0; d < num_total_derivatives; ++d) {
                        temp_deriv_upw(d, s) += 0.5 * (g_minus_bottom.deriv(d, l * NSTOKES + s) +
                                                        g_minus_top.deriv(d, l * NSTOKES + s));
                        temp_deriv_downw(d, s) += 0.5 * (g_plus_bottom.deriv(d, l * NSTOKES + s) +
                                                           g_plus_top.deriv(d, l * NSTOKES + s));
                    }



                    if (num_total_derivatives > 0) {
                        storage.radiances.deriv(source_index_downw + s, Eigen::all)
                            .setZero();
                        storage.radiances.deriv(source_index_upw + s, Eigen::all)
                            .setZero();
                    }
                    for (int k = 0; k < num_total_derivatives; ++k) {
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

                            storage.radiances.deriv(
                                source_index_downw + s, group_fraction.first) +=
                                group_fraction.second * temp_deriv_downw(k, s) *
                                extinction;

                            storage.radiances.deriv(
                                source_index_upw + s, group_fraction.first) +=
                                group_fraction.second * temp_deriv_upw(k, s) *
                                extinction;
                        }
                    }
                }
                // Translate temporary DO level derivative to input atmosphere level derivative


            }

            if(p == optical_layer.numLayers() - 1) {
                // ground point, so add in the ground incoming radiance
                // l = QUADRATURE ANGLE INDEX
                for (int l = 0; l < nstr / 2; ++l) {
                    temp_deriv_downw.setZero();
                    temp_deriv_upw.setZero();
                    auto [upwelling_index, downwelling_index] = do_to_grid_zenith_index(l);

                    int source_index_downw = m * NSTOKES +
                                        downwelling_index * num_azi * NSTOKES +
                                        szaidx * num_azi * NSTOKES * nstr *
                                            m_sza_grid.grid().size() + m_ground_start;

                    // j = SOLUTION INDEX
                    for (int j = 0; j < nstr / 2 * NSTOKES; ++j) {
                        double L = dual_L.value(j);
                        double M = dual_M.value(j);

                        double exp_f = std::exp(-layer.dual_thickness().value * eigval.value(j));

                        // s = STOKES INDEX
                        for (int s = 0; s < NSTOKES; ++s) {
                            int h_lidx = l * NSTOKES + s;
                            double neg =
                            sasktran_disco::stokes_negation_factor<NSTOKES>(s);
                            
                            storage.radiances
                                .value[source_index_downw + s] +=
                                (L * homog_plus_matrix(h_lidx, j) * exp_f + M * homog_minus_matrix(h_lidx, j)) *
                                neg;

                            // Derivatives
                            // homog_plus, homog_minus, thickness, and eigval all have local derivatives
                            for (int d = 0; d < num_layer_derivatives; ++d) {
                                temp_deriv_downw(d + layer_d_start, s) += (d_homog_plus(d, h_lidx, j) * L * exp_f + d_homog_minus(d, h_lidx, j) * M) * neg;

                                temp_deriv_downw(d + layer_d_start, s) += -exp_f * L * homog_plus_matrix(h_lidx, j) * (
                                    layer.dual_thickness().deriv(d) * eigval.value(j) + eigval.deriv(d, j) * layer.dual_thickness().value) * neg;
                            }

                            // and L and M have total derivatives
                            for(int d = 0; d < num_total_derivatives; ++d) {
                                temp_deriv_downw(d, s) += dual_L.deriv(d, j) * homog_plus_matrix(h_lidx, j) * exp_f * neg +
                                    dual_M.deriv(d, j) * homog_minus_matrix(h_lidx, j) * neg;
                            }
                        }
                    }

                    // s = STOKES INDEX
                    for (int s = 0; s < NSTOKES; ++s) {
                        storage.radiances.value[source_index_downw + s] +=
                            g_minus_bottom.value(l * NSTOKES + s);

                        // Easy derivatives
                        for(int d = 0; d < num_total_derivatives; ++d) {
                            temp_deriv_downw(d, s) +=
                                g_minus_bottom.deriv(d, l * NSTOKES + s);
                        }

                        if (num_total_derivatives > 0) {
                            storage.radiances.deriv(source_index_downw + s, Eigen::all)
                                .setZero();
                        }
                        for (int k = 0; k < num_total_derivatives; ++k) {
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

                                storage.radiances.deriv(
                                    source_index_downw + s, group_fraction.first) +=
                                    group_fraction.second * temp_deriv_downw(k, s) *
                                    extinction;
                            }
                        }
                    }
                }
            }
        }
    }
