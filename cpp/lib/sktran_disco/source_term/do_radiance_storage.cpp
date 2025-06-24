#include "sasktran2/do_source.h"
#include "sasktran2/geometry.h"
#include "sktran_disco/sktran_do.h"
#include "sktran_disco/sktran_do_polarization_types.h"

namespace sasktran2 {
    template <int NSTOKES, int CNSTR>
    DORadianceStorage<NSTOKES, CNSTR>::DORadianceStorage(
        const sasktran_disco::GeometryLayerArray<NSTOKES, CNSTR>&
            layer_geometry,
        const sasktran_disco::PersistentConfiguration<NSTOKES, CNSTR>&
            do_config,
        const sasktran2::grids::Grid& sza_grid, const Config& config,
        const sasktran2::Geometry1D& geometry)
        : m_config(config), m_geometry(geometry), m_sza_grid(sza_grid) {
        // Create an altitude grid using the layer geometry
        Eigen::VectorXd altitude_mid =
            (layer_geometry.layer_ceiling() + layer_geometry.layer_floor()) /
            2.0;

        // Layers are defined from TOA downards, so we have to reverse it
        altitude_mid.reverseInPlace();

        m_altitude_grid = std::make_unique<sasktran2::grids::AltitudeGrid>(
            std::move(altitude_mid), sasktran2::grids::gridspacing::variable,
            sasktran2::grids::outofbounds::extend,
            sasktran2::grids::interpolation::linear);

        m_num_azi = do_config.nstr();

        m_storage.resize(config.num_wavelength_threads());

        // Construct the cos_angle grid, which is the same as the quadrature grid
        auto& quadrature_angles = *do_config.quadrature_cos_angle();

        Eigen::VectorXd cos_angles;
        cos_angles.resize(quadrature_angles.size());
        // Ordering of the quadrature angles is a little weird,
        for(int i = 0; i < m_num_azi/2; ++i) {
            cos_angles(i) = quadrature_angles[quadrature_angles.size() - 1 - i];
            cos_angles(i + m_num_azi/2) = quadrature_angles[i];
        }

        m_cos_angle_grid = std::make_unique<sasktran2::grids::Grid>(
            std::move(cos_angles), sasktran2::grids::gridspacing::variable,
            sasktran2::grids::outofbounds::extend,
            sasktran2::grids::interpolation::linear);

        // We have alt_grid * cos_angle_grid * num_azi * sza_grid source points
        int num_source_points = (int)m_altitude_grid->grid().size() *
                                (int)m_sza_grid.grid().size() *
                                (int)m_cos_angle_grid->grid().size() *
                                m_config.num_do_streams() * NSTOKES;

        // We include the full cos_angle in storage, even though half of them are 0 at the ground point,
        // makes the interpolation easier
        int num_ground_points = (int)m_cos_angle_grid->grid().size() *
                                m_config.num_do_streams() *
                                (int)m_sza_grid.grid().size() * NSTOKES;
        
        m_ground_start = num_source_points;

        for (auto& storage : m_storage) {
            storage.radiances.resize(num_source_points + num_ground_points, 0, true);
        }
    }

    template <int NSTOKES, int CNSTR>
    void DORadianceStorage<NSTOKES, CNSTR>::accumulate_sources(
        sasktran_disco::OpticalLayerArray<NSTOKES, CNSTR>& optical_layer,
        sasktran_disco::AEOrder m,
        sasktran2::DOSourceThreadStorage<NSTOKES, CNSTR>& thread_storage,
        int szaidx, int thread_idx) {
        using MatrixView = Eigen::Map<Eigen::MatrixXd>;
        using ConstMatrixView = Eigen::Map<const Eigen::MatrixXd>;

        auto& storage = m_storage[thread_idx];

        const int nstr = m_config.num_do_streams();
        const int num_azi = nstr;

        
        int downwelling_offset = nstr/2 * num_azi * NSTOKES;

        for (int lidx = 0; lidx < m_altitude_grid->grid().size(); ++lidx) {
            // Pull out the layer and where we are inside the layer
            auto& layer = optical_layer.layer(
                (int)m_altitude_grid->grid().size() - lidx - 1);
            sasktran_disco::LayerIndex p = layer.index();
            double altitude = m_altitude_grid->grid()(lidx);

            double layer_fraction =
                (layer.altitude(sasktran_disco::Location::CEILING) -
                 altitude) /
                (layer.altitude(sasktran_disco::Location::CEILING) -
                 layer.altitude(sasktran_disco::Location::FLOOR));
            double x = layer_fraction * layer.dual_thickness().value;

            // For the layer, we need the homogeneous solutions and the BVP coefficients
            const auto& solution = layer.solution(m);

            const auto& h_minus = solution.value.dual_homog_minus();
            const auto& h_plus = solution.value.dual_homog_plus();


            ConstMatrixView homog_plus_matrix(
                h_plus.value.data(),
                nstr / 2 * NSTOKES,
                nstr / 2 * NSTOKES);
            ConstMatrixView homog_minus_matrix(
                h_minus.value.data(),
                nstr / 2 * NSTOKES,
                nstr / 2 * NSTOKES);

            const auto& eigval = solution.value.dual_eigval();

            const auto& dual_L = solution.boundary.L_coeffs;
            const auto& dual_M = solution.boundary.M_coeffs;

            // Loop over the homogeneous solutions
            for(int j = 0; j < nstr / 2 * NSTOKES; ++j) {
                double L = dual_L.value(j);
                double M = dual_M.value(j);

                double exp_f1 = std::exp(-x * eigval.value(j));
                double exp_f2 = std::exp(- (1 - layer_fraction) * eigval.value(j));

                for(int l = 0; l < nstr / 2; ++l) {
                    for(int s = 0; s < NSTOKES; ++s) {
                        int lidx = l*NSTOKES + s;

                        double neg = sasktran_disco::stokes_negation_factor<NSTOKES>(s);

                        int source_index = s +
                            m * NSTOKES + 
                            l * num_azi * NSTOKES + 
                            szaidx * num_azi * NSTOKES * nstr + 
                            j * num_azi * NSTOKES * nstr * m_sza_grid.grid().size();

                            storage.radiances.value[source_index] += 
                                (L * homog_plus_matrix(lidx, j) * exp_f1 + M * homog_minus_matrix(lidx, j) * exp_f2) 
                            ;

                            storage.radiances.value[source_index + downwelling_offset] += 
                                (L * homog_minus_matrix(lidx, j) * exp_f1 + M * homog_plus_matrix(lidx, j) * exp_f2) * neg;

                    }
                    
                }
            }
        }

    }

    template <int NSTOKES, int CNSTR>
    void DORadianceStorage<NSTOKES, CNSTR>::create_location_radiance_interpolator(
        const std::vector<Eigen::Vector3d>& locations,
        const std::vector<Eigen::Vector3d>& directions,
        const std::vector<bool>& ground_hit_flag,
        Eigen::SparseMatrix<double, Eigen::RowMajor>& interpolator) {

            interpolator.resize(locations.size() * NSTOKES,
                                m_storage[0].radiances.value_size());

            typedef Eigen::Triplet<double> T;
            std::vector<T> tripletList;

            std::array<int, 2> alt_index, angle_index, sza_index;
            std::array<double, 2> alt_weight, angle_weight, sza_weight;
            int num_alt_contrib, num_angle_contrib, num_sza_contrib;
            double csz, saa;

            sasktran2::Location temp;

            int num_azi = m_config.num_do_streams();

            double earth_radius = m_geometry.coordinates().earth_radius();

            for (int i = 0; i < locations.size(); ++i) {
                const auto& location = locations[i];
                const auto& direction = directions[i];

                temp.position = location;

                sasktran2::raytracing::calculate_csz_saz(
                    m_geometry.coordinates().sun_unit(), temp, direction, csz, saa);

                double cos_angle, altitude;

                if (m_geometry.coordinates().geometry_type() ==
                    sasktran2::geometrytype::spherical) {
                    cos_angle = temp.cos_zenith_angle(-1 * direction);
                    altitude = temp.radius() - earth_radius;
                } else {
                    cos_angle = direction.z();
                    altitude = temp.position.z() - earth_radius;
                }

                m_sza_grid.calculate_interpolation_weights(
                    csz, sza_index, sza_weight, num_sza_contrib);
                if (!ground_hit_flag[i]) {
                    // Interior point, interpolate in angle and altitude
                    m_altitude_grid->calculate_interpolation_weights(
                        altitude, alt_index, alt_weight, num_alt_contrib);
                    m_cos_angle_grid->calculate_interpolation_weights(
                        cos_angle, angle_index, angle_weight, num_angle_contrib);

                    for (int szaidx = 0; szaidx < num_sza_contrib; ++szaidx) {
                        for (int altidx = 0; altidx < num_alt_contrib; ++altidx) {
                            for (int angleidx = 0; angleidx < num_angle_contrib;
                                ++angleidx) {
                                double weight = alt_weight[altidx] *
                                                angle_weight[angleidx] *
                                                sza_weight[szaidx];

                                for (int k = 0; k < num_azi; ++k) {
                                    double azi_factor = cos(k * (EIGEN_PI - saa));
                                    int index = linear_storage_index(
                                        angle_index[angleidx], alt_index[altidx],
                                        sza_index[szaidx], k);

                                    if constexpr (NSTOKES == 1) {
                                        tripletList.emplace_back(
                                            T(i, index, azi_factor * weight));
                                    } else if constexpr (NSTOKES == 3) {
                                        double sin_azi_factor =
                                            sin(k * (EIGEN_PI - saa));

                                        tripletList.emplace_back(
                                            T(i * NSTOKES, index * NSTOKES,
                                            azi_factor * weight));
                                        tripletList.emplace_back(
                                            T(i * NSTOKES + 1, index * NSTOKES + 1,
                                            azi_factor * weight));
                                        tripletList.emplace_back(
                                            T(i * NSTOKES + 2, index * NSTOKES + 2,
                                            sin_azi_factor * weight));
                                    }
                                }
                            }
                        }
                    }
                } else {
                    // Ground point, just have to interpolate in SZA and angle

                    double mu;
                    if (m_geometry.coordinates().geometry_type() ==
                        sasktran2::geometrytype::spherical) {
                        mu = location.normalized().dot(-direction);
                    } else {
                        mu = -direction.z();
                    }

                    m_cos_angle_grid->calculate_interpolation_weights(
                        mu, angle_index, angle_weight, num_angle_contrib);

                    // If we have a cos_angle < 0 move the weight into the next one
                    if (m_cos_angle_grid->grid()(angle_index[0]) < 0) {
                        angle_index[0] = angle_index[1];
                    }

                    for (int szaidx = 0; szaidx < num_sza_contrib; ++szaidx) {
                        for (int angleidx = 0; angleidx < num_angle_contrib;
                            ++angleidx) {
                            double weight =
                                sza_weight[szaidx] * angle_weight[angleidx];

                            for (int m = 0; m < num_azi; ++m) {
                                int index = ground_storage_index(
                                    angle_index[angleidx], sza_index[szaidx], m);

                                double azi_factor = cos(m * (EIGEN_PI - saa));

                                if constexpr (NSTOKES == 1) {
                                    tripletList.emplace_back(
                                        T(i, index, weight * azi_factor));
                                } else if constexpr (NSTOKES == 3) {
                                    double sin_azi_factor =
                                        sin(m * (EIGEN_PI - saa));
                                    tripletList.emplace_back(
                                        T(i * NSTOKES, index * NSTOKES,
                                        weight * azi_factor));
                                    tripletList.emplace_back(
                                        T(i * NSTOKES + 1, index * NSTOKES + 1,
                                        weight * azi_factor));
                                    tripletList.emplace_back(
                                        T(i * NSTOKES + 2, index * NSTOKES + 2,
                                        weight * sin_azi_factor));
                                }
                            }
                        }
                    }
                }
            }
            interpolator.setFromTriplets(tripletList.begin(), tripletList.end());
        }

    template <int NSTOKES, int CNSTR>
    int DORadianceStorage<NSTOKES, CNSTR>::linear_storage_index(
        int angleidx, int layeridx, int szaidx, int aziidx) const {

        return aziidx + m_num_azi * angleidx + 
                m_num_azi * (int)m_cos_angle_grid->grid().size() * szaidx + 
                m_num_azi * (int)m_cos_angle_grid->grid().size() * (int)m_sza_grid.grid().size() * layeridx;
    }

    template <int NSTOKES, int CNSTR>
    void DORadianceStorage<NSTOKES, CNSTR>::initialize_atmosphere(
        const sasktran2::atmosphere::Atmosphere<NSTOKES>& atmo) {
        // TODO: This seems like the wrong place to do this, but I'm not sure
        // where else it should go
        m_atmosphere = &atmo;

        int numderiv = atmo.num_deriv();

        for (auto& storage : m_storage) {
            storage.radiances.deriv.resize(
                storage.radiances.value_size(), numderiv);
        }
    }

    template <int NSTOKES, int CNSTR>
    int DORadianceStorage<NSTOKES, CNSTR>::ground_storage_index(
        int angleidx, int szaidx, int aziidx) const {

        return aziidx + m_num_azi * angleidx + 
                m_num_azi * (int)m_cos_angle_grid->grid().size() * szaidx +
                + m_ground_start;
    }

    SASKTRAN_DISCO_INSTANTIATE_TEMPLATE(DORadianceStorage);
}