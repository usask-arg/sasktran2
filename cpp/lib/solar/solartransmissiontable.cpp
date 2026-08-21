#include "sasktran2/geometry.h"
#include <sasktran2/solartransmission.h>

namespace sasktran2::solartransmission {

    void SolarTableInterpolation::clear() {
        std::vector<std::uint32_t>().swap(m_outer);
        std::vector<std::uint32_t>().swap(m_inner);
        std::vector<double>().swap(m_values);
        m_rows = 0;
        m_cols = 0;
        m_next_row = 0;
    }

    void SolarTableInterpolation::initialize(Eigen::Index rows,
                                             Eigen::Index cols,
                                             Eigen::Index capacity) {
        if (rows < 0 || cols < 0 || capacity < 0 ||
            cols > static_cast<Eigen::Index>(
                       std::numeric_limits<std::uint32_t>::max())) {
            throw std::invalid_argument(
                "Invalid compact solar interpolation dimensions");
        }
        m_rows = rows;
        m_cols = cols;
        m_next_row = 0;
        m_outer.assign(static_cast<std::size_t>(rows) + 1, 0);
        m_inner.clear();
        m_values.clear();
        m_inner.reserve(static_cast<std::size_t>(capacity));
        m_values.reserve(static_cast<std::size_t>(capacity));
    }

    void SolarTableInterpolation::append_row(
        const std::vector<std::pair<int, double>>& interpolation_weights) {
        if (m_next_row >= m_rows) {
            throw std::logic_error(
                "Too many rows appended to compact solar interpolation");
        }
        m_outer[static_cast<std::size_t>(m_next_row)] =
            static_cast<std::uint32_t>(m_values.size());
        for (const auto& [column, value] : interpolation_weights) {
            if (column < 0 || column >= m_cols || !std::isfinite(value)) {
                throw std::invalid_argument(
                    "Invalid compact solar interpolation entry");
            }
            if (value == 0.0) {
                continue;
            }
            if (m_values.size() == std::numeric_limits<std::uint32_t>::max()) {
                throw std::length_error(
                    "Compact solar interpolation has too many entries");
            }
            m_inner.push_back(static_cast<std::uint32_t>(column));
            m_values.push_back(value);
        }
        ++m_next_row;
    }

    void SolarTableInterpolation::finalize() {
        if (m_next_row != m_rows) {
            throw std::logic_error(
                "Compact solar interpolation has missing rows");
        }
        m_outer[static_cast<std::size_t>(m_rows)] =
            static_cast<std::uint32_t>(m_values.size());
        m_inner.shrink_to_fit();
        m_values.shrink_to_fit();
    }

    void SolarTableInterpolation::apply(
        Eigen::Ref<const Eigen::VectorXd> table_values,
        Eigen::Ref<Eigen::VectorXd> endpoint_values) const {
        if (table_values.size() != m_cols || endpoint_values.size() != m_rows) {
            throw std::invalid_argument(
                "Invalid compact solar interpolation product dimensions");
        }
        for (Eigen::Index row = 0; row < m_rows; ++row) {
            double result = 0.0;
            const auto begin = m_outer[static_cast<std::size_t>(row)];
            const auto end = m_outer[static_cast<std::size_t>(row) + 1];
            for (std::uint32_t entry = begin; entry < end; ++entry) {
                result += m_values[entry] * table_values(m_inner[entry]);
            }
            endpoint_values(row) = result;
        }
    }

    void SolarTableInterpolation::apply_transpose(
        Eigen::Ref<const Eigen::VectorXd> endpoint_values,
        Eigen::Ref<Eigen::VectorXd> table_values) const {
        if (endpoint_values.size() != m_rows || table_values.size() != m_cols) {
            throw std::invalid_argument(
                "Invalid compact solar interpolation transpose dimensions");
        }
        table_values.setZero();
        for (Eigen::Index row = 0; row < m_rows; ++row) {
            const double value = endpoint_values(row);
            const auto begin = m_outer[static_cast<std::size_t>(row)];
            const auto end = m_outer[static_cast<std::size_t>(row) + 1];
            for (std::uint32_t entry = begin; entry < end; ++entry) {
                table_values(m_inner[entry]) += m_values[entry] * value;
            }
        }
    }

    std::size_t SolarTableInterpolation::storage_bytes() const {
        return m_outer.capacity() * sizeof(std::uint32_t) +
               m_inner.capacity() * sizeof(std::uint32_t) +
               m_values.capacity() * sizeof(double);
    }

    void SolarTransmissionTable::initialize_geometry(
        const std::vector<sasktran2::raytracing::TracedRay>& integration_rays) {
        // find the min/max SZA from the LOS rays and generate the cos_sza_grid
        std::pair<double, double> min_max_cos_sza =
            sasktran2::raytracing::min_max_cos_sza_of_all_rays(
                integration_rays);

        Eigen::VectorXd cos_sza_grid_values;

        if (m_geometry.coordinates().geometry_type() ==
            sasktran2::geometrytype::spherical) {
            // TODO: configure this resolution
            cos_sza_grid_values.setLinSpaced(100, min_max_cos_sza.first,
                                             min_max_cos_sza.second);
        } else {
            // TODO: Can we handle pseudo-spherical here?
            cos_sza_grid_values.resize(1);
            cos_sza_grid_values(0) = min_max_cos_sza.first;
        }

        Eigen::VectorXd alt_values = m_geometry_1d->altitude_grid().grid();

        // create the location interpolator
        m_location_interpolator = std::make_unique<
            sasktran2::grids::AltitudeSZASourceLocationInterpolator>(
            sasktran2::grids::AltitudeGrid(
                std::move(alt_values), sasktran2::grids::gridspacing::constant,
                sasktran2::grids::outofbounds::extend,
                sasktran2::grids::interpolation::linear),
            sasktran2::grids::Grid(std::move(cos_sza_grid_values),
                                   sasktran2::grids::gridspacing::constant,
                                   sasktran2::grids::outofbounds::extend,
                                   sasktran2::grids::interpolation::linear));

        // Construct the matrix that calculates OD on the solar transmission
        // table locations i.e. solar_od_on_grid = matrix @ extinction
        m_geometry_matrix.resize(m_location_interpolator->num_interior_points(),
                                 m_geometry.size());
        m_geometry_matrix.setZero();

        m_ground_hit_flag.resize(
            m_location_interpolator->num_interior_points());

        sasktran2::viewinggeometry::ViewingRay ray_to_sun;

        ray_to_sun.look_away = m_geometry.coordinates().sun_unit();

        raytracing::TracedRay traced_ray;

        for (int i = 0; i < m_location_interpolator->num_interior_points();
             ++i) {
            ray_to_sun.observer.position =
                m_location_interpolator->grid_location(m_geometry.coordinates(),
                                                       i);

            // This method specifically does not allow for refraction
            m_raytracer->trace_ray(ray_to_sun, traced_ray, false);

            if (!traced_ray.ground_is_hit) {
                assign_dense_matrix_column(i, traced_ray, m_geometry_matrix);
                m_ground_hit_flag[i] = false;
            } else {
                m_ground_hit_flag[i] = true;
            }
        }
    }

    void SolarTransmissionTable::generate_interpolation_matrix(
        const std::vector<sasktran2::raytracing::TracedRay>& rays,
        Eigen::SparseMatrix<double, Eigen::RowMajor>& interpolator,
        std::vector<bool>& ground_hit_flag) const {
        // First calculate the number of points we need to create the matrix for
        // We calculate solar transmission at the boundaries of layers, so it is
        // nlayer+1 for each ray
        int numpoints = 0;
        for (const auto& ray : rays) {
            numpoints += (int)ray.layers.size() + 1;
        }

        // od matrix is such that matrix @ extinction = od
        interpolator.resize(numpoints,
                            m_location_interpolator->num_interior_points());

        // Have to handle rays that hit the ground separately since they have no
        // solar transmission
        ground_hit_flag.resize(numpoints, false);

        typedef Eigen::Triplet<double> T;
        std::vector<T> tripletList;

        std::vector<std::pair<int, double>> interpolator_storage;
        int num_interp;

        int row = 0;
        for (int i = 0; i < rays.size(); ++i) {
            const auto& ray = rays[i];
            for (int j = 0; j < ray.layers.size(); ++j) {
                const auto& layer = ray.layers[j];

                if (j == 0) {
                    // End layer at TOA, need to use layer exit
                    m_location_interpolator->interior_interpolation_weights(
                        m_geometry.coordinates(), layer.exit,
                        interpolator_storage, num_interp);

                    for (int k = 0; k < num_interp; ++k) {
                        tripletList.emplace_back(
                            T(row, interpolator_storage[k].first,
                              interpolator_storage[k].second));
                    }
                    ++row;
                }

                m_location_interpolator->interior_interpolation_weights(
                    m_geometry.coordinates(), layer.entrance,
                    interpolator_storage, num_interp);
                for (int k = 0; k < num_interp; ++k) {
                    tripletList.emplace_back(T(row,
                                               interpolator_storage[k].first,
                                               interpolator_storage[k].second));
                }
                ++row;
            }
        }
        interpolator.setFromTriplets(tripletList.begin(), tripletList.end());
    }

    void SolarTransmissionTable::generate_interpolation(
        const std::vector<sasktran2::raytracing::TracedRay>& rays,
        SolarTableInterpolation& interpolator,
        std::vector<bool>& ground_hit_flag,
        std::vector<Eigen::Vector3d>* solar_propagation_directions) const {
        Eigen::Index numpoints = 0;
        for (const auto& ray : rays) {
            numpoints += static_cast<Eigen::Index>(ray.layers.size()) + 1;
        }
        interpolator.initialize(numpoints,
                                m_location_interpolator->num_interior_points(),
                                numpoints * 4);
        ground_hit_flag.assign(static_cast<std::size_t>(numpoints), false);
        if (solar_propagation_directions != nullptr) {
            solar_propagation_directions->assign(
                static_cast<std::size_t>(numpoints),
                -m_geometry.coordinates().sun_unit());
        }

        std::vector<std::pair<int, double>> weights;
        int num_weights = 0;
        for (const auto& ray : rays) {
            if (ray.layers.empty()) {
                weights.clear();
                interpolator.append_row(weights);
                continue;
            }
            for (int layer_index = 0; layer_index < ray.layers.size();
                 ++layer_index) {
                const auto& layer = ray.layers[layer_index];
                if (layer_index == 0) {
                    m_location_interpolator->interior_interpolation_weights(
                        m_geometry.coordinates(), layer.exit, weights,
                        num_weights);
                    weights.resize(num_weights);
                    interpolator.append_row(weights);
                }
                m_location_interpolator->interior_interpolation_weights(
                    m_geometry.coordinates(), layer.entrance, weights,
                    num_weights);
                weights.resize(num_weights);
                interpolator.append_row(weights);
            }
        }
        interpolator.finalize();
    }

    void
    SolarTransmissionTable::apply(Eigen::Ref<const Eigen::VectorXd> extinction,
                                  Eigen::Ref<Eigen::VectorXd> table_od) const {
        if (extinction.size() != m_geometry_matrix.cols() ||
            table_od.size() != m_geometry_matrix.rows()) {
            throw std::invalid_argument(
                "Invalid 1D solar-table product dimensions");
        }
        table_od.noalias() = m_geometry_matrix * extinction;
    }

    void SolarTransmissionTable::accumulate_transpose(
        Eigen::Ref<const Eigen::VectorXd> table_cotangent,
        Eigen::Ref<Eigen::VectorXd> extinction_cotangent, double scale) const {
        if (table_cotangent.size() != m_geometry_matrix.rows() ||
            extinction_cotangent.size() != m_geometry_matrix.cols()) {
            throw std::invalid_argument(
                "Invalid 1D solar-table transpose dimensions");
        }
        extinction_cotangent.noalias() +=
            scale * m_geometry_matrix.transpose() * table_cotangent;
    }

    std::size_t SolarTransmissionTable::storage_bytes() const {
        return static_cast<std::size_t>(m_geometry_matrix.size()) *
                   sizeof(double) +
               m_ground_hit_flag.capacity() * sizeof(bool);
    }
} // namespace sasktran2::solartransmission
