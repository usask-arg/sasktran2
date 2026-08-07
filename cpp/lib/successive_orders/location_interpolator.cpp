#include <sasktran2/geometry.h>
#include "location_interpolator.h"

namespace sasktran2::successive_orders {
    AltitudeHorizontalSourceLocationInterpolator::
        AltitudeHorizontalSourceLocationInterpolator(
            sasktran2::grids::AltitudeGrid&& altitude_grid,
            Eigen::VectorXd&& horizontal_angle_grid,
            const sasktran2::Geometry2D& geometry)
        : sasktran2::grids::SourceLocationInterpolator(
              std::move(altitude_grid)),
          m_geometry(geometry),
          m_horizontal_grid(std::move(horizontal_angle_grid),
                            sasktran2::grids::gridspacing::automatic,
                            sasktran2::grids::outofbounds::extend,
                            sasktran2::grids::interpolation::linear),
          m_ground_altitude(geometry.altitude_grid().grid()[0]) {}

    int
    AltitudeHorizontalSourceLocationInterpolator::num_ground_points() const {
        return static_cast<int>(m_horizontal_grid.grid().size());
    }

    int
    AltitudeHorizontalSourceLocationInterpolator::num_interior_points() const {
        return static_cast<int>(m_altitude_grid.grid().size()) *
               num_ground_points();
    }

    int AltitudeHorizontalSourceLocationInterpolator::interior_linear_index(
        int alt_index, int horizontal_index) const {
        return alt_index + horizontal_index *
                               static_cast<int>(m_altitude_grid.grid().size());
    }

    int AltitudeHorizontalSourceLocationInterpolator::ground_linear_index(
        int horizontal_index) const {
        return num_interior_points() + horizontal_index;
    }

    Eigen::Vector3d AltitudeHorizontalSourceLocationInterpolator::grid_location(
        const sasktran2::Coordinates&, int location_index) const {
        const int num_altitudes =
            static_cast<int>(m_altitude_grid.grid().size());
        const int alt_index = location_index % num_altitudes;
        const int horizontal_index = location_index / num_altitudes;
        const double radius = m_geometry.coordinates().earth_radius() +
                              m_altitude_grid.grid()[alt_index];
        return radius * m_geometry.coordinates().unit_vector_from_angles(
                            m_horizontal_grid.grid()[horizontal_index], 0.0);
    }

    Eigen::Vector3d
    AltitudeHorizontalSourceLocationInterpolator::ground_location(
        const sasktran2::Coordinates&, int ground_index) const {
        const double radius =
            m_geometry.coordinates().earth_radius() + m_ground_altitude;
        return radius * m_geometry.coordinates().unit_vector_from_angles(
                            m_horizontal_grid.grid()[ground_index], 0.0);
    }

    void AltitudeHorizontalSourceLocationInterpolator::
        interior_interpolation_weights(
            const sasktran2::Coordinates&, const sasktran2::Location& location,
            std::vector<std::pair<int, double>>& weights, int& num_interp) {
        std::array<int, 2> alt_index;
        std::array<int, 2> horizontal_index;
        std::array<double, 2> alt_weight;
        std::array<double, 2> horizontal_weight;
        int num_altitude;
        int num_horizontal;

        m_altitude_grid.calculate_interpolation_weights(
            m_geometry.altitude_at(location), alt_index, alt_weight,
            num_altitude);
        m_horizontal_grid.calculate_interpolation_weights(
            m_geometry.horizontal_angle_at(location), horizontal_index,
            horizontal_weight, num_horizontal);

        num_interp = num_altitude * num_horizontal;
        if (weights.size() < static_cast<std::size_t>(num_interp)) {
            weights.resize(num_interp);
        }
        for (int horizontal = 0; horizontal < num_horizontal; ++horizontal) {
            for (int altitude = 0; altitude < num_altitude; ++altitude) {
                const int output = altitude + horizontal * num_altitude;
                weights[output] = {
                    interior_linear_index(alt_index[altitude],
                                          horizontal_index[horizontal]),
                    alt_weight[altitude] * horizontal_weight[horizontal]};
            }
        }
    }

    void
    AltitudeHorizontalSourceLocationInterpolator::ground_interpolation_weights(
        const sasktran2::Coordinates&, const sasktran2::Location& location,
        std::vector<std::pair<int, double>>& weights, int& num_interp) const {
        std::array<int, 2> horizontal_index;
        std::array<double, 2> horizontal_weight;
        m_horizontal_grid.calculate_interpolation_weights(
            m_geometry.horizontal_angle_at(location), horizontal_index,
            horizontal_weight, num_interp);

        if (weights.size() < static_cast<std::size_t>(num_interp)) {
            weights.resize(num_interp);
        }
        for (int horizontal = 0; horizontal < num_interp; ++horizontal) {
            weights[horizontal] = {
                ground_linear_index(horizontal_index[horizontal]),
                horizontal_weight[horizontal]};
        }
    }
} // namespace sasktran2::successive_orders
