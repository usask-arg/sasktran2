#include <sasktran2/geometry.h>
#include <sasktran2/validation/validation.h>

#include <climits>
#include <limits>
#include <stdexcept>

namespace {
    constexpr double altitude_tolerance_m = 1e-8;
    constexpr double angular_tolerance_rad = 1e-8;
    constexpr double pi_rad = static_cast<double>(EIGEN_PI);
    constexpr double pi_boundary_tolerance_rad =
        16.0 * std::numeric_limits<double>::epsilon() * pi_rad;

    struct AxisStencil {
        std::array<int, 2> indices{};
        std::array<double, 2> weights{};
        int size = 0;
    };

    int exact_axis_index(const Eigen::VectorXd& values, double value,
                         double tolerance) {
        const auto lower =
            std::lower_bound(values.begin(), values.end(), value);
        if (lower != values.end() && std::abs(*lower - value) <= tolerance) {
            return static_cast<int>(lower - values.begin());
        }
        if (lower != values.begin()) {
            const auto previous = lower - 1;
            if (std::abs(*previous - value) <= tolerance) {
                return static_cast<int>(previous - values.begin());
            }
        }
        return -1;
    }

    AxisStencil
    axis_interpolation_weights(const Eigen::VectorXd& values, double value,
                               sasktran2::grids::interpolation interpolation,
                               double tolerance) {
        AxisStencil result;

        const int exact = exact_axis_index(values, value, tolerance);
        if (exact >= 0) {
            result.indices[0] = exact;
            result.weights[0] = 1.0;
            result.size = 1;
            return result;
        }

        if (value <= values[0]) {
            result.indices[0] = 0;
            result.weights[0] = 1.0;
            result.size = 1;
            return result;
        }

        const int last = static_cast<int>(values.size()) - 1;
        if (value >= values[last]) {
            result.indices[0] = last;
            result.weights[0] = 1.0;
            result.size = 1;
            return result;
        }

        const auto upper_iterator =
            std::upper_bound(values.begin(), values.end(), value);
        const int upper = static_cast<int>(upper_iterator - values.begin());
        const int lower = upper - 1;

        if (interpolation == sasktran2::grids::interpolation::lower) {
            result.indices[0] = lower;
            result.weights[0] = 1.0;
            result.size = 1;
        } else {
            result.indices[0] = lower;
            result.indices[1] = upper;
            result.size = 2;

            if (interpolation == sasktran2::grids::interpolation::shell) {
                result.weights[0] = 0.5;
                result.weights[1] = 0.5;
            } else {
                result.weights[1] =
                    (value - values[lower]) / (values[upper] - values[lower]);
                result.weights[0] = 1.0 - result.weights[1];
            }
        }

        return result;
    }

    std::optional<int> finite_axis_cell(const Eigen::VectorXd& values,
                                        double value, double tolerance) {
        if (value < values[0] - tolerance ||
            value > values[values.size() - 1] + tolerance) {
            return std::nullopt;
        }

        const int exact = exact_axis_index(values, value, tolerance);
        if (exact >= 0) {
            return std::min(exact, static_cast<int>(values.size()) - 2);
        }

        const auto upper =
            std::upper_bound(values.begin(), values.end(), value);
        if (upper == values.begin() || upper == values.end()) {
            return std::nullopt;
        }
        return static_cast<int>(upper - values.begin()) - 1;
    }

    int extended_axis_cell(const Eigen::VectorXd& values, double value,
                           double tolerance) {
        const int exact = exact_axis_index(values, value, tolerance);
        if (exact >= 0) {
            return std::min(exact, static_cast<int>(values.size()) - 2);
        }

        if (value <= values[0]) {
            return 0;
        }
        if (value >= values[values.size() - 1]) {
            return static_cast<int>(values.size()) - 2;
        }

        return static_cast<int>(
                   std::upper_bound(values.begin(), values.end(), value) -
                   values.begin()) -
               1;
    }

    double unwrap_angle_near(double angle, double center) {
        constexpr double two_pi = 2.0 * pi_rad;
        while (angle - center > pi_rad) {
            angle -= two_pi;
        }
        while (angle - center < -pi_rad) {
            angle += two_pi;
        }
        return angle;
    }
} // namespace

namespace sasktran2 {
    Geometry2D::Geometry2D(double cos_sza, double saa, double earth_radius,
                           Eigen::VectorXd&& altitude_grid,
                           Eigen::VectorXd&& horizontal_angle_grid,
                           grids::interpolation altitude_interpolation)
        : Geometry2D(Coordinates(cos_sza, saa, earth_radius,
                                 geometrytype::spherical, false),
                     std::move(altitude_grid), std::move(horizontal_angle_grid),
                     altitude_interpolation) {}

    Geometry2D::Geometry2D(Coordinates&& coordinates,
                           Eigen::VectorXd&& altitude_grid,
                           Eigen::VectorXd&& horizontal_angle_grid,
                           grids::interpolation altitude_interpolation)
        : Geometry(std::move(coordinates)),
          m_alt_grid(std::move(altitude_grid), grids::gridspacing::automatic,
                     grids::outofbounds::extend, altitude_interpolation),
          m_horizontal_angles(std::move(horizontal_angle_grid)) {
        validate();
    }

    int Geometry2D::size() const {
        return num_altitudes() * num_horizontal_locations();
    }

    int Geometry2D::num_altitudes() const {
        return static_cast<int>(m_alt_grid.grid().size());
    }

    int Geometry2D::num_horizontal_locations() const {
        return static_cast<int>(m_horizontal_angles.size());
    }

    int Geometry2D::num_cells() const {
        return (num_altitudes() - 1) * (num_horizontal_locations() - 1);
    }

    std::pair<int, int> Geometry2D::location_shape() const {
        return {num_horizontal_locations(), num_altitudes()};
    }

    std::pair<int, int> Geometry2D::cell_shape() const {
        return {num_horizontal_locations() - 1, num_altitudes() - 1};
    }

    int Geometry2D::location_index(int altitude_index,
                                   int horizontal_index) const {
        if (altitude_index < 0 || altitude_index >= num_altitudes() ||
            horizontal_index < 0 ||
            horizontal_index >= num_horizontal_locations()) {
            throw std::out_of_range("Geometry2D location index out of range");
        }
        return horizontal_index * num_altitudes() + altitude_index;
    }

    std::pair<int, int>
    Geometry2D::location_indices(int flattened_index) const {
        if (flattened_index < 0 || flattened_index >= size()) {
            throw std::out_of_range(
                "Geometry2D flattened location index out of range");
        }
        return {flattened_index % num_altitudes(),
                flattened_index / num_altitudes()};
    }

    int Geometry2D::cell_index(int altitude_cell, int horizontal_cell) const {
        const int num_altitude_cells = num_altitudes() - 1;
        const int num_horizontal_cells = num_horizontal_locations() - 1;
        if (altitude_cell < 0 || altitude_cell >= num_altitude_cells ||
            horizontal_cell < 0 || horizontal_cell >= num_horizontal_cells) {
            throw std::out_of_range("Geometry2D cell index out of range");
        }
        return horizontal_cell * num_altitude_cells + altitude_cell;
    }

    std::pair<int, int> Geometry2D::cell_indices(int flattened_index) const {
        if (flattened_index < 0 || flattened_index >= num_cells()) {
            throw std::out_of_range(
                "Geometry2D flattened cell index out of range");
        }
        const int num_altitude_cells = num_altitudes() - 1;
        return {flattened_index % num_altitude_cells,
                flattened_index / num_altitude_cells};
    }

    Eigen::Vector3d Geometry2D::grid_location(int altitude_index,
                                              int horizontal_index) const {
        // Validate both indices before accessing either grid.
        location_index(altitude_index, horizontal_index);
        const double radius =
            coordinates().earth_radius() + m_alt_grid.grid()[altitude_index];
        return radius * coordinates().unit_vector_from_angles(
                            m_horizontal_angles[horizontal_index], 0.0);
    }

    double Geometry2D::altitude_at(const Location& location) const {
        if (!location.position.allFinite()) {
            spdlog::critical("Cannot evaluate a non-finite 2D location");
            validation::throw_configuration_error();
        }
        return location.radius() - coordinates().earth_radius();
    }

    double Geometry2D::horizontal_angle_at(const Location& location) const {
        if (!location.position.allFinite()) {
            spdlog::critical("Cannot evaluate a non-finite 2D location");
            validation::throw_configuration_error();
        }

        const double x = location.position.dot(coordinates().reference_x());
        const double z = location.position.dot(coordinates().reference_z());
        const double projected_norm = std::hypot(x, z);
        const double scale = std::max(1.0, location.position.norm());
        double raw_angle = 0.0;
        if (projected_norm >
            64.0 * std::numeric_limits<double>::epsilon() * scale) {
            raw_angle = std::atan2(x, z);
        }

        const double center =
            (m_horizontal_angles[0] +
             m_horizontal_angles[m_horizontal_angles.size() - 1]) /
            2.0;
        return unwrap_angle_near(raw_angle, center);
    }

    std::pair<double, double>
    Geometry2D::cell_interpolation_coordinates(const Location& location,
                                               int altitude_cell,
                                               int horizontal_cell) const {
        if (altitude_cell < 0 ||
            altitude_cell >= m_alt_grid.grid().size() - 1 ||
            horizontal_cell < 0 ||
            horizontal_cell >= m_horizontal_angles.size() - 1) {
            throw std::out_of_range("Invalid Geometry2D cell index");
        }
        double altitude_upper;
        switch (m_alt_grid.interpolation_method()) {
        case grids::interpolation::lower:
            altitude_upper = 0.0;
            break;
        case grids::interpolation::shell:
            altitude_upper = 0.5;
            break;
        case grids::interpolation::linear: {
            const double lower = m_alt_grid.grid()[altitude_cell];
            const double upper = m_alt_grid.grid()[altitude_cell + 1];
            altitude_upper = std::clamp(
                (altitude_at(location) - lower) / (upper - lower), 0.0, 1.0);
            break;
        }
        default:
            altitude_upper = std::numeric_limits<double>::quiet_NaN();
        }

        const double horizontal_lower = m_horizontal_angles[horizontal_cell];
        const double horizontal_upper =
            m_horizontal_angles[horizontal_cell + 1];
        const double horizontal_fraction =
            std::clamp((horizontal_angle_at(location) - horizontal_lower) /
                           (horizontal_upper - horizontal_lower),
                       0.0, 1.0);
        return {altitude_upper, horizontal_fraction};
    }

    std::optional<std::pair<int, int>>
    Geometry2D::cell_indices(const Location& location) const {
        const double altitude = altitude_at(location);
        const auto altitude_cell =
            finite_axis_cell(m_alt_grid.grid(), altitude, altitude_tolerance_m);
        if (!altitude_cell.has_value()) {
            return std::nullopt;
        }

        const int horizontal_cell = extended_axis_cell(
            m_horizontal_angles, horizontal_angle_at(location),
            angular_tolerance_rad);
        return std::make_pair(*altitude_cell, horizontal_cell);
    }

    void Geometry2D::assign_interpolation_weights(
        const Location& location,
        std::vector<std::pair<int, double>>& index_weights) const {
        const auto altitude_weights = axis_interpolation_weights(
            m_alt_grid.grid(), altitude_at(location),
            m_alt_grid.interpolation_method(), altitude_tolerance_m);
        const auto horizontal_weights = axis_interpolation_weights(
            m_horizontal_angles, horizontal_angle_at(location),
            grids::interpolation::linear, angular_tolerance_rad);

        index_weights.resize(altitude_weights.size * horizontal_weights.size);
        int output_index = 0;
        for (int horizontal = 0; horizontal < horizontal_weights.size;
             ++horizontal) {
            for (int altitude = 0; altitude < altitude_weights.size;
                 ++altitude) {
                index_weights[output_index] = {
                    location_index(altitude_weights.indices[altitude],
                                   horizontal_weights.indices[horizontal]),
                    altitude_weights.weights[altitude] *
                        horizontal_weights.weights[horizontal]};
                ++output_index;
            }
        }
    }

    void Geometry2D::validate() const {
        if (coordinates().geometry_type() != geometrytype::spherical) {
            spdlog::critical("Geometry2D only supports spherical coordinates");
            validation::throw_configuration_error();
        }

        const auto& altitudes = m_alt_grid.grid();
        if (altitudes.size() < 2) {
            spdlog::critical(
                "Invalid 2D altitude grid size: {}, must be at least 2",
                altitudes.size());
            validation::throw_configuration_error();
        }
        if (!altitudes.allFinite()) {
            spdlog::critical("Invalid 2D altitude grid: values must be finite");
            validation::throw_configuration_error();
        }
        for (Eigen::Index index = 1; index < altitudes.size(); ++index) {
            if (altitudes[index] <= altitudes[index - 1]) {
                spdlog::critical(
                    "Invalid 2D altitude grid: must be strictly increasing");
                validation::throw_configuration_error();
            }
        }
        if (coordinates().earth_radius() + altitudes[0] <= 0.0) {
            spdlog::critical(
                "Invalid 2D altitude grid: radii must be positive");
            validation::throw_configuration_error();
        }

        switch (m_alt_grid.interpolation_method()) {
        case grids::interpolation::linear:
        case grids::interpolation::shell:
        case grids::interpolation::lower:
            break;
        default:
            spdlog::critical("Invalid 2D altitude interpolation mode");
            validation::throw_configuration_error();
        }

        if (m_horizontal_angles.size() < 2) {
            spdlog::critical(
                "Invalid horizontal angle grid size: {}, must be at least 2",
                m_horizontal_angles.size());
            validation::throw_configuration_error();
        }
        if (!m_horizontal_angles.allFinite()) {
            spdlog::critical(
                "Invalid horizontal angle grid: values must be finite");
            validation::throw_configuration_error();
        }
        for (Eigen::Index index = 1; index < m_horizontal_angles.size();
             ++index) {
            if (m_horizontal_angles[index] <= m_horizontal_angles[index - 1]) {
                spdlog::critical(
                    "Invalid horizontal angle grid: must be strictly "
                    "increasing");
                validation::throw_configuration_error();
            }
        }
        if (m_horizontal_angles[m_horizontal_angles.size() - 1] -
                m_horizontal_angles[0] >=
            pi_rad - pi_boundary_tolerance_rad) {
            spdlog::critical(
                "Invalid horizontal angle grid: span must be less than pi");
            validation::throw_configuration_error();
        }

        if (altitudes.size() > INT_MAX ||
            m_horizontal_angles.size() > INT_MAX ||
            altitudes.size() >
                INT_MAX / static_cast<int>(m_horizontal_angles.size())) {
            spdlog::critical(
                "Geometry2D location count exceeds the public index range");
            validation::throw_configuration_error();
        }
    }
} // namespace sasktran2
