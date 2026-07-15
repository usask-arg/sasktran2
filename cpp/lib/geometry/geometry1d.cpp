#include <sasktran2/geometry.h>
#include <sasktran2/validation/validation.h>

namespace sasktran2 {

    Geometry1D::Geometry1D(double cos_sza, double saa, double earth_radius,
                           Eigen::VectorXd&& grid_values,
                           sasktran2::grids::interpolation interp_method,
                           geometrytype geotype)
        : Geometry(Coordinates(cos_sza, saa, earth_radius, geotype, false)),
          m_alt_grid(std::move(grid_values),
                     sasktran2::grids::gridspacing::automatic,
                     sasktran2::grids::outofbounds::extend, interp_method) {

        m_refractive_index.resize(m_alt_grid.grid().size());
        m_refractive_index.setConstant(1.0);
        validate();
    }

    void Geometry1D::assign_interpolation_weights(
        const Location& loc,
        std::vector<std::pair<int, double>>& index_weights) const {
        double alt;

        if (coordinates().geometry_type() ==
            sasktran2::geometrytype::planeparallel) {
            alt = loc.position.z() - coordinates().earth_radius();
        } else {
            alt = loc.radius() - coordinates().earth_radius();
        }

        if (!std::isfinite(alt)) {
            spdlog::critical("Cannot interpolate a non-finite location");
            sasktran2::validation::throw_configuration_error();
        }

        // Shell intersections carry an authoritative grid index, but their
        // reconstructed Cartesian radius can differ slightly from the shell
        // radius.  Preserve the historical one-metre diagnostic tolerance so
        // harmless coordinate roundoff does not turn an exact endpoint into
        // an interpolated one, while still rejecting genuinely stale indices.
        constexpr double exact_altitude_tolerance_m = 1.0;
        const bool valid_exact_index =
            loc.lower_alt_index >= 0 &&
            loc.lower_alt_index < m_alt_grid.grid().size();
        // Pseudospherical geometry uses plane-parallel LOS rays but spherical
        // solar rays. An exact grid index is therefore the only unambiguous
        // altitude coordinate for locations produced by either tracer.
        const bool exact_index_is_authoritative =
            coordinates().geometry_type() ==
            sasktran2::geometrytype::pseudospherical;
        const bool exact_altitude_matches =
            valid_exact_index &&
            (exact_index_is_authoritative ||
             std::abs(alt - m_alt_grid.grid()[loc.lower_alt_index]) <=
                 exact_altitude_tolerance_m);

        if (loc.on_exact_altitude && exact_altitude_matches) {
            index_weights.resize(1);
            index_weights[0].first = loc.lower_alt_index;
            index_weights[0].second = 1;
            return;
        }

        std::array<double, 2> weight;
        std::array<int, 2> index;

        int num_contrib;

        m_alt_grid.calculate_interpolation_weights(alt, index, weight,
                                                   num_contrib);

        index_weights.resize(num_contrib);

        for (int i = 0; i < num_contrib; ++i) {
            index_weights[i].first = index[i];
            index_weights[i].second = weight[i];
        }
    }

    void Geometry1D::validate() const {
        // Check that the altitude grid is valid
        if (m_alt_grid.grid().size() < 2) {
            spdlog::critical(
                "Invalid altitude grid size: {}, must be at least 2",
                m_alt_grid.grid().size());

            sasktran2::validation::throw_configuration_error();
        }

        if (!m_alt_grid.grid().allFinite()) {
            spdlog::critical("Invalid altitude grid: values must be finite");
            sasktran2::validation::throw_configuration_error();
        }

        for (int i = 1; i < m_alt_grid.grid().size(); ++i) {
            if (m_alt_grid.grid()[i] <= m_alt_grid.grid()[i - 1]) {
                spdlog::critical(
                    "Invalid altitude grid: must be strictly increasing");
                sasktran2::validation::throw_configuration_error();
            }
        }

        const auto geometry_type = coordinates().geometry_type();
        if (geometry_type != sasktran2::geometrytype::planeparallel &&
            geometry_type != sasktran2::geometrytype::pseudospherical &&
            coordinates().earth_radius() + m_alt_grid.grid()[0] <= 0.0) {
            spdlog::critical(
                "Invalid spherical altitude grid: radii must be positive");
            sasktran2::validation::throw_configuration_error();
        }

        if (m_refractive_index.size() != m_alt_grid.grid().size()) {
            spdlog::critical("Invalid refractive index size: {}, expected {}",
                             m_refractive_index.size(),
                             m_alt_grid.grid().size());
            sasktran2::validation::throw_configuration_error();
        }

        if (!m_refractive_index.allFinite() ||
            (m_refractive_index.array() <= 0.0).any()) {
            spdlog::critical(
                "Invalid refractive index: values must be finite and positive");
            sasktran2::validation::throw_configuration_error();
        }
    }

} // namespace sasktran2
