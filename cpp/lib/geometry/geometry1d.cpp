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
        if (loc.on_exact_altitude && loc.lower_alt_index >= 0) {
            index_weights.resize(1);
            index_weights[0].first = loc.lower_alt_index;
            index_weights[0].second = 1;
        }

        std::array<double, 2> weight;
        std::array<int, 2> index;

        int num_contrib;

        m_alt_grid.calculate_interpolation_weights(
            loc.radius() - coordinates().earth_radius(), index, weight,
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

        for (int i = 1; i < m_alt_grid.grid().size(); ++i) {
            if (m_alt_grid.grid()[i] <= m_alt_grid.grid()[i - 1]) {
                spdlog::critical(
                    "Invalid altitude grid: must be strictly increasing");
                sasktran2::validation::throw_configuration_error();
            }
        }
    }

} // namespace sasktran2
