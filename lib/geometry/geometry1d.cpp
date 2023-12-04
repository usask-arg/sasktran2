#include <sasktran2/geometry.h>

namespace sasktran2 {

    Geometry1D::Geometry1D(double cos_sza, double saa, double earth_radius,
                           Eigen::VectorXd&& grid_values,
                           sasktran2::grids::interpolation interp_method,
                           geometrytype geotype)
        : Geometry(Coordinates(cos_sza, saa, earth_radius, geotype, false)),
          m_alt_grid(std::move(grid_values),
                     sasktran2::grids::gridspacing::automatic,
                     sasktran2::grids::outofbounds::extend, interp_method) {}

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
} // namespace sasktran2
