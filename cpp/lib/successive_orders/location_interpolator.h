#pragma once

#include <sasktran2/grids.h>

namespace sasktran2 {
    class Geometry2D;
}

namespace sasktran2::successive_orders {
    /** Source grid on the native horizontal coordinate of a Geometry2D.
     *
     * Interior points are ordered with altitude as the fastest-varying axis,
     * matching Geometry2D atmosphere storage. The source altitude grid may
     * differ from the atmosphere altitude grid.
     */
    class AltitudeHorizontalSourceLocationInterpolator final
        : public sasktran2::grids::SourceLocationInterpolator {
      private:
        const sasktran2::Geometry2D& m_geometry;
        const sasktran2::grids::Grid m_horizontal_grid;
        const double m_ground_altitude;

        int interior_linear_index(int alt_index, int horizontal_index) const;
        int ground_linear_index(int horizontal_index) const;

      public:
        AltitudeHorizontalSourceLocationInterpolator(
            sasktran2::grids::AltitudeGrid&& altitude_grid,
            Eigen::VectorXd&& horizontal_angle_grid,
            const sasktran2::Geometry2D& geometry);

        int num_interior_points() const override;
        int num_ground_points() const override;

        Eigen::Vector3d grid_location(const sasktran2::Coordinates& coords,
                                      int location_index) const override;
        Eigen::Vector3d ground_location(const sasktran2::Coordinates& coords,
                                        int ground_index) const override;

        void interior_interpolation_weights(
            const sasktran2::Coordinates& coords,
            const sasktran2::Location& location,
            std::vector<std::pair<int, double>>& weights,
            int& num_interp) override;
        void ground_interpolation_weights(
            const sasktran2::Coordinates& coords,
            const sasktran2::Location& location,
            std::vector<std::pair<int, double>>& weights,
            int& num_interp) const override;
    };
} // namespace sasktran2::successive_orders
