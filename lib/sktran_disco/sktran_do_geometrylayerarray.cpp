#include "sasktran2/geometry.h"
#include "sasktran2/refraction.h"
#include "sktran_disco/sktran_do.h"
#include "sktran_disco/sktran_do_geometrylayerarray.h"

template <int NSTOKES, int CNSTR>
sasktran_disco::GeometryLayerArray<NSTOKES, CNSTR>::GeometryLayerArray(
    const PersistentConfiguration<NSTOKES, CNSTR>& config,
    const sasktran2::Geometry1D& geometry)
    : m_config(config), GeometryLayerArrayROP<NSTOKES>(config) {
    // Initialize the chapman factor storage
    m_chapman_factors.resize(this->M_NLYR, this->M_NLYR);
    m_chapman_factors.setZero();

    // Calculate the ceiling heights and floor heights of each layer
    m_ceiling_h.resize(this->M_NLYR);
    m_floor_h.resize(this->M_NLYR);

    for (int p = 0; p < (int)this->M_NLYR; p++) {
        m_ceiling_h(p) = geometry.altitude_grid().grid().reverse()(p);

        if (p == this->M_NLYR - 1) {
            m_floor_h(p) = geometry.altitude_grid().grid()(0);
        } else {
            m_floor_h(p) = geometry.altitude_grid().grid().reverse()(p + 1);
        }
    }

    m_no_interp = geometry.altitude_grid().interpolation_method() ==
                  sasktran2::grids::interpolation::lower;

    // Construct the optical interpolating matrix
    m_optical_interpolator.resize(this->M_NLYR,
                                  geometry.altitude_grid().grid().size());
    m_optical_interpolator.setZero();

    if (this->M_NLYR != geometry.altitude_grid().grid().size() - 1) {
        throw std::runtime_error(
            "Number of layers does not match the number of grid points");
    }

    std::array<int, 2> index;
    std::array<double, 2> weight;
    int num_contributing;
    for (int p = 0; p < this->M_NLYR; p++) {
        double central_altitude = (m_ceiling_h(p) + m_floor_h(p)) / 2.0;

        geometry.altitude_grid().calculate_interpolation_weights(
            central_altitude, index, weight, num_contributing);

        for (int q = 0; q < num_contributing; q++) {
            m_optical_interpolator(p, index[q]) = weight[q];
        }
    }

    if (geometry.coordinates().geometry_type() ==
        sasktran2::geometrytype::planeparallel) {
        m_chapman_factors.triangularView<Eigen::Lower>().setConstant(
            1 / this->M_CSZ);
        m_chapman_factors.diagonal().setConstant(1 / this->M_CSZ);
    } else {
        // Now we have the layer locations, we can calculate the chapman factors
        calculate_chapman_factors(geometry.coordinates().earth_radius(),
                                  geometry);
    }
}

template <int NSTOKES, int CNSTR>
void sasktran_disco::GeometryLayerArray<NSTOKES, CNSTR>::
    calculate_chapman_factors(double earth_rad,
                              const sasktran2::Geometry1D& geometry) {
    m_chapman_factors.setZero();
    // Calculate the chapman factors for the layer
    double sinthetasq = 1 - this->M_CSZ * this->M_CSZ;
    std::vector<std::pair<int, double>> index_weights;

    for (sasktran_disco::LayerIndex p = 0; p < this->M_NLYR; ++p) {
        double rp = earth_rad + m_floor_h(p);

        // Straight line rt = rp sin(theta)
        double rt = rp * sqrt(sinthetasq);

        if (this->M_SOLAR_REFRACTION) {
            rt = sasktran2::raytracing::refraction::tangent_radius(
                geometry, rt, index_weights);
        }
        double nt =
            sasktran2::raytracing::refraction::refractive_index_at_altitude(
                geometry, rt - earth_rad, index_weights);

        if (rt > rp) {
            // Tangent viewing ray
            spdlog::warn(
                "Tangent viewing ray at layer {}, results may not be accurate",
                p);
            continue;
        }

        for (sasktran_disco::LayerIndex q = 0; q <= p; ++q) {
            double rfloor = earth_rad + m_floor_h(q);
            double rceil = earth_rad + m_ceiling_h(q);

            // If we are refracting do the full integral
            if (this->M_SOLAR_REFRACTION) {
                std::pair<double, double> refraction_result =
                    sasktran2::raytracing::refraction::integrate_path(
                        geometry, rt, nt, rfloor, rceil, index_weights);
                m_chapman_factors(p, q) =
                    refraction_result.first / (rceil - rfloor);
            } else {
                // Just use the straight line equation
                m_chapman_factors(p, q) =
                    (sqrt(rceil * rceil - rp * rp * sinthetasq) -
                     sqrt(rfloor * rfloor - rp * rp * sinthetasq)) /
                    (rceil - rfloor);
            }
        }
    }
}

SASKTRAN_DISCO_INSTANTIATE_TEMPLATE(sasktran_disco::GeometryLayerArray);
