#include "sasktran2/geometry.h"
#include "sasktran2/raytracing.h"
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
        calculate_chapman_factors_raytracer(
            geometry.coordinates().earth_radius(), geometry);
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

template <int NSTOKES, int CNSTR>
void sasktran_disco::GeometryLayerArray<NSTOKES, CNSTR>::
    calculate_chapman_factors_raytracer(double earth_rad,
                                        const sasktran2::Geometry1D& geometry) {

    double earth_radius = geometry.coordinates().earth_radius();

    m_chapman_factors.setZero();
    auto ray_tracer = sasktran2::raytracing::SphericalShellRayTracer(geometry);

    double csz = this->M_CSZ;
    bool refraction = this->M_SOLAR_REFRACTION;

    sasktran2::viewinggeometry::ViewingRay ray;
    ray.look_away = geometry.coordinates().sun_unit();

    sasktran2::raytracing::TracedRay result;

    // For every layer, we have to construct a ray from the sun to the bottom of
    // the layer
    for (sasktran_disco::LayerIndex p = 0; p < this->M_NLYR; ++p) {

        ray.observer.position = geometry.coordinates().solar_coordinate_vector(
            csz, 0.0, m_floor_h(p));

        ray_tracer.trace_ray(ray, result, refraction);

        if (result.ground_is_hit) {
            // No transmission, but this means the chapman factors are actually
            // really large, so we set it negative and detect it later
            m_chapman_factors.row(p).setConstant(-1.0);
            continue;
        }

        // Now we have the traced ray, we can extract the chapman factors
        for (auto& layer : result.layers) {
            double layer_distance =
                layer.layer_distance * layer.curvature_factor;

            // Now we have to figure out which layer this corresponds to,
            double average_alttude =
                (layer.entrance.radius() + layer.exit.radius()) / 2.0 -
                earth_radius;

            // Start at the top and work down until we find the layer where the
            // floor radius is <
            sasktran_disco::LayerIndex q = 0;
            for (; q < this->M_NLYR; ++q) {
                if (average_alttude >= m_floor_h(q) &&
                    average_alttude <= m_ceiling_h(q)) {
                    break;
                }
            }
            if (q == this->M_NLYR) {
                // We will just assume we are at the bottom layer
                q = this->M_NLYR - 1;
            }
            double vertical_layer_distance = (m_ceiling_h(q) - m_floor_h(q));
            m_chapman_factors(p, q) += layer_distance / vertical_layer_distance;
        }
    }
}

SASKTRAN_DISCO_INSTANTIATE_TEMPLATE(sasktran_disco::GeometryLayerArray);
