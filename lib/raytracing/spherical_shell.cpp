#include <sasktran2/raytracing.h>
#include <sasktran2/refraction.h>

namespace sasktran2::raytracing {
    void SphericalShellRayTracer::trace_ray(
        const sasktran2::viewinggeometry::ViewingRay& ray, TracedRay& result,
        bool include_refraction) const {
        // Set the ray to 0
        result.reset();

        result.is_straight = !include_refraction ||
                             (fabs(ray.cos_viewing()) > NADIR_VIEWING_CUTOFF);

        // Calculate the tangent point details in a straight geometry
        // Ensure we don't get negative values from rounding errors
        double rt =
            ray.observer.radius() *
            sqrt(std::max(0.0, 1 - ray.cos_viewing() * ray.cos_viewing()));

        if (include_refraction) {
            // If including refraction, adjust the tangent radius
            rt = refraction::tangent_radius(m_geometry, rt,
                                            result.interpolation_index_weights);
        }
        result.tangent_radius = rt;

        double tangent_altitude = rt - m_earth_radius;

        // If the tangent altitude is greater than the TOA altitude then we have
        // an empty ray
        if (tangent_altitude >= m_alt_grid.grid()(Eigen::last)) {
            // Empty ray
            result.observer_and_look = ray;
            result.ground_is_hit = false;
        }

        if (ray.observer.radius() - m_earth_radius >=
            m_alt_grid.grid()(Eigen::last)) {
            // Outside atmosphere, probably
            if (ray.cos_viewing() > 0) {
                // We are looking up, so this is just an empty ray
                result.observer_and_look = ray;
                result.ground_is_hit = false;
                return;
            }

            // There are two cases, limb viewing and ground viewing
            if (tangent_altitude > m_alt_grid.grid()(0)) {
                // Limb viewing
                trace_ray_observer_outside_limb_viewing(ray, result);
            } else {
                // Ground hitting
                trace_ray_observer_outside_ground_viewing(ray, result);
            }
        } else {
            // We are inside the atmosphere
            if (ray.cos_viewing() > 0) {
                // Looking up, not limb viewing
                trace_ray_observer_inside_looking_up(ray, result);
            } else {
                // Looking downwards, could either be limb viewing or ground
                // hitting
                if (tangent_altitude > m_alt_grid.grid()(0)) {
                    // Limb viewing
                    trace_ray_observer_inside_looking_limb(ray, result);
                } else {
                    // Ground hitting
                    trace_ray_observer_inside_looking_ground(ray, result);
                }
            }
        }
        // Now that we have found all the layers, we go through and add in
        // all of the necessary geometry information for each layer
        finalize_ray_geometry(result);
    }

    /**
     * Goes through all of the layers, computes path lengths, and finds the
     * locations of the layers within the atmosphere.  Then we add in the
     * optical depth quadrature factors and compute solar angles.
     *
     * @param result
     */
    void
    SphericalShellRayTracer::finalize_ray_geometry(TracedRay& result) const {
        // Iterate through all of the rays, starting from the observer
        double total_deflection = 0.0;
        std::vector<std::pair<int, double>> index_weights;

        double costh = result.observer_and_look.cos_viewing();

        double rt_refracted = result.tangent_radius;

        double nt = refraction::refractive_index_at_altitude(
            m_geometry, rt_refracted - m_geometry.coordinates().earth_radius(),
            index_weights);

        Eigen::Vector3d x_basis, y_basis;

        for (int i = 0; i < result.layers.size(); ++i) {
            auto& layer = result.layers[result.layers.size() - i - 1];

            if (i == 0) {
                // First layer, set the position of the entrance layer

                // If the observer is inside the atmosphere
                if (result.observer_and_look.observer.radius() -
                        m_earth_radius <
                    m_alt_grid.grid()(Eigen::last)) {
                    // The observer is inside the atmosphere
                    // And so our ray tracing starts with the observer
                    layer.entrance.position =
                        result.observer_and_look.observer.position;
                } else {
                    // The observer is outside the atmosphere
                    // and so our raytracing starts at TOA
                    layer.entrance.position =
                        result.observer_and_look.observer.position +
                        result.observer_and_look.look_away *
                            distance_to_altitude(result.observer_and_look,
                                                 m_alt_grid.grid()(Eigen::last),
                                                 ViewingDirection::down,
                                                 TangentSide::nearside);
                }
                // Get the basis vectors
                x_basis = layer.entrance.position.normalized();
                y_basis = (x_basis.cross(result.observer_and_look.look_away)
                               .normalized())
                              .cross(x_basis)
                              .normalized();
            } else {
                // Have to assign entrance to be the exit of the previous layer
                layer.entrance = result.layers[result.layers.size() - i].exit;
            }

            // Now we have to assign the curvature factors and deflection angles
            if (!result.is_straight) {
                std::pair<double, double> refraction_result =
                    refraction::integrate_path(
                        m_geometry, rt_refracted, nt, layer.r_entrance,
                        layer.r_exit,
                        result
                            .interpolation_index_weights); // (layer distance,
                                                           // deflection angle)

                // Set the layer deflection angle
                total_deflection += refraction_result.second;
                // Can use the deflection to calculate the exit point
                layer.exit.position =
                    layer.r_exit * (x_basis * cos(total_deflection) +
                                    y_basis * sin(total_deflection));

                // And then get the average look vector in the layer
                layer.average_look_away =
                    (layer.exit.position - layer.entrance.position)
                        .normalized();

                // And the straight line layer distance
                layer.layer_distance =
                    (layer.exit.position - layer.entrance.position).norm();

                // Which gives us the curvature
                layer.curvature_factor =
                    refraction_result.first / layer.layer_distance;

            } else {
                // Layer is straight
                layer.curvature_factor = 1;

                // Straight line distance
                layer.layer_distance =
                    std::abs(sqrt(fmax(layer.r_entrance * layer.r_entrance -
                                           rt_refracted * rt_refracted,
                                       0)) -
                             sqrt(fmax(layer.r_exit * layer.r_exit -
                                           rt_refracted * rt_refracted,
                                       0.0)));

                // Same look vector as observer
                layer.average_look_away = result.observer_and_look.look_away;

                // Exit can be calculated from the look vector and straight line
                // distance
                layer.exit.position =
                    layer.entrance.position +
                    layer.average_look_away * layer.layer_distance;
            }

            add_od_quadrature(layer);
            add_interpolation_weights(layer, m_geometry);
            add_solar_parameters(m_geometry.coordinates().sun_unit(), layer);
        }
    }

    void SphericalShellRayTracer::trace_ray_observer_outside_ground_viewing(
        const sasktran2::viewinggeometry::ViewingRay& ray,
        TracedRay& tracedray) const {
        // We have len(altitudes) - 1 spherical layers, starting from the ground
        auto& result = tracedray;
        result.observer_and_look = ray;
        result.ground_is_hit = true;

        result.layers.resize(m_alt_grid.grid().size() - 1);

        for (int i = 0; i < m_alt_grid.grid().size() - 1; ++i) {
            complete_layer(result.layers[i], ray, i, ViewingDirection::down,
                           TangentSide::nearside);
        }
    }

    void SphericalShellRayTracer::trace_ray_observer_outside_limb_viewing(
        const sasktran2::viewinggeometry::ViewingRay& ray,
        TracedRay& tracedray) const {
        double rt = tracedray.tangent_radius;

        double tangent_altitude = rt - m_earth_radius;

        // Find the index to the first altitude ABOVE the tangent altitude
        auto it = std::upper_bound(m_alt_grid.grid().begin(),
                                   m_alt_grid.grid().cend(), tangent_altitude);
        size_t above_tangent_idx = std::distance(m_alt_grid.grid().begin(), it);

        auto& result = tracedray;
        result.observer_and_look = ray;
        result.ground_is_hit = false;

        size_t numlayer = 2 * (m_alt_grid.grid().size() - above_tangent_idx);
        result.layers.resize(numlayer);

        if (numlayer == 0) {
            // Empty ray
            return;
        }

        size_t layer_c = 0;
        // We have complete layers from the TOA down to the tangent layer
        for (size_t i = m_alt_grid.grid().size() - 1; i != above_tangent_idx;
             --i) {
            complete_layer(result.layers[layer_c], ray, i, ViewingDirection::up,
                           TangentSide::farside);
            ++layer_c;
        }

        tangent_layer(result.layers[layer_c], ray, above_tangent_idx,
                      tangent_altitude, ViewingDirection::up,
                      TangentSide::farside);
        ++layer_c;
        tangent_layer(result.layers[layer_c], ray, above_tangent_idx,
                      tangent_altitude, ViewingDirection::down,
                      TangentSide::nearside);
        ++layer_c;

        // And the complete layers from the tangent layer up to the instrument
        for (int i = (int)above_tangent_idx;
             i < (int)m_alt_grid.grid().size() - 1; ++i) {
            complete_layer(result.layers[layer_c], ray, i,
                           ViewingDirection::down, TangentSide::nearside);
            ++layer_c;
        }
        assert(layer_c == result.layers.size());
    }

    void SphericalShellRayTracer::complete_layer(
        SphericalLayer& layer,
        const sasktran2::viewinggeometry::ViewingRay& ray, size_t exit_index,
        ViewingDirection direction, TangentSide side) const {
        layer.type = LayerType::complete;

        double entrance_altitude = m_alt_grid.grid()(exit_index + direction);
        double exit_altitude = m_alt_grid.grid()(exit_index);

        layer.r_entrance = entrance_altitude + m_earth_radius;
        layer.r_exit = exit_altitude + m_earth_radius;

        layer.entrance.on_exact_altitude = true;
        layer.entrance.lower_alt_index = int(exit_index + direction);

        layer.exit.on_exact_altitude = true;
        layer.exit.lower_alt_index = int(exit_index);
    }

    void SphericalShellRayTracer::partial_layer(
        SphericalLayer& layer,
        const sasktran2::viewinggeometry::ViewingRay& ray, size_t start_index,
        ViewingDirection direction, TangentSide side) const {
        layer.type = LayerType::partial;

        double entrance_altitude = ray.observer.radius() - m_earth_radius;
        double exit_altitude = m_alt_grid.grid()(start_index);

        layer.r_entrance = entrance_altitude + m_earth_radius;
        layer.r_exit = exit_altitude + m_earth_radius;

        layer.exit.on_exact_altitude = true;
        layer.exit.lower_alt_index = int(start_index);

        layer.entrance.on_exact_altitude = false;
        layer.entrance.lower_alt_index =
            direction < 0 ? int(start_index + direction) : int(start_index);
    }

    void SphericalShellRayTracer::tangent_layer(
        SphericalLayer& layer,
        const sasktran2::viewinggeometry::ViewingRay& ray, size_t upper_index,
        double tangent_altitude, ViewingDirection direction,
        TangentSide side) const {
        double entrance_altitude, exit_altitude;

        layer.type = LayerType::tangent;

        if (direction == ViewingDirection::up) {
            entrance_altitude = tangent_altitude;
            exit_altitude = m_alt_grid.grid()(upper_index);

            layer.exit.on_exact_altitude = true;
            layer.exit.lower_alt_index = int(upper_index);

            layer.entrance.on_exact_altitude = false;
            layer.entrance.lower_alt_index = int(upper_index - 1);
        } else {
            exit_altitude = tangent_altitude;
            entrance_altitude = m_alt_grid.grid()(upper_index);

            layer.entrance.on_exact_altitude = true;
            layer.entrance.lower_alt_index = int(upper_index);

            layer.exit.on_exact_altitude = false;
            layer.exit.lower_alt_index = int(upper_index - 1);
        }

        layer.r_entrance = entrance_altitude + m_earth_radius;
        layer.r_exit = exit_altitude + m_earth_radius;
    }

    void SphericalShellRayTracer::partial_tangent_layer(
        SphericalLayer& layer,
        const sasktran2::viewinggeometry::ViewingRay& ray, size_t start_index,
        double tangent_altitude, ViewingDirection direction,
        TangentSide side) const {
        double entrance_altitude, exit_altitude;

        layer.type = LayerType::tangent;

        if (direction == ViewingDirection::up) {
            spdlog::error("Trying to construct a partial tangent layer looking "
                          "up, this shouldn't be a thing");
            throw std::runtime_error("critical raytracing error");
        } else {
            exit_altitude = tangent_altitude;
            entrance_altitude = ray.observer.radius() - m_earth_radius;

            layer.entrance.on_exact_altitude = false;
            layer.entrance.lower_alt_index = int(start_index - 1);

            layer.exit.on_exact_altitude = false;
            layer.exit.lower_alt_index = int(start_index - 1);

            layer.r_exit = exit_altitude + m_earth_radius;
            layer.r_entrance = entrance_altitude + m_earth_radius;
        }
    }

    void SphericalShellRayTracer::trace_ray_observer_inside_looking_up(
        const sasktran2::viewinggeometry::ViewingRay& ray,
        TracedRay& tracedray) const {
        // Find the index to the first altitude ABOVE the observer
        auto it =
            std::upper_bound(m_alt_grid.grid().begin(), m_alt_grid.grid().end(),
                             ray.observer.radius() - m_earth_radius);
        size_t start_index = std::distance(m_alt_grid.grid().begin(), it);

        auto& result = tracedray;
        result.observer_and_look = ray;
        result.ground_is_hit = false;

        result.layers.resize(m_alt_grid.grid().size() - start_index);

        int layer_c = 0;
        for (size_t i = m_alt_grid.grid().size() - 1; i != start_index; --i) {
            auto& layer = result.layers[layer_c];
            complete_layer(layer, ray, i, ViewingDirection::up,
                           TangentSide::nearside);
            ++layer_c;
        }
        assert(layer_c == result.layers.size() - 1);
        partial_layer(result.layers[layer_c], ray, start_index,
                      ViewingDirection::up, TangentSide::nearside);

        for (int i = 0; i < result.layers.size() - 1; ++i) {
            assert(result.layers[i + 1].exit.radius() ==
                   result.layers[i].entrance.radius());
        }

        assert(abs(result.layers[0].exit.radius() -
                   m_alt_grid.grid()(Eigen::last) - m_earth_radius) < 1e-8);
        assert(abs(result.layers[result.layers.size() - 1].entrance.radius() -
                   ray.observer.radius()) < 1e-8);
    }

    void SphericalShellRayTracer::trace_ray_observer_inside_looking_ground(
        const sasktran2::viewinggeometry::ViewingRay& ray,
        TracedRay& tracedray) const {
        auto& result = tracedray;

        result.observer_and_look = ray;
        result.ground_is_hit = true;

        // Find the index to the first altitude ABOVE the observer
        auto it =
            std::upper_bound(m_alt_grid.grid().begin(), m_alt_grid.grid().end(),
                             ray.observer.radius() - m_earth_radius);
        size_t start_index = std::distance(m_alt_grid.grid().begin(), it);

        if (start_index == 0) {
            // Have 0 layers, most likely rounding error
            return;
        }

        result.layers.resize(start_index);

        // Complete layers from the ground
        for (size_t i = 0; i < start_index - 1; ++i) {
            auto& layer = result.layers[i];
            complete_layer(layer, ray, i, ViewingDirection::down,
                           TangentSide::nearside);
        }
        partial_layer(result.layers[start_index - 1], ray, start_index - 1,
                      ViewingDirection::down, TangentSide::nearside);
    }
    void SphericalShellRayTracer::trace_ray_observer_inside_looking_limb(
        const sasktran2::viewinggeometry::ViewingRay& ray,
        TracedRay& tracedray) const {
        auto& result = tracedray;

        result.observer_and_look = ray;
        result.ground_is_hit = false;

        // Find the index to the first altitude ABOVE the observer
        auto it =
            std::upper_bound(m_alt_grid.grid().begin(), m_alt_grid.grid().end(),
                             ray.observer.radius() - m_earth_radius);
        size_t observer_idx = std::distance(m_alt_grid.grid().begin(), it);

        double rt = ray.observer.radius() *
                    sqrt(1 - ray.cos_viewing() * ray.cos_viewing());
        double tangent_altitude = rt - m_earth_radius;

        // Find the index to the first altitude ABOVE the tangent altitude
        auto it_tangent =
            std::upper_bound(m_alt_grid.grid().begin(),
                             m_alt_grid.grid().cend(), tangent_altitude);
        size_t above_tangent_idx =
            std::distance(m_alt_grid.grid().begin(), it_tangent);

        // We have (len(alt) - above_tangent) full layers
        // +2 tangent layers
        // + (observer - above_tangent) layers on the opposite side

        // Above the tangent, far side
        int num_layers =
            int(m_alt_grid.grid().size() - 1) - int(above_tangent_idx);

        // Tangent layers
        num_layers += 2;

        // Layers on the far side
        if (observer_idx > above_tangent_idx) {
            num_layers += int(observer_idx - above_tangent_idx);
        }

        result.layers.resize(num_layers);

        size_t layer_c = 0;
        // We have complete layers from the TOA down to the tangent layer
        for (size_t i = m_alt_grid.grid().size() - 1; i != above_tangent_idx;
             --i) {
            complete_layer(result.layers[layer_c], ray, i, ViewingDirection::up,
                           TangentSide::farside);
            ++layer_c;
        }

        // Then always one tangent layer
        tangent_layer(result.layers[layer_c], ray, above_tangent_idx,
                      tangent_altitude, ViewingDirection::up,
                      TangentSide::farside);
        ++layer_c;

        // We have a special case if above_tangent_idx == observer_idx, then we
        // have a partial tangent layer
        if (above_tangent_idx == observer_idx) {
            partial_tangent_layer(result.layers[layer_c], ray,
                                  above_tangent_idx, tangent_altitude,
                                  ViewingDirection::down,
                                  TangentSide::nearside);
            ++layer_c;
        } else {
            tangent_layer(result.layers[layer_c], ray, above_tangent_idx,
                          tangent_altitude, ViewingDirection::down,
                          TangentSide::nearside);
            ++layer_c;
        }

        // Complete layers from the tangent to the observer
        for (size_t i = above_tangent_idx; i < observer_idx - 1; ++i) {
            auto& layer = result.layers[layer_c];
            complete_layer(layer, ray, i, ViewingDirection::down,
                           TangentSide::nearside);
            ++layer_c;
        }

        // Final partial layer, only if we didn't have a tangent partial layer
        if (observer_idx > above_tangent_idx) {
            partial_layer(result.layers[layer_c], ray, observer_idx - 1,
                          ViewingDirection::down, TangentSide::nearside);
            ++layer_c;
        }

        assert(layer_c == result.layers.size());
    }

} // namespace sasktran2::raytracing
