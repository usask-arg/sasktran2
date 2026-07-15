#pragma once

#include <sasktran2/raytracing.h>

namespace sasktran2::viewinggeometry {
    struct InternalViewingGeometry {
        std::vector<sasktran2::raytracing::TracedRay> traced_rays;

        std::vector<sasktran2::viewinggeometry::FluxObserver> flux_observers;

        std::size_t num_rays() const { return traced_rays.size(); }

        const sasktran2::viewinggeometry::ViewingRay&
        viewing_ray(std::size_t index) const {
            return traced_rays.at(index).observer_and_look;
        }
    };
} // namespace sasktran2::viewinggeometry
