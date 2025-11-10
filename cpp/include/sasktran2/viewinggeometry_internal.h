#pragma once

#include <sasktran2/raytracing.h>

namespace sasktran2::viewinggeometry {
    struct InternalViewingGeometry {
        std::vector<sasktran2::raytracing::TracedRay> traced_rays;

        std::vector<sasktran2::viewinggeometry::FluxObserver> flux_observers;
    };
} // namespace sasktran2::viewinggeometry
