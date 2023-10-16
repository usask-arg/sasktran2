#include <sasktran2.h>
#include "sasktran2/viewinggeometry.h"

namespace sasktran2::viewinggeometry {
    void ViewingGeometryContainer::add_ray(
        const sasktran2::viewinggeometry::ViewingGeometryBase& ray) {

        if (dynamic_cast<
                const sasktran2::viewinggeometry::TangentAltitudeSolar*>(
                &ray)) {
            m_observer_rays.push_back(
                std::make_unique<
                    sasktran2::viewinggeometry::TangentAltitudeSolar>(
                    *dynamic_cast<const sasktran2::viewinggeometry::
                                      TangentAltitudeSolar*>(&ray)));
        } else if (dynamic_cast<
                       const sasktran2::viewinggeometry::GroundViewingSolar*>(
                       &ray)) {
            m_observer_rays.push_back(
                std::make_unique<
                    sasktran2::viewinggeometry::GroundViewingSolar>(
                    *dynamic_cast<
                        const sasktran2::viewinggeometry::GroundViewingSolar*>(
                        &ray)));
        } else {
            spdlog::error("Unsupported viewing ray type");
        }
    }
} // namespace sasktran2::viewinggeometry
