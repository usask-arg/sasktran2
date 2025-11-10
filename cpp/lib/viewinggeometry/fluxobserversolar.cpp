#include <sasktran2.h>
#include "sasktran2/viewinggeometry.h"

namespace sasktran2::viewinggeometry {
    FluxObserverSolar::FluxObserverSolar(double cos_sza,
                                         double observeraltitude)
        : m_cos_sza(cos_sza), m_observer_altitude(observeraltitude) {}

    FluxObserver FluxObserverSolar::construct_flux_observer(
        const sasktran2::Coordinates& geometry) {
        FluxObserver flux_observer;

        flux_observer.observer.position = geometry.solar_coordinate_vector(
            m_cos_sza, 0.0, m_observer_altitude);

        return flux_observer;
    }

    std::string FluxObserverSolar::to_string() const {
        return fmt::format(
            "FluxObserverSolar(cos_sza={}, observer_altitude={})", m_cos_sza,
            m_observer_altitude);
    }
} // namespace sasktran2::viewinggeometry
