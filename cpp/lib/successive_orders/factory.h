#pragma once

#include <sasktran2/source_interface.h>

#include <memory>

namespace sasktran2 {
    class Geometry1D;
    namespace raytracing {
        class RayTracerBase;
    }
} // namespace sasktran2

namespace sasktran2::successive_orders {
    /** Constructs the private C++ successive-orders implementation. */
    template <int NSTOKES>
    std::unique_ptr<SourceTermInterface<NSTOKES>>
    make_successive_orders_source(
        const sasktran2::raytracing::RayTracerBase& raytracer,
        const sasktran2::Geometry1D& geometry);

    extern template std::unique_ptr<SourceTermInterface<1>>
    make_successive_orders_source<1>(
        const sasktran2::raytracing::RayTracerBase&,
        const sasktran2::Geometry1D&);
    extern template std::unique_ptr<SourceTermInterface<3>>
    make_successive_orders_source<3>(
        const sasktran2::raytracing::RayTracerBase&,
        const sasktran2::Geometry1D&);
} // namespace sasktran2::successive_orders
