#pragma once

#ifdef SKTRAN_RUST_SUPPORT

#include <sasktran2/hr/diffuse_source.h>
#include <stdexcept>

namespace sasktran2::successive_orders {
    /**
     * A separate successive-orders source whose fixed-point operator and
     * convergence logic are implemented in Rust.
     *
     * The first implementation deliberately reuses the established C++ 1-D
     * diffuse-grid construction, traced-ray integration, and line-of-sight
     * projection. This keeps the language boundary at immutable operator
     * structure plus wavelength-dependent coefficients. The Rust operator
     * types themselves do not assume a 1-D grid.
     */
    template <int NSTOKES>
    class RustSource final : public sasktran2::hr::DiffuseTable<NSTOKES> {
      public:
        RustSource(const sasktran2::raytracing::RayTracerBase& ray_tracer,
                   const sasktran2::Geometry1D& geometry)
            : sasktran2::hr::DiffuseTable<NSTOKES>(ray_tracer, geometry, true) {
        }

        void initialize_geometry(
            const sasktran2::viewinggeometry::InternalViewingGeometry&
                internal_viewing) override {
            if (!internal_viewing.flux_observers.empty()) {
                throw std::invalid_argument(
                    "The Rust successive-orders source currently supports "
                    "radiances only, not flux observers");
            }
            sasktran2::hr::DiffuseTable<NSTOKES>::initialize_geometry(
                internal_viewing);
        }

        void initialize_atmosphere(
            const sasktran2::atmosphere::Atmosphere<NSTOKES>& atmosphere)
            override {
            if (atmosphere.num_deriv() != 0) {
                throw std::invalid_argument(
                    "The Rust successive-orders source does not yet support "
                    "weighting functions");
            }
            sasktran2::hr::DiffuseTable<NSTOKES>::initialize_atmosphere(
                atmosphere);
        }
    };
} // namespace sasktran2::successive_orders

#endif
