#pragma once

#ifdef SKTRAN_RUST_SUPPORT

#include <sasktran2/hr/diffuse_source.h>
#include <stdexcept>

namespace sasktran2::successive_orders {
    /**
     * A separate successive-orders source whose fixed-point operator and
     * convergence logic are implemented in Rust.
     *
     * C++ constructs the dimension-dependent diffuse grid and traces its rays.
     * Rust owns the immutable fixed-point operator, wavelength-dependent solve,
     * and integration of the converged source over the C++-traced LOS. The
     * Rust operator types do not assume an atmospheric dimensionality.
     */
    template <int NSTOKES>
    class RustSource final : public sasktran2::hr::DiffuseTable<NSTOKES> {
      public:
        RustSource(const sasktran2::raytracing::RayTracerBase& ray_tracer,
                   const sasktran2::Geometry1D& geometry)
            : sasktran2::hr::DiffuseTable<NSTOKES>(ray_tracer, geometry, true) {
        }

        RustSource(const sasktran2::raytracing::RustRayTracer2D& ray_tracer,
                   const sasktran2::Geometry2D& geometry)
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

        bool supports_linearization(
            sasktran2::LinearizationMode mode) const override {
            if (mode == sasktran2::LinearizationMode::Jacobian) {
                return false;
            }
            return this->native_products_available();
        }

        bool requires_integration() const override { return false; }

        bool has_interior_source() const override { return false; }
    };
} // namespace sasktran2::successive_orders

#endif
