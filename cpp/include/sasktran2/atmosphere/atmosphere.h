#pragma once

#include "../geometry.h"
#include <atomic>
#include <cstdint>
#include <sasktran2/config.h>
#include <sasktran2/atmosphere/grid_storage.h>
#include <sasktran2/atmosphere/surface.h>

namespace sasktran2::atmosphere {
    /** Essentially void base class for the Atmosphere to remove the NSTOKES
     * parameter for SWIG.
     *
     */
    class AtmosphereInterface {
      public:
        virtual ~AtmosphereInterface() {}
    };

    /** Stores all of the atmosphere information for SASKTRAN2.  Essentially
     * this is extinction/single scatter albedo/phase information on a grid that
     * matches the global geometry object, as well as surface parameters.
     * Eventually terms needed for emission sources likely will be added here as
     * well.
     *
     * @tparam NSTOKES
     */
    template <int NSTOKES> class Atmosphere : public AtmosphereInterface {
      private:
        inline static std::atomic<std::uint64_t> s_next_instance_id{1};

        std::shared_ptr<AtmosphereGridStorageFull<NSTOKES>>
            m_storage_holder; /** The internal storage object */
        std::shared_ptr<Surface<NSTOKES>> m_surface_holder; /** The surface */

        AtmosphereGridStorageFull<NSTOKES>& m_storage; /** The internal storage
                                                         object */
        Surface<NSTOKES>& m_surface;                   /** The surface */
        bool m_calculate_derivatives; /** True if we are going to be calculating
                                         derivatives */
        bool m_include_emission_derivatives; /** True if we are going to include
                                                emission derivatives */
        std::uint64_t m_revision = 0; /** Monotonic revision of the built native
                                         atmosphere state. */
        std::uint64_t m_volume_revision =
            0; /** Monotonic revision of volume optical and emission state.
                   Surface-only changes deliberately leave this unchanged so
                   ray optical-depth and solar-transmission caches can remain
                   valid. */
        const std::uint64_t m_instance_id = s_next_instance_id.fetch_add(
            1, std::memory_order_relaxed); /** Stable identity that cannot be
                                               confused by allocator address
                                               reuse. */

      public:
        /** Directly constructs the atmosphere from it's base objects, taking
         * ownership
         *
         * @param storage
         * @param surface
         * @param calculate_derivatives
         */
        Atmosphere(AtmosphereGridStorageFull<NSTOKES>&& storage,
                   Surface<NSTOKES>&& surface,
                   bool calculate_derivatives = false,
                   bool include_emission_derivatives = false);

        /** Directly constructs the atmosphere from it's base objects, sharing
         * ownership
         *
         * @param storage
         * @param surface
         * @param calculate_derivatives
         */
        Atmosphere(AtmosphereGridStorageFull<NSTOKES>& storage,
                   Surface<NSTOKES>& surface,
                   bool calculate_derivatives = false,
                   bool include_emission_derivatives = false);

        /** Constructs an empty atmosphere that we can then modify afterwards
         *
         * @param nwavel
         * @param geometry
         * @param config
         * @param calculate_derivatives
         */
        Atmosphere(int nwavel, const sasktran2::Geometry& geometry,
                   const sasktran2::Config& config,
                   bool calculate_derivatives = false);

        /** Copies the atmosphere view while assigning the new object a unique
         * lifetime identity. Owned storage remains shared, matching the
         * previous implicit-copy behavior. */
        Atmosphere(const Atmosphere& other)
            : m_storage_holder(other.m_storage_holder),
              m_surface_holder(other.m_surface_holder),
              m_storage(other.m_storage), m_surface(other.m_surface),
              m_calculate_derivatives(other.m_calculate_derivatives),
              m_include_emission_derivatives(
                  other.m_include_emission_derivatives),
              m_revision(other.m_revision),
              m_volume_revision(other.m_volume_revision) {}

        /** Moves the atmosphere view while assigning the new object a unique
         * lifetime identity. */
        Atmosphere(Atmosphere&& other) noexcept
            : m_storage_holder(std::move(other.m_storage_holder)),
              m_surface_holder(std::move(other.m_surface_holder)),
              m_storage(other.m_storage), m_surface(other.m_surface),
              m_calculate_derivatives(other.m_calculate_derivatives),
              m_include_emission_derivatives(
                  other.m_include_emission_derivatives),
              m_revision(other.m_revision),
              m_volume_revision(other.m_volume_revision) {}

        Atmosphere& operator=(const Atmosphere&) = delete;
        Atmosphere& operator=(Atmosphere&&) = delete;

        virtual ~Atmosphere() {}

        /** Applies delta_m scaling of a specific order to the internal storage
         * object, overwriting it. Note this is a "half" delta-m scaling.  We
         * scale the extinction/ssa by the regular scaling factors, and then
         * scale the phase function by 1-f.  This completes the TMS single
         * scatter correction.  For multiple scatter it is still necessary to
         * further scale the legendre coefficients.
         *
         * @param order
         */
        void apply_delta_m_scaling(int order);

        const AtmosphereGridStorageFull<NSTOKES>& storage() const {
            return m_storage;
        };
        AtmosphereGridStorageFull<NSTOKES>& storage() { return m_storage; }
        int num_wavel() const { return (int)m_storage.total_extinction.cols(); }

        Surface<NSTOKES>& surface() { return m_surface; }
        const Surface<NSTOKES>& surface() const { return m_surface; }

        // TODO: refactor the below functions into a derivative handler class of
        // some kind
        int ssa_deriv_start_index() const {
            return (int)m_storage.total_extinction.rows();
        }
        int scat_deriv_start_index() const {
            return (int)m_storage.total_extinction.rows() * 2;
        }
        int surface_deriv_start_index() const {
            return scat_deriv_start_index() +
                   num_source_deriv_groups() *
                       m_storage.total_extinction.rows();
        }
        int surface_emission_deriv_start_index() const {
            return surface_deriv_start_index() + m_surface.num_deriv();
        }

        int emission_deriv_start_index() const {
            return scat_deriv_start_index() +
                   num_scattering_deriv_groups() *
                       m_storage.total_extinction.rows();
        }

        /**
         *  Number of internal derivative values (extinction, ssa, phase,
         * albedo, etc)
         */
        int num_deriv() const;

        int num_scattering_deriv_groups() const {
            return m_storage.numscatderiv;
        }

        int num_source_deriv_groups() const {
            return num_scattering_deriv_groups() +
                   int(m_include_emission_derivatives);
        }

        /**
         *  Number of output derivative values, computed from the
         * DerivativeMappings
         */
        int num_output_deriv() const;

        bool include_emission_derivatives() const {
            return m_include_emission_derivatives;
        }

        /** Marks volume and surface native atmosphere state as changed. */
        void mark_changed() {
            ++m_revision;
            ++m_volume_revision;
        }

        /** Marks only native surface state as changed. */
        void mark_surface_changed() { ++m_revision; }

        /** Revision of the native atmosphere state. */
        std::uint64_t revision() const { return m_revision; }

        /** Revision of volume optical and emission state. */
        std::uint64_t volume_revision() const { return m_volume_revision; }

        /** Identity of this native atmosphere lifetime. */
        std::uint64_t instance_id() const { return m_instance_id; }
    };
} // namespace sasktran2::atmosphere
