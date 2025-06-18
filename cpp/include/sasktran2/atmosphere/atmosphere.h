#pragma once

#include "../geometry.h"
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
                   bool calculate_derivatives = false);

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
        Atmosphere(int nwavel, const sasktran2::Geometry1D& geometry,
                   const sasktran2::Config& config,
                   bool calculate_derivatives = false);

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
    };
} // namespace sasktran2::atmosphere
