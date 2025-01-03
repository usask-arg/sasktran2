#pragma once
#include <optional>
#include <sasktran2/internal_common.h>

namespace sasktran2 {

    /**
     *  A mapping of the derivatives on the native grid coordinates (i.e. the
     * atmosphere geometry levels) from the internal derivative quantities
     * (d_od, d_ssa, d_legendre, d_brdf, etc.) to a resultant quantity.
     *
     *  Roughly we are performing the calculation:
     *
     *  dI / dx = dI/d_od * d_od/ dx + dI/d_ssa * d_ssa / dx + ...
     *
     *
     */
    struct NativeDerivativeMapping {
      public:
        /**  d_k / dx
         *   Dimension: [geometry, wavel]
         */
        std::optional<Eigen::MatrixXd> d_extinction;

        /**  d_ssa / dx
         *   Dimension: [geometry, wavel]
         */
        std::optional<Eigen::MatrixXd> d_ssa;

        /**  d_legendre / dx
         *   Dimension: [legendre, geometry, wavel]
         */
        std::optional<Eigen::Tensor<double, 3>> d_legendre;

        /**
         *  The scattering index that this mapping corresponds to
         */
        int scat_deriv_index;

        /**
         *  Factors applieds to the scattering derivative.
         *  Dimension: [geometry, wavel]
         */
        std::optional<Eigen::MatrixXd> scat_factor;
    };

    /**
     *  Same as NativeDerivativeMapping but for surface quantities.
     */
    struct NativeSurfaceDerivativeMapping {
      public:
        /**  d_brdf / dx
         *   Dimension: [wavel, brdf_coeff]
         */
        std::optional<Eigen::MatrixXd> d_brdf;
    };

    /**
     *  A full mapping of the internal derivatives to the final output
     * derivative.  Typically this involves a mapping on the native grid
     * coordinates first, and then an optional interpolation step.
     */
    class DerivativeMapping {
      private:
        NativeDerivativeMapping
            m_native_mapping; /** Quantities on the internal atmosphere grid */

        std::optional<Eigen::MatrixXd>
            m_interpolator; /** Optional interpolation matrix from the native
                               grid  */

        int m_nwavel;            /** Number of wavelengths in the storage */
        int m_ninternallocation; /** Number of locations in the storage */
        int m_nlegendre; /** Number of legendre coefficients in the storage */

        std::string m_interp_dim;  /** Constituent name */
        bool m_log_radiance_space; /** Whether the radiance space is in log
                                      space */

        std::string m_assign_name; /** Name of the assign variable */

      public:
        DerivativeMapping(int nwavel, int ninternallocation, int nlegendre);

        NativeDerivativeMapping& native_mapping() { return m_native_mapping; }

        void allocate_extinction_derivatives();
        void allocate_ssa_derivatives();
        void allocate_legendre_derivatives();

        int get_scattering_index() { return m_native_mapping.scat_deriv_index; }
        void set_scattering_index(int index) {
            m_native_mapping.scat_deriv_index = index;
        }

        bool is_scattering_derivative() {
            return m_native_mapping.d_legendre.has_value();
        }

        std::optional<Eigen::MatrixXd>& get_interpolator() {
            return m_interpolator;
        }
        void set_interpolator(const Eigen::MatrixXd& interpolator) {
            m_interpolator = interpolator;
        }

        const std::string& get_interp_dim() { return m_interp_dim; }
        void set_interp_dim(const std::string& name) { m_interp_dim = name; }

        const std::string& get_assign_name() { return m_assign_name; }
        void set_assign_name(const std::string& name) { m_assign_name = name; }

        bool log_radiance_space() { return m_log_radiance_space; }
        void set_log_radiance_space(bool log) { m_log_radiance_space = log; }

        void set_zero();
    };

    class SurfaceDerivativeMapping {
      private:
        NativeSurfaceDerivativeMapping
            m_native_surface_mapping; /** Quantities on the surface grid */

        std::optional<Eigen::MatrixXd>
            m_interpolator; /** Optional interpolation matrix from native
                               wavelength to other */

        int m_nwavel;             /** Number of wavelengths in the storage */
        int m_nbrdf_args;         /** Number of BRDF args in the storage */
        std::string m_interp_dim; /** Interpolating Dimension name */
      public:
        SurfaceDerivativeMapping(int nwavel, int nbrdf_args);

        NativeSurfaceDerivativeMapping& native_surface_mapping() {
            return m_native_surface_mapping;
        }

        void allocate_brdf_derivatives();

        std::optional<Eigen::MatrixXd>& get_interpolator() {
            return m_interpolator;
        }
        void set_interpolator(const Eigen::MatrixXd& interpolator) {
            m_interpolator = interpolator;
        }

        const std::string& get_interp_dim() { return m_interp_dim; }
        void set_interp_dim(const std::string& name) { m_interp_dim = name; }

        void set_zero();
    };

}; // namespace sasktran2
