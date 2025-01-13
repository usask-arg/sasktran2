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
         *
         *   Note, this quantitiy is only used for storage to set the atmosphere
         * storage d_legendre parameters later on
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

        /**
         *  The emission derivative
         *  Dimension: [geometry, wavel]
         *
         *   Note, this quantity is only used when the Emission source term is
         * included in the engine
         */
        std::optional<Eigen::MatrixXd> d_emission;
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

        /**
         * d_emission / dx
         * Dimension: [wavel, 1]
         *
         */
        std::optional<Eigen::MatrixXd> d_emission;
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
        /**
         * @brief Construct a new Derivative Mapping object
         *
         * @param nwavel Number of wavelengths
         * @param ninternallocation Number of internal geometry locations
         * @param nlegendre Number of legendre coefficients (in the atmosphere
         * storage)
         */
        DerivativeMapping(int nwavel, int ninternallocation, int nlegendre);

        /**
         * @brief Returns the native derivative mapping
         *
         * @return NativeDerivativeMapping&
         */
        NativeDerivativeMapping& native_mapping() { return m_native_mapping; }

        /**
         * @brief Returns the native derivative mapping
         *
         * @return const NativeDerivativeMapping&
         */
        const NativeDerivativeMapping& native_mapping() const {
            return m_native_mapping;
        }

        /**
         * @brief Allocates the extinction derivative object
         *
         */
        void allocate_extinction_derivatives();

        /**
         * @brief Allocates the ssa derivative object
         *
         */
        void allocate_ssa_derivatives();

        /**
         * @brief Allocates the legendre derivative object
         *
         */
        void allocate_legendre_derivatives();

        /**
         * @brief Allocates the emission derivative object
         *
         */
        void allocate_emission_derivatives();

        /**
         * @brief Get the scattering index
         *
         * @return int
         */
        int get_scattering_index() const {
            return m_native_mapping.scat_deriv_index;
        }

        /**
         * @brief Set the scattering index
         *
         * @param index
         */
        void set_scattering_index(int index) {
            m_native_mapping.scat_deriv_index = index;
        }

        /**
         * @return true If the derivative mapping is a scattering derivative
         * @return false If the derivative mapping is not a scattering
         * derivative
         */
        bool is_scattering_derivative() const {
            return m_native_mapping.d_legendre.has_value();
        }

        /**
         * @brief Get the interpolator object
         *
         * @return std::optional<Eigen::MatrixXd>&
         */
        std::optional<Eigen::MatrixXd>& get_interpolator() {
            return m_interpolator;
        }

        /**
         * @brief Get the interpolator const object
         *
         * @return const std::optional<Eigen::MatrixXd>&
         */
        const std::optional<Eigen::MatrixXd>& get_interpolator_const() const {
            return m_interpolator;
        }

        /**
         * @brief Set the interpolator object
         *
         * @param interpolator
         */
        void set_interpolator(const Eigen::MatrixXd& interpolator) {
            m_interpolator = interpolator;
        }

        /**
         * @brief Get the interp dim
         *
         * The interpolating dimension name is only used by the python interface
         *
         * @return const std::string&
         */
        const std::string& get_interp_dim() { return m_interp_dim; }

        /**
         * @brief Set the interp dim
         *
         * @param name
         */
        void set_interp_dim(const std::string& name) { m_interp_dim = name; }

        /**
         * @brief Get the assign name
         *
         * The assign name is only used by the python interface
         *
         * @return const std::string&
         */
        const std::string& get_assign_name() { return m_assign_name; }

        /**
         * @brief Set the assign name
         *
         * @param name
         */
        void set_assign_name(const std::string& name) { m_assign_name = name; }

        /**
         * @brief True if the derivative should be considered in log radiance
         * space (divide by the radiance at the end)
         *
         * @return true
         * @return false
         */
        bool log_radiance_space() const { return m_log_radiance_space; }

        /**
         * @brief Set the log radiance space
         *
         * @param log
         */
        void set_log_radiance_space(bool log) { m_log_radiance_space = log; }

        /**
         * @brief Resets all internal elements to 0. (Does not delete them if
         * they exist)
         *
         */
        void set_zero();

        /**
         * @brief The number of output values in the final derivative (per
         * wavelength/los)
         *
         * @return int
         */
        int num_output() const {
            if (m_interpolator.has_value()) {
                return m_interpolator.value().cols();
            } else {
                return m_ninternallocation;
            }
        }
    };

    /**
     * @brief maps derivatives of the surface quantities
     *
     * Unlike the DerivativeMapping class this only deals with native
     * coordinates.  An interpolation object is stored to be used by the Python
     * code, but not used internally
     *
     */
    class SurfaceDerivativeMapping {
      private:
        NativeSurfaceDerivativeMapping
            m_native_surface_mapping; /** Quantities on the surface grid */

        std::optional<Eigen::MatrixXd>
            m_interpolator; /** Optional interpolation matrix from native
                               wavelength to other, only used in the Python code
                             */

        int m_nwavel;             /** Number of wavelengths in the storage */
        int m_nbrdf_args;         /** Number of BRDF args in the storage */
        std::string m_interp_dim; /** Interpolating Dimension name */
      public:
        /**
         * @brief Construct a new Surface Derivative Mapping
         *
         * @param nwavel Number of wavelengths
         * @param nbrdf_args Number of arguments the BRDF takes
         */
        SurfaceDerivativeMapping(int nwavel, int nbrdf_args);

        /**
         * @brief The Native mapping object
         *
         * @return NativeSurfaceDerivativeMapping&
         */
        NativeSurfaceDerivativeMapping& native_surface_mapping() {
            return m_native_surface_mapping;
        }

        /**
         * @brief The native mapping object
         *
         * @return const NativeSurfaceDerivativeMapping&
         */
        const NativeSurfaceDerivativeMapping& native_surface_mapping() const {
            return m_native_surface_mapping;
        }

        /**
         * @brief Allocates derivatives with respect to the BRDF
         *
         */
        void allocate_brdf_derivatives();

        /**
         * @brief Allocates derivatives with respect to surface emission
         *
         */
        void allocate_emission_derivatives();

        /**
         * @brief Get the interpolator object
         *
         * @return std::optional<Eigen::MatrixXd>&
         */
        std::optional<Eigen::MatrixXd>& get_interpolator() {
            return m_interpolator;
        }

        /**
         * @brief Set the interpolator object
         *
         * @param interpolator
         */
        void set_interpolator(const Eigen::MatrixXd& interpolator) {
            m_interpolator = interpolator;
        }

        /**
         * @brief Get the interp dim object
         *
         * Only used by the Python code
         *
         * @return const std::string&
         */
        const std::string& get_interp_dim() { return m_interp_dim; }

        /**
         * @brief Set the interp dim object
         *
         * @param name
         */
        void set_interp_dim(const std::string& name) { m_interp_dim = name; }

        /**
         * @brief Sets all internal storage to 0
         *
         */
        void set_zero();
    };

}; // namespace sasktran2
