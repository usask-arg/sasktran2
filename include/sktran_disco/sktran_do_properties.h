#pragma once
#include "sktran_disco/sktran_do.h"
#include "sktran_do_specs.h"

namespace sasktran_disco {
    /**
     * An object which inherits the requested read-only properties and safely
     * copies them to this object. This object should be inherited by an
     * object which would like access to read-only properties.
     *
     * @tparam Ts
     */
    template <class... Ts> class ReadOnlyProperties : public Ts... {
      protected:
        ReadOnlyProperties() = delete;

        // Only constructible from a super property.
        template <class SuperProperty>
        ReadOnlyProperties(const SuperProperty& super)
            : Ts(static_cast<const Ts*>(&super))... {}
    };

#pragma region "Read only properties"
    // Requirements for a properties:
    // Public: Copy constructor from const pointer
    // Protected: Data members. Default constructor.

    /**
     * @brief Basic properties
     *
     * @tparam NSTOKES
     * @tparam CNSTR
     */
    template <int NSTOKES, int CNSTR = -1> class BasicProperties {
      public:
        BasicProperties(const BasicProperties* other)
            : M_NSTR(other->M_NSTR), M_NLYR(other->M_NLYR), M_MU(other->M_MU),
              M_WT(other->M_WT), M_LP_MU(other->M_LP_MU),
              M_BACKPROP_BVP(other->M_BACKPROP_BVP),
              M_SS_ONLY(other->M_SS_ONLY) {
            // empty
        }

      protected:
        BasicProperties()
            : M_NSTR(0), M_NLYR(0), M_MU(nullptr), M_WT(nullptr),
              M_LP_MU(nullptr), M_BACKPROP_BVP(false), M_SS_ONLY(false) {
            // empty
        }
        const uint M_NSTR; /** Number of streams*/

        const uint M_NLYR; /** Number of atmospheric layers*/

        const bool M_BACKPROP_BVP; /** Whether to enable backprop to calculate
                                      the BVP derivatives */

        const bool M_SS_ONLY; /** Debug option to only include SS source terms*/

        const VectorDim1<double>*
            M_MU; /** Pointer to a vector of stream cos angles */

        const VectorDim1<double>*
            M_WT; /** Pointer to a vector of quadrature weights*/

        // Legendre polynomials evaluated at stream angles. Access is [m][i][l]
        // where m is AEOrder, i is stream index, l is polynomial order.
        const VectorDim3<LegendrePhaseContainer<NSTOKES>>* M_LP_MU;
    };

    /**
     * @brief Properties related to the solar beam
     *
     * @tparam NSTOKES
     * @tparam CNSTR
     */
    template <int NSTOKES, int CNSTR = -1> class SolarProperties {
      public:
        SolarProperties(const SolarProperties* other)
            : M_CSZ(other->M_CSZ), M_SAZ(other->M_SAZ),
              M_SOLAR_REFRACTION(other->M_SOLAR_REFRACTION),
              M_LP_CSZ(other->M_LP_CSZ) {
            // empty
        }

      protected:
        SolarProperties()
            : M_CSZ(std::nan("1")), M_SAZ(std::nan("1")),
              M_SOLAR_REFRACTION(false), M_LP_CSZ(nullptr) {
            // empty
        }

        const bool M_SOLAR_REFRACTION; /** Whether to include solar refraction
                                         in the calculations */
        const double M_CSZ;            /** Cosine solar zenith angle */
        const double M_SAZ;            /** Solar azimuth angle in radians*/
        // Legendre polynomials evaluated at M_CSZ. Accessed by [m][l] where
        // m is the AEOrder and l is the polynomial order.
        VectorDim2<LegendrePhaseContainer<NSTOKES>>* M_LP_CSZ;
    };

    /**
     * @brief Stores properties related to the user spec
     *
     */
    class UserSpecProperties {
      public:
        UserSpecProperties(const UserSpecProperties* other)
            : m_userspec(other->m_userspec) {
            // empty
        }

      protected:
        UserSpecProperties() : m_userspec(nullptr) {
            // empty
        }
        // Pointer to user specifications.
        const SKTRAN_DO_UserSpec* m_userspec;
    };

#pragma endregion
} // namespace sasktran_disco
