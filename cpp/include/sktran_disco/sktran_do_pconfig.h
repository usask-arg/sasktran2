#pragma once
#include "sktran_disco/sktran_do.h"
#include "sktran_disco/sktran_do_polarization_types.h"
#include "sktran_do_specs.h"
#include <sasktran2/config.h>
#include <sasktran2/raytracing.h>
#include <sasktran2/geometry.h>

namespace sasktran_disco {

    /**
     * @brief The main DO configuration object
     *
     * TODO: Deprecate this entirely in favor of the sasktran2 config?
     *
     *
     * @tparam NSTOKES
     * @tparam CNSTR
     */
    template <int NSTOKES, int CNSTR = -1>
    class PersistentConfiguration : public BasicProperties<NSTOKES>,
                                    public SolarProperties<NSTOKES>,
                                    public UserSpecProperties {
      public:
        /**
         * @brief Construct a new Persistent Configuration object
         *
         */
        PersistentConfiguration() {}
        ~PersistentConfiguration() {}

      public:
        /**
         * @brief Sets the configuration
         *
         * @param userspecmemory
         * @param config
         * @param cos_sza
         * @param nlyr
         * @param traced_rays
         */
        void configure(
            SKTRAN_DO_UserSpec& userspecmemory, const sasktran2::Config& config,
            double cos_sza, int nlyr,
            const std::vector<sasktran2::raytracing::TracedRay>& traced_rays);

      public: // Configuration getters
        /**
         * @brief The number of streams in the full space
         *
         * @return const uint
         */
        inline const uint nstr() const { return this->M_NSTR; }

        /**
         * @brief The number of lyaers
         *
         * @return const uint
         */
        inline const uint nlyr() const { return this->M_NLYR; }

        /**
         * @brief User spec object
         *
         * @return const SKTRAN_DO_UserSpec*
         */
        inline const SKTRAN_DO_UserSpec* userSpec() const { return m_userspec; }

        /**
         * @brief True if we should perform backpropogation for the BVP
         * derivatives
         *
         * @return true
         * @return false
         */
        inline bool backprop_bvp() const { return this->M_BACKPROP_BVP; }

        /**
         * @brief True if we are in single scatter only debug mode
         *
         * @return true
         * @return false
         */
        inline bool ss_only() const { return this->M_SS_ONLY; }

        /**
         * @brief Solar azimuth angle
         *
         * @return const double
         */
        inline const double saz() const { return this->M_SAZ; }

        /**
         * @brief Cosine of solar zenith angle
         *
         * @return const double
         */
        inline const double csz() const { return this->M_CSZ; }

        /**
         * @brief Memory thread pool
         *
         * @return MemoryPool<NSTOKES, CNSTR>&
         */
        inline MemoryPool<NSTOKES, CNSTR>& pool() const { return m_pool; }

        /**
         * @brief Quadrature weights
         *
         * @return const VectorDim1<double>*
         */
        inline const VectorDim1<double>* quadrature_weights() const {
            return this->M_WT;
        }

        /**
         * @brief Quadrature cos angles
         *
         * @return const VectorDim1<double>*
         */
        inline const VectorDim1<double>* quadrature_cos_angle() const {
            return this->M_MU;
        }

        inline const VectorDim3<LegendrePhaseContainer<NSTOKES>>*
        legendre_streams() const {
            return this->M_LP_MU;
        }

      protected: // Private configuration functions
        void configureModelSpecs(const SKTRAN_DO_UserSpec* userspec);
        void configureLP(const SKTRAN_DO_UserSpec* userspec);

        mutable MemoryPool<NSTOKES, CNSTR> m_pool;
        int m_poolindex;

        /**
         * @brief Legendre polynomials for the cosine of the solar zenith angle
         *
         */
        VectorDim2<LegendrePhaseContainer<NSTOKES>> m_lp_csz_storage;
    };

} // namespace sasktran_disco
