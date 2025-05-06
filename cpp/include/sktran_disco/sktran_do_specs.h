#pragma once

#include "sktran_disco/sktran_do.h"

namespace sasktran_disco {
    /*!
     * @class SKTRAN_DO_UserSpec
     * @ingroup SASKTRAN-DO
     *
     * @brief An object which stores the users settings for a
     * SKTRAN_DO_Engine instance.
     *
     * @throws An invalid configuration will likely throw std::invalid_argument.
     *
     */
    class SKTRAN_DO_UserSpec {
      public: // Constructors and destructor
        /*!
         * Default construct. Configure must be called before this object can be
         * used.
         */
        SKTRAN_DO_UserSpec();

        /*!
         * Default destructor.
         */
        virtual ~SKTRAN_DO_UserSpec() {}

      public: // Setters - typically called by user
        /*!
         * Minimum required user configuration.
         *
         * @param[in] num_streams The number of streams used in the calculation.
         * \p num_streams must be an even number greater than or equal to 4.
         * @param[in] num_atmo_layers The number of homogeneous layers used to
         * approximate the atmosphere. \p num_atmo_layers Must be greater than
         * 0.
         * @param[in] solar_position The solar position vector.
         * @param[in] solar_direct_intensity The intensity of the sun at the
         * TOA. The default value is 1.0.
         *
         * @pre None.
         * @post This specification now meets the absolute minimum required
         * configuration for a SKTRAN_DO_Engine instance. Additional
         * configuration after is fine.
         * @warning Repeated calls can diminish performance. This is because a
         * call to this function generates some cached data.
         */
        void configure(unsigned int num_streams, unsigned int num_atmo_layers,
                       double solar_direct_intensity = 1.0);

        /*!
         * Sets the number of streams used in the calculation. The number of
         * streams corresponds to the number of quadrature angles used to
         * calculate the numeric integral of the source function in the RTE. An
         * increased number of streams can improve precision but is also
         * computationally more expensive.
         *
         * @param[in] num_streams The number of streams used in the calculation.
         * \p num_streams must be an even number greater than or equal to 4 and
         * less than or equal to 40. If you would like to use more than 40
         * streams you must supply your own. See setNumberOfStreams.
         *
         * @warning Using an insufficient number of streams can cause the
         * calculation to fail. In this event SKTRAN_DO_Engine will safely exit
         * and report the error.
         */
        SKTRAN_DO_UserSpec* setNumberOfStreams(unsigned int num_streams);

        /*!
         * Sets the number of homogeneous layers used to approximate the
         * atmosphere. An increased number of layers will improve precision but
         * is also more computationally expensive.
         *
         * @param[in] num_atmo_layers The number of homogeneous layers used to
         * approximate the atmosphere. \p num_atmo_layers Must be greater than
         * 0.
         */
        SKTRAN_DO_UserSpec* setNumberOfLayers(unsigned int num_atmo_layers);

        /*!
         * Sets the incident intensity at the top of the atmosphere.
         *
         * @param[in] solar_direct_intensity The intensity of the sun at the
         * TOA. default value is 0.0.
         */
        SKTRAN_DO_UserSpec* setTOAIntensities(double direct = 1);

        /*!
         * Sets a debug option to only include single scatter source terms.
         *
         * @param[in] ss_only If true, only single scatter source terms will be
         * included
         */
        SKTRAN_DO_UserSpec* setSSOnly(bool ss_only = false);

        /*!
         * Sets the single-scatter dither amount in the event a purely
         * scattering atmospheric layer is encountered. This is a special case
         * which discrete- ordinate algorithms must handle by dithering the SSA.
         *
         * @param[in] dither The amount to dither the SSA by.
         */
        SKTRAN_DO_UserSpec* setSSAEqual1Dither(double dither = 1e-9);

      public:
      public: // Getters - typically called by the internals of SASKTRAN-DO
        uint getNumberOfStreams() const;
        uint getNumberOfLayers() const;
        const VectorDim1<double>* getStreamAbscissae() const;
        const VectorDim1<double>* getStreamWeights() const;
        const VectorDim3<LegendrePhaseContainer<4>>*
        getAbscissaeLegendreP4() const;
        const VectorDim3<LegendrePhaseContainer<1>>*
        getAbscissaeLegendreP1() const;
        double getTopDirectIntensity() const;
        double getSSAEqual1Dither() const;
        void configureDefaultDetails();
        void cacheLPOfStreamAngles();
        bool getSSOnly() const { return m_ss_only; }

      private: // Members
        // Gaussian Quadrature
        // Just store the full stokes vector components, shouldn't be that much
        // extra calculation
        VectorDim3<LegendrePhaseContainer<4>> m_lp_abscissae4;
        VectorDim3<LegendrePhaseContainer<1>> m_lp_abscissae1;

        VectorDim1<double> m_abscissae;
        VectorDim1<double> m_weights;

        // Incident Intensity Specs
        double m_itop_direct;

        // Discrete-Ordinate Method Configuration
        uint m_nstr;
        uint m_nlyr;

        bool m_ss_only;

        // Dither for handling single-scatter albedo = 1 special case
        double m_ssalb1_dither;
    };
} // namespace sasktran_disco
