#pragma once
#include "sktran_disco/sktran_do.h"
#include "sktran_disco/sktran_do_linearization_types.h"

namespace sasktran_disco {

    /**
     * @brief A representation of a single layer in the RTE solution
     *
     * @tparam NSTOKES
     * @tparam CNSTR
     */
    template <int NSTOKES, int CNSTR = -1>
    using OpticalLayerROP =
        ReadOnlyProperties<BasicProperties<NSTOKES>, SolarProperties<NSTOKES>,
                           UserSpecProperties>;
    template <int NSTOKES, int CNSTR = -1>
    class OpticalLayer : public OpticalLayerROP<NSTOKES> {
        using HomogType = typename std::conditional<NSTOKES != 5, double,
                                                    std::complex<double>>::type;

      public:
        /**
         * @brief Construct a new Optical Layer object without setting the
         * optical properties
         *
         * @param config
         * @param index
         * @param altitude_ceiling
         * @param altitude_floor
         * @param input_derivs
         * @param thread_data
         */
        OpticalLayer(const PersistentConfiguration<NSTOKES, CNSTR>& config,
                     LayerIndex index, double altitude_ceiling,
                     double altitude_floor,
                     const InputDerivatives<NSTOKES>& input_derivs,
                     const ThreadData<NSTOKES, CNSTR>& thread_data);

        /**
         * @brief Set the optical properties
         *
         * @param scat_ext
         * @param tot_ext
         * @param phasef_expansion
         * @param optical_depth_ceiling
         * @param optical_depth_floor
         */
        void set_optical(
            double scat_ext, double tot_ext,
            const VectorDim1<sasktran_disco::LegendreCoefficient<NSTOKES>>&
                phasef_expansion,
            double optical_depth_ceiling, double optical_depth_floor);

        /**
         * @brief Construct a new Optical Layer object completely
         *
         * @param config
         * @param index
         * @param scat_ext
         * @param tot_ext
         * @param phasef_expansion
         * @param optical_depth_ceiling
         * @param optical_depth_floor
         * @param altitude_ceiling
         * @param altitude_floor
         * @param input_derivs
         */
        OpticalLayer(
            const PersistentConfiguration<NSTOKES, CNSTR>& config,
            LayerIndex index, double scat_ext, double tot_ext,
            std::unique_ptr<
                VectorDim1<sasktran_disco::LegendreCoefficient<NSTOKES>>>
                phasef_expansion,
            double optical_depth_ceiling, double optical_depth_floor,
            double altitude_ceiling, double altitude_floor,
            const InputDerivatives<NSTOKES>& input_derivs);

        /**
         *  Called after the optical properties are set, does things like
         * constructs duals of the SSA and Thickness objects
         *
         */
        void configureDerivative();

        /**
         * @brief The single scatter albedo for the layer
         *
         * @return double
         */
        inline double ssa() const { return M_SSA; }

        /**
         * If location is set to FLOOR or CEILING, returns the optical depth at
         * the bottom or middle of layer.  If set to INSIDE then return the
         * optical thickness of the layer
         *
         * @param loc
         * @return double
         */
        inline double opticalDepth(Location loc) const {
            switch (loc) {
            case Location::FLOOR:
                return M_OPTICALDEPTH_FLOOR;
                break;
            case Location::CEILING:
                return M_OPTICALDEPTH_CEILING;
                break;
            case Location::INSIDE:
                return M_OPTICAL_THICKNESS;
                break;
            default:
                abort();
            }
        }

        /**
         * Returns either the floor or ceiling altitude of the layer
         *
         * @param loc
         * @return double
         */
        inline double altitude(Location loc) const {
            switch (loc) {
            case Location::FLOOR:
                return M_ALT_FLOOR;
                break;
            case Location::CEILING:
                return M_ALT_CEILING;
                break;
            default:
                abort();
            }
        }

        /**
         * TODO: REFACTOR
         *
         * @param loc
         * @param x
         * @return double
         */
        inline double beamTransmittance(Location loc, double x = -1) const {
            switch (loc) {
            case Location::FLOOR:
                return floor_beam_transmittance().value;
                break;
            case Location::CEILING:
                return ceiling_beam_transmittanc().value;
                break;
            default:
                abort();
            }
        }

        /**
         * TODO: REFACTOR
         *
         * @param loc
         * @param deriv
         * @param derivindex
         * @param x
         * @return double
         */
        inline double
        d_beamTransmittance(Location loc,
                            const LayerInputDerivative<NSTOKES>& deriv,
                            uint derivindex, double x = -1) const {
            LayerIndex d_idx = deriv.layer_index;
            LayerIndex cur_layer = M_INDEX;

            switch (loc) {
            case Location::FLOOR:
                return floor_beam_transmittance().deriv(derivindex);
                break;
            case Location::CEILING:
                return ceiling_beam_transmittanc().deriv(derivindex);
                break;
            default:
                abort();
            }
        }

        /**
         * TODO: REFACTOR
         *
         * @param loc
         * @param deriv
         * @param x
         * @return Dual<double>
         */
        inline Dual<double>
        dual_beamTransmittance(Location loc,
                               const InputDerivatives<NSTOKES>& deriv,
                               double x = -1) const {
            Dual<double> result(deriv.numDerivative());

            result.value = beamTransmittance(loc, x);

            for (uint i = 0; i < deriv.numDerivative(); ++i) {
                result.deriv[i] =
                    d_beamTransmittance(loc, deriv.layerDerivatives()[i], i, x);
            }

            return result;
        }

        /**
         * TODO: REFACTOR
         *
         * @param loc
         * @param m
         * @param j
         * @return HomogType
         */
        inline HomogType streamTransmittance(Location loc, AEOrder m,
                                             SolutionIndex j) const {
            if (loc != Location::INSIDE) {
                abort();
            }
            return exp(-std::abs(m_solutions[m].value.eigval(j)) *
                       opticalDepth(Location::INSIDE));
        }

        /**
         * TODO: REFACTOR
         *
         * @param loc
         * @param m
         * @param j
         * @param derivIndex
         * @param deriv
         * @return HomogType
         */
        inline HomogType d_streamTransmittance(
            Location loc, AEOrder m, SolutionIndex j, uint derivIndex,
            const LayerInputDerivative<NSTOKES>& deriv) const {
            if (loc != Location::INSIDE) {
                abort();
            }
            double d_od = deriv.d_optical_depth;

            return -(m_solutions[m].value.dual_eigval().deriv(derivIndex, j) *
                         M_OPTICAL_THICKNESS +
                     m_solutions[m].value.eigval(j) * d_od) *
                   streamTransmittance(loc, m, j);
        }

        /**
         * TODO: REFACTOR
         *
         * @param loc
         * @param m
         * @param j
         * @param deriv
         * @return Dual<HomogType>
         */
        inline Dual<HomogType>
        dual_streamTransmittance(Location loc, AEOrder m, SolutionIndex j,
                                 const InputDerivatives<NSTOKES>& deriv) const {
            // TODO: Layerdual instead?

            size_t derivStart = deriv.layerStartIndex(M_INDEX);
            if (loc != Location::INSIDE) {
                abort();
            }
            Dual<HomogType> result(deriv.numDerivative());
            result.value = streamTransmittance(loc, m, j);
            // Have no cross derivatives
            for (uint i = 0; i < deriv.numDerivativeLayer(M_INDEX); ++i) {
                result.deriv[i + derivStart] = d_streamTransmittance(
                    loc, m, j, i, deriv.layerDerivatives()[derivStart + i]);
            }

            return result;
        }

        /**
         * @brief Index of the layer, 0 = TOA
         *
         * @return LayerIndex
         */
        inline LayerIndex index() const { return M_INDEX; }

        /**
         * @brief The solutions for azimuth order m
         *
         * @param m
         * @return const LayerSolution<NSTOKES, CNSTR>&
         */
        inline const LayerSolution<NSTOKES, CNSTR>& solution(AEOrder m) const {
            return m_solutions[m];
        }

        /**
         * @brief The solutions for azimuth order M
         *
         * @param m
         * @return LayerSolution<NSTOKES, CNSTR>&
         */
        inline LayerSolution<NSTOKES, CNSTR>& solution(AEOrder m) {
            return m_solutions[m];
        }

        /**
         * @brief Scattering extinction for the layer
         *
         * @return double
         */
        inline double scatExt() const { return M_SCAT_EXT; }

        /**
         * @brief Total extinction for the layer
         *
         * @return double
         */
        inline double totalExt() const { return M_TOT_EXT; }

        /**
         * @brief Dual of the layer optical thickness
         *
         * @return const LayerDual<double>&
         */
        inline const LayerDual<double>& dual_thickness() const {
            return m_dual_thickness;
        }

        /**
         * @brief Dual of the layer single scatter albedo
         *
         * @return const LayerDual<double>&
         */
        inline const LayerDual<double>& dual_ssa() const { return m_dual_ssa; }

        /**
         * @brief Dual of the average secant in the layer
         *
         * @return const Dual<double>&
         */
        inline const Dual<double>& dual_average_secant() const {
            return m_average_secant;
        }

        /**
         * @brief Dual of the beam transittance at the top of the layer
         *
         * @return const Dual<double>&
         */
        inline const Dual<double>& ceiling_beam_transmittanc() const {
            return *m_bt_ceiling;
        }

        /**
         * @brief Dual of the beam transmittance at the top of the layer
         *
         * @return Dual<double>&
         */
        inline Dual<double>& ceiling_beam_transmittance() const {
            return *m_bt_ceiling;
        }

        /**
         * @brief Dual of the beam transmittance at the bottom of the layer
         *
         * @return Dual<double>&
         */
        inline Dual<double>& floor_beam_transmittance() const {
            return *m_bt_floor;
        }

        inline void set_transmittances(Dual<double>& bt_ceiling,
                                       Dual<double>& bt_floor) {
            m_bt_ceiling = &bt_ceiling;
            m_bt_floor = &bt_floor;
        }

        /**
         * @brief Dual of the average secant in the layer
         *
         * @return Dual<double>&
         */
        inline Dual<double>& dual_average_secant() { return m_average_secant; }

        /**
         * @brief Legendre coefficients for the phase function expansion
         *
         * @return const std::vector<LegendreCoefficient<NSTOKES>>&
         */
        inline const std::vector<LegendreCoefficient<NSTOKES>>&
        legendre_coeff() const {
            return *m_lephasef;
        }

        /**
         * TODO: Move to standalone function
         *
         * @param m
         * @param mu
         * @param obsod
         * @param lp_mu
         * @param result
         * @param reverse_linearization_trace
         * @param manual_ss_source
         * @param include_ss
         */
        void integrate_source(
            AEOrder m, double mu, double obsod,
            const std::vector<LegendrePhaseContainer<NSTOKES>>& lp_mu,
            Radiance<NSTOKES>& result,
            ReverseLinearizationTrace<NSTOKES>& reverse_linearization_trace,
            const sasktran_disco::Radiance<NSTOKES>* manual_ss_source = nullptr,
            bool include_ss = true) const;

        /**
         * @brief TODO: Move to standalone function
         *
         * @param m
         * @param mu
         * @param x
         * @param thickness
         * @param j
         * @param xform
         */
        void h_plus(AEOrder m, double mu, double x, double thickness,
                    SolutionIndex j, LayerDual<HomogType>& xform) const;

        /**
         * @brief TODO: Move to standalone function
         *
         * @param m
         * @param mu
         * @param x
         * @param thickness
         * @param j
         * @param xform
         */
        void h_minus(AEOrder m, double mu, double x, double thickness,
                     SolutionIndex j, LayerDual<HomogType>& xform) const;

        /**
         * @brief TODO: Move to standalone function
         *
         * @param mu
         * @param x
         * @param thickness
         * @param transmission
         * @param xform
         */
        void E(double mu, double x, double thickness,
               const Dual<double>& transmission, Dual<double>& xform) const;

      protected:
        double M_SSA;                /** Layer single scatter albedo */
        double M_SCAT_EXT;           /** Layer scattering extinction*/
        double M_TOT_EXT;            /** Layer total extinction */
        double M_OPTICALDEPTH_FLOOR; /** Optical depth at the bottom of the
                                        layer */
        double
            M_OPTICALDEPTH_CEILING; /** Optical depth at the top of the layre */
        double M_OPTICAL_THICKNESS; /** Optical thickness of the layer */
        const double M_ALT_CEILING; /** Altitude at the top of the layer */
        const double M_ALT_FLOOR;   /** Altitude at the bottom of the layer */
        std::unique_ptr<
            VectorDim1<sasktran_disco::LegendreCoefficient<NSTOKES>>>
            m_lephasef;           /** Legendre phase expansion */
        const LayerIndex M_INDEX; /** Layer index */

        LayerCache<NSTOKES>& m_layercache; /** Layer cache */

        PostProcessingCache<NSTOKES>&
            m_postprocessing_cache; /** Cache for post processing */

        std::vector<LayerSolution<NSTOKES, CNSTR>>&
            m_solutions; /** Layer solutions for each azimuth */
        const InputDerivatives<NSTOKES>&
            m_input_derivs;                  /** User input derivatives */
        LayerDual<double>& m_dual_thickness; /** Dual of the optical thickness*/
        LayerDual<double>& m_dual_ssa;  /** Dual of the single scatter albedo */
        Dual<double>& m_average_secant; /** Dual of the average secant */

        Dual<double>* m_bt_ceiling; /** Dual of transmission at the top of
                                             the layer */
        Dual<double>* m_bt_floor;   /** Dual of transmission at the bottom of
                                             the layer */
    };

} // namespace sasktran_disco
