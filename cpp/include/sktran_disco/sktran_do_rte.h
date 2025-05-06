#pragma once
#include "sktran_disco/sktran_do.h"
#include "sktran_disco/sktran_do_linearization_types.h"

namespace sasktran_disco {

    // RTESolver is the heart of SASKATRAN-DO. This class is a calculation
    // object which calculated the discrete-ordinate solutions to the RTE.
    // Note that this class uses the boost::pool_allocator to improve
    // allocation/deallocation performance.
    template <int NSTOKES, int CNSTR = -1>
    using RTESProperties =
        ReadOnlyProperties<BasicProperties<NSTOKES>, SolarProperties<NSTOKES>,
                           UserSpecProperties>;
    template <int NSTOKES, int CNSTR = -1>
    class RTESolver : public RTESProperties<NSTOKES> {
        using HomogType = typename std::conditional<NSTOKES != 5, double,
                                                    std::complex<double>>::type;

        using Matrix = typename std::conditional<
            CNSTR != -1,
            Eigen::Matrix<double, CNSTR / 2 * NSTOKES, CNSTR / 2 * NSTOKES>,
            Eigen::MatrixXd>::type;
        using MatrixH = typename std::conditional<
            CNSTR != -1,
            Eigen::Matrix<HomogType, CNSTR / 2 * NSTOKES, CNSTR / 2 * NSTOKES>,
            Eigen::MatrixXd>::type;
        using MatrixView = typename Eigen::Map<Matrix>;
        using MatrixViewH = typename Eigen::Map<MatrixH>;

        using Vector = typename std::conditional<
            CNSTR != -1, Eigen::Vector<double, CNSTR / 2 * NSTOKES>,
            Eigen::VectorXd>::type;
        using VectorH = typename std::conditional<
            CNSTR != -1, Eigen::Vector<HomogType, CNSTR / 2 * NSTOKES>,
            Eigen::VectorXd>::type;

        using VectorViewH = typename Eigen::Map<VectorH>;

        using MatrixHLHS = typename std::conditional<
            CNSTR != -1,
            Eigen::Matrix<double, CNSTR / 2 * NSTOKES + 1,
                          CNSTR / 2 * NSTOKES + 1>,
            Eigen::MatrixXd>::type;
        using MatrixHRHS = typename std::conditional<
            CNSTR != -1, Eigen::Matrix<double, CNSTR / 2 * NSTOKES + 1, -1>,
            Eigen::MatrixXd>::type;

      public:
        // Only RTESolver constructor. Enforcing this constructor means that
        // an RTESolver object is always in a valid state and thus calls to
        // solve the orders of the azimuth expansion are always valid.
        RTESolver(const PersistentConfiguration<NSTOKES, CNSTR>& config,
                  OpticalLayerArray<NSTOKES, CNSTR>& layers);

        // Solves order m of the azimuth expansion. After this call, solutions
        // stored in layers are valid for azimuth expansion order m.
        void solve(AEOrder m);

        void backprop(AEOrder m, ReverseLinearizationTrace<NSTOKES>& trace,
                      sasktran_disco::Radiance<NSTOKES>& component);

        // Frees all pooled memory.
        ~RTESolver() {}

      private: // Helper functions
        // Calculates the homogeneous solution in layer p for azimuth order m.
        void solveHomogeneous(AEOrder m, OpticalLayer<NSTOKES, CNSTR>& layer);
        void linearizeHomogeneous(AEOrder m,
                                  OpticalLayer<NSTOKES, CNSTR>& layer);

        void assignHomogenousSplusMinus(AEOrder m,
                                        OpticalLayer<NSTOKES, CNSTR>& layer);

        void assignParticularQ(AEOrder m,
                               const OpticalLayer<NSTOKES, CNSTR>& layer,
                               VectorLayerDual<double>& Qplus,
                               VectorLayerDual<double>& Qminus);

        // Calculates the particular solution in layer p for azimuth order m.
        void solveParticularGreen(AEOrder m,
                                  OpticalLayer<NSTOKES, CNSTR>& layer);

        // Calculates the homogeneous coefficients which satisfy boundary
        // conditions.
        void solveBVP(AEOrder m);

        // Configures the cache
        void configureCache();

        using BVPMatrix = la::BVPMatrix<NSTOKES>;

        // Boundary Conditions for the coefficient matrix
        void
        bvpTOACondition(AEOrder m, BoundaryIndex p, BVPMatrix& A,
                        VectorDim1<BVPMatrixDenseBlock<NSTOKES>>& d_A) const;
        void bvpContinuityCondition(
            AEOrder m, BoundaryIndex p, BVPMatrix& A,
            VectorDim1<BVPMatrixDenseBlock<NSTOKES>>& d_A) const;
        void
        bvpGroundCondition(AEOrder m, BoundaryIndex p, BVPMatrix& A,
                           VectorDim1<BVPMatrixDenseBlock<NSTOKES>>& d_A) const;

        // Boundary conditions for the RHS of the boundary value problem
        void bvpCouplingCondition_BC1(AEOrder m, BoundaryIndex p, uint& loc,
                                      Eigen::VectorXd& b, Eigen::MatrixXd& d_b);
        void bvpCouplingCondition_BC2(AEOrder m, BoundaryIndex p, uint& loc,
                                      Eigen::VectorXd& b, Eigen::MatrixXd& d_b);
        void bvpCouplingCondition_BC3(AEOrder m, BoundaryIndex p, uint& loc,
                                      Eigen::VectorXd& b, Eigen::MatrixXd& d_b);

        inline HomogType v_plus(AEOrder m,
                                const OpticalLayer<NSTOKES, CNSTR>& layer,
                                StreamIndex j, SolutionIndex a) const {
            // TODO: Polarized surface
            int s1 = j % NSTOKES;

            HomogType sum = layer.solution(m).value.homog_plus(j, a);
            if ((m >= m_layers.surface().sk2_surface().max_azimuthal_order()) ||
                s1 != 0)
                return sum; // rho will be zero
            auto& rho = m_layers.surface().storage().brdf.stream_stream;
            for (StreamIndex i = 0; i < this->M_NSTR / 2; ++i) {
                sum -= (1 + kronDelta(m, 0)) * rho(j / NSTOKES, i) *
                       (*this->M_WT)[i] * (*this->M_MU)[i] *
                       layer.solution(m).value.homog_minus(i * NSTOKES, a);
            }
            return sum;
        }

        inline HomogType v_minus(AEOrder m,
                                 const OpticalLayer<NSTOKES, CNSTR>& layer,
                                 StreamIndex j, SolutionIndex a) const {
            // TODO: Polarized surface
            int s1 = j % NSTOKES;

            HomogType sum = layer.solution(m).value.homog_minus(j, a);
            if ((m >= m_layers.surface().sk2_surface().max_azimuthal_order()) ||
                s1 != 0)
                return sum; // rho will be zero
            const auto& rho = m_layers.surface().storage().brdf.stream_stream;
            for (StreamIndex i = 0; i < this->M_NSTR / 2; ++i) {
                sum -= (1 + kronDelta(m, 0)) * rho(j / NSTOKES, i) *
                       (*this->M_WT)[i] * (*this->M_MU)[i] *
                       layer.solution(m).value.homog_plus(i * NSTOKES, a);
            }
            return sum;
        }

        inline HomogType
        d_v_plus(AEOrder m, const OpticalLayer<NSTOKES, CNSTR>& layer,
                 StreamIndex j, SolutionIndex a, uint derivindex,
                 const LayerInputDerivative<NSTOKES>& deriv) const {
            // TODO: Polarized surface
            int s1 = j % NSTOKES;

            HomogType d_sum = layer.solution(m).value.dual_homog_plus().deriv(
                derivindex,
                j + (layer.solution(m).value.nstr() * NSTOKES / 2) * a);
            if (m >= m_layers.surface().sk2_surface().max_azimuthal_order() ||
                s1 != 0)
                return d_sum; // rho will be zero
            const auto& rho = m_layers.surface().storage().brdf.stream_stream;
            const auto& d_rho = m_layers.surface()
                                    .storage()
                                    .d_brdf[deriv.surface_deriv_index]
                                    .stream_stream;
            for (StreamIndex i = 0; i < this->M_NSTR / 2; ++i) {
                HomogType d_minus =
                    layer.solution(m).value.dual_homog_minus().deriv(
                        derivindex,
                        i * NSTOKES +
                            (layer.solution(m).value.nstr() * NSTOKES / 2) * a);
                d_sum -= (1 + kronDelta(m, 0)) * rho(j / NSTOKES, i) *
                         (*this->M_WT)[i] * (*this->M_MU)[i] * d_minus;
                d_sum -= (1 + kronDelta(m, 0)) * deriv.d_albedo *
                         d_rho(j / NSTOKES, i) * (*this->M_WT)[i] *
                         (*this->M_MU)[i] *
                         layer.solution(m).value.homog_minus(i * NSTOKES, a);
            }
            return d_sum;
        }

        inline HomogType
        d_v_minus(AEOrder m, const OpticalLayer<NSTOKES, CNSTR>& layer,
                  StreamIndex j, SolutionIndex a, uint derivindex,
                  const LayerInputDerivative<NSTOKES>& deriv) const {
            // TODO: Polarized surface
            int s1 = j % NSTOKES;

            HomogType d_sum = layer.solution(m).value.dual_homog_minus().deriv(
                derivindex,
                j + (layer.solution(m).value.nstr() * NSTOKES / 2) * a);
            if ((m >= m_layers.surface().sk2_surface().max_azimuthal_order()) ||
                s1 != 0)
                return d_sum; // rho will be zero
            const auto& rho = m_layers.surface().storage().brdf.stream_stream;
            const auto& d_rho = m_layers.surface()
                                    .storage()
                                    .d_brdf[deriv.surface_deriv_index]
                                    .stream_stream;
            for (StreamIndex i = 0; i < this->M_NSTR / 2; ++i) {
                HomogType d_plus =
                    layer.solution(m).value.dual_homog_plus().deriv(
                        derivindex,
                        i * NSTOKES +
                            (layer.solution(m).value.nstr() * NSTOKES / 2) * a);
                d_sum -= (1 + kronDelta(m, 0)) * rho(j / NSTOKES, i) *
                         (*this->M_WT)[i] * (*this->M_MU)[i] * d_plus;
                d_sum -= (1 + kronDelta(m, 0)) * deriv.d_albedo *
                         d_rho(j / NSTOKES, i) * (*this->M_WT)[i] *
                         (*this->M_MU)[i] *
                         layer.solution(m).value.homog_plus(i * NSTOKES, a);
            }
            return d_sum;
        }

        inline double
        ground_direct_sun(AEOrder m, const OpticalLayer<NSTOKES, CNSTR>& layer,
                          StreamIndex out) const {
            // TODO: Polarized surface
            int s1 = out % NSTOKES;

            if (m_layers.surface().sk2_surface().max_azimuthal_order() <= m ||
                s1 != 0) {
                return 0;
            } else {
                return this->M_CSZ *
                       m_layers.surface().storage().brdf.stream_solar(out /
                                                                      NSTOKES) /
                       PI * layer.beamTransmittance(Location::FLOOR);
            }
        }

        inline double d_ground_direct_sun(
            AEOrder m, const OpticalLayer<NSTOKES, CNSTR>& layer,
            StreamIndex out, const LayerInputDerivative<NSTOKES>& deriv,
            uint derivindex) const {
            // TODO: Polarized surface
            int s1 = out % NSTOKES;

            if (m_layers.surface().sk2_surface().max_azimuthal_order() <= m ||
                s1 != 0) {
                return 0;
            } else {
                double result = this->M_CSZ *
                                m_layers.surface().storage().brdf.stream_solar(
                                    out / NSTOKES) /
                                PI *
                                layer.d_beamTransmittance(Location::FLOOR,
                                                          deriv, derivindex);

                double d_albedo =
                    deriv.d_albedo * m_layers.surface()
                                         .storage()
                                         .d_brdf[deriv.surface_deriv_index]
                                         .stream_solar(out / NSTOKES);

                result += this->M_CSZ * d_albedo / PI *
                          layer.beamTransmittance(Location::FLOOR);
                return result;
            }
        }

        inline double u_minus(AEOrder m,
                              const OpticalLayer<NSTOKES, CNSTR>& layer,
                              StreamIndex j) const {
            // TODO: Polarized surface
            int s1 = j % NSTOKES;

            const Eigen::VectorXd& part_soln_minus =
                layer.solution(m).value.dual_Gminus_bottom().value;
            double psi = part_soln_minus(j);

            if ((m < m_layers.surface().sk2_surface().max_azimuthal_order()) &&
                s1 == 0) {
                const auto& rho =
                    m_layers.surface().storage().brdf.stream_stream;
                const Eigen::VectorXd& part_soln_plus =
                    layer.solution(m).value.dual_Gplus_bottom().value;
                for (StreamIndex i = 0; i < this->M_NSTR / 2; ++i) {
                    psi -= (1 + kronDelta(m, 0)) * rho(j / NSTOKES, i) *
                           (*this->M_WT)[i] * (*this->M_MU)[i] *
                           part_soln_plus(i * NSTOKES);
                }
            }
            return psi;
        }

        inline double
        d_u_minus(AEOrder m, const OpticalLayer<NSTOKES, CNSTR>& layer,
                  StreamIndex j, uint derivindex,
                  const LayerInputDerivative<NSTOKES>& deriv) const {
            // TODO: Polarized surface
            int s1 = j % NSTOKES;

            // We don't have to calculate psi/d_psi separately since the
            // greens solution already includes the optical depth dependence
            double d_psi = layer.solution(m).value.dual_Gminus_bottom().deriv(
                derivindex, j);
            if ((m < m_layers.surface().sk2_surface().max_azimuthal_order()) &&
                s1 == 0) {
                const auto& rho =
                    m_layers.surface().storage().brdf.stream_stream;
                const auto& d_rho = m_layers.surface()
                                        .storage()
                                        .d_brdf[deriv.surface_deriv_index]
                                        .stream_stream;
                for (StreamIndex i = 0; i < this->M_NSTR / 2; ++i) {
                    d_psi -= (1 + kronDelta(m, 0)) * rho(j / NSTOKES, i) *
                             (*this->M_WT)[i] * (*this->M_MU)[i] *
                             layer.solution(m).value.dual_Gplus_bottom().deriv(
                                 derivindex, i * NSTOKES);
                    d_psi -= (1 + kronDelta(m, 0)) * deriv.d_albedo *
                             d_rho(j / NSTOKES, i) * (*this->M_WT)[i] *
                             (*this->M_MU)[i] *
                             layer.solution(m).value.dual_Gplus_bottom().value(
                                 i * NSTOKES);
                }
            }
            return d_psi;
        }

      private:
        // Reference to the layer array on which this object is operating.
        OpticalLayerArray<NSTOKES, CNSTR>& m_layers;
        // Flags indicating whether or not the order of the azimuth expansion
        // has been solved.
        std::vector<bool> m_is_solved;

        // Cached memory so we don't have to realloc for every layer/azimuth
        // direction
        RTEMemoryCache<NSTOKES, CNSTR>& m_cache;
    };

} // namespace sasktran_disco
