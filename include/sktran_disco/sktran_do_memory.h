#pragma once
#include "sktran_disco/sktran_do.h"
#include "sktran_disco/sktran_do_surface.h"

#ifdef SKTRAN_OPENMP_SUPPORT
#include "omp.h"
#endif

namespace sasktran_disco {
    /**
     * @brief temporary storage needed for layer postprocessing
     *
     * Contains all of the memory needed for the postprocessing of the radiative
     * transfer solution through the layers.
     *
     * @tparam NSTOKES
     * @tparam CNSTR
     */
    template <int NSTOKES, int CNSTR = -1> struct PostProcessingCache {
        VectorLayerDual<double> dual_lpsum_plus,
            dual_lpsum_minus; /** Duals of the legendre sums */
        VectorLayerDual<double> Y_plus,
            Y_minus; /** Duals of the homogeneous interpolated solutions */
        Radiance<NSTOKES> V, Q,
            J; /** To store the individual radiance sources */

        std::vector<LayerDual<double>> hp, hm; /** Homogenous multipliers */

        std::vector<Dual<double>> Dm, Dp,
            Eform; /** Greens function multipliers*/

        LayerIndex
            cached_layer; /** Index of the layer this cache corresponds to */

        PostProcessingCache() { cached_layer = -1; }

        void resize(uint NSTR, LayerIndex p, uint numlayerderiv,
                    uint layerstart, uint numtotalderiv) {
            if (p == cached_layer) {
                // Already cached
                return;
            }
            // Vector Layer Duals
            dual_lpsum_plus.resize(NSTR / 2 * NSTOKES * NSTOKES, numlayerderiv,
                                   p, layerstart);
            dual_lpsum_minus.resize(NSTR / 2 * NSTOKES * NSTOKES, numlayerderiv,
                                    p, layerstart);
            Y_plus.resize(NSTR / 2 * NSTOKES * NSTOKES, numlayerderiv, p,
                          layerstart);
            Y_minus.resize(NSTR / 2 * NSTOKES * NSTOKES, numlayerderiv, p,
                           layerstart);

            // Radiance containers
            V.resize(numtotalderiv, false);
            Q.resize(numtotalderiv, false);
            J.resize(numtotalderiv, false);

            hp.resize(NSTR / 2 * NSTOKES);
            hm.resize(NSTR / 2 * NSTOKES);
            Dm.resize(NSTR / 2 * NSTOKES);
            Dp.resize(NSTR / 2 * NSTOKES);
            Eform.resize(NSTR / 2 * NSTOKES);

            for (int i = 0; i < NSTR / 2 * NSTOKES; ++i) {
                // Layer Duals
                hp[i].resize(numlayerderiv);
                hm[i].resize(numlayerderiv);

                // Full Duals
                Dm[i].resize(numtotalderiv);
                Dp[i].resize(numtotalderiv);
                Eform[i].resize(numtotalderiv);
            }

            cached_layer = p;
        }
    };

    /**
     * @brief Temporaries that are used inside every OpticalLayer
     *
     * @tparam NSTOKES
     * @tparam CNSTR
     */
    template <int NSTOKES, int CNSTR = -1> struct LayerCache {
        LayerDual<double> dual_thickness; /** The thickness of the layer */
        LayerDual<double> dual_ssa;       /** Single scatter albedo */
        Dual<double> average_secant;      /** Average secant */

        LayerCache(uint NSTR) {}
    };

    /**
     * @brief Memory for the RTE solution
     *
     * @tparam NSTOKES
     * @tparam CNSTR
     */
    template <int NSTOKES, int CNSTR> struct RTEMemoryCache {
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

        Matrix h_eigmtx_destroy;
        Matrix h_MX_plus;
        Matrix h_MX_minus;
        Vector h_eigvalsq;
        Vector h_reigval_imag;
        Matrix h_identity;
        MatrixHLHS h_lhs;
        std::vector<MatrixHRHS> h_rhs;

        Eigen::PartialPivLU<MatrixHLHS> h_partiallu;
        Eigen::FullPivLU<MatrixHLHS> h_fullpivlu;

        std::vector<MatrixHRHS> h_d_X_d_k;

        std::vector<VectorLayerDual<double>> p_Qplus, p_Qminus;
        Dual<double> p_Cplus, p_Cminus;
        std::vector<LayerDual<double>> p_norm;

        sasktran_disco::TripleProductDerivativeHolder<NSTOKES> h_l_upwelling;
        sasktran_disco::TripleProductDerivativeHolder<NSTOKES> h_l_downwelling;

        // BVP things
        VectorDim1<BVPMatrixDenseBlock<NSTOKES>> d_mat;
        Eigen::MatrixXd d_b;

        std::unique_ptr<la::BVPMatrix<NSTOKES>> bvp_mat;
        Eigen::VectorXd bvp_b;
        Eigen::VectorXd bvp_temp;
        Eigen::VectorXi ipiv;

        // 2 Stream pentadiagonal cache, only allocated if in 2stream mode
        Eigen::VectorXd bvp_pd_alpha;
        Eigen::VectorXd bvp_pd_beta;
        Eigen::MatrixXd bvp_pd_d_z;
        Eigen::MatrixXd bvp_pd_z;
        Eigen::VectorXd bvp_pd_gamma;
        Eigen::VectorXd bvp_pd_mu;

        // Reverse Linearization Traces
        Eigen::MatrixXd m_Cplus_to_b;
        Eigen::MatrixXd m_Cminus_to_b;

        Eigen::MatrixXd m_trans_to_Cplus;
        Eigen::MatrixXd m_trans_to_Cminus;

        Eigen::MatrixXd m_secant_to_Cplus;
        Eigen::MatrixXd m_secant_to_Cminus;

        Eigen::MatrixXd m_bvp_backprop_z;

        Eigen::Matrix<double, NSTOKES, -1> m_secant_weights;
        Eigen::Matrix<double, NSTOKES, -1> m_trans_weights;

#ifdef SKTRAN_USE_ACCELERATE
        Eigen::VectorXd homog_work;
#endif

        bool has_been_configured_by_rte_solver;
    };

    /**
     * @brief Data that each thread will need, reused across wavelengths.
     *
     * @tparam NSTOKES
     * @tparam CNSTR
     */
    template <int NSTOKES, int CNSTR = -1> class ThreadData {
      private:
        mutable VectorDim2<LayerSolution<NSTOKES, CNSTR>>
            m_rte_solution; /** RTE Solutions*/
        mutable VectorDim1<PostProcessingCache<NSTOKES>>
            m_postprocessing_cache; /** Post processing caches*/
        mutable VectorDim1<LayerCache<NSTOKES>>
            m_layer_cache; /** Layer caches*/
        mutable InputDerivatives<NSTOKES>
            m_input_derivatives; /** Input derivatives*/
        mutable std::vector<Dual<double>> m_transmission;   /** Transmission
                                                                    cache*/
        mutable RTEMemoryCache<NSTOKES, CNSTR> m_rte_cache; /** RTE Cache*/

        mutable SurfaceStorage<NSTOKES, CNSTR>
            m_surface_storage; /** Storage for the expansion coefficients of the
                                  surface BRDF */

      public:
        /**
         * @brief Construct a new Thread Data object
         *
         * @param NLYR
         * @param NSTR
         */
        ThreadData(uint NLYR, uint NSTR) {
            m_rte_solution.resize(NLYR);
            for (auto& soln : m_rte_solution) {
                soln.resize(NSTR);
            }
            m_postprocessing_cache.resize(NLYR);
            m_layer_cache.resize(NLYR, NSTR);
            m_transmission.resize(NLYR + 1);

            m_rte_cache.has_been_configured_by_rte_solver = false;
        }

        /**
         * @brief Get the RTE solution for a given layer
         *
         * @param layerindex
         * @return std::vector<LayerSolution<NSTOKES, CNSTR>>&
         */
        std::vector<LayerSolution<NSTOKES, CNSTR>>&
        rte_solution(uint layerindex) const {
            return m_rte_solution[layerindex];
        }

        /**
         * @brief Get the post processing cache for a given layer
         *
         * @param layerindex
         * @return PostProcessingCache<NSTOKES>&
         */
        PostProcessingCache<NSTOKES>&
        postprocessing_cache(uint layerindex) const {
            return m_postprocessing_cache[layerindex];
        }

        /**
         * @brief Get the layer cache
         *
         * @param layerindex
         * @return LayerCache<NSTOKES>&
         */
        LayerCache<NSTOKES>& layer_cache(uint layerindex) const {
            return m_layer_cache[layerindex];
        }

        /**
         * @brief Get the input derivatives
         *
         * @return InputDerivatives<NSTOKES>&
         */
        InputDerivatives<NSTOKES>& input_derivatives() const {
            return m_input_derivatives;
        }

        /**
         * @brief Layer transmissions
         *
         * @return std::vector<Dual<double>>&
         */
        std::vector<Dual<double>>& transmission() const {
            return m_transmission;
        }

        /**
         * @brief Get the RTE cache
         *
         * @return RTEMemoryCache<NSTOKES, CNSTR>&
         */
        RTEMemoryCache<NSTOKES, CNSTR>& rte_cache() const {
            return m_rte_cache;
        }

        /**
         * @brief Get the storage for the surface expansion coefficients
         *
         * @return SurfaceStorage<NSTOKES, CNSTR>&
         */
        SurfaceStorage<NSTOKES, CNSTR>& surface_storage() const {
            return m_surface_storage;
        }
    };

    /**
     * Pool of memory for a single Engine instance.  One ThreadData object is
     * instantiated for each thread that is then intended to be reused across
     * wavelengths.  Currently the engine does not know exactly how many threads
     * will be used at the time of calling init, and so we instantiate the
     * maximum number of available threads, but this should be fixed
     *
     * @tparam NSTOKES
     * @tparam CNSTR
     */
    template <int NSTOKES, int CNSTR = -1> class MemoryPool {
      private:
        // Use a map because we emplace a new object and we don't want reallocs
        // triggering a move
        mutable std::map<int, ThreadData<NSTOKES, CNSTR>> m_threaddata;
        uint m_nlyr;
        uint m_nstr;

      public:
        MemoryPool() {}

        MemoryPool(uint NLYR, uint NSTR) { init(NLYR, NSTR); }

        void init(uint NLYR, uint NSTR) {
#ifdef SKTRAN_OPENMP_SUPPORT
            int num_threads = omp_get_max_threads();
#else
            int num_threads = 1;
#endif

            m_nlyr = NLYR;
            m_nstr = NSTR;

            for (int k = 0; k < num_threads; ++k) {
                m_threaddata.emplace(k, ThreadData<NSTOKES, CNSTR>(NLYR, NSTR));
            }
        }

        ThreadData<NSTOKES, CNSTR>& thread_data() const {
#ifdef SKTRAN_OPENMP_SUPPORT
            int thread_num = omp_get_thread_num();
#else
            int thread_num = 0;
#endif

            auto it = m_threaddata.find(thread_num);
            if (it == m_threaddata.end()) {
                m_threaddata.emplace(
                    thread_num, ThreadData<NSTOKES, CNSTR>(m_nlyr, m_nstr));
            }

            return m_threaddata.at(thread_num);
        }

        ThreadData<NSTOKES, CNSTR>& thread_data(int thread_num) const {
            auto it = m_threaddata.find(thread_num);
            if (it == m_threaddata.end()) {
                m_threaddata.emplace(
                    thread_num, ThreadData<NSTOKES, CNSTR>(m_nlyr, m_nstr));
            }

            return m_threaddata.at(thread_num);
        }
    };

} // namespace sasktran_disco
