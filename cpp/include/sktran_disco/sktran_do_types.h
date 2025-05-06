#pragma once
#include "sktran_disco/sktran_do.h"
#include "sktran_disco/sktran_do_linearization_types.h"

#ifdef _WIN32
const double M_PI = 3.14159265358979323846264338327950288419716939937510582;
#endif

// Absolute difference between the solar secant and an eigenvalue before a
// taylor expansion is used in the solution
#define SKTRAN_DO_GREENS_EPS 1e-4

namespace sasktran_disco {
#pragma region "Index aliases"

    // General purpose unsigned integer
    typedef unsigned int uint;

    // Specifies a specific order of the azimuth expansion
    typedef unsigned int AEOrder;

    // Specifies a polynomial order
    typedef unsigned int LPOrder;

    // Specifies a solution index
    typedef unsigned int SolutionIndex;

    // Specifies a stream index
    typedef unsigned int StreamIndex;

    // Specifies an atmospheric layer index. TOA layer index is 0, ground
    // layer index is NLYR-1.
    typedef unsigned int LayerIndex;

    // Specifies the boundary between two atmospheric layer. A boundary
    // index, p, has a top layer with LayerIndex=p-1 and bottom layer with
    // LayerIndex=p.
    typedef unsigned int BoundaryIndex;
#pragma endregion

#pragma region "Type aliases"
    // Alias for 1 dimensional std::vector
    template <class StoredType> using VectorDim1 = std::vector<StoredType>;

    // Alias for 2 dimensional std::vector
    template <class StoredType>
    using VectorDim2 = std::vector<VectorDim1<StoredType>>;

    // Alias for 3 dimensional std::vector
    template <class StoredType>
    using VectorDim3 = std::vector<VectorDim2<StoredType>>;
#pragma endregion

#pragma region "Common enumerations"
    // Specifies a location in an optical layer
    enum struct Location { CEILING, INSIDE, FLOOR };

    // Specifies a propagating direction
    enum class Propagating { UP, DOWN };

#pragma endregion

#pragma region "Linearization"
    /**
     * Derivative of the layer quantities with respect to a parameter.  These
     * form the input parameters to the model to calculate the derivative.  The
     * core LIDORT algorithm needs to know the changes in SSA, Optical Depth,
     * and Legendre Coefficients.
     *
     * @tparam NSTOKES
     * @tparam CNSTR
     */
    template <int NSTOKES, int CNSTR> class LayerInputDerivative {
      public:
        /**
         * @brief Construct a new Layer Input Derivative object. Need to know
         * the layer and the number of streams
         *
         * @param nstr
         * @param p
         */
        LayerInputDerivative(uint nstr, LayerIndex p) {
            d_legendre_coeff.resize(nstr);
            layer_index = p;
            surface_deriv_index = 0;

            setZero();
        }

        std::vector<LegendreCoefficient<NSTOKES>>
            d_legendre_coeff; /** Derivatives of the legendre coefficients with
                                 respect to the parameter */
        double d_optical_depth; /** Derivative of optical depth with respect to
                                   the parameter */
        double d_SSA;    /** Derivative of single scatter albedo with respect to
                            the parameter */
        double d_albedo; /** Derivative of the surface albedo with respect to
                            the parameter */
        int surface_deriv_index;

        LayerIndex layer_index; /** Layer index that is varying */
        std::vector<std::pair<uint, double>>
            group_and_triangle_fraction; /** Mapping between engine weightig
                                            functions and internal weighting
                                            functions */
        std::vector<double>
            extinctions; /** Multipliers for each group and triangle fraction */

      private:
        /**
         * @brief Sets all values to 0
         *
         */
        void setZero();
    };

    /**
     * @brief Class which stores all of the user input derivatives
     *
     * @tparam NSTOKES
     * @tparam CNSTR
     */
    template <int NSTOKES, int CNSTR> class InputDerivatives {
      public:
        /**
         * @brief Construct a new Input Derivatives object
         *
         */
        InputDerivatives() { m_geometry_configured = false; };

        /**
         * @brief The full storage container of all derivatives
         *
         * @return const VectorDim1<LayerInputDerivative<NSTOKES>>&
         */
        const VectorDim1<LayerInputDerivative<NSTOKES>>&
        layerDerivatives() const {
            return m_layerderivs;
        }

        /**
         * @brief The full storage container of all derivatives
         *
         * @return VectorDim1<LayerInputDerivative<NSTOKES>>&
         */
        VectorDim1<LayerInputDerivative<NSTOKES>>& layerDerivatives() {
            return m_layerderivs;
        }

        /**
         * @brief Get the Layer Derivative for an index
         *
         * @param i
         * @return const LayerInputDerivative<NSTOKES>&
         */
        const LayerInputDerivative<NSTOKES>& operator[](int i) const {
            return m_layerderivs[i];
        }

        /**
         * @brief Adds a derivative to the container
         *
         * @param deriv
         */
        void addDerivative(LayerInputDerivative<NSTOKES>&& deriv) {
            m_layerderivs.emplace_back(deriv);
        }

        /**
         * @brief Adds an empty derivative to the container and returns back a
         * reference to it
         *
         * @param nstr
         * @param p
         * @return LayerInputDerivative<NSTOKES>&
         */
        inline LayerInputDerivative<NSTOKES>& addDerivative(uint nstr,
                                                            LayerIndex p) {
            m_layerderivs.emplace_back(nstr, p);

            return m_layerderivs.back();
        }

        /**
         * @brief Total number of derivatives in the container
         *
         * @return size_t
         */
        inline size_t numDerivative() const { return m_layerderivs.size(); }

        /**
         * @brief Total number of derivatives for a given layer
         *
         * @param p
         * @return size_t
         */
        inline size_t numDerivativeLayer(LayerIndex p) const {
            if (m_layerderivs.size() == 0) {
                // We are not calculating derivatives at all
                return 0;
            }
            return m_numderivlayer[p];
        }

        /**
         * @brief The index that the layer derivatives for index p start at
         *
         * @param p
         * @return size_t
         */
        inline size_t layerStartIndex(LayerIndex p) const {
            if (m_layerderivs.size() == 0) {
                // Not calculating derivatives
                return 0;
            }
            return m_layerstartindex[p];
        }

        /**
         * @brief Sorts all the derivatives so they are ordered by layer
         *
         * @param nlyr
         */
        void sort(uint nlyr) {
            std::sort(std::begin(m_layerderivs), std::end(m_layerderivs),
                      [](const LayerInputDerivative<NSTOKES>& a,
                         const LayerInputDerivative<NSTOKES>& b) -> bool {
                          return a.layer_index < b.layer_index;
                      });
            // Construct the layer start indicies
            m_layerstartindex.resize(nlyr);
            m_numderivlayer.resize(nlyr, 0);
            for (uint i = 0; i < m_layerderivs.size(); ++i) {
                m_numderivlayer[m_layerderivs[i].layer_index]++;
            }

            m_layerstartindex[0] = 0;
            for (uint i = 0; i < m_numderivlayer.size() - 1; i++) {
                m_layerstartindex[i + 1] =
                    m_layerstartindex[i] + m_numderivlayer[i];
            }
        }

        /**
         *
         * @return true If the geometry is configured
         * @return false If the geometry has not been configured
         */
        bool is_geometry_configured() const { return m_geometry_configured; }

        /**
         * @brief Sets the geometry as configured
         *
         */
        void set_geometry_configured() { m_geometry_configured = true; }

        /**
         * @brief Reverse linearization trace for a given los index
         *
         * @param i
         * @return ReverseLinearizationTrace<NSTOKES>&
         */
        ReverseLinearizationTrace<NSTOKES>& reverse_trace(size_t i) {
            return m_reverse_linearization_traces[i];
        }

        /**
         * @brief Sets the number of lines of sight to initialize the reverse
         * traces
         *
         * @param numlos
         * @param nstr
         * @param nlyr
         */
        void set_num_los(size_t numlos, uint nstr, uint nlyr) {
            m_reverse_linearization_traces.resize(numlos);

            for (auto& trace : m_reverse_linearization_traces) {
                trace.resize(nstr, nlyr);
            }
        }

        /**
         * @brief Sets all of the reverse linearization traces to be 0
         *
         */
        void set_zero_traces() {
            for (auto& trace : m_reverse_linearization_traces) {
                trace.set_zero();
            }
        }

      private:
        VectorDim1<LayerInputDerivative<NSTOKES>>
            m_layerderivs; /** Storage of the derivatives */
        std::vector<size_t>
            m_layerstartindex; /** Start indexes for each layer */
        std::vector<size_t>
            m_numderivlayer; /** Number of derivatives in each layer */

        std::vector<ReverseLinearizationTrace<NSTOKES>>
            m_reverse_linearization_traces; /** Revers linearization traces */

        bool m_geometry_configured; /** True if the geometry has been configured
                                     */
    };

    /**
     * @brief Class to store forward derivatives in the DO model
     *
     * A dual is a combination of a value, and the derivatives of value with
     * respect to ALL quantities
     *
     * @tparam T
     */
    template <typename T> struct Dual {
        T value;
        Eigen::VectorX<T> deriv;

        Dual(size_t numderiv) { resize(numderiv); }
        Dual() { value = 0.0; };

        void resize(size_t numderiv, bool setzero = true) {
            deriv.resize(numderiv);
            if (setzero) {
                deriv.setZero();
                value = 0.0;
            }
        }
    };

    /**
     * A layer dual is a Dual where the derivative is 0 except for inside a
     * single layer.  We define a separate class to represent this since
     * operations involving LayerDuals are significantly faster
     *
     * @tparam T
     */
    template <typename T> struct LayerDual {
        T value;
        uint layer_start;
        LayerIndex layer_index;
        Eigen::VectorX<T> deriv;

        LayerDual(size_t numderiv, LayerIndex l, uint layerstartindex) {
            resize(numderiv);
            layer_index = l;
            layer_start = layerstartindex;
            value = 0.0;
        }

        LayerDual() { value = 0.0; };

        LayerDual<T>& operator=(const LayerDual<T>& other) {
            this->value = other.value;
            this->deriv.noalias() = other.deriv;

            return *this;
        }

        void resize(size_t numderiv) {
            deriv.resize(numderiv);
            deriv.setZero();
        }
    };

#pragma endregion

#pragma region "Math helpers"
    // Kronecker delta function.
    // Returns : 1.0 if i == j, else return 0.0
    inline double kronDelta(uint i, uint j) { return (i == j) ? 1.0 : 0.0; }

    static const double PI = M_PI;

    /**
     * @brief Efficiently calculates the triple product of vectors of Legendre
     * polynomials
     *
     * @tparam NSTOKES
     * @tparam CNSTR
     */
    template <int NSTOKES, int CNSTR = -1> class LPTripleProduct {
      public:
        // Calculates an auxilary result for the given Legendre polynomial
        // vectors.
        LPTripleProduct(
            AEOrder m, const std::vector<LegendreCoefficient<NSTOKES>>& lephase,
            const std::vector<LegendrePhaseContainer<NSTOKES>>& lp1,
            const std::vector<LegendrePhaseContainer<NSTOKES>>& lp2)
            : m_association_order(m),
              m_aux(std::piecewise_construct,
                    std::forward_as_tuple((uint)lephase.size()),
                    std::forward_as_tuple((uint)lephase.size())),
              m_nstr((uint)lephase.size()) {
            calculate(lephase, lp1, lp2);
        }

        LPTripleProduct(AEOrder m, uint size)
            : m_association_order(m),
              m_aux(std::piecewise_construct, std::forward_as_tuple(size),
                    std::forward_as_tuple(size)),
              m_nstr(size) {}

        LPTripleProduct(uint size)
            : m_aux(std::piecewise_construct, std::forward_as_tuple(size),
                    std::forward_as_tuple(size)),
              m_nstr(size) {}

        void
        calculate(AEOrder m,
                  const std::vector<LegendreCoefficient<NSTOKES>>& lephase,
                  const std::vector<LegendrePhaseContainer<NSTOKES>>& lp1,
                  const std::vector<LegendrePhaseContainer<NSTOKES>>& lp2) {
            m_association_order = m;
            calculate(lephase, lp1, lp2);
        }

        void calculate(const std::vector<LegendreCoefficient<NSTOKES>>& lephase,
                       const std::vector<LegendrePhaseContainer<NSTOKES>>& lp1,
                       const std::vector<LegendrePhaseContainer<NSTOKES>>& lp2);

        void calculate_and_emplace(
            AEOrder m, const std::vector<LegendreCoefficient<NSTOKES>>& lephase,
            const std::vector<LegendrePhaseContainer<NSTOKES>>& lp1,
            const std::vector<LegendrePhaseContainer<NSTOKES>>& lp2,
            sasktran_disco::TripleProductDerivativeHolder<NSTOKES>&
                holder_0_negation,
            sasktran_disco::TripleProductDerivativeHolder<NSTOKES>&
                holder_1_negation,
            double ssa);

        void negations_derivative_emplace(
            uint num,
            sasktran_disco::TripleProductDerivativeHolder<NSTOKES>& holder);

      private:
        // A pair of the result of the triple product of even and odd order
        // polynomials. Triple product of even polynomials are stored in first,
        // the triple product of odd polynomials are stored in second.
        using AuxilaryResult =
            std::pair<sasktran_disco::TripleProductDerivativeHolder<NSTOKES>,
                      sasktran_disco::TripleProductDerivativeHolder<NSTOKES>>;

        AuxilaryResult m_aux;
        AEOrder m_association_order;
        uint m_nstr;
    };

#pragma endregion

#pragma region "Useful DO classes"
    /**
     * @brief Storage for DO postprocessing lines of sight
     *
     */
    struct LineOfSight {
        double coszenith;
        double azimuth;
        double cos_scattering_angle;
        uint unsorted_index;
        double observeraltitude;
    };

    /**
     * @brief The RTE solution for a given layer and azimuth order
     *
     * @tparam NSTOKES
     * @tparam CNSTR
     */
    template <int NSTOKES, int CNSTR = -1> class RTEGeneralSolution {

        // If NSTOKES == 1 all of the eigenvalues are real, but for NSTOKES > 1
        // we can have complex eigenvalues and complex homogeneous solution
        // vectors
        using HomogType = typename std::conditional<NSTOKES != 5, double,
                                                    std::complex<double>>::type;

        using HomogVector = typename std::conditional<
            CNSTR == -1, Eigen::VectorX<HomogType>,
            Eigen::Vector<HomogType, CNSTR / 2 * NSTOKES>>::type;
        using HomogMatrix = typename std::conditional<
            CNSTR == -1, Eigen::VectorX<HomogType>,
            Eigen::Vector<HomogType,
                          CNSTR / 2 * NSTOKES * CNSTR / 2 * NSTOKES>>::type;

        using Vector = typename std::conditional<
            CNSTR == -1, Eigen::VectorX<double>,
            Eigen::Vector<double, CNSTR / 2 * NSTOKES>>::type;

        using HomogVectorLayerDual = typename std::conditional<
            CNSTR == -1, VectorLayerDual<HomogType>,
            VectorLayerDual<HomogType,
                            CNSTR / 2 * NSTOKES * CNSTR / 2 * NSTOKES>>::type;
        using EigVectorLayerDual = typename std::conditional<
            CNSTR == -1, VectorLayerDual<HomogType>,
            VectorLayerDual<HomogType, CNSTR / 2 * NSTOKES>>::type;

        using VectorLayerDualType = typename std::conditional<
            CNSTR == -1, VectorLayerDual<double>,
            VectorLayerDual<double, CNSTR / 2 * NSTOKES>>::type;
        using VectorDualType = typename std::conditional<
            CNSTR == -1, VectorDual<double>,
            VectorDual<double, CNSTR / 2 * NSTOKES>>::type;

      public:
        RTEGeneralSolution() : M_NSTR(0) {}

        /**
         * @brief Resizes the solution storage
         *
         * @param nstr
         * @param nderivlayer
         * @param l
         * @param layer_start
         * @param nderivtotal
         */
        void resize(size_t nstr, size_t nderivlayer, LayerIndex l,
                    uint layer_start, uint nderivtotal) {
            const_cast<uint&>(M_NSTR) = (uint)nstr;
            m_eigval.resize(nstr / 2 * NSTOKES, nderivlayer, l, layer_start);

            m_homog_plus.resize((nstr / 2) * (nstr / 2) * NSTOKES * NSTOKES,
                                nderivlayer, l, layer_start);
            m_homog_minus.resize((nstr / 2) * (nstr / 2) * NSTOKES * NSTOKES,
                                 nderivlayer, l, layer_start);

            m_green_A_minus.resize(nstr / 2 * NSTOKES, nderivlayer, l,
                                   layer_start);
            m_green_A_plus.resize(nstr / 2 * NSTOKES, nderivlayer, l,
                                  layer_start);

            m_Gplus_top.resize(nstr / 2 * NSTOKES, nderivtotal);
            m_Gplus_bottom.resize(nstr / 2 * NSTOKES, nderivtotal);
            m_Gminus_top.resize(nstr / 2 * NSTOKES, nderivtotal);
            m_Gminus_bottom.resize(nstr / 2 * NSTOKES, nderivtotal);
        }

        inline HomogType eigval(SolutionIndex a) const {
            return m_eigval.value(a);
        }

        inline HomogType homog_plus(StreamIndex j, SolutionIndex a) const {
            return m_homog_plus.value(j + (M_NSTR * NSTOKES / 2) * a);
        }

        inline HomogType homog_minus(StreamIndex j, SolutionIndex a) const {
            return m_homog_minus.value(j + (M_NSTR * NSTOKES / 2) * a);
        }

        inline HomogType d_homog_minus(StreamIndex j, SolutionIndex a,
                                       uint deriv) const {
            return m_homog_minus.deriv(deriv, j + (M_NSTR * NSTOKES / 2) * a);
        }

        inline HomogType d_homog_plus(StreamIndex j, SolutionIndex a,
                                      uint deriv) const {
            return m_homog_plus.deriv(deriv, j + (M_NSTR * NSTOKES / 2) * a);
        }

        inline HomogVector& eigval() { return m_eigval.value; }
        inline const HomogVector& eigval() const { return m_eigval.value; }

        inline EigVectorLayerDual& dual_eigval() { return m_eigval; }
        inline const EigVectorLayerDual& dual_eigval() const {
            return m_eigval;
        }

        inline HomogMatrix& homog_plus() { return m_homog_plus.value; }
        inline const HomogMatrix& homog_plus() const {
            return m_homog_plus.value;
        }

        inline HomogVectorLayerDual& dual_homog_plus() { return m_homog_plus; }
        inline const HomogVectorLayerDual& dual_homog_plus() const {
            return m_homog_plus;
        }

        inline HomogMatrix& homog_minus() { return m_homog_minus.value; }
        inline const HomogMatrix& homog_minus() const {
            return m_homog_minus.value;
        }

        inline HomogVectorLayerDual& dual_homog_minus() {
            return m_homog_minus;
        }
        inline const HomogVectorLayerDual& dual_homog_minus() const {
            return m_homog_minus;
        }

        inline VectorLayerDualType& dual_green_A_minus() {
            return m_green_A_minus;
        }

        inline const VectorLayerDualType& dual_green_A_minus() const {
            return m_green_A_minus;
        }

        inline VectorLayerDualType& dual_green_A_plus() {
            return m_green_A_plus;
        }

        inline const VectorLayerDualType& dual_green_A_plus() const {
            return m_green_A_plus;
        }

        inline VectorDualType& dual_Gplus_top() { return m_Gplus_top; }

        inline const VectorDualType& dual_Gplus_top() const {
            return m_Gplus_top;
        }

        inline VectorDualType& dual_Gplus_bottom() { return m_Gplus_bottom; }

        inline const VectorDualType& dual_Gplus_bottom() const {
            return m_Gplus_bottom;
        }

        inline VectorDualType& dual_Gminus_top() { return m_Gminus_top; }

        inline const VectorDualType& dual_Gminus_top() const {
            return m_Gminus_top;
        }

        inline VectorDualType& dual_Gminus_bottom() { return m_Gminus_bottom; }

        inline const VectorDualType& dual_Gminus_bottom() const {
            return m_Gminus_bottom;
        }

        inline const uint nstr() const { return M_NSTR; }

      private:
        const uint M_NSTR;

        EigVectorLayerDual m_eigval; // k,  See eq (16)

        HomogVectorLayerDual m_homog_plus;  // W+, See eq (17)
        HomogVectorLayerDual m_homog_minus; // W-, See eq (17)

        VectorLayerDualType m_green_A_plus;
        VectorLayerDualType m_green_A_minus;

        VectorDualType m_Gplus_top;
        VectorDualType m_Gplus_bottom;
        VectorDualType m_Gminus_top;
        VectorDualType m_Gminus_bottom;
    };

    /**
     * Stores some cached quantities that are needed between different steps of
     * the RTE solution
     *
     * @tparam NSTOKES
     * @tparam CNSTR
     */
    template <int NSTOKES, int CNSTR = -1> class RTESolutionCache {
        using Matrix =
            typename std::conditional<CNSTR == -1, Eigen::MatrixXd,
                                      Eigen::Matrix<double, CNSTR / 2 * NSTOKES,
                                                    CNSTR / 2 * NSTOKES>>::type;

      public:
        RTESolutionCache() : M_NSTR(0) {
            // empty
        }
        RTESolutionCache(uint nstr)
            : M_NSTR(nstr), m_s_plus(nstr / 2 * NSTOKES, nstr / 2 * NSTOKES),
              m_s_minus(nstr / 2 * NSTOKES, nstr / 2 * NSTOKES),
              m_eigmtx(nstr / 2 * NSTOKES, nstr / 2 * NSTOKES) {
            // empty
        }
        void resize(uint nstr) {
            const_cast<uint&>(M_NSTR) = nstr;
            m_s_plus.resize(nstr / 2 * NSTOKES, nstr / 2 * NSTOKES);
            m_s_minus.resize(nstr / 2 * NSTOKES, nstr / 2 * NSTOKES);
            m_eigmtx.resize(nstr / 2 * NSTOKES, nstr / 2 * NSTOKES);
        }

        inline Matrix& s_plus() { return m_s_plus; }
        inline const Matrix& s_plus() const { return m_s_plus; }

        inline Matrix& s_minus() { return m_s_minus; }
        inline const Matrix& s_minus() const { return m_s_minus; }

        inline Matrix& eigmtx() { return m_eigmtx; }
        inline const Matrix& eigmtx() const { return m_eigmtx; }

      private:
        const uint M_NSTR;

        Matrix m_s_plus;
        Matrix m_s_minus;
        Matrix m_eigmtx;
    };

    /**
     * Stores the homogeneous coefficients, L, M, for a single solution
     *
     */
    struct RTEBoundarySolution {
        VectorDual<double> L_coeffs;
        VectorDual<double> M_coeffs;
    };

    /**
     * @brief The full solution for a single layer
     *
     * @tparam NSTOKES
     * @tparam CNSTR
     */
    template <int NSTOKES, int CNSTR = -1> struct LayerSolution {
        using HomogType = typename std::conditional<NSTOKES != 5, double,
                                                    std::complex<double>>::type;

        /**
         * @brief Sets up the solution memory
         *
         * @param nstr
         * @param l
         * @param input_deriv
         */
        void configure(size_t nstr, LayerIndex l,
                       const InputDerivatives<NSTOKES>& input_deriv) {
            // Core solution contains no cross derivatives
            value.resize(nstr, input_deriv.numDerivativeLayer(l), l,
                         (uint)input_deriv.layerStartIndex(l),
                         (uint)input_deriv.numDerivative());

            // Cache uses a different format and is a little messy
            cache.resize((uint)nstr);

            boundary.L_coeffs.resize(nstr / 2 * NSTOKES,
                                     input_deriv.numDerivative());
            boundary.M_coeffs.resize(nstr / 2 * NSTOKES,
                                     input_deriv.numDerivative());

            layer_index = l;

            configureDerivatives(nstr, layer_index, input_deriv);
        }

        /**
         * @brief Sets up any memory needed for derivative propagation of the
         * cache
         *
         * @param nstr
         * @param layer_index
         * @param input_deriv
         */
        void
        configureDerivatives(size_t nstr, LayerIndex layer_index,
                             const InputDerivatives<NSTOKES>& input_deriv) {
            if (input_deriv.numDerivative() == 0) {
                // We have no derivatives, just return
                return;
            }

            // Terms that don't contain cross derivatives
            auto numderiv = input_deriv.numDerivativeLayer(layer_index);
            d_cache.resize(numderiv);
            for (uint i = 0; i < numderiv; ++i) {
                d_cache[i].resize((uint)nstr);
            }
        }

        LayerIndex layer_index;
        RTEGeneralSolution<NSTOKES, CNSTR> value;
        RTESolutionCache<NSTOKES, CNSTR> cache;
        std::vector<RTESolutionCache<NSTOKES, CNSTR>> d_cache;
        RTEBoundarySolution boundary;
    };

#pragma endregion
    class InvalidConfiguration : public std::exception {
      public:
        explicit InvalidConfiguration(const char* message) : m_msg(message) {}
        explicit InvalidConfiguration(const std::string& message)
            : m_msg(message) {}
        virtual ~InvalidConfiguration() throw() {}
        virtual const char* what() const throw() { return m_msg.c_str(); }

      private:
        std::string m_msg;
    };

    class InternalError : public std::exception {
      public:
        explicit InternalError(const char* message = "<none>") {
            m_msg = "An unexpected internal exception was thrown. This is "
                    "likely a bug! "
                    "Please submit a issue at: "
                    "https://arggit.usask.ca/ARGPackages/SasktranDO. "
                    "The following error message was given: " +
                    std::string(message);
        }
        explicit InternalError(const std::string& message = "<none>") {
            m_msg = "An unexpected internal exception was thrown. This is "
                    "likely a bug! "
                    "Please submit a issue at: "
                    "https://arggit.usask.ca/ARGPackages/SasktranDO. "
                    "The following error message was given: " +
                    std::string(message);
        }
        virtual ~InternalError() throw() {}
        virtual const char* what() const throw() { return m_msg.c_str(); }

      protected:
        std::string m_msg;
    };

    class InternalRuntimeError : public InternalError {
      public:
        explicit InternalRuntimeError(const char* message)
            : InternalError(message) {
            m_msg = "An internal runtime exception has occured. Likey due to "
                    "insufficient precision. ERROR MESSAGE:" +
                    std::string(message);
        }
        explicit InternalRuntimeError(const std::string& message)
            : InternalError(message) {
            m_msg = "An internal runtime exception has occured. Likey due to "
                    "insufficient precision. ERROR MESSAGE:" +
                    std::string(message);
        }
        virtual ~InternalRuntimeError() throw() {}
        virtual const char* what() const throw() { return m_msg.c_str(); }
    };
} // namespace sasktran_disco
