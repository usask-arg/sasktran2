#pragma once
// #include "sktran_disco/sktran_do.h"

namespace sasktran_disco {
    /**
     * Base storage for the value/derivative pairs, used for terms that contain
     * cross derivatives.  Unfortunately some terms we want to store are
     * complex.
     *
     * @tparam T
     * @tparam CSIZE
     */
    template <typename T, int CSIZE = -1> struct VectorDual {
        Eigen::Vector<T, CSIZE> value;
        Eigen::Matrix<T, -1, CSIZE> deriv;

        /**
         * @brief Construct a new Vector Dual object
         *
         * @param numvalue
         * @param numderiv
         */
        VectorDual(size_t numvalue, size_t numderiv) {
            resize(numvalue, numderiv);
        }

        VectorDual(){};

        /**
         * @brief Resize the object
         *
         * @param numvalue
         * @param numderiv
         */
        void resize(size_t numvalue, size_t numderiv) {
            if (numvalue != value.size()) {
                value.resize(numvalue);
            }
            if (numvalue != deriv.cols() || numderiv != deriv.rows()) {
                deriv.resize(numderiv, numvalue);
            }
        }
    };

    /**
     * Base storage for value derivative pairs that do not contain cross
     * derivatives. Essentially the same storage but also contains the layer
     * indicies.
     *
     * @tparam T
     * @tparam CSIZE
     */
    template <typename T, int CSIZE = -1>
    struct VectorLayerDual : public VectorDual<T, CSIZE> {
        unsigned int layer_index;
        unsigned int layer_start;

        /**
         * @brief Construct a new Vector Layer Dual object
         *
         * @param numvalue
         * @param numderiv
         * @param l
         * @param lstart
         */
        VectorLayerDual(size_t numvalue, size_t numderiv, unsigned int l,
                        unsigned int lstart)
            : VectorDual<T>(numvalue, numderiv), layer_index(l),
              layer_start(lstart) {}

        VectorLayerDual(const VectorLayerDual& other) {
            resize(other.value.size(), other.deriv.rows(), other.layer_index,
                   other.layer_start);
        }

        VectorLayerDual(){};

        /**
         * @brief Resize the object
         *
         * @param numvalue
         * @param numderiv
         * @param l
         * @param lstart
         */
        void resize(size_t numvalue, size_t numderiv, unsigned int l,
                    unsigned int lstart) {
            VectorDual<T, CSIZE>::resize(numvalue, numderiv);

            layer_index = l;
            layer_start = lstart;
        }
    };

    /**
     *  Trace that is passed through the calculation in order to store extra
     * parameters that are necessary for the back-propagation of the
     * derivatives.
     *
     * @tparam NSTOKES
     */
    template <int NSTOKES> class ReverseLinearizationTrace {
      private:
        // The weights for the BVP coefficients, these are the values so that
        // the contributions from the BVP coeffs to the radiance are (weights,
        // bvp_coeffs)
        Eigen::Matrix<double, -1, NSTOKES> m_bvp_coeff_weights;

      public:
        ReverseLinearizationTrace() {}
        ~ReverseLinearizationTrace() {}

        /**
         * @brief Weights for the BVP coefficients
         *
         * @return Eigen::Matrix<double, -1, NSTOKES>&
         */
        Eigen::Matrix<double, -1, NSTOKES>& bvp_coeff_weights() {
            return m_bvp_coeff_weights;
        }

        /**
         * @brief Multiplies all the weights by a constant
         *
         * @param x
         */
        void multiply_by_constant(double x) { m_bvp_coeff_weights *= x; }

        /**
         * @brief Resizes the object
         *
         * @param nstr
         * @param nlyr
         */
        void resize(int nstr, int nlyr) {
            int N = nstr * nlyr * NSTOKES;
            m_bvp_coeff_weights.resize(N, NSTOKES);
        }

        /**
         * @brief Set the object to zero
         *
         */
        void set_zero() { m_bvp_coeff_weights.setZero(); }
    };

} // namespace sasktran_disco
