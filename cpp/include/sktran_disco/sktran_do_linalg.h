#pragma once
#include "sktran_disco/sktran_do.h"

namespace sasktran_disco {
    // Linear algebra namespace.
    namespace la {
        /**
         * @brief Class used to store the boundary value problem matrix that is
         * band diagonal
         *
         * @tparam NSTOKES
         * @tparam CNSTR
         */
        template <int NSTOKES, int CNSTR = -1> class BVPMatrix {
          public:
            /**
             * @brief Construct a new BVPMatrix object
             *
             * @param NSTR Number of streams
             * @param NLYR Number of layers
             */
            BVPMatrix(uint NSTR, uint NLYR)
                : M_NCD(3 * NSTR * NSTOKES / 2 - 1),
                  M_NCOLS(NSTR * NSTOKES * NLYR),
                  M_NROWS(3 * (3 * NSTR * NSTOKES / 2 - 1) + 1),
                  M_NRM1(3 * (3 * NSTR * NSTOKES / 2 - 1)), M_NSTR(NSTR),
                  M_NLYR(NLYR) {
                m_data =
                    std::unique_ptr<double[]>(new double[M_NROWS * M_NCOLS]);
            }

            ~BVPMatrix(){};

          public: // Interface
            inline double& operator()(uint i, uint j) {
                return m_data[2 * M_NCD + i + j * M_NRM1];
            }
            inline double operator()(uint i, uint j) const {
                return m_data[2 * M_NCD + i + j * M_NRM1];
            }

            /**
             * @brief Set boundary matrix to zero
             *
             */
            inline void setZero() const {
                std::fill_n(m_data.get(), M_NCOLS * M_NROWS, 0.0);
            }

            /**
             * @brief Class which represents a block in the boundary value
             * matrix
             *
             */
            class Block {
              public:
                /**
                 * @brief Construct a new Block object
                 *
                 * @param mat
                 * @param b
                 * @param nstr
                 * @param nlyr
                 */
                Block(BVPMatrix& mat, BoundaryIndex b, uint nstr, uint nlyr)
                    : m_mat(mat) {
                    m_row_offset = (b == 0) ? 0
                                            : nstr / 2 * NSTOKES +
                                                  (b - 1) * nstr * NSTOKES;
                    m_col_offset = (b == 0) ? 0
                                   : (b == nlyr)
                                       ? m_mat.M_NCOLS - nstr * NSTOKES
                                       : (b - 1) * nstr * NSTOKES;
                }

              public:
                inline double& operator()(StreamIndex i, SolutionIndex j) {
                    return m_mat(m_row_offset + i, m_col_offset + j);
                }
                inline double operator()(StreamIndex i, SolutionIndex j) const {
                    return m_mat(m_row_offset + i, m_col_offset + j);
                }

              private:
                uint m_row_offset;
                uint m_col_offset;
                BVPMatrix& m_mat;
            };

            /**
             * @brief Get the Block object for a given boundary index
             *
             * @param b
             * @return Block
             */
            inline Block getBlock(BoundaryIndex b) {
                return Block(*this, b, M_NSTR, M_NLYR);
            }

            /**
             * @brief Pointer to the underlying data
             *
             * @return double*
             */
            inline double* data() { return m_data.get(); }

            /**
             * @brief Number of elements in data
             *
             * @return size_t
             */
            inline size_t size() { return M_NSTR * NSTOKES * M_NLYR; }

            /**
             * @brief Number of elements above/below the diagonal
             *
             * @return lapack_int
             */
            inline lapack_int NCD() const { return M_NCD; }

            /**
             * @brief Number of columns (band storage)
             *
             * @return lapack_int
             */
            inline lapack_int N() const { return M_NCOLS; }

            /**
             * @brief Leading dimension of the matrix (band storage)
             *
             * @return lapack_int
             */
            inline lapack_int LD() const { return M_NROWS; }

          private: // Members
            // Number of elements above/below the diagonal
            const uint M_NCD; /** Number of elements above/below the diagonal */

            const uint M_NRM1; /** Number of rows minus 1 */

            const uint M_NROWS; /** Number of rows in band-storage scheme*/

            const uint M_NCOLS; /** Number of columns in the original matrix
                                   (also in band-storage scheme)*/

            const uint M_NSTR;
            const uint M_NLYR;
            std::unique_ptr<double[]> m_data;
        };

        /**
         * Equivalent to the LAPACK routine dgbsv, but with the additional
         * constraint that the matrix is pentadiagonal this is used for the
         * special 2-stream case
         *
         * @param N
         * @param NRHS
         * @param AB
         * @param B
         * @param LDB
         * @param alpha
         * @param beta
         * @param z
         * @param gamma
         * @param mu
         * @param transpose
         * @return int
         */
        int dgbsv_pentadiagonal(int N, int NRHS, double* AB, double* B, int LDB,
                                Eigen::VectorXd& alpha, Eigen::VectorXd& beta,
                                Eigen::MatrixXd& z, Eigen::VectorXd& gamma,
                                Eigen::VectorXd& mu, bool transpose);

    } // namespace la

    /**
     * @brief Dense blocks of the BVP matrix used to propagate derivatives
     *
     * @tparam NSTOKES
     * @tparam CNSTR
     */
    template <int NSTOKES, int CNSTR = -1> class BVPMatrixDenseBlock {
      public:
        /**
         * @brief Construct a new BVPMatrixDenseBlock object
         *
         * @param b
         * @param nstr
         * @param nlyr
         */
        BVPMatrixDenseBlock(BoundaryIndex b, uint nstr, uint nlyr) {
            m_b = b;
            if (b + 1 == nlyr) {
                m_upper_data.resize(nstr / 2 * NSTOKES, nstr * NSTOKES);
            } else {
                m_upper_data.resize(nstr * NSTOKES, 2 * nstr * NSTOKES);
            }
            m_upper_row_offset = nstr / 2 * NSTOKES + (b)*nstr * NSTOKES;
            m_upper_col_offset = (b + 1 == nlyr)
                                     ? (nlyr * nstr * NSTOKES) - nstr * NSTOKES
                                     : (b)*nstr * NSTOKES;
            if (m_b == 0) {
                m_data.resize(nstr * NSTOKES / 2, nstr * NSTOKES);
            } else {
                m_data.resize(nstr * NSTOKES, 2 * nstr * NSTOKES);
            }
            m_row_offset =
                (b == 0) ? 0 : nstr / 2 * NSTOKES + (b - 1) * nstr * NSTOKES;
            m_col_offset = (b == 0) ? 0
                           : (b == nlyr)
                               ? (nlyr * nstr * NSTOKES) - nstr * NSTOKES
                               : (b - 1) * nstr * NSTOKES;
        }

        /**
         * @brief Upper storage
         *
         * @return Eigen::MatrixXd&
         */
        Eigen::MatrixXd& upper() { return m_upper_data; }

        /**
         * @brief Layer storage
         *
         * @return Eigen::MatrixXd&
         */
        Eigen::MatrixXd& layer() { return m_data; }

        /**
         * @brief Set the block to zero
         *
         */
        inline void setZero() {
            m_data.setZero();
            m_upper_data.setZero();
        }

        /**
         * @brief Assigns the right hand side of the boundary value problem
         *
         * Basically takes d_b and subtracts off this block * b
         *
         * @param colindex
         * @param d_b
         * @param b
         */
        inline void assign_rhs_d_bvp(int colindex, Eigen::MatrixXd& d_b,
                                     const Eigen::VectorXd& b) {

            // Then subtract off the multiplication terms
            for (uint i = 0; i < m_data.rows(); ++i) {
                for (uint j = 0; j < m_data.cols(); ++j) {
                    d_b(i + m_row_offset, colindex) -=
                        b[j + m_col_offset] * m_data(i, j);
                }
            }

            for (uint i = 0; i < m_upper_data.rows(); ++i) {
                for (uint j = 0; j < m_upper_data.cols(); ++j) {
                    d_b(i + m_upper_row_offset, colindex) -=
                        b[j + m_upper_col_offset] * m_upper_data(i, j);
                }
            }
        }

      private:
        Eigen::MatrixXd m_upper_data;
        Eigen::MatrixXd m_data;
        BoundaryIndex m_b;
        uint m_row_offset;
        uint m_col_offset;

        uint m_upper_row_offset;
        uint m_upper_col_offset;
    };

} // namespace sasktran_disco
