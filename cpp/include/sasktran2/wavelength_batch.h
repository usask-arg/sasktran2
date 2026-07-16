#pragma once

#include <sasktran2/dual.h>
#include <sasktran2/internal_common.h>

namespace sasktran2 {
    struct WavelengthBatch {
        int start = 0;
        int count = 0;

        int end() const { return start + count; }
        int wavelength(int lane) const { return start + lane; }
    };

    template <int NSTOKES> class WavelengthBatchDual {
      public:
        using ValueStorage =
            Eigen::Matrix<double, NSTOKES, Eigen::Dynamic, Eigen::RowMajor>;
        using DerivativeStorage =
            Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic,
                          Eigen::RowMajor>;

        ValueStorage value;
        DerivativeStorage deriv;

        void resize(int batch_capacity, int num_derivatives,
                    bool set_zero = false) {
            m_batch_capacity = batch_capacity;
            m_num_derivatives = num_derivatives;
            value.resize(NSTOKES, batch_capacity);
            deriv.resize(NSTOKES * num_derivatives, batch_capacity);
            if (set_zero) {
                value.setZero();
                deriv.setZero();
            }
        }

        int batch_capacity() const { return m_batch_capacity; }
        int derivative_size() const { return m_num_derivatives; }

        void set_zero(int count) {
            value.leftCols(count).setZero();
            deriv.leftCols(count).setZero();
        }

        auto derivative(int derivative_index, int count) {
            return deriv.block(derivative_index * NSTOKES, 0, NSTOKES, count);
        }

        auto derivative(int derivative_index, int count) const {
            return deriv.block(derivative_index * NSTOKES, 0, NSTOKES, count);
        }

        auto derivative_range(int derivative_start, int derivative_count,
                              int wavelength_count) {
            return deriv.block(derivative_start * NSTOKES, 0,
                               derivative_count * NSTOKES, wavelength_count);
        }

        auto derivative_range(int derivative_start, int derivative_count,
                              int wavelength_count) const {
            return deriv.block(derivative_start * NSTOKES, 0,
                               derivative_count * NSTOKES, wavelength_count);
        }

        void copy_lane_to(int lane,
                          sasktran2::Dual<double, sasktran2::dualstorage::dense,
                                          NSTOKES>& target) const {
            target.value = value.col(lane);
            for (int derivative_index = 0; derivative_index < m_num_derivatives;
                 ++derivative_index) {
                target.deriv.col(derivative_index) =
                    derivative(derivative_index, m_batch_capacity).col(lane);
            }
        }

      private:
        int m_batch_capacity = 0;
        int m_num_derivatives = 0;
    };

    struct WavelengthBatchODView {
        const double* od_data;
        const double* exp_minus_od_data;
        int count;
        Eigen::SparseMatrix<double, Eigen::RowMajor>::InnerIterator deriv_iter;

        WavelengthBatchODView(
            const double* od, const double* exp_minus_od, int wavelength_count,
            const Eigen::SparseMatrix<double, Eigen::RowMajor>& deriv_matrix,
            int row)
            : od_data(od), exp_minus_od_data(exp_minus_od),
              count(wavelength_count), deriv_iter(deriv_matrix, row) {}

        Eigen::Map<const Eigen::RowVectorXd> od() const {
            return {od_data, count};
        }

        Eigen::Map<const Eigen::RowVectorXd> exp_minus_od() const {
            return {exp_minus_od_data, count};
        }
    };
} // namespace sasktran2
