#pragma once

#include <sasktran2/dual.h>
#include <sasktran2/internal_common.h>

namespace sasktran2 {
    /** A contiguous range of wavelengths processed together. */
    struct WavelengthBlock {
        int start = 0;
        int count = 0;

        int end() const { return start + count; }
        int wavelength(int lane) const { return start + lane; }
    };

    template <int NSTOKES> class WavelengthBlockDual {
      public:
        using ValueStorage =
            Eigen::Matrix<double, NSTOKES, Eigen::Dynamic, Eigen::RowMajor>;
        using DerivativeStorage =
            Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic,
                          Eigen::RowMajor>;

        ValueStorage value;
        DerivativeStorage deriv;

        void resize(int block_capacity, int num_derivatives,
                    bool set_zero = false) {
            m_block_capacity = block_capacity;
            m_num_derivatives = num_derivatives;
            value.resize(NSTOKES, block_capacity);
            deriv.resize(NSTOKES * num_derivatives, block_capacity);
            if (set_zero) {
                value.setZero();
                deriv.setZero();
            }
        }

        int block_capacity() const { return m_block_capacity; }
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

      private:
        int m_block_capacity = 0;
        int m_num_derivatives = 0;
    };

    /** Mutable Dual-like view of one lane in wavelength-block storage. */
    template <int NSTOKES> class WavelengthBlockLaneDualView {
      private:
        using ValueStride = Eigen::InnerStride<Eigen::Dynamic>;
        using DerivativeStride = Eigen::Stride<Eigen::Dynamic, Eigen::Dynamic>;

      public:
        using ValueMap =
            Eigen::Map<Eigen::Vector<double, NSTOKES>, 0, ValueStride>;
        using DerivativeMap =
            Eigen::Map<Eigen::Matrix<double, NSTOKES, Eigen::Dynamic>, 0,
                       DerivativeStride>;

        WavelengthBlockLaneDualView(WavelengthBlockDual<NSTOKES>& block,
                                    int lane)
            : value(block.value.data() + lane, NSTOKES,
                    ValueStride(block.block_capacity())),
              deriv(block.deriv.data() + lane, NSTOKES, block.derivative_size(),
                    DerivativeStride(NSTOKES * block.block_capacity(),
                                     block.block_capacity())) {}

        ValueMap value;
        DerivativeMap deriv;

        int derivative_size() const { return deriv.cols(); }

        auto d_extinction(int n) {
            return deriv(Eigen::placeholders::all, Eigen::seq(0, n - 1));
        }

        auto d_ssa(int n) {
            return deriv(Eigen::placeholders::all, Eigen::seq(n, 2 * n - 1));
        }

        auto d_scatterer(int n, int scatterer_index) {
            return deriv(Eigen::placeholders::all,
                         Eigen::seq((2 + scatterer_index) * n,
                                    (3 + scatterer_index) * n - 1));
        }

        auto d_emission(int n, int num_scattering_deriv_groups) {
            return deriv(Eigen::placeholders::all,
                         Eigen::seq((2 + num_scattering_deriv_groups) * n,
                                    (3 + num_scattering_deriv_groups) * n - 1));
        }

        auto d_brdf(int n, int num_source_groups, int num_brdf_args) {
            return deriv(
                Eigen::placeholders::all,
                Eigen::seq((2 + num_source_groups) * n,
                           (2 + num_source_groups) * n + num_brdf_args - 1));
        }
    };

    /** Read-only Dual-like view of one lane in wavelength-block storage. */
    template <int NSTOKES> class WavelengthBlockConstLaneDualView {
      private:
        using ValueStride = Eigen::InnerStride<Eigen::Dynamic>;
        using DerivativeStride = Eigen::Stride<Eigen::Dynamic, Eigen::Dynamic>;

      public:
        using ValueMap =
            Eigen::Map<const Eigen::Vector<double, NSTOKES>, 0, ValueStride>;
        using DerivativeMap =
            Eigen::Map<const Eigen::Matrix<double, NSTOKES, Eigen::Dynamic>, 0,
                       DerivativeStride>;

        WavelengthBlockConstLaneDualView(
            const WavelengthBlockDual<NSTOKES>& block, int lane)
            : value(block.value.data() + lane, NSTOKES,
                    ValueStride(block.block_capacity())),
              deriv(block.deriv.data() + lane, NSTOKES, block.derivative_size(),
                    DerivativeStride(NSTOKES * block.block_capacity(),
                                     block.block_capacity())) {}

        ValueMap value;
        DerivativeMap deriv;

        auto d_extinction(int n) const {
            return deriv(Eigen::placeholders::all, Eigen::seq(0, n - 1));
        }

        auto d_ssa(int n) const {
            return deriv(Eigen::placeholders::all, Eigen::seq(n, 2 * n - 1));
        }

        auto d_scatterer(int n, int scatterer_index) const {
            return deriv(Eigen::placeholders::all,
                         Eigen::seq((2 + scatterer_index) * n,
                                    (3 + scatterer_index) * n - 1));
        }

        auto d_emission(int n, int num_scattering_deriv_groups) const {
            return deriv(Eigen::placeholders::all,
                         Eigen::seq((2 + num_scattering_deriv_groups) * n,
                                    (3 + num_scattering_deriv_groups) * n - 1));
        }

        auto d_brdf(int n, int num_source_groups, int num_brdf_args) const {
            return deriv(
                Eigen::placeholders::all,
                Eigen::seq((2 + num_source_groups) * n,
                           (2 + num_source_groups) * n + num_brdf_args - 1));
        }
    };

    class WavelengthBlockODView {
      public:
        WavelengthBlockODView(
            const double* od, const double* exp_minus_od, int wavelength_count,
            const Eigen::SparseMatrix<double, Eigen::RowMajor>& deriv_matrix,
            int row)
            : m_od_data(od), m_exp_minus_od_data(exp_minus_od),
              m_count(wavelength_count), m_deriv_iter(deriv_matrix, row) {}

        Eigen::Map<const Eigen::RowVectorXd> od() const {
            return {m_od_data, m_count};
        }

        Eigen::Map<const Eigen::RowVectorXd> exp_minus_od() const {
            return {m_exp_minus_od_data, m_count};
        }

        double od(int lane) const { return m_od_data[lane]; }
        double exp_minus_od(int lane) const {
            return m_exp_minus_od_data[lane];
        }

        auto derivative_iterator() const { return m_deriv_iter; }

        int count() const { return m_count; }

      private:
        const double* m_od_data = nullptr;
        const double* m_exp_minus_od_data = nullptr;
        int m_count = 0;
        Eigen::SparseMatrix<double, Eigen::RowMajor>::InnerIterator
            m_deriv_iter;
    };
} // namespace sasktran2
