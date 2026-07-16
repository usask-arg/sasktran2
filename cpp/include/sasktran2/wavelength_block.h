#pragma once

#include <sasktran2/dual.h>
#include <sasktran2/internal_common.h>

namespace sasktran2 {
    /** A contiguous wavelength range and the storage reserved for it.
     *
     * count is the number of active lanes. capacity remains the configured
     * block width for a partial tail, allowing a one-lane tail to continue
     * using block storage. A capacity of one selects the scalar storage path.
     */
    struct WavelengthBlock {
        int start = 0;
        int count = 0;
        int capacity = 1;

        int end() const { return start + count; }
        int wavelength(int lane) const { return start + lane; }
        bool is_scalar() const { return capacity == 1; }
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

    template <int NSTOKES> class WavelengthBlockLaneDualView {
      private:
        using ValueStride = Eigen::InnerStride<Eigen::Dynamic>;
        using DerivativeStride = Eigen::Stride<Eigen::Dynamic, Eigen::Dynamic>;

      public:
        using ValueMap =
            Eigen::Map<const Eigen::Vector<double, NSTOKES>, 0, ValueStride>;
        using DerivativeMap =
            Eigen::Map<const Eigen::Matrix<double, NSTOKES, Eigen::Dynamic>, 0,
                       DerivativeStride>;

        WavelengthBlockLaneDualView(const WavelengthBlockDual<NSTOKES>& block,
                                    int lane)
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

    template <int NSTOKES> class WavelengthBlockDualView {
      public:
        using ScalarDual =
            sasktran2::Dual<double, sasktran2::dualstorage::dense, NSTOKES>;

        explicit WavelengthBlockDualView(ScalarDual& scalar)
            : m_scalar(&scalar) {}

        explicit WavelengthBlockDualView(WavelengthBlockDual<NSTOKES>& block)
            : m_block(&block) {}

        bool is_scalar() const { return m_scalar != nullptr; }

        ScalarDual& scalar() const {
            if (m_scalar == nullptr) {
                throw std::logic_error(
                    "Wavelength block does not use scalar storage");
            }
            return *m_scalar;
        }

        WavelengthBlockDual<NSTOKES>& block() const {
            if (m_block == nullptr) {
                throw std::logic_error(
                    "Wavelength block does not use block storage");
            }
            return *m_block;
        }

        WavelengthBlockLaneDualView<NSTOKES> lane(int lane_index) const {
            return {block(), lane_index};
        }

      private:
        ScalarDual* m_scalar = nullptr;
        WavelengthBlockDual<NSTOKES>* m_block = nullptr;
    };

    class WavelengthBlockODView {
      public:
        explicit WavelengthBlockODView(const SparseODDualView& scalar)
            : m_scalar(&scalar), m_od_data(&scalar.od),
              m_exp_minus_od_data(&scalar.exp_minus_od), m_count(1) {}

        WavelengthBlockODView(
            const double* od, const double* exp_minus_od, int wavelength_count,
            const Eigen::SparseMatrix<double, Eigen::RowMajor>& deriv_matrix,
            int row)
            : m_od_data(od), m_exp_minus_od_data(exp_minus_od),
              m_count(wavelength_count),
              m_deriv_iter(std::in_place, deriv_matrix, row) {}

        bool is_scalar() const { return m_scalar != nullptr; }

        const SparseODDualView& scalar() const {
            if (m_scalar == nullptr) {
                throw std::logic_error(
                    "Wavelength block optical depth is not scalar");
            }
            return *m_scalar;
        }

        Eigen::Map<const Eigen::RowVectorXd> od() const {
            return {m_od_data, m_count};
        }

        Eigen::Map<const Eigen::RowVectorXd> exp_minus_od() const {
            return {m_exp_minus_od_data, m_count};
        }

        double od(int lane) const { return m_od_data[lane]; }

        auto derivative_iterator() const {
            return is_scalar() ? m_scalar->deriv_iter : *m_deriv_iter;
        }

        int count() const { return m_count; }

      private:
        const SparseODDualView* m_scalar = nullptr;
        const double* m_od_data = nullptr;
        const double* m_exp_minus_od_data = nullptr;
        int m_count = 0;
        std::optional<
            Eigen::SparseMatrix<double, Eigen::RowMajor>::InnerIterator>
            m_deriv_iter;
    };
} // namespace sasktran2
