#pragma once

#include <sasktran2/dual.h>
#include <sasktran2/internal_common.h>

namespace sasktran2 {
    /** A contiguous range of wavelengths processed together.
     *
     * A fixed N lets hot kernels expose their width to Eigen and the compiler.
     * Eigen::Dynamic is the type-erased form used at virtual/API boundaries.
     */
    template <int N = Eigen::Dynamic> struct WavelengthBlock {
        static_assert(N == Eigen::Dynamic || N > 0,
                      "Wavelength block size must be positive or dynamic");

        static constexpr int static_size = N;

        int start = 0;
        int count = N == Eigen::Dynamic ? 0 : N;

        WavelengthBlock() = default;

        WavelengthBlock(int block_start, int block_count)
            : start(block_start), count(block_count) {
            if constexpr (N != Eigen::Dynamic) {
                if (block_count != N) {
                    throw std::invalid_argument(
                        "Runtime wavelength block size does not match its "
                        "static size");
                }
            }
        }

        template <int M,
                  std::enable_if_t<N == Eigen::Dynamic && M != Eigen::Dynamic,
                                   int> = 0>
        WavelengthBlock(const WavelengthBlock<M>& other)
            : start(other.start), count(other.count) {}

        int end() const { return start + count; }
        int wavelength(int lane) const { return start + lane; }

        WavelengthBlock<> dynamic() const { return {start, count}; }
    };

    /** Dispatches supported fixed widths and preserves an arbitrary-width
     * dynamic fallback. Dispatch is based on the active count, so tails also
     * use fixed-width kernels when possible. */
    template <typename Function>
    decltype(auto) dispatch_wavelength_block(const WavelengthBlock<>& block,
                                             Function&& function) {
        switch (block.count) {
        case 1:
            return std::forward<Function>(function)(
                WavelengthBlock<1>{block.start, block.count});
        case 4:
            return std::forward<Function>(function)(
                WavelengthBlock<4>{block.start, block.count});
        // Larger fixed Eigen expressions increased register pressure in the
        // single-scatter derivative kernels. Keep them on the dynamic path;
        // fixed 1- and 4-wide tails still avoid scalar redispatch overhead.
        default:
            return std::forward<Function>(function)(block);
        }
    }

    template <int N, typename Derived>
    auto wavelength_head(Derived&& vector, const WavelengthBlock<N>& block) {
        if constexpr (N == Eigen::Dynamic) {
            return std::forward<Derived>(vector).head(block.count);
        } else {
            return std::forward<Derived>(vector).template head<N>();
        }
    }

    template <int N, typename Derived>
    auto wavelength_left_cols(Derived&& matrix,
                              const WavelengthBlock<N>& block) {
        if constexpr (N == Eigen::Dynamic) {
            return std::forward<Derived>(matrix).leftCols(block.count);
        } else {
            return std::forward<Derived>(matrix).template leftCols<N>();
        }
    }

    template <int N, typename Derived>
    auto wavelength_middle_cols(Derived&& matrix,
                                const WavelengthBlock<N>& block) {
        if constexpr (N == Eigen::Dynamic) {
            return std::forward<Derived>(matrix).middleCols(block.start,
                                                            block.count);
        } else {
            return std::forward<Derived>(matrix).template middleCols<N>(
                block.start);
        }
    }

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

        template <int N> void set_zero(const WavelengthBlock<N>& block) {
            wavelength_left_cols(value, block).setZero();
            wavelength_left_cols(deriv, block).setZero();
        }

        auto derivative(int derivative_index, int count) {
            return deriv.block(derivative_index * NSTOKES, 0, NSTOKES, count);
        }

        auto derivative(int derivative_index, int count) const {
            return deriv.block(derivative_index * NSTOKES, 0, NSTOKES, count);
        }

        template <int N>
        auto derivative(int derivative_index, const WavelengthBlock<N>& block) {
            if constexpr (N == Eigen::Dynamic) {
                return deriv.block(derivative_index * NSTOKES, 0, NSTOKES,
                                   block.count);
            } else {
                return deriv.template block<NSTOKES, N>(
                    derivative_index * NSTOKES, 0);
            }
        }

        template <int N>
        auto derivative(int derivative_index,
                        const WavelengthBlock<N>& block) const {
            if constexpr (N == Eigen::Dynamic) {
                return deriv.block(derivative_index * NSTOKES, 0, NSTOKES,
                                   block.count);
            } else {
                return deriv.template block<NSTOKES, N>(
                    derivative_index * NSTOKES, 0);
            }
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

        template <int N>
        auto derivative_range(int derivative_start, int derivative_count,
                              const WavelengthBlock<N>& block) {
            if constexpr (N == Eigen::Dynamic) {
                return deriv.block(derivative_start * NSTOKES, 0,
                                   derivative_count * NSTOKES, block.count);
            } else {
                return deriv.template block<Eigen::Dynamic, N>(
                    derivative_start * NSTOKES, 0, derivative_count * NSTOKES,
                    N);
            }
        }

        template <int N>
        auto derivative_range(int derivative_start, int derivative_count,
                              const WavelengthBlock<N>& block) const {
            if constexpr (N == Eigen::Dynamic) {
                return deriv.block(derivative_start * NSTOKES, 0,
                                   derivative_count * NSTOKES, block.count);
            } else {
                return deriv.template block<Eigen::Dynamic, N>(
                    derivative_start * NSTOKES, 0, derivative_count * NSTOKES,
                    N);
            }
        }

        using ConstDerivativeStokesMap =
            Eigen::Map<const Eigen::Matrix<double, Eigen::Dynamic,
                                           Eigen::Dynamic, Eigen::RowMajor>,
                       0, Eigen::Stride<Eigen::Dynamic, Eigen::Dynamic>>;

        ConstDerivativeStokesMap derivative_stokes(int derivative_start,
                                                   int derivative_count,
                                                   int stokes,
                                                   int wavelength_count) const {
            static const double dummy = 0.0;
            const double* derivative_data =
                m_num_derivatives == 0
                    ? &dummy
                    : deriv.data() + (derivative_start * NSTOKES + stokes) *
                                         m_block_capacity;
            return {derivative_data, derivative_count, wavelength_count,
                    Eigen::Stride<Eigen::Dynamic, Eigen::Dynamic>(
                        NSTOKES * m_block_capacity, 1)};
        }

      private:
        int m_block_capacity = 0;
        int m_num_derivatives = 0;
    };

    /** Mutable Dual-like view of one lane in wavelength-block storage. */
    template <int NSTOKES, int BlockCapacity = Eigen::Dynamic>
    class WavelengthBlockLaneDualView {
      private:
        static_assert(BlockCapacity == Eigen::Dynamic || BlockCapacity > 0,
                      "Block capacity must be positive or dynamic");
        static constexpr int DerivativeOuterStride =
            BlockCapacity == Eigen::Dynamic ? Eigen::Dynamic
                                            : NSTOKES * BlockCapacity;
        using ValueStride = Eigen::InnerStride<BlockCapacity>;
        using DerivativeStride =
            Eigen::Stride<DerivativeOuterStride, BlockCapacity>;

        static double* value_data(WavelengthBlockDual<NSTOKES>& block,
                                  int lane) {
            if constexpr (BlockCapacity != Eigen::Dynamic) {
                if (block.block_capacity() != BlockCapacity) {
                    throw std::invalid_argument(
                        "Wavelength block storage capacity does not match "
                        "its static lane view");
                }
            }
            return block.value.data() + lane;
        }

        static double* derivative_data(WavelengthBlockDual<NSTOKES>& block,
                                       int lane) {
            static double dummy = 0.0;
            return block.derivative_size() == 0 ? &dummy
                                                : block.deriv.data() + lane;
        }

      public:
        using ValueMap =
            Eigen::Map<Eigen::Vector<double, NSTOKES>, 0, ValueStride>;
        using DerivativeMap =
            Eigen::Map<Eigen::Matrix<double, NSTOKES, Eigen::Dynamic>, 0,
                       DerivativeStride>;

        WavelengthBlockLaneDualView(WavelengthBlockDual<NSTOKES>& block,
                                    int lane)
            : value(value_data(block, lane), NSTOKES,
                    ValueStride(block.block_capacity())),
              deriv(derivative_data(block, lane), NSTOKES,
                    block.derivative_size(),
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

        static const double*
        derivative_data(const WavelengthBlockDual<NSTOKES>& block, int lane) {
            static const double dummy = 0.0;
            return block.derivative_size() == 0 ? &dummy
                                                : block.deriv.data() + lane;
        }

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
              deriv(derivative_data(block, lane), NSTOKES,
                    block.derivative_size(),
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
