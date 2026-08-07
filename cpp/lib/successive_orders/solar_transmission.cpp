#include "solar_transmission.h"

#ifdef SKTRAN_RUST_SUPPORT
#include "sasktran2-core/src/successive_orders/cxx.rs.h"
#endif

#include <algorithm>
#include <limits>
#include <stdexcept>

namespace {
#ifdef SKTRAN_RUST_SUPPORT
    template <typename T>
    ::rust::Slice<const T> as_rust_slice(const std::vector<T>& values) {
        return {values.data(), values.size()};
    }

    template <typename T>
    ::rust::Slice<T> as_rust_mut_slice(std::vector<T>& values) {
        return {values.data(), values.size()};
    }

    ::rust::Slice<const double>
    as_rust_slice(Eigen::Ref<const Eigen::VectorXd> values) {
        return {values.data(), static_cast<std::size_t>(values.size())};
    }

    ::rust::Slice<double> as_rust_slice(Eigen::Ref<Eigen::VectorXd> values) {
        return {values.data(), static_cast<std::size_t>(values.size())};
    }

    template <typename Matrix>
    ::rust::Slice<const double> as_rust_matrix_slice(const Matrix& values) {
        return {values.data(), static_cast<std::size_t>(values.size())};
    }

    template <typename Matrix>
    ::rust::Slice<double> as_rust_mut_matrix_slice(Matrix& values) {
        return {values.data(), static_cast<std::size_t>(values.size())};
    }
#endif

    void pack_dense(const Eigen::MatrixXd& matrix,
                    std::vector<double>& values) {
        values.resize(static_cast<std::size_t>(matrix.rows()) *
                      static_cast<std::size_t>(matrix.cols()));
        for (Eigen::Index row = 0; row < matrix.rows(); ++row) {
            for (Eigen::Index column = 0; column < matrix.cols(); ++column) {
                values[static_cast<std::size_t>(row * matrix.cols() + column)] =
                    matrix(row, column);
            }
        }
    }

    void pack_sparse(
        const sasktran2::successive_orders::SolarTransmission::SparseMatrix&
            input,
        std::vector<std::uint32_t>& offsets,
        std::vector<std::uint32_t>& indices, std::vector<double>& values) {
        auto matrix = input;
        matrix.makeCompressed();
        if (matrix.nonZeros() >
                static_cast<Eigen::Index>(
                    std::numeric_limits<std::uint32_t>::max()) ||
            matrix.cols() > static_cast<Eigen::Index>(
                                std::numeric_limits<std::uint32_t>::max())) {
            throw std::length_error(
                "Solar transmission geometry exceeds packed Rust index "
                "range");
        }
        offsets.resize(static_cast<std::size_t>(matrix.rows()) + 1);
        for (Eigen::Index row = 0; row <= matrix.rows(); ++row) {
            offsets[static_cast<std::size_t>(row)] =
                static_cast<std::uint32_t>(matrix.outerIndexPtr()[row]);
        }
        indices.resize(static_cast<std::size_t>(matrix.nonZeros()));
        values.resize(static_cast<std::size_t>(matrix.nonZeros()));
        for (Eigen::Index index = 0; index < matrix.nonZeros(); ++index) {
            indices[static_cast<std::size_t>(index)] =
                static_cast<std::uint32_t>(matrix.innerIndexPtr()[index]);
            values[static_cast<std::size_t>(index)] = matrix.valuePtr()[index];
        }
    }
} // namespace

namespace sasktran2::successive_orders {
    struct SolarTransmission::Impl {
#ifdef SKTRAN_RUST_SUPPORT
        ::rust::Box<
            sasktran2::rust::successive_orders::RustSolarTransmissionOperator>
            rust_operator;
        mutable std::vector<std::vector<double>> worker_scratch;

        Impl(::rust::Box<sasktran2::rust::successive_orders::
                             RustSolarTransmissionOperator>&& rust_operator,
             int wavelength_workers)
            : rust_operator(std::move(rust_operator)),
              worker_scratch(static_cast<std::size_t>(wavelength_workers)) {
            const std::size_t scratch_size = sasktran2::rust::
                successive_orders::solar_transmission_forward_scratch_size(
                    *this->rust_operator);
            for (auto& scratch : worker_scratch) {
                scratch.resize(scratch_size);
            }
        }
#else
        explicit Impl(int) {}
#endif
    };

    namespace {
        std::shared_ptr<SolarTransmission>
        create(const Eigen::MatrixXd& dense_path,
               const SolarTransmission::SparseMatrix& sparse_path,
               const SolarTransmission::SparseMatrix* interpolation,
               const std::vector<bool>& ground_hit, std::size_t input_size,
               int wavelength_workers) {
            if (wavelength_workers < 1) {
                throw std::invalid_argument(
                    "Solar transmission requires a wavelength worker");
            }
#ifdef SKTRAN_RUST_SUPPORT
            std::size_t path_rows = 0;
            std::size_t path_columns = 0;
            std::vector<double> dense_path_values;
            std::vector<std::uint32_t> path_row_offsets;
            std::vector<std::uint32_t> path_column_indices;
            std::vector<double> sparse_path_values;
            if (dense_path.size() > 0) {
                path_rows = static_cast<std::size_t>(dense_path.rows());
                path_columns = static_cast<std::size_t>(dense_path.cols());
                pack_dense(dense_path, dense_path_values);
            } else {
                path_rows = static_cast<std::size_t>(sparse_path.rows());
                path_columns = static_cast<std::size_t>(sparse_path.cols());
                pack_sparse(sparse_path, path_row_offsets, path_column_indices,
                            sparse_path_values);
            }

            std::size_t interpolation_rows = 0;
            std::size_t interpolation_columns = 0;
            std::vector<std::uint32_t> interpolation_row_offsets;
            std::vector<std::uint32_t> interpolation_column_indices;
            std::vector<double> interpolation_values;
            if (interpolation != nullptr) {
                interpolation_rows =
                    static_cast<std::size_t>(interpolation->rows());
                interpolation_columns =
                    static_cast<std::size_t>(interpolation->cols());
                pack_sparse(*interpolation, interpolation_row_offsets,
                            interpolation_column_indices, interpolation_values);
            }

            std::vector<std::uint8_t> packed_ground_hit(ground_hit.size());
            std::transform(
                ground_hit.begin(), ground_hit.end(), packed_ground_hit.begin(),
                [](bool hit) { return static_cast<std::uint8_t>(hit); });
            auto rust_operator = sasktran2::rust::successive_orders::
                new_solar_transmission_operator(
                    path_rows, path_columns, as_rust_slice(dense_path_values),
                    as_rust_slice(path_row_offsets),
                    as_rust_slice(path_column_indices),
                    as_rust_slice(sparse_path_values), interpolation_rows,
                    interpolation_columns,
                    as_rust_slice(interpolation_row_offsets),
                    as_rust_slice(interpolation_column_indices),
                    as_rust_slice(interpolation_values),
                    as_rust_slice(packed_ground_hit));
            const std::size_t output_size = sasktran2::rust::successive_orders::
                solar_transmission_output_size(*rust_operator);
            if (sasktran2::rust::successive_orders::
                        solar_transmission_input_size(*rust_operator) !=
                    input_size ||
                output_size != ground_hit.size()) {
                throw std::logic_error(
                    "Packed Rust solar transmission dimensions do not match "
                    "the source geometry");
            }
            return std::shared_ptr<SolarTransmission>(
                new SolarTransmission(std::make_unique<SolarTransmission::Impl>(
                    std::move(rust_operator), wavelength_workers)));
#else
            (void)dense_path;
            (void)sparse_path;
            (void)interpolation;
            (void)ground_hit;
            (void)input_size;
            throw std::logic_error(
                "Rust solar transmission requires Rust support");
#endif
        }
    } // namespace

    SolarTransmission::SolarTransmission(std::unique_ptr<Impl> impl)
        : m_impl(std::move(impl)) {}

    SolarTransmission::~SolarTransmission() = default;
    SolarTransmission::SolarTransmission(SolarTransmission&&) noexcept =
        default;
    SolarTransmission&
    SolarTransmission::operator=(SolarTransmission&&) noexcept = default;

    std::shared_ptr<SolarTransmission> SolarTransmission::create_exact(
        const Eigen::MatrixXd& dense_path, const SparseMatrix& sparse_path,
        const std::vector<bool>& ground_hit, std::size_t input_size,
        int wavelength_workers) {
        return create(dense_path, sparse_path, nullptr, ground_hit, input_size,
                      wavelength_workers);
    }

    std::shared_ptr<SolarTransmission> SolarTransmission::create_table(
        const Eigen::MatrixXd& dense_path, const SparseMatrix& interpolation,
        const std::vector<bool>& ground_hit, std::size_t input_size,
        int wavelength_workers) {
        const SparseMatrix unused_path;
        return create(dense_path, unused_path, &interpolation, ground_hit,
                      input_size, wavelength_workers);
    }

    std::size_t SolarTransmission::input_size() const {
#ifdef SKTRAN_RUST_SUPPORT
        return sasktran2::rust::successive_orders::
            solar_transmission_input_size(*m_impl->rust_operator);
#else
        return 0;
#endif
    }

    std::size_t SolarTransmission::output_size() const {
#ifdef SKTRAN_RUST_SUPPORT
        return sasktran2::rust::successive_orders::
            solar_transmission_output_size(*m_impl->rust_operator);
#else
        return 0;
#endif
    }

    std::size_t SolarTransmission::storage_bytes() const {
#ifdef SKTRAN_RUST_SUPPORT
        return sasktran2::rust::successive_orders::
            solar_transmission_storage_bytes(*m_impl->rust_operator);
#else
        return 0;
#endif
    }

    std::size_t SolarTransmission::forward_scratch_bytes() const {
#ifdef SKTRAN_RUST_SUPPORT
        constexpr std::size_t lanes = 4;
        return lanes *
               sasktran2::rust::successive_orders::
                   solar_transmission_forward_scratch_size(
                       *m_impl->rust_operator) *
               m_impl->worker_scratch.size() * sizeof(double);
#else
        return 0;
#endif
    }

    std::size_t SolarTransmission::vjp_scratch_bytes() const {
#ifdef SKTRAN_RUST_SUPPORT
        constexpr std::size_t lanes = 4;
        const std::size_t scalar_size =
            sasktran2::rust::successive_orders::solar_transmission_scratch_size(
                *m_impl->rust_operator);
        const std::size_t packed_size =
            lanes *
            sasktran2::rust::successive_orders::
                solar_transmission_forward_scratch_size(*m_impl->rust_operator);
        return std::max(scalar_size, packed_size) *
               m_impl->worker_scratch.size() * sizeof(double);
#else
        return 0;
#endif
    }

    void SolarTransmission::calculate(
        int worker, Eigen::Ref<const Eigen::VectorXd> extinction,
        double solar_irradiance, Eigen::Ref<Eigen::VectorXd> transmission) {
#ifdef SKTRAN_RUST_SUPPORT
        auto& scratch = m_impl->worker_scratch.at(worker);
        scratch.resize(sasktran2::rust::successive_orders::
                           solar_transmission_forward_scratch_size(
                               *m_impl->rust_operator));
        sasktran2::rust::successive_orders::calculate_solar_transmission(
            *m_impl->rust_operator, as_rust_slice(extinction), solar_irradiance,
            as_rust_slice(transmission), as_rust_mut_slice(scratch));
#else
        (void)worker;
        (void)extinction;
        (void)solar_irradiance;
        (void)transmission;
#endif
    }

    void SolarTransmission::calculate_batch4(
        int worker, const double* extinction, std::size_t extinction_size,
        const double* solar_irradiance, std::size_t irradiance_size,
        double* transmission, std::size_t transmission_size) {
#ifdef SKTRAN_RUST_SUPPORT
        auto& scratch = m_impl->worker_scratch.at(worker);
        scratch.resize(sasktran2::rust::successive_orders::
                           solar_transmission_forward_scratch_size(
                               *m_impl->rust_operator) *
                       4);
        sasktran2::rust::successive_orders::calculate_solar_transmission_batch4(
            *m_impl->rust_operator, {extinction, extinction_size},
            {solar_irradiance, irradiance_size},
            {transmission, transmission_size}, as_rust_mut_slice(scratch));
#else
        (void)worker;
        (void)extinction;
        (void)extinction_size;
        (void)solar_irradiance;
        (void)irradiance_size;
        (void)transmission;
        (void)transmission_size;
#endif
    }

    void SolarTransmission::calculate_jvp(
        int worker, Eigen::Ref<const Eigen::VectorXd> extinction_tangent,
        Eigen::Ref<const Eigen::VectorXd> transmission,
        Eigen::Ref<Eigen::VectorXd> transmission_tangent) {
#ifdef SKTRAN_RUST_SUPPORT
        auto& scratch = m_impl->worker_scratch.at(worker);
        scratch.resize(sasktran2::rust::successive_orders::
                           solar_transmission_forward_scratch_size(
                               *m_impl->rust_operator));
        sasktran2::rust::successive_orders::calculate_solar_transmission_jvp(
            *m_impl->rust_operator, as_rust_slice(extinction_tangent),
            as_rust_slice(transmission), as_rust_slice(transmission_tangent),
            as_rust_mut_slice(scratch));
#else
        (void)worker;
        (void)extinction_tangent;
        (void)transmission;
        (void)transmission_tangent;
#endif
    }

    void SolarTransmission::calculate_jvp_batch4(
        int worker, Eigen::Ref<const RowMajorMatrix> extinction_tangent,
        Eigen::Ref<const RowMajorMatrix> transmission,
        Eigen::Ref<RowMajorMatrix> transmission_tangent) {
#ifdef SKTRAN_RUST_SUPPORT
        auto& scratch = m_impl->worker_scratch.at(worker);
        scratch.resize(sasktran2::rust::successive_orders::
                           solar_transmission_forward_scratch_size(
                               *m_impl->rust_operator) *
                       4);
        sasktran2::rust::successive_orders::
            calculate_solar_transmission_jvp_batch4(
                *m_impl->rust_operator,
                as_rust_matrix_slice(extinction_tangent),
                as_rust_matrix_slice(transmission),
                as_rust_mut_matrix_slice(transmission_tangent),
                as_rust_mut_slice(scratch));
#else
        (void)worker;
        (void)extinction_tangent;
        (void)transmission;
        (void)transmission_tangent;
#endif
    }

    void SolarTransmission::accumulate_vjp(
        int worker, Eigen::Ref<const Eigen::VectorXd> transmission,
        Eigen::Ref<const Eigen::VectorXd> transmission_cotangent,
        Eigen::Ref<Eigen::VectorXd> extinction_gradient) {
#ifdef SKTRAN_RUST_SUPPORT
        auto& scratch = m_impl->worker_scratch.at(worker);
        scratch.resize(
            sasktran2::rust::successive_orders::solar_transmission_scratch_size(
                *m_impl->rust_operator));
        sasktran2::rust::successive_orders::accumulate_solar_transmission_vjp(
            *m_impl->rust_operator, as_rust_slice(transmission),
            as_rust_slice(transmission_cotangent),
            as_rust_slice(extinction_gradient), as_rust_mut_slice(scratch));
        scratch.resize(sasktran2::rust::successive_orders::
                           solar_transmission_forward_scratch_size(
                               *m_impl->rust_operator));
#else
        (void)worker;
        (void)transmission;
        (void)transmission_cotangent;
        (void)extinction_gradient;
#endif
    }

    void SolarTransmission::accumulate_vjp_batch4(
        int worker, Eigen::Ref<const RowMajorMatrix> transmission,
        Eigen::Ref<RowMajorMatrix> transmission_cotangent,
        Eigen::Ref<RowMajorMatrix> extinction_gradient) {
#ifdef SKTRAN_RUST_SUPPORT
        auto& scratch = m_impl->worker_scratch.at(worker);
        scratch.resize(sasktran2::rust::successive_orders::
                           solar_transmission_forward_scratch_size(
                               *m_impl->rust_operator) *
                       4);
        sasktran2::rust::successive_orders::
            accumulate_solar_transmission_vjp_batch4(
                *m_impl->rust_operator, as_rust_matrix_slice(transmission),
                as_rust_mut_matrix_slice(transmission_cotangent),
                as_rust_mut_matrix_slice(extinction_gradient),
                as_rust_mut_slice(scratch));
#else
        (void)worker;
        (void)transmission;
        (void)transmission_cotangent;
        (void)extinction_gradient;
#endif
    }
} // namespace sasktran2::successive_orders
