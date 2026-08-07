#pragma once

#include <sasktran2/solar_transmission_operator.h>
#include <cstddef>
#include <cstdint>
#include <memory>
#include <vector>

namespace sasktran2::successive_orders {
    /** Opaque C++ owner for the Rust solar-transmission operator used while
     * constructing the successive-orders first-order forcing.
     *
     * The generic single-scatter source only sees this interface. Rust/CXX
     * types, SIMD scratch, and packed operator storage remain private to the
     * successive-orders module.
     */
    class SolarTransmission final : public sasktran2::solartransmission::
                                        detail::SolarTransmissionOperator {
      public:
        using SparseMatrix = Eigen::SparseMatrix<double, Eigen::RowMajor>;
        using RowMajorMatrix = Eigen::Matrix<double, Eigen::Dynamic,
                                             Eigen::Dynamic, Eigen::RowMajor>;

        static std::shared_ptr<SolarTransmission>
        create_exact(const Eigen::MatrixXd& dense_path,
                     const SparseMatrix& sparse_path,
                     const std::vector<bool>& ground_hit,
                     std::size_t input_size, int wavelength_workers);

        static std::shared_ptr<SolarTransmission>
        create_table(const Eigen::MatrixXd& dense_path,
                     const SparseMatrix& interpolation,
                     const std::vector<bool>& ground_hit,
                     std::size_t input_size, int wavelength_workers);

        ~SolarTransmission() override;
        SolarTransmission(SolarTransmission&&) noexcept;
        SolarTransmission& operator=(SolarTransmission&&) noexcept;
        SolarTransmission(const SolarTransmission&) = delete;
        SolarTransmission& operator=(const SolarTransmission&) = delete;

        std::size_t input_size() const override;
        std::size_t output_size() const override;
        std::size_t storage_bytes() const override;
        std::size_t forward_scratch_bytes() const override;
        std::size_t vjp_scratch_bytes() const override;

        void calculate(int worker, Eigen::Ref<const Eigen::VectorXd> extinction,
                       double solar_irradiance,
                       Eigen::Ref<Eigen::VectorXd> transmission) override;

        void calculate_batch4(int worker, const double* extinction,
                              std::size_t extinction_size,
                              const double* solar_irradiance,
                              std::size_t irradiance_size, double* transmission,
                              std::size_t transmission_size) override;

        void calculate_jvp(
            int worker, Eigen::Ref<const Eigen::VectorXd> extinction_tangent,
            Eigen::Ref<const Eigen::VectorXd> transmission,
            Eigen::Ref<Eigen::VectorXd> transmission_tangent) override;

        void calculate_jvp_batch4(
            int worker, Eigen::Ref<const RowMajorMatrix> extinction_tangent,
            Eigen::Ref<const RowMajorMatrix> transmission,
            Eigen::Ref<RowMajorMatrix> transmission_tangent) override;

        void accumulate_vjp(
            int worker, Eigen::Ref<const Eigen::VectorXd> transmission,
            Eigen::Ref<const Eigen::VectorXd> transmission_cotangent,
            Eigen::Ref<Eigen::VectorXd> extinction_gradient) override;

        void accumulate_vjp_batch4(
            int worker, Eigen::Ref<const RowMajorMatrix> transmission,
            Eigen::Ref<RowMajorMatrix> transmission_cotangent,
            Eigen::Ref<RowMajorMatrix> extinction_gradient) override;

        // Internal construction surface used by the module-local factories.
        // Impl remains incomplete outside the implementation file.
        struct Impl;
        explicit SolarTransmission(std::unique_ptr<Impl> impl);

      private:
        std::unique_ptr<Impl> m_impl;
    };
} // namespace sasktran2::successive_orders
