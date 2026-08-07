#pragma once

#include <Eigen/Core>
#include <Eigen/SparseCore>

#include <cstddef>
#include <functional>
#include <memory>
#include <vector>

namespace sasktran2::solartransmission::detail {
    /** Backend-neutral solar-path operator used by delegated source terms.
     *
     * This narrow interface prevents generic single-scatter code from
     * depending on a concrete Rust source or any CXX bridge types.
     */
    class SolarTransmissionOperator {
      public:
        using SparseMatrix = Eigen::SparseMatrix<double, Eigen::RowMajor>;
        using RowMajorMatrix = Eigen::Matrix<double, Eigen::Dynamic,
                                             Eigen::Dynamic, Eigen::RowMajor>;

        virtual ~SolarTransmissionOperator() = default;

        virtual std::size_t input_size() const = 0;
        virtual std::size_t output_size() const = 0;
        virtual std::size_t storage_bytes() const = 0;
        virtual std::size_t forward_scratch_bytes() const = 0;
        virtual std::size_t vjp_scratch_bytes() const = 0;

        virtual void calculate(int worker,
                               Eigen::Ref<const Eigen::VectorXd> extinction,
                               double solar_irradiance,
                               Eigen::Ref<Eigen::VectorXd> transmission) = 0;
        virtual void calculate_batch4(int worker, const double* extinction,
                                      std::size_t extinction_size,
                                      const double* solar_irradiance,
                                      std::size_t irradiance_size,
                                      double* transmission,
                                      std::size_t transmission_size) = 0;
        virtual void
        calculate_jvp(int worker,
                      Eigen::Ref<const Eigen::VectorXd> extinction_tangent,
                      Eigen::Ref<const Eigen::VectorXd> transmission,
                      Eigen::Ref<Eigen::VectorXd> transmission_tangent) = 0;
        virtual void calculate_jvp_batch4(
            int worker, Eigen::Ref<const RowMajorMatrix> extinction_tangent,
            Eigen::Ref<const RowMajorMatrix> transmission,
            Eigen::Ref<RowMajorMatrix> transmission_tangent) = 0;
        virtual void
        accumulate_vjp(int worker,
                       Eigen::Ref<const Eigen::VectorXd> transmission,
                       Eigen::Ref<const Eigen::VectorXd> transmission_cotangent,
                       Eigen::Ref<Eigen::VectorXd> extinction_gradient) = 0;
        virtual void accumulate_vjp_batch4(
            int worker, Eigen::Ref<const RowMajorMatrix> transmission,
            Eigen::Ref<RowMajorMatrix> transmission_cotangent,
            Eigen::Ref<RowMajorMatrix> extinction_gradient) = 0;
    };

    using SolarTransmissionOperatorFactory =
        std::function<std::shared_ptr<SolarTransmissionOperator>(
            const Eigen::MatrixXd&,
            const SolarTransmissionOperator::SparseMatrix&,
            const std::vector<bool>&, std::size_t, int)>;
} // namespace sasktran2::solartransmission::detail
