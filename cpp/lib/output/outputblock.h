#pragma once

#include <sasktran2/output.h>

namespace sasktran2::detail {
    template <typename Output>
    auto strided_output_vector(Output& output, int first, int count,
                               int stride) {
        using Scalar = typename std::decay_t<Output>::Scalar;
        using StridedMap =
            Eigen::Map<Eigen::Matrix<Scalar, Eigen::Dynamic, 1>,
                       Eigen::Unaligned, Eigen::InnerStride<Eigen::Dynamic>>;
        return StridedMap(output.data() + first, count,
                          Eigen::InnerStride<Eigen::Dynamic>(stride));
    }

    template <typename Output>
    auto strided_output_matrix(Output& output, int first, int rows, int columns,
                               int inner_stride) {
        using Scalar = typename std::decay_t<Output>::Scalar;
        using StridedMap =
            Eigen::Map<Eigen::Matrix<Scalar, Eigen::Dynamic, Eigen::Dynamic>,
                       Eigen::Unaligned,
                       Eigen::Stride<Eigen::Dynamic, Eigen::Dynamic>>;
        return StridedMap(output.data() + first, rows, columns,
                          Eigen::Stride<Eigen::Dynamic, Eigen::Dynamic>(
                              output.outerStride(), inner_stride));
    }

    template <int NSTOKES, typename RadianceVector,
              typename DerivativeCollection,
              typename SurfaceDerivativeCollection>
    void assign_mapped_block(
        const sasktran2::WavelengthBlock<>& block,
        const sasktran2::WavelengthBlockDual<NSTOKES>& radiance, int losidx,
        int nlos, int ngeometry, double stokes_c, double stokes_s,
        bool include_emission_derivatives,
        const sasktran2::atmosphere::Atmosphere<NSTOKES>& atmosphere,
        RadianceVector& output, DerivativeCollection& derivatives,
        SurfaceDerivativeCollection& surface_derivatives,
        Eigen::MatrixXd& native_storage, Eigen::MatrixXd& mapped_storage) {
        const int count = block.count;
        const int capacity = native_storage.rows() / NSTOKES;
        const int output_stride = NSTOKES * nlos;
        const int first_output = output_stride * block.start + NSTOKES * losidx;

        const auto output_value = [&](int stokes) {
            return strided_output_vector(output, first_output + stokes, count,
                                         output_stride);
        };
        output_value(0) = radiance.value.row(0).head(count).transpose();
        if constexpr (NSTOKES >= 3) {
            output_value(1) = (stokes_c * radiance.value.row(1).head(count) -
                               stokes_s * radiance.value.row(2).head(count))
                                  .transpose();
            output_value(2) = (stokes_s * radiance.value.row(1).head(count) +
                               stokes_c * radiance.value.row(2).head(count))
                                  .transpose();
        }
        if constexpr (NSTOKES == 4) {
            output_value(3) = radiance.value.row(3).head(count).transpose();
        }

        const auto mapping_columns = Eigen::seqN(block.start, count);
        for (auto& [name, derivative_output] : derivatives) {
            const auto& mapping =
                atmosphere.storage().derivative_mappings_const().at(name);
            const auto& d_ssa = mapping.native_mapping().d_ssa.value();
            const auto& d_extinction =
                mapping.native_mapping().d_extinction.value();

            for (int stokes = 0; stokes < NSTOKES; ++stokes) {
                auto native =
                    native_storage.middleRows(stokes * capacity, count);
                const auto radiance_ssa = radiance.derivative_stokes(
                    ngeometry, ngeometry, stokes, count);
                const auto radiance_extinction =
                    radiance.derivative_stokes(0, ngeometry, stokes, count);
                native.array() =
                    d_ssa(Eigen::placeholders::all, mapping_columns)
                            .transpose()
                            .array() *
                        radiance_ssa.transpose().array() +
                    d_extinction(Eigen::placeholders::all, mapping_columns)
                            .transpose()
                            .array() *
                        radiance_extinction.transpose().array();

                if (mapping.is_scattering_derivative()) {
                    const auto& scatter_factor =
                        mapping.native_mapping().scat_factor.value();
                    const auto radiance_scatter = radiance.derivative_stokes(
                        (2 + mapping.get_scattering_index()) * ngeometry,
                        ngeometry, stokes, count);
                    native.array() += scatter_factor(Eigen::placeholders::all,
                                                     mapping_columns)
                                          .transpose()
                                          .array() *
                                      radiance_scatter.transpose().array();
                }

                if (stokes == 0 &&
                    mapping.native_mapping().d_emission.has_value() &&
                    include_emission_derivatives) {
                    const auto& d_emission =
                        mapping.native_mapping().d_emission.value();
                    const auto radiance_emission = radiance.derivative_stokes(
                        (2 + atmosphere.num_scattering_deriv_groups()) *
                            ngeometry,
                        ngeometry, 0, count);
                    native.array() +=
                        d_emission(Eigen::placeholders::all, mapping_columns)
                            .transpose()
                            .array() *
                        radiance_emission.transpose().array();
                }
            }

            const int num_output = derivative_output.cols();
            auto& derivative_result = derivative_output;
            const auto assign_stokes = [&](int stokes, const auto& values) {
                strided_output_matrix(derivative_result, first_output + stokes,
                                      count, num_output, output_stride) =
                    values;
            };
            if (mapping.get_interpolator_const().has_value()) {
                const auto& interpolator =
                    mapping.get_interpolator_const().value();
                for (int stokes = 0; stokes < NSTOKES; ++stokes) {
                    auto mapped =
                        mapped_storage.middleRows(stokes * capacity, count)
                            .leftCols(num_output);
                    mapped.noalias() =
                        native_storage.middleRows(stokes * capacity, count) *
                        interpolator;
                }
                assign_stokes(0, mapped_storage.middleRows(0, count).leftCols(
                                     num_output));
                if constexpr (NSTOKES >= 3) {
                    assign_stokes(
                        1, stokes_c * mapped_storage.middleRows(capacity, count)
                                          .leftCols(num_output) -
                               stokes_s * mapped_storage
                                              .middleRows(2 * capacity, count)
                                              .leftCols(num_output));
                    assign_stokes(
                        2, stokes_s * mapped_storage.middleRows(capacity, count)
                                          .leftCols(num_output) +
                               stokes_c * mapped_storage
                                              .middleRows(2 * capacity, count)
                                              .leftCols(num_output));
                }
            } else {
                assign_stokes(0, native_storage.middleRows(0, count).leftCols(
                                     num_output));
                if constexpr (NSTOKES >= 3) {
                    assign_stokes(
                        1, stokes_c * native_storage.middleRows(capacity, count)
                                          .leftCols(num_output) -
                               stokes_s * native_storage
                                              .middleRows(2 * capacity, count)
                                              .leftCols(num_output));
                    assign_stokes(
                        2, stokes_s * native_storage.middleRows(capacity, count)
                                          .leftCols(num_output) +
                               stokes_c * native_storage
                                              .middleRows(2 * capacity, count)
                                              .leftCols(num_output));
                }
            }

            if (mapping.log_radiance_space()) {
                const auto radiance_i =
                    radiance.value.row(0).head(count).transpose();
                for (int stokes = 0; stokes < NSTOKES; ++stokes) {
                    auto result = strided_output_matrix(
                        derivative_output, first_output + stokes, count,
                        num_output, output_stride);
                    result.array().colwise() /= radiance_i.array();
                }
            }
        }

        const int surface_derivative_count = atmosphere.surface().num_deriv();
        const int surface_derivative_start =
            (2 + atmosphere.num_scattering_deriv_groups()) * ngeometry;
        for (auto& [name, derivative_output] : surface_derivatives) {
            const auto& mapping =
                atmosphere.surface().derivative_mappings().at(name);
            if (mapping.native_surface_mapping().d_brdf.has_value()) {
                const auto brdf_rows =
                    mapping.native_surface_mapping().d_brdf.value().middleRows(
                        block.start, count);
                for (int stokes = 0; stokes < NSTOKES; ++stokes) {
                    const auto radiance_surface = radiance.derivative_stokes(
                        surface_derivative_start, surface_derivative_count,
                        stokes, count);
                    auto mapped =
                        mapped_storage.middleRows(stokes * capacity, count)
                            .leftCols(1);
                    mapped.col(0) = (radiance_surface.transpose().array() *
                                     brdf_rows.array())
                                        .rowwise()
                                        .sum();
                }
                strided_output_matrix(derivative_output, first_output, count, 1,
                                      output_stride) =
                    mapped_storage.middleRows(0, count).col(0);
                if constexpr (NSTOKES >= 3) {
                    strided_output_matrix(derivative_output, first_output + 1,
                                          count, 1, output_stride) =
                        stokes_c *
                            mapped_storage.middleRows(capacity, count).col(0) -
                        stokes_s *
                            mapped_storage.middleRows(2 * capacity, count)
                                .col(0);
                    strided_output_matrix(derivative_output, first_output + 2,
                                          count, 1, output_stride) =
                        stokes_s *
                            mapped_storage.middleRows(capacity, count).col(0) +
                        stokes_c *
                            mapped_storage.middleRows(2 * capacity, count)
                                .col(0);
                }
            }

            if (mapping.native_surface_mapping().d_emission.has_value() &&
                include_emission_derivatives) {
                const auto radiance_emission = radiance.derivative_stokes(
                    atmosphere.surface_emission_deriv_start_index(), 1, 0,
                    count);
                strided_output_matrix(derivative_output, first_output, count, 1,
                                      output_stride)
                    .array() = radiance_emission.transpose().array().col(0) *
                               mapping.native_surface_mapping()
                                   .d_emission.value()
                                   .middleRows(block.start, count)
                                   .col(0)
                                   .array();
            }
        }
    }
} // namespace sasktran2::detail
