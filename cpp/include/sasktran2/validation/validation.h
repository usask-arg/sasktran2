#pragma once

#include <sasktran2/internal_common.h>

namespace sasktran2::validation {

    inline void throw_configuration_error() {
        throw std::runtime_error(
            "Invalid input. Check log for more information");
    }

    inline void verify_finite(const Eigen::MatrixXd& mat,
                              const std::string& name) {
        if (!mat.allFinite()) {
            spdlog::critical("{} contains non-finite values", name);

            for (int loc = 0; loc < mat.rows(); ++loc) {
                for (int wav = 0; wav < mat.cols(); ++wav) {
                    if (!std::isfinite(mat(loc, wav))) {
                        spdlog::critical("{} contains non-finite values at "
                                         "location: {} and wavelength: {}",
                                         name, loc, wav);
                    }
                }
            }
            throw_configuration_error();
        }
    }

    inline void verify_greater_than(const Eigen::MatrixXd& mat,
                                    const std::string& name, double value) {
        if (mat.minCoeff() < value) {
            spdlog::critical("{} contains values less than {}", name, value);

            for (int loc = 0; loc < mat.rows(); ++loc) {
                for (int wav = 0; wav < mat.cols(); ++wav) {
                    if (mat(loc, wav) < value) {
                        spdlog::critical("{} contains values less than {} at "
                                         "location: {} and wavelength: {}",
                                         name, value, loc, wav);
                    }
                }
            }
            throw_configuration_error();
        }
    }

    inline void verify_less_than(const Eigen::MatrixXd& mat,
                                 const std::string& name, double value) {
        if (mat.maxCoeff() > value) {
            spdlog::critical("{} contains values greater than {}", name, value);

            for (int loc = 0; loc < mat.rows(); ++loc) {
                for (int wav = 0; wav < mat.cols(); ++wav) {
                    if (mat(loc, wav) > value) {
                        spdlog::critical("{} contains values greater than {} "
                                         "at location: {} and wavelength: {}",
                                         name, value, loc, wav);
                    }
                }
            }
            throw_configuration_error();
        }
    }

} // namespace sasktran2::validation
