#pragma once

#include <sasktran2/internal_common.h>
#include <sasktran2/dual.h>

namespace sasktran2::sourcealgo {

    template <int NSTOKES, typename SourceType>
    void add_integrated_constant_source(
        const sasktran2::SparseODDualView& shell_od,
        const SourceType& start_source, const SourceType& end_source,
        double start_weight, double end_weight,
        sasktran2::Dual<double, sasktran2::dualstorage::dense, NSTOKES>&
            result_source);

    template <int NSTOKES, typename SourceType>
    void add_integrated_exponential_source(
        const sasktran2::SparseODDualView& shell_od,
        const SourceType& start_source, const SourceType& end_source,
        sasktran2::Dual<double, sasktran2::dualstorage::dense, NSTOKES>&
            result_source);

} // namespace sasktran2::sourcealgo
