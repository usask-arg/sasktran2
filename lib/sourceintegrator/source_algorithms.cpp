
#include "sasktran2/dual.h"
#include <sasktran2/source_algorithms.h>

namespace sasktran2::sourcealgo {
    template <int NSTOKES, typename SourceType>
    void add_integrated_constant_source(
        const sasktran2::SparseODDualView& shell_od,
        const SourceType& start_source, const SourceType& end_source,
        double start_weight, double end_weight,
        sasktran2::Dual<double, sasktran2::dualstorage::dense, NSTOKES>&
            result_source) {

        double source_factor1 = (1 - shell_od.exp_minus_od);
        // Note dsource_factor = d_od * (1 - source_factor)

        result_source.value.array() +=
            source_factor1 * (start_source.value.array() * start_weight +
                              end_source.value.array() * end_weight);
    }

    template <int NSTOKES, typename SourceType>
    void add_integrated_exponential_source(
        const sasktran2::SparseODDualView& shell_od,
        const SourceType& start_source, const SourceType& end_source,
        sasktran2::Dual<double, sasktran2::dualstorage::dense, NSTOKES>&
            result_source) {
        // Here we assume that the source varies exponentially in optical depth
        // across the layer.  I.e. J(x) = A exp(b x)

        // We can calculate A, and B. At the start of the layer, x=0, so A =
        // start_source

        // end_source / start_source = exp(b * OD)
        // b = log(end_source / start_source) / OD

        // Since we are trying to integrate J(x) exp(-x) from x=0 to OD, this
        // integral becomes A / (1 - b) * (1 - exp(-(1-b) * OD))

        // If end_source = start_source the source is constant and this is the
        // standard form source * (1 - exp(-OD))

        // Also, exp(-b) = exp(-OD) - end_source / start_source
        // Which avoids one exponentiation?

        auto A = start_source.value;
        auto b = log(end_source.value.array() / start_source.value.array()) /
                 shell_od.od;

        result_source.value.array() +=
            A.array() / (1 - b.array()) *
            (1 - exp(-(1 - b.array()) * shell_od.od));

        // db =
    }

    template void add_integrated_exponential_source<
        1, sasktran2::Dual<double, sasktran2::dualstorage::dense, 1>>(
        const sasktran2::SparseODDualView& shell_od,
        const sasktran2::Dual<double, sasktran2::dualstorage::dense, 1>&
            start_source,
        const sasktran2::Dual<double, sasktran2::dualstorage::dense, 1>&
            end_source,
        sasktran2::Dual<double, sasktran2::dualstorage::dense, 1>&
            result_source);

    template void add_integrated_exponential_source<
        3, sasktran2::Dual<double, sasktran2::dualstorage::dense, 3>>(
        const sasktran2::SparseODDualView& shell_od,
        const sasktran2::Dual<double, sasktran2::dualstorage::dense, 3>&
            start_source,
        const sasktran2::Dual<double, sasktran2::dualstorage::dense, 3>&
            end_source,
        sasktran2::Dual<double, sasktran2::dualstorage::dense, 3>&
            result_source);

    template void add_integrated_constant_source<
        1, sasktran2::Dual<double, sasktran2::dualstorage::dense, 1>>(
        const sasktran2::SparseODDualView& shell_od,
        const sasktran2::Dual<double, sasktran2::dualstorage::dense, 1>&
            start_source,
        const sasktran2::Dual<double, sasktran2::dualstorage::dense, 1>&
            end_source,
        double start_weight, double end_weight,
        sasktran2::Dual<double, sasktran2::dualstorage::dense, 1>&
            result_source);

    template void add_integrated_constant_source<
        3, sasktran2::Dual<double, sasktran2::dualstorage::dense, 3>>(
        const sasktran2::SparseODDualView& shell_od,
        const sasktran2::Dual<double, sasktran2::dualstorage::dense, 3>&
            start_source,
        const sasktran2::Dual<double, sasktran2::dualstorage::dense, 3>&
            end_source,
        double start_weight, double end_weight,
        sasktran2::Dual<double, sasktran2::dualstorage::dense, 3>&
            result_source);

} // namespace sasktran2::sourcealgo
