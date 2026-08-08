#pragma once

#include <sasktran2/config.h>

namespace sasktran2::detail {

    /** Resolves machine-dependent implementation choices once, while an
     * engine-owned Config snapshot is being constructed. */
    class RuntimeBackendTuner {
      public:
        /** Resolve all applicable runtime backends for a one-dimensional
         * engine with the supplied number of atmospheric layers. */
        static void resolve(Config& config, int num_layers);

        /** Apply the DISCO timing policy. LAPACK wins ties and cases where the
         * unblocked implementation is less than 10% faster. */
        static BandedLUBackend
        select_banded_lu_backend(double lapack_nanoseconds,
                                 double unblocked_nanoseconds);

        /** Inspect a resolved value in tests and diagnostic tooling without
         * exposing mutable backend policy on Config. */
        static BandedLUBackend resolved_banded_lu_backend(const Config& config);
    };

} // namespace sasktran2::detail
