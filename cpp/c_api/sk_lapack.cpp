#include "sk_lapack.h"

#include <limits>
#include <vector>

#include <sasktran2/internal_common.h>

namespace {
int convert_lapack_int(int64_t in, lapack_int& out) {
    constexpr int64_t min_v = static_cast<int64_t>(std::numeric_limits<lapack_int>::min());
    constexpr int64_t max_v = static_cast<int64_t>(std::numeric_limits<lapack_int>::max());

    if (in < min_v || in > max_v) {
        return -1;
    }

    out = static_cast<lapack_int>(in);
    return 0;
}
} // namespace

extern "C" {

int64_t sk_lapack_dgesv(int64_t n, int64_t nrhs, double* a, int64_t lda,
                        int64_t* ipiv, double* b, int64_t ldb) {
    if (a == nullptr || ipiv == nullptr || b == nullptr) {
        return -2;
    }

    lapack_int n_lapack;
    lapack_int nrhs_lapack;
    lapack_int lda_lapack;
    lapack_int ldb_lapack;

    if (convert_lapack_int(n, n_lapack) != 0 ||
        convert_lapack_int(nrhs, nrhs_lapack) != 0 ||
        convert_lapack_int(lda, lda_lapack) != 0 ||
        convert_lapack_int(ldb, ldb_lapack) != 0) {
        return -3;
    }

    std::vector<lapack_int> ipiv_local(static_cast<size_t>(n_lapack));
    lapack_int info = 0;

        dgesv_(&n_lapack, &nrhs_lapack, a, &lda_lapack, ipiv_local.data(), b,
            &ldb_lapack, &info);

    for (lapack_int i = 0; i < n_lapack; ++i) {
        ipiv[i] = static_cast<int64_t>(ipiv_local[static_cast<size_t>(i)]);
    }

    return static_cast<int64_t>(info);
}

} // extern "C"
