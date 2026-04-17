#pragma once

#include <stdint.h>

#ifdef __cplusplus
extern "C" {
#endif

/**
 * C API wrapper around the configured LAPACK dgesv routine.
 *
 * All integer parameters are int64_t for FFI stability. Internally these are
 * converted to the build's LAPACK integer type before calling dgesv_.
 *
 * Returns LAPACK's INFO value. Negative values indicate an illegal argument in
 * this wrapper conversion/validation step.
 */
int64_t sk_lapack_dgesv(int64_t n, int64_t nrhs, double* a, int64_t lda,
                        int64_t* ipiv, double* b, int64_t ldb);

#ifdef __cplusplus
}
#endif
