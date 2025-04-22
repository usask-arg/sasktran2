#pragma once

#ifdef __cplusplus
extern "C" {
#endif

typedef struct OutputC OutputC;

OutputC* sk_output_create(double* radiance, int nrad, int nstokes);

void sk_output_destroy(OutputC* config);

#ifdef __cplusplus
}
#endif
