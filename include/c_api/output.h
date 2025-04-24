#pragma once

#ifdef __cplusplus
extern "C" {
#endif

typedef struct OutputC OutputC;

OutputC* sk_output_create(double* radiance, int nrad, int nstokes);
void sk_output_destroy(OutputC* config);

int sk_output_assign_derivative_memory(
    OutputC* output, const char* name, double* derivative_mapping, int nrad, int nstokes, int nderiv);

#ifdef __cplusplus
}
#endif
