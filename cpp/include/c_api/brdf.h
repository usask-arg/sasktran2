#pragma once

#ifdef __cplusplus
extern "C" {
#endif

typedef struct BRDF BRDF;

BRDF* sk_brdf_create_lambertian(int nstokes);
BRDF* sk_brdf_create_kokhanovsky(int nstokes);
BRDF* sk_brdf_create_modis(int nstokes);

int sk_brdf_get_num_deriv(BRDF* config, int* num_deriv);
int sk_brdf_get_num_args(BRDF* config, int* num_args);

void sk_brdf_destroy(BRDF* config);

#ifdef __cplusplus
}
#endif
