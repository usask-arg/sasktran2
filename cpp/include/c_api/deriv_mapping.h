#pragma once

#ifdef __cplusplus
extern "C" {
#endif

typedef struct DerivativeMapping DerivativeMapping;
typedef struct SurfaceDerivativeMapping SurfaceDerivativeMapping;

int sk_deriv_mapping_destroy(DerivativeMapping* mapping);
int sk_deriv_mapping_set_zero(DerivativeMapping* mapping);
int sk_deriv_mapping_get_d_ssa(DerivativeMapping* mapping, double** ssa);
int sk_deriv_mapping_get_d_extinction(DerivativeMapping* mapping,
                                      double** extinction);
int sk_deriv_mapping_get_scat_factor(DerivativeMapping* mapping,
                                     double** scat_factor);
int sk_deriv_mapping_get_d_legendre(DerivativeMapping* mapping,
                                    double** d_legendre);
int sk_deriv_mapping_get_d_emission(DerivativeMapping* mapping,
                                    double** d_emission);
int sk_deriv_mapping_get_scat_deriv_index(DerivativeMapping* mapping,
                                          int* scat_deriv_index);
int sk_deriv_mapping_set_scat_deriv_index(DerivativeMapping* mapping,
                                          int scat_deriv_index);

int sk_deriv_mapping_get_num_location(DerivativeMapping* mapping,
                                      int* num_location);
int sk_deriv_mapping_get_num_wavel(DerivativeMapping* mapping, int* num_wavel);
int sk_deriv_mapping_get_num_legendre(DerivativeMapping* mapping,
                                      int* num_legendre);

int sk_deriv_mapping_set_interp_dim(DerivativeMapping* mapping,
                                    const char* name);
int sk_deriv_mapping_set_assign_name(DerivativeMapping* mapping,
                                     const char* name);
int sk_deriv_mapping_set_log_radiance_space(DerivativeMapping* mapping,
                                            int log_radiance_space);
int sk_deriv_mapping_is_scattering_derivative(DerivativeMapping* mapping,
                                              int* is_scattering_derivative);

int sk_deriv_mapping_get_num_output(DerivativeMapping* mapping,
                                    int* num_output);

int sk_deriv_mapping_get_assign_name(DerivativeMapping* mapping,
                                     const char** name);
int sk_deriv_mapping_get_interp_dim(DerivativeMapping* mapping,
                                    const char** name);
int sk_deriv_mapping_set_interpolator(DerivativeMapping* mapping,
                                      double* interpolator, int dim1, int dim2);
int sk_deriv_mapping_get_interpolator(DerivativeMapping* mapping,
                                      double** interpolator, int* dim1,
                                      int* dim2);

// Surface mapping methods
int sk_surface_deriv_mapping_get_d_emission(SurfaceDerivativeMapping* mapping,
                                            double** emission);

int sk_surface_deriv_mapping_get_num_wavel(SurfaceDerivativeMapping* mapping,
                                           int* num_wavel);

int sk_surface_deriv_mapping_get_num_brdf_args(
    SurfaceDerivativeMapping* mapping, int* num_brdf_args);

int sk_surface_deriv_mapping_get_d_brdf(SurfaceDerivativeMapping* mapping,
                                        double** brdf);

int sk_surface_deriv_mapping_get_interpolator(SurfaceDerivativeMapping* mapping,
                                              double** interpolator, int* dim1,
                                              int* dim2);

int sk_surface_deriv_mapping_set_interpolator(SurfaceDerivativeMapping* mapping,
                                              double* interpolator, int dim1,
                                              int dim2);

int sk_surface_deriv_mapping_get_interp_dim(SurfaceDerivativeMapping* mapping,
                                            const char** name);

int sk_surface_deriv_mapping_set_interp_dim(SurfaceDerivativeMapping* mapping,
                                            const char* name);

int sk_surface_deriv_mapping_set_zero(SurfaceDerivativeMapping* mapping);

int sk_surface_deriv_mapping_destroy(SurfaceDerivativeMapping* mapping);

#ifdef __cplusplus
}
#endif
