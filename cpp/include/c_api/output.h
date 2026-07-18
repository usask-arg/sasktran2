#pragma once

#ifdef __cplusplus
extern "C" {
#endif

typedef struct OutputC OutputC;
typedef struct OutputJVP OutputJVP;
typedef struct OutputVJP OutputVJP;

OutputC* sk_output_create(double* radiance, int nrad, int nstokes, double* flux,
                          int nflux);
void sk_output_destroy(OutputC* config);

int sk_output_assign_derivative_memory(OutputC* output, const char* name,
                                       double* derivative_mapping, int nrad,
                                       int nstokes, int nderiv);

int sk_output_assign_surface_derivative_memory(OutputC* output,
                                               const char* name,
                                               double* derivative_mapping,
                                               int nrad, int nstokes);

int sk_output_assign_flux_derivative_memory(OutputC* output, const char* name,
                                            double* derivative_mapping,
                                            int nrad, int nderiv);

int sk_output_assign_surface_flux_derivative_memory(OutputC* output,
                                                    const char* name,
                                                    double* derivative_mapping,
                                                    int nrad);

int sk_output_get_los_optical_depth(OutputC* output, double** od);

OutputJVP* sk_output_jvp_create(double* radiance, double* jvp, int nrad,
                                int nstokes);
void sk_output_jvp_destroy(OutputJVP* output);
int sk_output_jvp_assign_derivative_tangent(OutputJVP* output, const char* name,
                                            const double* tangent, int nparam);
int sk_output_jvp_assign_surface_tangent(OutputJVP* output, const char* name,
                                         const double* tangent, int nparam);

OutputVJP* sk_output_vjp_create(double* radiance, const double* cotangent,
                                int nrad, int nstokes);
void sk_output_vjp_destroy(OutputVJP* output);
int sk_output_vjp_assign_derivative_gradient(OutputVJP* output,
                                             const char* name, double* gradient,
                                             int nparam);
int sk_output_vjp_assign_surface_gradient(OutputVJP* output, const char* name,
                                          double* gradient, int nparam);
int sk_output_vjp_finalize(OutputVJP* output);

#ifdef __cplusplus
}
#endif
