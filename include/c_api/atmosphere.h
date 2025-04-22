#pragma once

#ifdef __cplusplus
extern "C" {
#endif
typedef struct AtmosphereStorage AtmosphereStorage;
typedef struct Atmosphere Atmosphere;
typedef struct Surface Surface;

// STORAGE METHODS
AtmosphereStorage*
sk_atmosphere_storage_create(int nlocation, int nwavel, int nphase_moments,
                             int nstokes, int nderiv, double* ssa,
                             double* total_extinction, double* emission_source,
                             double* f, double* leg_coeff, double* d_leg_coeff,
                             double* d_f, double* solar_irradiance);

void sk_atmosphere_storage_destroy(AtmosphereStorage* storage);

// ATMOSPHERE METHODS
Atmosphere* sk_atmosphere_create(AtmosphereStorage* storage, Surface* surface,
                                 int calculate_derivatives);

void sk_atmosphere_destroy(Atmosphere* atmosphere);

// SURFACE METHODS
Surface* sk_surface_create(int nwavel, int nstokes);

void sk_surface_destroy(Surface* config);

#ifdef __cplusplus
}
#endif
