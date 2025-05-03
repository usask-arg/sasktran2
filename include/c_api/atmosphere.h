/**
 * @file atmosphere.h
 * @brief Functions and types for managing atmospheric storage in radiative
 * transfer.
 *
 * This header defines the C API for creating, destroying, and accessing data
 * in the atmosphere storage structure used for radiative transfer calculations.
 */

#pragma once

#include "brdf.h"
#include "deriv_mapping.h"

#ifdef __cplusplus
extern "C" {
#endif

/**
 * @defgroup AtmosphereStorage Atmosphere Storage
 * @brief Creation and manipulation of atmospheric storage containers.
 * @{
 */

/**
 * @brief Opaque handle to atmospheric storage.
 *
 * An AtmosphereStorage object contains all of the optical data necessary to
 * perform the radiative transfer calculation at the given grid points.
 * Typically this is things like single scatter albedo, extinction, legendre
 * coefficients at a set of geometry grid points.
 */
typedef struct AtmosphereStorage AtmosphereStorage;

/**
 * @brief Opaque handle to a derivative mapping object.
 *
 * A DerivativeMapping object is the input quantities necessary to map
 * derivatives calculated by the model internally (with respect to SSA,
 * total_extinction, scattering quantities) to user defined quantities.
 */
typedef struct DerivativeMapping DerivativeMapping;

typedef struct SurfaceDerivativeMapping SurfaceDerivativeMapping;

/**
 * @brief Opaque handle to an atmosphere object.
 *
 * An atmosphere object is all of the geometry independent information required
 * to perform the radiative transfer calculation.  This is a combination of an
 * AtmosphereStorage object that contains layer/grid quantities, and a Surface
 * object.
 */
typedef struct Atmosphere Atmosphere;

/**
 * @brief Opaque handle to a surface model.
 *
 * The surface object contains information of the Earth's surface, including the
 * BRDF and emission properties.
 */
typedef struct Surface Surface;

typedef struct BRDF BRDF;

// ---------------------
// STORAGE METHODS
// ---------------------

/**
 * @brief Creates a new AtmosphereStorage object.
 *
 * @param nlocation Number of spatial locations.
 * @param nwavel Number of wavelengths.
 * @param nphase_moments Number of phase function moments.
 * @param nstokes Number of Stokes parameters (usually 1 or 3).
 * @param ssa Pointer to single scattering albedo data. [nlocation, nwavel]
 * @param total_extinction Pointer to total extinction coefficient data.
 * [nlocation, nwavel]
 * @param emission_source Pointer to emission source data. [nlocation, nwavel]
 * @param leg_coeff Pointer to Legendre coefficients. [nphase_moments,
 * nlocation, nwavel]
 * @param solar_irradiance Pointer to top-of-atmosphere solar irradiance.
 * @return Pointer to a new AtmosphereStorage object.
 */
AtmosphereStorage*
sk_atmosphere_storage_create(int nlocation, int nwavel, int nphase_moments,
                             int nstokes, double* ssa, double* total_extinction,
                             double* emission_source, double* leg_coeff,
                             double* solar_irradiance);

/**
 * @brief Destroys a previously created AtmosphereStorage object.
 *
 * @param storage The storage object to destroy.
 */
void sk_atmosphere_storage_destroy(AtmosphereStorage* storage);

/**
 * @brief Gets a pointer to a derivative mapping for a named quantity.
 *
 * @param storage The atmosphere storage to query.
 * @param name The name of the quantity.
 * @param mapping Output pointer to the derivative mapping object.
 * @return 0 on success, non-zero on error.
 */
int sk_atmosphere_storage_get_derivative_mapping(AtmosphereStorage* storage,
                                                 const char* name,
                                                 DerivativeMapping** mapping);

int sk_atmosphere_storage_get_derivative_mapping_by_index(
    AtmosphereStorage* storage, int index, DerivativeMapping** mapping);

int sk_atmosphere_storage_get_num_derivative_mappings(
    AtmosphereStorage* storage, int* num_mappings);

int sk_atmosphere_storage_get_derivative_mapping_name(
    AtmosphereStorage* storage, int index, const char** name);

int sk_atmosphere_storage_finalize_scattering_derivatives(
    AtmosphereStorage* storage);

int sk_atmosphere_storage_set_zero(AtmosphereStorage* storage);

// ---------------------
// ATMOSPHERE METHODS
// ---------------------

/**
 * @brief Creates a new Atmosphere object.
 *
 * @param storage Pointer to atmosphere storage.
 * @param surface Pointer to surface model.
 * @param calculate_derivatives Whether to compute derivatives (non-zero = yes).
 * @param calculate_emission_derivatives Whether to compute emission
 * derivatives, (non-zero = yes)
 * @return Pointer to a new Atmosphere object.
 */
Atmosphere* sk_atmosphere_create(AtmosphereStorage* storage, Surface* surface,
                                 int calculate_derivatives,
                                 int calculate_emission_derivatives);

/**
 * @brief Destroys a previously created Atmosphere object.
 *
 * @param atmosphere The atmosphere object to destroy.
 */
void sk_atmosphere_destroy(Atmosphere* atmosphere);

int sk_atmosphere_apply_delta_m_scaling(Atmosphere* atmosphere, int order);

// ---------------------
// SURFACE METHODS
// ---------------------

/**
 * @brief Creates a new Surface object.
 *
 * @param nwavel Number of wavelengths.
 * @param nstokes Number of Stokes parameters.
 * @param emission Pointer to emission data. [nwavel]
 * @return Pointer to a new Surface object.
 */
Surface* sk_surface_create(int nwavel, int nstokes, double* emission);

/**
 * @brief Destroys a previously created Surface object.
 *
 * @param surface The surface object to destroy.
 */
void sk_surface_destroy(Surface* surface);

int sk_surface_set_brdf(Surface* surface, BRDF* brdf, double* brdf_args);

int sk_surface_get_derivative_mapping(Surface* storage, const char* name,
                                      SurfaceDerivativeMapping** mapping);

int sk_surface_get_num_derivative_mappings(Surface* storage, int* num_mappings);

int sk_surface_get_derivative_mapping_name(Surface* storage, int index,
                                           const char** name);

int sk_surface_set_zero(Surface* storage);

#ifdef __cplusplus
}
#endif

/** @} */ // end of AtmosphereStorage
