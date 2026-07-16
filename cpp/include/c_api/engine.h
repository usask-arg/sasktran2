#pragma once

#include "output.h"
#ifdef __cplusplus
extern "C" {
#endif
typedef struct Engine Engine;
typedef struct Config Config;
typedef struct Atmosphere Atmosphere;
typedef struct Geometry1D Geometry1D;
typedef struct Geometry2D Geometry2D;
typedef struct ViewingGeometry ViewingGeometry;
typedef struct OutputC OutputC;

Engine* sk_engine_create(Config* engine, Geometry1D* geometry,
                         ViewingGeometry* viewing_geometry);

Engine* sk_engine_create_2d(Config* engine, Geometry2D* geometry,
                            ViewingGeometry* viewing_geometry);

int sk_engine_calculate_radiance(Engine* engine, Atmosphere* atmosphere,
                                 OutputC* output, int only_initialize);
int sk_engine_calculate_radiance_thread(Engine* engine, Atmosphere* atmosphere,
                                        OutputC* output, int wavelength_idx,
                                        int thread_idx);
int sk_engine_effective_wavelength_batch_size(Engine* engine,
                                              int num_wavelengths);
int sk_engine_calculate_radiance_batch_thread(
    Engine* engine, Atmosphere* atmosphere, OutputC* output,
    int wavelength_start, int wavelength_count, int thread_idx);

void sk_engine_destroy(Engine* engine);

int sk_openmp_support_enabled();

#ifdef __cplusplus
}
#endif
