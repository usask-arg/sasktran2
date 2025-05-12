#pragma once

#include "output.h"
#ifdef __cplusplus
extern "C" {
#endif
typedef struct Engine Engine;
typedef struct Config Config;
typedef struct Atmosphere Atmosphere;
typedef struct Geometry1D Geometry1D;
typedef struct ViewingGeometry ViewingGeometry;
typedef struct OutputC OutputC;

Engine* sk_engine_create(Config* engine, Geometry1D* geometry,
                         ViewingGeometry* viewing_geometry);

int sk_engine_calculate_radiance(Engine* engine, Atmosphere* atmosphere,
                                 OutputC* output, int only_initialize);
int sk_engine_calculate_radiance_thread(Engine* engine, Atmosphere* atmosphere,
                                        OutputC* output, int wavelength_idx,
                                        int thread_idx);

void sk_engine_destroy(Engine* engine);

int sk_openmp_support_enabled();

#ifdef __cplusplus
}
#endif
