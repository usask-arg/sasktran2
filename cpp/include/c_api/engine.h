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
                                 OutputC* output);

void sk_engine_destroy(Engine* engine);

#ifdef __cplusplus
}
#endif
