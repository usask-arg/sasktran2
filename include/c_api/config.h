#pragma once

#ifdef __cplusplus
extern "C" {
#endif

typedef struct Config Config;

Config* sk_config_create();

void sk_config_destroy(Config* config);

int sk_config_get_num_stokes(Config* config, int* num_stokes);

#ifdef __cplusplus
}
#endif
