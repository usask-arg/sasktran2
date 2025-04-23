#pragma once

#ifdef __cplusplus
extern "C" {
#endif

typedef struct Config Config;

Config* sk_config_create();

void sk_config_destroy(Config* config);

int sk_config_get_num_stokes(Config* config, int* num_stokes);
int sk_config_set_num_stokes(Config* config, int num_stokes);

int sk_config_set_multiple_scatter_source(Config* config,
                                          int multiple_scatter_source);
int sk_config_set_single_scatter_source(Config* config,
                                        int single_scatter_source);

int sk_config_set_num_streams(Config* config, int num_streams);

#ifdef __cplusplus
}
#endif
