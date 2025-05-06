#pragma once

#ifdef __cplusplus
extern "C" {
#endif

typedef struct Config Config;

Config* sk_config_create();

void sk_config_destroy(Config* config);

int sk_config_get_num_stokes(Config* config, int* num_stokes);
int sk_config_set_num_stokes(Config* config, int num_stokes);

int sk_config_get_multiple_scatter_source(Config* config,
                                          int* multiple_scatter_source);
int sk_config_set_multiple_scatter_source(Config* config,
                                          int multiple_scatter_source);

int sk_config_get_single_scatter_source(Config* config,
                                        int* single_scatter_source);
int sk_config_set_single_scatter_source(Config* config,
                                        int single_scatter_source);

int sk_config_get_num_streams(Config* config, int* num_streams);
int sk_config_set_num_streams(Config* config, int num_streams);

// Thread-related methods
int sk_config_get_num_threads(Config* config, int* num_threads);
int sk_config_set_num_threads(Config* config, int num_threads);
int sk_config_get_threading_model(Config* config, int* threading_model);
int sk_config_set_threading_model(Config* config, int threading_model);

// Single scatter configuration methods
int sk_config_get_num_singlescatter_moments(Config* config, int* num_moments);
int sk_config_set_num_singlescatter_moments(Config* config, int num_moments);
int sk_config_get_singlescatter_phasemode(Config* config, int* phasemode);
int sk_config_set_singlescatter_phasemode(Config* config, int phasemode);

// Delta scaling configuration
int sk_config_get_apply_delta_scaling(Config* config, int* apply_delta_scaling);
int sk_config_set_apply_delta_scaling(Config* config, int apply_delta_scaling);

// DO configuration methods
int sk_config_get_num_do_sza(Config* config, int* num_sza);
int sk_config_set_num_do_sza(Config* config, int num_sza);
int sk_config_get_num_do_forced_azimuth(Config* config,
                                        int* num_forced_azimuth);
int sk_config_set_num_do_forced_azimuth(Config* config, int num_forced_azimuth);
int sk_config_get_do_backprop(Config* config, int* do_backprop);
int sk_config_set_do_backprop(Config* config, int do_backprop);
int sk_config_get_num_do_spherical_iterations(Config* config,
                                              int* num_iterations);
int sk_config_set_num_do_spherical_iterations(Config* config,
                                              int num_iterations);

// HR configuration methods
int sk_config_get_num_hr_spherical_iterations(Config* config,
                                              int* num_iterations);
int sk_config_set_num_hr_spherical_iterations(Config* config,
                                              int num_iterations);
int sk_config_get_num_hr_incoming(Config* config, int* num_incoming);
int sk_config_set_num_hr_incoming(Config* config, int num_incoming);
int sk_config_get_num_hr_outgoing(Config* config, int* num_outgoing);
int sk_config_set_num_hr_outgoing(Config* config, int num_outgoing);
int sk_config_get_num_hr_full_incoming_points(Config* config, int* num_points);
int sk_config_set_num_hr_full_incoming_points(Config* config, int num_points);
int sk_config_get_initialize_hr_with_do(Config* config, int* initialize);
int sk_config_set_initialize_hr_with_do(Config* config, int initialize);

// Source configuration methods
int sk_config_get_occultation_source(Config* config, int* occultation_source);
int sk_config_set_occultation_source(Config* config, int occultation_source);
int sk_config_get_emission_source(Config* config, int* emission_source);
int sk_config_set_emission_source(Config* config, int emission_source);

// Refraction configuration methods
int sk_config_get_los_refraction(Config* config, int* refraction);
int sk_config_set_los_refraction(Config* config, int refraction);
int sk_config_get_solar_refraction(Config* config, int* refraction);
int sk_config_set_solar_refraction(Config* config, int refraction);
int sk_config_get_multiple_scatter_refraction(Config* config, int* refraction);
int sk_config_set_multiple_scatter_refraction(Config* config, int refraction);

// Weighting function configuration methods
int sk_config_get_wf_enabled(Config* config, int* enabled);
int sk_config_set_wf_enabled(Config* config, int enabled);
int sk_config_get_wf_precision(Config* config, int* precision);
int sk_config_set_wf_precision(Config* config, int precision);

// Stokes configuration methods
int sk_config_get_stokes_basis(Config* config, int* basis);
int sk_config_set_stokes_basis(Config* config, int basis);

// Input validation configuration
int sk_config_get_input_validation_mode(Config* config, int* mode);
int sk_config_set_input_validation_mode(Config* config, int mode);

// Diagnostic output configuration
int sk_config_get_output_los_optical_depth(Config* config, int* output);
int sk_config_set_output_los_optical_depth(Config* config, int output);

#ifdef __cplusplus
}
#endif
