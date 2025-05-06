#include "config.h"
#include "sasktran2/config.h"
#include "internal_types.h"
#include <sasktran2.h>

extern "C" {
Config* sk_config_create() { return new Config(); }

void sk_config_destroy(Config* storage) { delete storage; }

int sk_config_get_num_stokes(Config* config, int* num_stokes) {
    if (config == nullptr || num_stokes == nullptr) {
        return -1; // Error: null pointer
    }
    *num_stokes = static_cast<int>(config->impl.num_stokes());
    return 0; // Success
}

int sk_config_set_num_stokes(Config* config, int num_stokes) {
    if (config == nullptr) {
        return -1; // Error: null pointer
    }
    if (num_stokes != 1 && num_stokes != 3) {
        return -2; // Error: invalid number of Stokes parameters
    }
    config->impl.set_num_stokes(num_stokes);
    return 0; // Success
}

int sk_config_get_multiple_scatter_source(Config* config,
                                          int* multiple_scatter_source) {
    if (config == nullptr || multiple_scatter_source == nullptr) {
        return -1; // Error: null pointer
    }
    *multiple_scatter_source =
        static_cast<int>(config->impl.multiple_scatter_source());
    return 0; // Success
}

int sk_config_set_multiple_scatter_source(Config* config,
                                          int multiple_scatter_source) {
    if (config == nullptr) {
        return -1; // Error: null pointer
    }
    config->impl.set_multiple_scatter_source(
        static_cast<sasktran2::Config::MultipleScatterSource>(
            multiple_scatter_source));
    return 0; // Success
}

int sk_config_get_single_scatter_source(Config* config,
                                        int* single_scatter_source) {
    if (config == nullptr || single_scatter_source == nullptr) {
        return -1; // Error: null pointer
    }
    *single_scatter_source =
        static_cast<int>(config->impl.single_scatter_source());
    return 0; // Success
}

int sk_config_set_single_scatter_source(Config* config,
                                        int single_scatter_source) {
    if (config == nullptr) {
        return -1; // Error: null pointer
    }
    config->impl.set_single_scatter_source(
        static_cast<sasktran2::Config::SingleScatterSource>(
            single_scatter_source));
    return 0; // Success
}

int sk_config_get_num_streams(Config* config, int* num_streams) {
    if (config == nullptr || num_streams == nullptr) {
        return -1; // Error: null pointer
    }
    *num_streams = config->impl.num_do_streams();
    return 0; // Success
}

int sk_config_set_num_streams(Config* config, int num_streams) {
    if (config == nullptr) {
        return -1; // Error: null pointer
    }
    config->impl.set_num_do_streams(num_streams);
    return 0; // Success
}

int sk_config_get_num_threads(Config* config, int* num_threads) {
    if (config == nullptr || num_threads == nullptr) {
        return -1; // Error: null pointer
    }
    *num_threads = config->impl.num_threads();
    return 0; // Success
}

int sk_config_set_num_threads(Config* config, int num_threads) {
    if (config == nullptr) {
        return -1; // Error: null pointer
    }
    config->impl.set_num_threads(num_threads);
    return 0; // Success
}

int sk_config_get_threading_model(Config* config, int* threading_model) {
    if (config == nullptr || threading_model == nullptr) {
        return -1; // Error: null pointer
    }
    *threading_model = static_cast<int>(config->impl.threading_model());
    return 0; // Success
}

int sk_config_set_threading_model(Config* config, int threading_model) {
    if (config == nullptr) {
        return -1; // Error: null pointer
    }
    config->impl.set_threading_model(
        static_cast<sasktran2::Config::ThreadingModel>(threading_model));
    return 0; // Success
}

int sk_config_get_num_singlescatter_moments(Config* config, int* num_moments) {
    if (config == nullptr || num_moments == nullptr) {
        return -1; // Error: null pointer
    }
    *num_moments = config->impl.num_singlescatter_moments();
    return 0; // Success
}

int sk_config_set_num_singlescatter_moments(Config* config, int num_moments) {
    if (config == nullptr) {
        return -1; // Error: null pointer
    }
    config->impl.set_num_singlescatter_moments(num_moments);
    return 0; // Success
}

int sk_config_get_singlescatter_phasemode(Config* config, int* phasemode) {
    if (config == nullptr || phasemode == nullptr) {
        return -1; // Error: null pointer
    }
    *phasemode = static_cast<int>(config->impl.singlescatter_phasemode());
    return 0; // Success
}

int sk_config_set_singlescatter_phasemode(Config* config, int phasemode) {
    if (config == nullptr) {
        return -1; // Error: null pointer
    }
    config->impl.set_singlescatter_phasemode(
        static_cast<sasktran2::Config::SingleScatterPhaseMode>(phasemode));
    return 0; // Success
}

int sk_config_get_apply_delta_scaling(Config* config,
                                      int* apply_delta_scaling) {
    if (config == nullptr || apply_delta_scaling == nullptr) {
        return -1; // Error: null pointer
    }
    *apply_delta_scaling = static_cast<int>(config->impl.apply_delta_scaling());
    return 0; // Success
}

int sk_config_set_apply_delta_scaling(Config* config, int apply_delta_scaling) {
    if (config == nullptr) {
        return -1; // Error: null pointer
    }
    config->impl.set_apply_delta_scaling(apply_delta_scaling != 0);
    return 0; // Success
}

int sk_config_get_num_do_sza(Config* config, int* num_sza) {
    if (config == nullptr || num_sza == nullptr) {
        return -1; // Error: null pointer
    }
    *num_sza = config->impl.num_do_sza();
    return 0; // Success
}

int sk_config_set_num_do_sza(Config* config, int num_sza) {
    if (config == nullptr) {
        return -1; // Error: null pointer
    }
    config->impl.set_num_do_sza(num_sza);
    return 0; // Success
}

int sk_config_get_num_do_forced_azimuth(Config* config,
                                        int* num_forced_azimuth) {
    if (config == nullptr || num_forced_azimuth == nullptr) {
        return -1; // Error: null pointer
    }
    *num_forced_azimuth = config->impl.num_do_forced_azimuth();
    return 0; // Success
}

int sk_config_set_num_do_forced_azimuth(Config* config,
                                        int num_forced_azimuth) {
    if (config == nullptr) {
        return -1; // Error: null pointer
    }
    config->impl.set_num_do_forced_azimuth(num_forced_azimuth);
    return 0; // Success
}

int sk_config_get_do_backprop(Config* config, int* do_backprop) {
    if (config == nullptr || do_backprop == nullptr) {
        return -1; // Error: null pointer
    }
    *do_backprop = static_cast<int>(config->impl.do_backprop());
    return 0; // Success
}

int sk_config_set_do_backprop(Config* config, int do_backprop) {
    if (config == nullptr) {
        return -1; // Error: null pointer
    }
    config->impl.set_do_backprop(do_backprop != 0);
    return 0; // Success
}

int sk_config_get_num_do_spherical_iterations(Config* config,
                                              int* num_iterations) {
    if (config == nullptr || num_iterations == nullptr) {
        return -1; // Error: null pointer
    }
    *num_iterations = config->impl.num_do_spherical_iterations();
    return 0; // Success
}

int sk_config_set_num_do_spherical_iterations(Config* config,
                                              int num_iterations) {
    if (config == nullptr) {
        return -1; // Error: null pointer
    }
    config->impl.set_num_do_spherical_iterations(num_iterations);
    return 0; // Success
}

int sk_config_get_num_hr_spherical_iterations(Config* config,
                                              int* num_iterations) {
    if (config == nullptr || num_iterations == nullptr) {
        return -1; // Error: null pointer
    }
    *num_iterations = config->impl.num_hr_spherical_iterations();
    return 0; // Success
}

int sk_config_set_num_hr_spherical_iterations(Config* config,
                                              int num_iterations) {
    if (config == nullptr) {
        return -1; // Error: null pointer
    }
    config->impl.set_num_hr_spherical_iterations(num_iterations);
    return 0; // Success
}

int sk_config_get_num_hr_incoming(Config* config, int* num_incoming) {
    if (config == nullptr || num_incoming == nullptr) {
        return -1; // Error: null pointer
    }
    *num_incoming = config->impl.num_hr_incoming();
    return 0; // Success
}

int sk_config_set_num_hr_incoming(Config* config, int num_incoming) {
    if (config == nullptr) {
        return -1; // Error: null pointer
    }
    config->impl.set_num_hr_incoming(num_incoming);
    return 0; // Success
}

int sk_config_get_num_hr_outgoing(Config* config, int* num_outgoing) {
    if (config == nullptr || num_outgoing == nullptr) {
        return -1; // Error: null pointer
    }
    *num_outgoing = config->impl.num_hr_outgoing();
    return 0; // Success
}

int sk_config_set_num_hr_outgoing(Config* config, int num_outgoing) {
    if (config == nullptr) {
        return -1; // Error: null pointer
    }
    config->impl.set_num_hr_outgoing(num_outgoing);
    return 0; // Success
}

int sk_config_get_num_hr_full_incoming_points(Config* config, int* num_points) {
    if (config == nullptr || num_points == nullptr) {
        return -1; // Error: null pointer
    }
    *num_points = config->impl.num_hr_full_incoming_points();
    return 0; // Success
}

int sk_config_set_num_hr_full_incoming_points(Config* config, int num_points) {
    if (config == nullptr) {
        return -1; // Error: null pointer
    }
    config->impl.set_num_hr_full_incoming_points(num_points);
    return 0; // Success
}

int sk_config_get_initialize_hr_with_do(Config* config, int* initialize) {
    if (config == nullptr || initialize == nullptr) {
        return -1; // Error: null pointer
    }
    *initialize = static_cast<int>(config->impl.initialize_hr_with_do());
    return 0; // Success
}

int sk_config_set_initialize_hr_with_do(Config* config, int initialize) {
    if (config == nullptr) {
        return -1; // Error: null pointer
    }
    config->impl.set_initialize_hr_with_do(initialize != 0);
    return 0; // Success
}

int sk_config_get_occultation_source(Config* config, int* occultation_source) {
    if (config == nullptr || occultation_source == nullptr) {
        return -1; // Error: null pointer
    }
    *occultation_source = static_cast<int>(config->impl.occultation_source());
    return 0; // Success
}

int sk_config_set_occultation_source(Config* config, int occultation_source) {
    if (config == nullptr) {
        return -1; // Error: null pointer
    }
    config->impl.set_occultation_source(
        static_cast<sasktran2::Config::OccultationSource>(occultation_source));
    return 0; // Success
}

int sk_config_get_emission_source(Config* config, int* emission_source) {
    if (config == nullptr || emission_source == nullptr) {
        return -1; // Error: null pointer
    }
    *emission_source = static_cast<int>(config->impl.emission_source());
    return 0; // Success
}

int sk_config_set_emission_source(Config* config, int emission_source) {
    if (config == nullptr) {
        return -1; // Error: null pointer
    }
    config->impl.set_emission_source(
        static_cast<sasktran2::Config::EmissionSource>(emission_source));
    return 0; // Success
}

int sk_config_get_los_refraction(Config* config, int* refraction) {
    if (config == nullptr || refraction == nullptr) {
        return -1; // Error: null pointer
    }
    *refraction = static_cast<int>(config->impl.los_refraction());
    return 0; // Success
}

int sk_config_set_los_refraction(Config* config, int refraction) {
    if (config == nullptr) {
        return -1; // Error: null pointer
    }
    config->impl.set_los_refraction(refraction != 0);
    return 0; // Success
}

int sk_config_get_solar_refraction(Config* config, int* refraction) {
    if (config == nullptr || refraction == nullptr) {
        return -1; // Error: null pointer
    }
    *refraction = static_cast<int>(config->impl.solar_refraction());
    return 0; // Success
}

int sk_config_set_solar_refraction(Config* config, int refraction) {
    if (config == nullptr) {
        return -1; // Error: null pointer
    }
    config->impl.set_solar_refraction(refraction != 0);
    return 0; // Success
}

int sk_config_get_multiple_scatter_refraction(Config* config, int* refraction) {
    if (config == nullptr || refraction == nullptr) {
        return -1; // Error: null pointer
    }
    *refraction = static_cast<int>(config->impl.multiple_scatter_refraction());
    return 0; // Success
}

int sk_config_set_multiple_scatter_refraction(Config* config, int refraction) {
    if (config == nullptr) {
        return -1; // Error: null pointer
    }
    config->impl.set_multiple_scatter_refraction(refraction != 0);
    return 0; // Success
}

int sk_config_get_wf_enabled(Config* config, int* enabled) {
    if (config == nullptr || enabled == nullptr) {
        return -1; // Error: null pointer
    }
    *enabled = static_cast<int>(config->impl.wf_enabled());
    return 0; // Success
}

int sk_config_set_wf_enabled(Config* config, int enabled) {
    if (config == nullptr) {
        return -1; // Error: null pointer
    }
    config->impl.set_wf_enabled(enabled != 0);
    return 0; // Success
}

int sk_config_get_wf_precision(Config* config, int* precision) {
    if (config == nullptr || precision == nullptr) {
        return -1; // Error: null pointer
    }
    *precision = static_cast<int>(config->impl.wf_precision());
    return 0; // Success
}

int sk_config_set_wf_precision(Config* config, int precision) {
    if (config == nullptr) {
        return -1; // Error: null pointer
    }
    config->impl.set_wf_precision(
        static_cast<sasktran2::Config::WeightingFunctionPrecision>(precision));
    return 0; // Success
}

int sk_config_get_stokes_basis(Config* config, int* basis) {
    if (config == nullptr || basis == nullptr) {
        return -1; // Error: null pointer
    }
    *basis = static_cast<int>(config->impl.stokes_basis());
    return 0; // Success
}

int sk_config_set_stokes_basis(Config* config, int basis) {
    if (config == nullptr) {
        return -1; // Error: null pointer
    }
    config->impl.set_stokes_basis(
        static_cast<sasktran2::Config::StokesBasis>(basis));
    return 0; // Success
}

int sk_config_get_input_validation_mode(Config* config, int* mode) {
    if (config == nullptr || mode == nullptr) {
        return -1; // Error: null pointer
    }
    *mode = static_cast<int>(config->impl.input_validation_mode());
    return 0; // Success
}

int sk_config_set_input_validation_mode(Config* config, int mode) {
    if (config == nullptr) {
        return -1; // Error: null pointer
    }
    config->impl.set_input_validation_mode(
        static_cast<sasktran2::Config::InputValidationMode>(mode));
    return 0; // Success
}

int sk_config_get_output_los_optical_depth(Config* config, int* output) {
    if (config == nullptr || output == nullptr) {
        return -1; // Error: null pointer
    }
    *output = static_cast<int>(config->impl.output_los_optical_depth());
    return 0; // Success
}

int sk_config_set_output_los_optical_depth(Config* config, int output) {
    if (config == nullptr) {
        return -1; // Error: null pointer
    }
    config->impl.set_output_los_optical_depth(output != 0);
    return 0; // Success
}
}
