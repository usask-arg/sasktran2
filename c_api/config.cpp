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

int sk_config_set_num_streams(Config* config, int num_streams) {
    if (config == nullptr) {
        return -1; // Error: null pointer
    }
    config->impl.set_num_do_streams(num_streams);
    return 0; // Success
}
}
