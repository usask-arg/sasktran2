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
}
