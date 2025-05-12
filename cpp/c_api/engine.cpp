#include "sasktran2/engine.h"
#include "engine.h"
#include "internal_types.h"
#include "output.h"
#include "sasktran2/output.h"
#include <sasktran2.h>

struct Engine {
    std::unique_ptr<Sasktran2Interface> impl;
    Config* _config;

    Engine(Config* config, Geometry1D* geometry,
           ViewingGeometry* viewing_geometry)
        : _config(config) {
        try {
            if (_config->impl.num_stokes() == 1) {
                impl = std::make_unique<Sasktran2<1>>(
                    config->impl, geometry->impl.get(), viewing_geometry->impl);
            } else if (config->impl.num_stokes() == 3) {
                impl = std::make_unique<Sasktran2<3>>(
                    config->impl, geometry->impl.get(), viewing_geometry->impl);
            } else {
                impl = nullptr;
            }
        } catch (const std::exception& e) {
            // Handle the exception, log it, etc.
            impl = nullptr;
        }
    }

    int calculate_radiance(Atmosphere* atmosphere, OutputC* output,
                           int only_initialize) {
        bool initialize = only_initialize != 0;
        try {
            if (impl) {
                if (_config->impl.num_stokes() == 1) {
                    Sasktran2<1>* impl1 =
                        dynamic_cast<Sasktran2<1>*>(impl.get());
                    impl1->calculate_radiance(
                        *static_cast<sasktran2::atmosphere::Atmosphere<1>*>(
                            atmosphere->impl.get()),
                        *static_cast<sasktran2::Output<1>*>(output->impl.get()),
                        initialize);
                    return 0;
                } else if (_config->impl.num_stokes() == 3) {
                    Sasktran2<3>* impl3 =
                        dynamic_cast<Sasktran2<3>*>(impl.get());
                    impl3->calculate_radiance(
                        *static_cast<sasktran2::atmosphere::Atmosphere<3>*>(
                            atmosphere->impl.get()),
                        *static_cast<sasktran2::Output<3>*>(output->impl.get()),
                        initialize);
                    return 0;
                } else {
                    return -2; // Error: invalid number of Stokes parameters
                }
            }
            return -1; // Error: impl is null
        } catch (const std::exception& e) {
            // Handle the exception, log it, etc.
            return -3; // Error: exception occurred
        }
    }
};

extern "C" {
Engine* sk_engine_create(Config* config, Geometry1D* geometry,
                         ViewingGeometry* viewing_geometry) {
    return new Engine(config, geometry, viewing_geometry);
}

void sk_engine_destroy(Engine* engine) { delete engine; }

int sk_engine_calculate_radiance(Engine* engine, Atmosphere* atmosphere,
                                 OutputC* output, int only_initialize) {
    return engine->calculate_radiance(atmosphere, output, only_initialize);
}

int sk_engine_calculate_radiance_thread(Engine* engine, Atmosphere* atmosphere,
                                        OutputC* output, int wavelength_idx,
                                        int thread_idx) {
    try {
        if (engine->impl) {
            if (engine->_config->impl.num_stokes() == 1) {
                Sasktran2<1>* impl1 =
                    dynamic_cast<Sasktran2<1>*>(engine->impl.get());
                impl1->calculate_radiance_thread(
                    *static_cast<sasktran2::atmosphere::Atmosphere<1>*>(
                        atmosphere->impl.get()),
                    *static_cast<sasktran2::Output<1>*>(output->impl.get()),
                    wavelength_idx, thread_idx);
                return 0;
            } else if (engine->_config->impl.num_stokes() == 3) {
                Sasktran2<3>* impl3 =
                    dynamic_cast<Sasktran2<3>*>(engine->impl.get());
                impl3->calculate_radiance_thread(
                    *static_cast<sasktran2::atmosphere::Atmosphere<3>*>(
                        atmosphere->impl.get()),
                    *static_cast<sasktran2::Output<3>*>(output->impl.get()),
                    wavelength_idx, thread_idx);
                return 0;
            } else {
                return -2; // Error: invalid number of Stokes parameters
            }
        }
        return -1; // Error: impl is null
    } catch (const std::exception& e) {
        // Handle the exception, log it, etc.
        return -3; // Error: exception occurred
    }
}

int sk_openmp_support_enabled() {
#ifdef SKTRAN_OPENMP_SUPPORT
    return 1; // OpenMP support is enabled
#else
    return 0; // OpenMP support is not enabled
#endif
}
}
