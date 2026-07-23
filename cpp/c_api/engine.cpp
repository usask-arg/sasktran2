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
            // Prinit the error message
            spdlog::error("Error initializing Engine: {}", e.what());
        }
    }

    Engine(Config* config, Geometry2D* geometry,
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
            impl = nullptr;
            spdlog::error("Error initializing Geometry2D Engine: {}", e.what());
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
            spdlog::error("Error calculating radiance: {}", e.what());
            return -3; // Error: exception occurred
        }
    }
};

extern "C" {
Engine* sk_engine_create(Config* config, Geometry1D* geometry,
                         ViewingGeometry* viewing_geometry) {

    // create a new engine instance
    Engine* engine = new Engine(config, geometry, viewing_geometry);

    if (engine->impl == nullptr) {
        delete engine;
        return nullptr;
    }

    return engine;
}

Engine* sk_engine_create_2d(Config* config, Geometry2D* geometry,
                            ViewingGeometry* viewing_geometry) {
    Engine* engine = new Engine(config, geometry, viewing_geometry);

    if (engine->impl == nullptr) {
        delete engine;
        return nullptr;
    }

    return engine;
}

void sk_engine_destroy(Engine* engine) { delete engine; }

int sk_engine_calculate_radiance(Engine* engine, Atmosphere* atmosphere,
                                 OutputC* output, int only_initialize) {
    return engine->calculate_radiance(atmosphere, output, only_initialize);
}

int sk_engine_effective_wavelength_batch_size(Engine* engine,
                                              int num_wavelengths) {
    try {
        if (engine == nullptr || !engine->impl) {
            return -1;
        }
        if (engine->_config->impl.num_stokes() == 1) {
            auto* impl = dynamic_cast<Sasktran2<1>*>(engine->impl.get());
            return impl->effective_wavelength_batch_size(num_wavelengths);
        }
        if (engine->_config->impl.num_stokes() == 3) {
            auto* impl = dynamic_cast<Sasktran2<3>*>(engine->impl.get());
            return impl->effective_wavelength_batch_size(num_wavelengths);
        }
        return -2;
    } catch (const std::exception&) {
        return -3;
    }
}

int sk_engine_supports_linearization(Engine* engine, int mode, int* supported) {
    try {
        if (engine == nullptr || !engine->impl || supported == nullptr) {
            return -1;
        }
        if (mode < static_cast<int>(sasktran2::LinearizationMode::Jacobian) ||
            mode > static_cast<int>(sasktran2::LinearizationMode::VJP)) {
            return -2;
        }

        const auto linearization_mode =
            static_cast<sasktran2::LinearizationMode>(mode);
        if (engine->_config->impl.num_stokes() == 1) {
            auto* impl = dynamic_cast<Sasktran2<1>*>(engine->impl.get());
            *supported =
                impl->supports_linearization(linearization_mode) ? 1 : 0;
            return 0;
        }
        if (engine->_config->impl.num_stokes() == 3) {
            auto* impl = dynamic_cast<Sasktran2<3>*>(engine->impl.get());
            *supported =
                impl->supports_linearization(linearization_mode) ? 1 : 0;
            return 0;
        }
        return -2;
    } catch (const std::exception&) {
        return -3;
    }
}

int sk_engine_linearization_backend(Engine* engine, int mode, int* backend) {
    try {
        if (engine == nullptr || !engine->impl || backend == nullptr) {
            return -1;
        }
        if (mode < static_cast<int>(sasktran2::LinearizationMode::Jacobian) ||
            mode > static_cast<int>(sasktran2::LinearizationMode::VJP)) {
            return -2;
        }
        const auto linearization_mode =
            static_cast<sasktran2::LinearizationMode>(mode);
        if (engine->_config->impl.num_stokes() == 1) {
            auto* impl = dynamic_cast<Sasktran2<1>*>(engine->impl.get());
            *backend = static_cast<int>(
                impl->linearization_backend(linearization_mode));
            return 0;
        }
        if (engine->_config->impl.num_stokes() == 3) {
            auto* impl = dynamic_cast<Sasktran2<3>*>(engine->impl.get());
            *backend = static_cast<int>(
                impl->linearization_backend(linearization_mode));
            return 0;
        }
        return -2;
    } catch (const std::exception&) {
        return -3;
    }
}

int sk_engine_calculate_jvp(Engine* engine, Atmosphere* atmosphere,
                            OutputJVP* output) {
    if (engine == nullptr || atmosphere == nullptr || output == nullptr) {
        return -1;
    }
    try {
        if (engine->_config->impl.num_stokes() == 1) {
            auto* impl = dynamic_cast<Sasktran2<1>*>(engine->impl.get());
            impl->calculate_radiance(
                *static_cast<sasktran2::atmosphere::Atmosphere<1>*>(
                    atmosphere->impl.get()),
                *static_cast<sasktran2::Output<1>*>(output->impl.get()));
            return 0;
        }
        if (engine->_config->impl.num_stokes() == 3) {
            auto* impl = dynamic_cast<Sasktran2<3>*>(engine->impl.get());
            impl->calculate_radiance(
                *static_cast<sasktran2::atmosphere::Atmosphere<3>*>(
                    atmosphere->impl.get()),
                *static_cast<sasktran2::Output<3>*>(output->impl.get()));
            return 0;
        }
        return -2;
    } catch (const std::exception&) {
        return -3;
    }
}

int sk_engine_calculate_vjp(Engine* engine, Atmosphere* atmosphere,
                            OutputVJP* output) {
    if (engine == nullptr || atmosphere == nullptr || output == nullptr) {
        return -1;
    }
    try {
        int result = -2;
        if (engine->_config->impl.num_stokes() == 1) {
            auto* impl = dynamic_cast<Sasktran2<1>*>(engine->impl.get());
            impl->calculate_radiance(
                *static_cast<sasktran2::atmosphere::Atmosphere<1>*>(
                    atmosphere->impl.get()),
                *static_cast<sasktran2::Output<1>*>(output->impl.get()));
            result = 0;
        } else if (engine->_config->impl.num_stokes() == 3) {
            auto* impl = dynamic_cast<Sasktran2<3>*>(engine->impl.get());
            impl->calculate_radiance(
                *static_cast<sasktran2::atmosphere::Atmosphere<3>*>(
                    atmosphere->impl.get()),
                *static_cast<sasktran2::Output<3>*>(output->impl.get()));
            result = 0;
        }
        if (result == 0) {
            result = output->finalize();
        }
        return result;
    } catch (const std::exception&) {
        return -3;
    }
}

int sk_engine_calculate_radiance_block_thread(Engine* engine, OutputC* output,
                                              int wavelength_start,
                                              int wavelength_count,
                                              int thread_idx) {
    try {
        if (engine == nullptr || !engine->impl || output == nullptr) {
            return -1;
        }
        if (wavelength_start < 0 || wavelength_count < 1 || thread_idx < 0) {
            return -2;
        }
        const sasktran2::WavelengthBlock block{wavelength_start,
                                               wavelength_count};
        if (engine->_config->impl.num_stokes() == 1) {
            auto* impl = dynamic_cast<Sasktran2<1>*>(engine->impl.get());
            impl->calculate_radiance_block_thread(
                *static_cast<sasktran2::Output<1>*>(output->impl.get()), block,
                thread_idx);
            return 0;
        }
        if (engine->_config->impl.num_stokes() == 3) {
            auto* impl = dynamic_cast<Sasktran2<3>*>(engine->impl.get());
            impl->calculate_radiance_block_thread(
                *static_cast<sasktran2::Output<3>*>(output->impl.get()), block,
                thread_idx);
            return 0;
        }
        return -2;
    } catch (const std::exception&) {
        return -3;
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
