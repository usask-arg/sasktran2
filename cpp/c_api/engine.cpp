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

int sk_engine_set_2d_refractive_profiles(Engine* engine, const double* profiles,
                                         int num_rays, int num_altitudes) {
    if (engine == nullptr || !engine->impl || profiles == nullptr ||
        num_rays < 0 || num_altitudes < 0) {
        return -1;
    }
    try {
        using RowMajorMatrix = Eigen::Matrix<double, Eigen::Dynamic,
                                             Eigen::Dynamic, Eigen::RowMajor>;
        const Eigen::Map<const RowMajorMatrix> mapped(profiles, num_rays,
                                                      num_altitudes);
        if (engine->_config->impl.num_stokes() == 1) {
            auto* impl = dynamic_cast<Sasktran2<1>*>(engine->impl.get());
            impl->set_2d_refractive_profiles(mapped);
            return 0;
        }
        if (engine->_config->impl.num_stokes() == 3) {
            auto* impl = dynamic_cast<Sasktran2<3>*>(engine->impl.get());
            impl->set_2d_refractive_profiles(mapped);
            return 0;
        }
        return -2;
    } catch (const std::exception& error) {
        spdlog::error("Error setting Geometry2D refractive profiles: {}",
                      error.what());
        return -3;
    }
}

int sk_engine_get_2d_surface_interpolation_weights(Engine* engine,
                                                   double* weights,
                                                   int num_rays,
                                                   int num_horizontal) {
    if (engine == nullptr || !engine->impl || weights == nullptr ||
        num_rays < 0 || num_horizontal < 0) {
        return -1;
    }
    try {
        using RowMajorMatrix = Eigen::Matrix<double, Eigen::Dynamic,
                                             Eigen::Dynamic, Eigen::RowMajor>;
        Eigen::Map<RowMajorMatrix> mapped(weights, num_rays, num_horizontal);
        if (engine->_config->impl.num_stokes() == 1) {
            auto* impl = dynamic_cast<Sasktran2<1>*>(engine->impl.get());
            impl->assign_2d_surface_interpolation_weights(mapped);
            return 0;
        }
        if (engine->_config->impl.num_stokes() == 3) {
            auto* impl = dynamic_cast<Sasktran2<3>*>(engine->impl.get());
            impl->assign_2d_surface_interpolation_weights(mapped);
            return 0;
        }
        return -2;
    } catch (const std::exception& error) {
        spdlog::error("Error obtaining Geometry2D surface weights: {}",
                      error.what());
        return -3;
    }
}

int sk_engine_get_2d_horizontal_edge_usage(Engine* engine, int* usage,
                                           int num_rays) {
    if (engine == nullptr || !engine->impl || usage == nullptr ||
        num_rays < 0) {
        return -1;
    }
    try {
        using RowMajorIntMatrix =
            Eigen::Matrix<int, Eigen::Dynamic, Eigen::Dynamic, Eigen::RowMajor>;
        Eigen::Map<RowMajorIntMatrix> mapped(usage, num_rays, 2);
        if (engine->_config->impl.num_stokes() == 1) {
            auto* impl = dynamic_cast<Sasktran2<1>*>(engine->impl.get());
            impl->assign_2d_horizontal_edge_usage(mapped);
            return 0;
        }
        if (engine->_config->impl.num_stokes() == 3) {
            auto* impl = dynamic_cast<Sasktran2<3>*>(engine->impl.get());
            impl->assign_2d_horizontal_edge_usage(mapped);
            return 0;
        }
        return -2;
    } catch (const std::exception& error) {
        spdlog::error("Error obtaining Geometry2D horizontal edge usage: {}",
                      error.what());
        return -3;
    }
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
    if (engine == nullptr || !engine->impl || engine->_config == nullptr ||
        atmosphere == nullptr || !atmosphere->impl || output == nullptr ||
        !output->impl) {
        return -1;
    }
    try {
        if (engine->_config->impl.num_stokes() == 1) {
            auto* impl = dynamic_cast<Sasktran2<1>*>(engine->impl.get());
            auto* native_atmosphere =
                dynamic_cast<sasktran2::atmosphere::Atmosphere<1>*>(
                    atmosphere->impl.get());
            auto* native_output =
                dynamic_cast<sasktran2::OutputJVP<1>*>(output->impl.get());
            if (impl == nullptr || native_atmosphere == nullptr ||
                native_output == nullptr) {
                return -2;
            }
            if (impl->supports_linearization(
                    sasktran2::LinearizationMode::JVP)) {
                impl->calculate_jvp(*native_atmosphere, *native_output);
            } else {
                impl->calculate_radiance(*native_atmosphere, *native_output);
            }
            return 0;
        }
        if (engine->_config->impl.num_stokes() == 3) {
            auto* impl = dynamic_cast<Sasktran2<3>*>(engine->impl.get());
            auto* native_atmosphere =
                dynamic_cast<sasktran2::atmosphere::Atmosphere<3>*>(
                    atmosphere->impl.get());
            auto* native_output =
                dynamic_cast<sasktran2::OutputJVP<3>*>(output->impl.get());
            if (impl == nullptr || native_atmosphere == nullptr ||
                native_output == nullptr) {
                return -2;
            }
            if (impl->supports_linearization(
                    sasktran2::LinearizationMode::JVP)) {
                impl->calculate_jvp(*native_atmosphere, *native_output);
            } else {
                impl->calculate_radiance(*native_atmosphere, *native_output);
            }
            return 0;
        }
        return -2;
    } catch (const std::exception& error) {
        spdlog::error("Error calculating JVP: {}", error.what());
        return -3;
    }
}

int sk_engine_initialize_jvp(Engine* engine, Atmosphere* atmosphere,
                             OutputJVP* output) {
    if (engine == nullptr || !engine->impl || engine->_config == nullptr ||
        atmosphere == nullptr || !atmosphere->impl || output == nullptr ||
        !output->impl) {
        return -1;
    }
    try {
        if (engine->_config->impl.num_stokes() == 1) {
            auto* impl = dynamic_cast<Sasktran2<1>*>(engine->impl.get());
            auto* native_atmosphere =
                dynamic_cast<sasktran2::atmosphere::Atmosphere<1>*>(
                    atmosphere->impl.get());
            auto* native_output =
                dynamic_cast<sasktran2::OutputJVP<1>*>(output->impl.get());
            if (impl == nullptr || native_atmosphere == nullptr ||
                native_output == nullptr) {
                return -2;
            }
            impl->initialize_jvp(*native_atmosphere, *native_output);
            return 0;
        }
        if (engine->_config->impl.num_stokes() == 3) {
            auto* impl = dynamic_cast<Sasktran2<3>*>(engine->impl.get());
            auto* native_atmosphere =
                dynamic_cast<sasktran2::atmosphere::Atmosphere<3>*>(
                    atmosphere->impl.get());
            auto* native_output =
                dynamic_cast<sasktran2::OutputJVP<3>*>(output->impl.get());
            if (impl == nullptr || native_atmosphere == nullptr ||
                native_output == nullptr) {
                return -2;
            }
            impl->initialize_jvp(*native_atmosphere, *native_output);
            return 0;
        }
        return -2;
    } catch (const std::exception& error) {
        spdlog::error("Error initializing JVP: {}", error.what());
        return -3;
    }
}

int sk_engine_calculate_jvp_wavelength_thread(Engine* engine, OutputJVP* output,
                                              int wavelength, int thread_idx) {
    if (engine == nullptr || !engine->impl || engine->_config == nullptr ||
        output == nullptr || !output->impl) {
        return -1;
    }
    try {
        if (engine->_config->impl.num_stokes() == 1) {
            auto* impl = dynamic_cast<Sasktran2<1>*>(engine->impl.get());
            auto* native_output =
                dynamic_cast<sasktran2::OutputJVP<1>*>(output->impl.get());
            if (impl == nullptr || native_output == nullptr) {
                return -2;
            }
            impl->calculate_jvp_wavelength_thread(*native_output, wavelength,
                                                  thread_idx);
            return 0;
        }
        if (engine->_config->impl.num_stokes() == 3) {
            auto* impl = dynamic_cast<Sasktran2<3>*>(engine->impl.get());
            auto* native_output =
                dynamic_cast<sasktran2::OutputJVP<3>*>(output->impl.get());
            if (impl == nullptr || native_output == nullptr) {
                return -2;
            }
            impl->calculate_jvp_wavelength_thread(*native_output, wavelength,
                                                  thread_idx);
            return 0;
        }
        return -2;
    } catch (const std::exception& error) {
        spdlog::error("Error calculating JVP wavelength: {}", error.what());
        return -3;
    }
}

int sk_engine_calculate_vjp(Engine* engine, Atmosphere* atmosphere,
                            OutputVJP* output) {
    if (engine == nullptr || !engine->impl || engine->_config == nullptr ||
        atmosphere == nullptr || !atmosphere->impl || output == nullptr ||
        !output->impl) {
        return -1;
    }
    try {
        if (engine->_config->impl.num_stokes() == 1) {
            auto* impl = dynamic_cast<Sasktran2<1>*>(engine->impl.get());
            auto* native_atmosphere =
                dynamic_cast<sasktran2::atmosphere::Atmosphere<1>*>(
                    atmosphere->impl.get());
            auto* native_output =
                dynamic_cast<sasktran2::OutputVJP<1>*>(output->impl.get());
            if (impl == nullptr || native_atmosphere == nullptr ||
                native_output == nullptr) {
                return -2;
            }
            if (impl->supports_linearization(
                    sasktran2::LinearizationMode::VJP)) {
                impl->calculate_vjp(*native_atmosphere, *native_output);
            } else {
                impl->calculate_radiance(*native_atmosphere, *native_output);
            }
            native_output->finalize();
            return 0;
        } else if (engine->_config->impl.num_stokes() == 3) {
            auto* impl = dynamic_cast<Sasktran2<3>*>(engine->impl.get());
            auto* native_atmosphere =
                dynamic_cast<sasktran2::atmosphere::Atmosphere<3>*>(
                    atmosphere->impl.get());
            auto* native_output =
                dynamic_cast<sasktran2::OutputVJP<3>*>(output->impl.get());
            if (impl == nullptr || native_atmosphere == nullptr ||
                native_output == nullptr) {
                return -2;
            }
            if (impl->supports_linearization(
                    sasktran2::LinearizationMode::VJP)) {
                impl->calculate_vjp(*native_atmosphere, *native_output);
            } else {
                impl->calculate_radiance(*native_atmosphere, *native_output);
            }
            native_output->finalize();
            return 0;
        }
        return -2;
    } catch (const std::exception& error) {
        spdlog::error("Error calculating VJP: {}", error.what());
        return -3;
    }
}

int sk_engine_initialize_vjp(Engine* engine, Atmosphere* atmosphere,
                             OutputVJP* output) {
    if (engine == nullptr || !engine->impl || engine->_config == nullptr ||
        atmosphere == nullptr || !atmosphere->impl || output == nullptr ||
        !output->impl) {
        return -1;
    }
    try {
        if (engine->_config->impl.num_stokes() == 1) {
            auto* impl = dynamic_cast<Sasktran2<1>*>(engine->impl.get());
            auto* native_atmosphere =
                dynamic_cast<sasktran2::atmosphere::Atmosphere<1>*>(
                    atmosphere->impl.get());
            auto* native_output =
                dynamic_cast<sasktran2::OutputVJP<1>*>(output->impl.get());
            if (impl == nullptr || native_atmosphere == nullptr ||
                native_output == nullptr) {
                return -2;
            }
            impl->initialize_vjp(*native_atmosphere, *native_output);
            return 0;
        }
        if (engine->_config->impl.num_stokes() == 3) {
            auto* impl = dynamic_cast<Sasktran2<3>*>(engine->impl.get());
            auto* native_atmosphere =
                dynamic_cast<sasktran2::atmosphere::Atmosphere<3>*>(
                    atmosphere->impl.get());
            auto* native_output =
                dynamic_cast<sasktran2::OutputVJP<3>*>(output->impl.get());
            if (impl == nullptr || native_atmosphere == nullptr ||
                native_output == nullptr) {
                return -2;
            }
            impl->initialize_vjp(*native_atmosphere, *native_output);
            return 0;
        }
        return -2;
    } catch (const std::exception& error) {
        spdlog::error("Error initializing VJP: {}", error.what());
        return -3;
    }
}

int sk_engine_calculate_vjp_block_thread(Engine* engine, OutputVJP* output,
                                         int wavelength_start,
                                         int wavelength_count, int thread_idx) {
    if (engine == nullptr || !engine->impl || engine->_config == nullptr ||
        output == nullptr || !output->impl || wavelength_start < 0 ||
        wavelength_count < 1 || thread_idx < 0) {
        return -1;
    }
    try {
        const sasktran2::WavelengthBlock block{wavelength_start,
                                               wavelength_count};
        if (engine->_config->impl.num_stokes() == 1) {
            auto* impl = dynamic_cast<Sasktran2<1>*>(engine->impl.get());
            auto* native_output =
                dynamic_cast<sasktran2::OutputVJP<1>*>(output->impl.get());
            if (impl == nullptr || native_output == nullptr) {
                return -2;
            }
            impl->calculate_vjp_block_thread(*native_output, block, thread_idx);
            return 0;
        }
        if (engine->_config->impl.num_stokes() == 3) {
            auto* impl = dynamic_cast<Sasktran2<3>*>(engine->impl.get());
            auto* native_output =
                dynamic_cast<sasktran2::OutputVJP<3>*>(output->impl.get());
            if (impl == nullptr || native_output == nullptr) {
                return -2;
            }
            impl->calculate_vjp_block_thread(*native_output, block, thread_idx);
            return 0;
        }
        return -2;
    } catch (const std::exception& error) {
        spdlog::error("Error calculating VJP wavelength block: {}",
                      error.what());
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
