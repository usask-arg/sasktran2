#include "atmosphere.h"
#include "sasktran2/atmosphere/atmosphere.h"
#include "sasktran2/atmosphere/grid_storage.h"
#include <sasktran2.h>
#include "internal_types.h"

AtmosphereStorage::AtmosphereStorage(int nlocation, int nwavel,
                                     int nphase_moments, int nstokes,
                                     double* ssa, double* total_extinction,
                                     double* emission_source, double* leg_coeff,
                                     double* solar_irradiance) {
    // Create the Eigen Maps
    Eigen::Map<Eigen::MatrixXd> ssa_map(ssa, nlocation, nwavel);
    Eigen::Map<Eigen::MatrixXd> total_extinction_map(total_extinction,
                                                     nlocation, nwavel);
    Eigen::Map<Eigen::MatrixXd> emission_source_map(emission_source, nlocation,
                                                    nwavel);
    Eigen::TensorMap<Eigen::Tensor<double, 3>> leg_coeff_map(
        leg_coeff, nphase_moments, nlocation, nwavel);
    Eigen::Map<Eigen::VectorXd> solar_irradiance_map(solar_irradiance, nwavel);

    if (nstokes == 1) {
        impl = std::make_unique<
            sasktran2::atmosphere::AtmosphereGridStorageFull<1>>(
            ssa_map, total_extinction_map, emission_source_map, leg_coeff_map,
            solar_irradiance_map);
    } else if (nstokes == 3) {
        impl = std::make_unique<
            sasktran2::atmosphere::AtmosphereGridStorageFull<3>>(
            ssa_map, total_extinction_map, emission_source_map, leg_coeff_map,
            solar_irradiance_map);
    } else {
        impl = nullptr;
    }
}

int AtmosphereStorage::get_derivative_mapping(const char* name,
                                              DerivativeMapping** mapping) {
    if (impl) {
        auto storage1 =
            dynamic_cast<sasktran2::atmosphere::AtmosphereGridStorageFull<1>*>(
                impl.get());
        auto storage3 =
            dynamic_cast<sasktran2::atmosphere::AtmosphereGridStorageFull<3>*>(
                impl.get());
        if (storage1) {
            auto& map = storage1->get_derivative_mapping(name);
            *mapping = new DerivativeMapping(&map);
            return 0;
        } else if (storage3) {
            auto& map = storage3->get_derivative_mapping(name);
            *mapping = new DerivativeMapping(&map);
            return 0;
        } else {
            return -3; // Error: storage implementation is not
                       // AtmosphereGridStorageFull
        }
    } else {
        return -2; // Error: storage implementation is null
    }
}

Surface::Surface(int nwavel, int nstokes, double* emission) {
    // Create the Eigen Map
    Eigen::Map<Eigen::VectorXd> emission_map(emission, nwavel);
    if (nstokes == 1) {
        impl = std::make_unique<sasktran2::atmosphere::Surface<1>>(
            nwavel, emission_map);
    } else if (nstokes == 3) {
        impl = std::make_unique<sasktran2::atmosphere::Surface<3>>(
            nwavel, emission_map);
    } else {
        impl = nullptr;
    }
}

Atmosphere::Atmosphere(AtmosphereStorage* storage, Surface* surface,
                       bool calculate_derivatives,
                       bool calculate_emission_derivatives) {
    if (storage->impl) {
        if (surface->impl) {
            sasktran2::atmosphere::Surface<1>* surface1 =
                dynamic_cast<sasktran2::atmosphere::Surface<1>*>(
                    surface->impl.get());

            if (surface1) {
                sasktran2::atmosphere::AtmosphereGridStorageFull<1>* storage1 =
                    dynamic_cast<
                        sasktran2::atmosphere::AtmosphereGridStorageFull<1>*>(
                        storage->impl.get());

                if (storage1) {
                    impl =
                        std::make_unique<sasktran2::atmosphere::Atmosphere<1>>(
                            *storage1, *surface1, calculate_derivatives,
                            calculate_emission_derivatives);
                } else {
                    impl = nullptr;
                }
            }

            sasktran2::atmosphere::Surface<3>* surface3 =
                dynamic_cast<sasktran2::atmosphere::Surface<3>*>(
                    surface->impl.get());

            if (surface3) {
                sasktran2::atmosphere::AtmosphereGridStorageFull<3>* storage3 =
                    dynamic_cast<
                        sasktran2::atmosphere::AtmosphereGridStorageFull<3>*>(
                        storage->impl.get());

                if (storage3) {
                    impl =
                        std::make_unique<sasktran2::atmosphere::Atmosphere<3>>(
                            *storage3, *surface3, calculate_derivatives,
                            calculate_emission_derivatives);
                } else {
                    impl = nullptr;
                }
            }

        } else {
            impl = nullptr;
        }
    } else {
        impl = nullptr;
    }
}

extern "C" {
AtmosphereStorage*
sk_atmosphere_storage_create(int nlocation, int nwavel, int nphase_moments,
                             int nstokes, double* ssa, double* total_extinction,
                             double* emission_source, double* leg_coeff,
                             double* solar_irradiance) {
    return new AtmosphereStorage(nlocation, nwavel, nphase_moments, nstokes,
                                 ssa, total_extinction, emission_source,
                                 leg_coeff, solar_irradiance);
}

void sk_atmosphere_storage_destroy(AtmosphereStorage* storage) {
    delete storage;
}

int sk_atmosphere_storage_get_derivative_mapping(AtmosphereStorage* storage,
                                                 const char* name,
                                                 DerivativeMapping** mapping) {
    if (storage == nullptr) {
        return -1; // Error: storage is null
    }

    return storage->get_derivative_mapping(name, mapping);
}

int sk_atmosphere_storage_get_derivative_mapping_by_index(
    AtmosphereStorage* storage, int index, DerivativeMapping** mapping) {
    if (storage == nullptr) {
        return -1; // Error: storage is null
    }
    if (storage->impl == nullptr) {
        return -2; // Error: storage implementation is null
    }

    auto* impl1 =
        dynamic_cast<sasktran2::atmosphere::AtmosphereGridStorageFull<1>*>(
            storage->impl.get());
    auto* impl3 =
        dynamic_cast<sasktran2::atmosphere::AtmosphereGridStorageFull<3>*>(
            storage->impl.get());

    if (impl1) {
        auto& mapping_list = impl1->derivative_mappings();
        auto it = mapping_list.begin();
        std::advance(it, index);

        if (it != mapping_list.end()) {
            auto& map = it->second;
            *mapping = new DerivativeMapping(&map);
            return 0;
        } else {
            return -3; // Error: index out of bounds
        }
    } else if (impl3) {
        auto& mapping_list = impl1->derivative_mappings();
        auto it = mapping_list.begin();
        std::advance(it, index);

        if (it != mapping_list.end()) {
            auto& map = it->second;
            *mapping = new DerivativeMapping(&map);
            return 0;
        } else {
            return -3; // Error: index out of bounds
        }
    } else {
        return -4; // Error: storage implementation is not
                   // AtmosphereGridStorageFull
    }
}

int sk_atmosphere_storage_get_num_derivative_mappings(
    AtmosphereStorage* storage, int* num_mappings) {
    if (storage == nullptr) {
        return -1; // Error: storage is null
    }
    if (storage->impl == nullptr) {
        return -2; // Error: storage implementation is null
    }

    auto* impl1 =
        dynamic_cast<sasktran2::atmosphere::AtmosphereGridStorageFull<1>*>(
            storage->impl.get());
    auto* impl3 =
        dynamic_cast<sasktran2::atmosphere::AtmosphereGridStorageFull<3>*>(
            storage->impl.get());

    if (impl1) {
        *num_mappings = impl1->derivative_mappings().size();
        return 0;
    } else if (impl3) {
        *num_mappings = impl3->derivative_mappings().size();
        return 0;
    } else {
        return -3; // Error: storage implementation is not
                   // AtmosphereGridStorageFull
    }
}

int sk_atmosphere_storage_get_derivative_mapping_name(
    AtmosphereStorage* storage, int index, const char** name) {
    if (storage == nullptr) {
        return -1; // Error: storage is null
    }
    if (storage->impl == nullptr) {
        return -2; // Error: storage implementation is null
    }

    auto* impl1 =
        dynamic_cast<sasktran2::atmosphere::AtmosphereGridStorageFull<1>*>(
            storage->impl.get());
    auto* impl3 =
        dynamic_cast<sasktran2::atmosphere::AtmosphereGridStorageFull<3>*>(
            storage->impl.get());

    if (impl1) {
        auto& mapping_list = impl1->derivative_mappings();
        auto it = mapping_list.begin();
        std::advance(it, index);

        if (it != mapping_list.end()) {
            auto& map = it->second;
            *name = it->first.c_str();
            return 0;
        } else {
            return -3; // Error: index out of bounds
        }
    } else if (impl3) {
        auto& mapping_list = impl3->derivative_mappings();
        auto it = mapping_list.begin();
        std::advance(it, index);

        if (it != mapping_list.end()) {
            auto& map = it->second;
            *name = it->first.c_str();
            return 0;
        } else {
            return -3; // Error: index out of bounds
        }
    } else {
        return -4; // Error: storage implementation is not
                   // AtmosphereGridStorageFull
    }
}

int sk_atmosphere_storage_set_zero(AtmosphereStorage* storage) {
    if (storage == nullptr) {
        return -1; // Error: storage is null
    }
    if (storage->impl == nullptr) {
        return -2; // Error: storage implementation is null
    }

    auto* impl1 =
        dynamic_cast<sasktran2::atmosphere::AtmosphereGridStorageFull<1>*>(
            storage->impl.get());
    auto* impl3 =
        dynamic_cast<sasktran2::atmosphere::AtmosphereGridStorageFull<3>*>(
            storage->impl.get());

    if (impl1) {
        impl1->set_zero();
    } else if (impl3) {
        impl3->set_zero();
    } else {
        return -3; // Error: storage implementation is not
                   // AtmosphereGridStorageFull
    }
    return 0;
}

Surface* sk_surface_create(int nwavel, int nstokes, double* emission) {
    if (emission == nullptr) {
        return nullptr; // Error: emission is null
    }
    return new Surface(nwavel, nstokes, emission);
}

void sk_surface_destroy(Surface* storage) { delete storage; }

int sk_surface_set_brdf(Surface* surface, BRDF* brdf, double* brdf_args) {
    if (surface == nullptr || brdf == nullptr) {
        return -1; // Error: surface or brdf is null
    }
    if (surface->impl == nullptr || brdf->impl == nullptr) {
        return -2; // Error: surface or brdf implementation is null
    }

    auto* impl1 =
        dynamic_cast<sasktran2::atmosphere::Surface<1>*>(surface->impl.get());
    auto* impl3 =
        dynamic_cast<sasktran2::atmosphere::Surface<3>*>(surface->impl.get());

    if (impl1) {
        auto brdf1 =
            std::dynamic_pointer_cast<sasktran2::atmosphere::brdf::BRDF<1>>(
                brdf->impl);
        if (!brdf1) {
            return -4; // Error: brdf implementation is not BRDF<1>
        }
        int num_brdf_args = brdf1->num_args();
        int nwavel = impl1->emission().size();
        Eigen::Map<Eigen::MatrixXd> brdf_args_map(brdf_args, num_brdf_args,
                                                  nwavel);

        impl1->set_brdf_object_with_memory(brdf1, brdf_args_map);
    } else if (impl3) {
        auto brdf3 =
            std::dynamic_pointer_cast<sasktran2::atmosphere::brdf::BRDF<3>>(
                brdf->impl);
        if (!brdf3) {
            return -4; // Error: brdf implementation is not BRDF<1>
        }
        int num_brdf_args = brdf3->num_args();
        int nwavel = impl3->emission().size();
        Eigen::Map<Eigen::MatrixXd> brdf_args_map(brdf_args, num_brdf_args,
                                                  nwavel);

        impl3->set_brdf_object_with_memory(brdf3, brdf_args_map);
    } else {
        return -3; // Error: surface implementation is not Surface
    }
    return 0;
}

int sk_surface_get_derivative_mapping(Surface* storage, const char* name,
                                      SurfaceDerivativeMapping** mapping) {
    if (storage == nullptr) {
        return -1; // Error: storage is null
    }
    if (storage->impl == nullptr) {
        return -2; // Error: storage implementation is null
    }

    auto* impl1 =
        dynamic_cast<sasktran2::atmosphere::Surface<1>*>(storage->impl.get());
    auto* impl3 =
        dynamic_cast<sasktran2::atmosphere::Surface<3>*>(storage->impl.get());

    if (impl1) {
        auto& map = impl1->get_derivative_mapping(name);
        *mapping = new SurfaceDerivativeMapping(&map);
    } else if (impl3) {
        auto& map = impl3->get_derivative_mapping(name);
        *mapping = new SurfaceDerivativeMapping(&map);
    } else {
        return -3; // Error: storage implementation is not Surface
    }
    return 0;
}

int sk_surface_get_num_derivative_mappings(Surface* storage,
                                           int* num_mappings) {
    if (storage == nullptr) {
        return -1; // Error: storage is null
    }
    if (storage->impl == nullptr) {
        return -2; // Error: storage implementation is null
    }

    auto* impl1 =
        dynamic_cast<sasktran2::atmosphere::Surface<1>*>(storage->impl.get());
    auto* impl3 =
        dynamic_cast<sasktran2::atmosphere::Surface<3>*>(storage->impl.get());

    if (impl1) {
        *num_mappings = impl1->derivative_mappings().size();
    } else if (impl3) {
        *num_mappings = impl3->derivative_mappings().size();
    } else {
        return -3; // Error: storage implementation is not Surface
    }
    return 0;
}

int sk_surface_get_derivative_mapping_name(Surface* storage, int index,
                                           const char** name) {
    if (storage == nullptr) {
        return -1; // Error: storage is null
    }
    if (storage->impl == nullptr) {
        return -2; // Error: storage implementation is null
    }

    auto* impl1 =
        dynamic_cast<sasktran2::atmosphere::Surface<1>*>(storage->impl.get());
    auto* impl3 =
        dynamic_cast<sasktran2::atmosphere::Surface<3>*>(storage->impl.get());

    if (impl1) {
        auto& mapping_list = impl1->derivative_mappings();
        auto it = mapping_list.begin();
        std::advance(it, index);

        if (it != mapping_list.end()) {
            *name = it->first.c_str();
            return 0;
        } else {
            return -3; // Error: index out of bounds
        }
    } else if (impl3) {
        auto& mapping_list = impl3->derivative_mappings();
        auto it = mapping_list.begin();
        std::advance(it, index);

        if (it != mapping_list.end()) {
            *name = it->first.c_str();
            return 0;
        } else {
            return -3; // Error: index out of bounds
        }
    } else {
        return -4; // Error: storage implementation is not Surface
    }
}

Atmosphere* sk_atmosphere_create(AtmosphereStorage* storage, Surface* surface,
                                 int calculate_derivatives,
                                 int calculate_emission_derivatives) {
    return new Atmosphere(storage, surface, calculate_derivatives == 1,
                          calculate_emission_derivatives == 1);
}

int sk_atmosphere_apply_delta_m_scaling(Atmosphere* atmosphere, int order) {
    if (atmosphere->impl) {
        auto* impl1 = dynamic_cast<sasktran2::atmosphere::Atmosphere<1>*>(
            atmosphere->impl.get());
        auto* impl3 = dynamic_cast<sasktran2::atmosphere::Atmosphere<3>*>(
            atmosphere->impl.get());

        if (impl1) {
            impl1->apply_delta_m_scaling(order);
        } else if (impl3) {
            impl3->apply_delta_m_scaling(order);
        } else {
            return -2; // Error: atmosphere implementation is not Atmosphere
        }
    } else {
        return -1; // Error: atmosphere implementation is null
    }
    return 0;
}

int sk_atmosphere_storage_finalize_scattering_derivatives(
    AtmosphereStorage* storage) {
    if (storage == nullptr) {
        return -1; // Error: storage is null
    }
    if (storage->impl == nullptr) {
        return -2; // Error: storage implementation is null
    }

    auto* impl1 =
        dynamic_cast<sasktran2::atmosphere::AtmosphereGridStorageFull<1>*>(
            storage->impl.get());
    auto* impl3 =
        dynamic_cast<sasktran2::atmosphere::AtmosphereGridStorageFull<3>*>(
            storage->impl.get());

    if (impl1) {
        impl1->finalize_scattering_derivatives(0);
    } else if (impl3) {
        impl3->finalize_scattering_derivatives(0);
    } else {
        return -3; // Error: storage implementation is not
                   // AtmosphereGridStorageFull
    }
    return 0;
}

int sk_surface_set_zero(Surface* storage) {
    if (storage == nullptr) {
        return -1; // Error: storage is null
    }
    if (storage->impl == nullptr) {
        return -2; // Error: storage implementation is null
    }

    auto* impl1 =
        dynamic_cast<sasktran2::atmosphere::Surface<1>*>(storage->impl.get());
    auto* impl3 =
        dynamic_cast<sasktran2::atmosphere::Surface<3>*>(storage->impl.get());

    if (impl1) {
        impl1->set_zero();
    } else if (impl3) {
        impl3->set_zero();
    } else {
        return -3; // Error: storage implementation is not Surface
    }
    return 0;
}

void sk_atmosphere_destroy(Atmosphere* atmosphere) { delete atmosphere; }
}
