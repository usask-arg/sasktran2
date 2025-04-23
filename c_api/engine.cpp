#include "sasktran2/engine.h"
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
        if (_config->impl.num_stokes() == 1) {
            impl = std::make_unique<Sasktran2<1>>(
                config->impl, geometry->impl.get(), viewing_geometry->impl);
        } else if (config->impl.num_stokes() == 3) {
            impl = std::make_unique<Sasktran2<3>>(
                config->impl, geometry->impl.get(), viewing_geometry->impl);
        } else {
            impl = nullptr;
        }
    }

    int calculate_radiance(Atmosphere* atmosphere, OutputC* output) {
        if (impl) {
            if (_config->impl.num_stokes() == 1) {
                Sasktran2<1>* impl1 = dynamic_cast<Sasktran2<1>*>(impl.get());
                impl1->calculate_radiance(
                    *static_cast<sasktran2::atmosphere::Atmosphere<1>*>(
                        atmosphere->impl.get()),
                    *static_cast<sasktran2::Output<1>*>(output->impl.get()));
                return 0;
            } else if (_config->impl.num_stokes() == 3) {
                Sasktran2<3>* impl3 = dynamic_cast<Sasktran2<3>*>(impl.get());
                impl3->calculate_radiance(
                    *static_cast<sasktran2::atmosphere::Atmosphere<3>*>(
                        atmosphere->impl.get()),
                    *static_cast<sasktran2::Output<3>*>(output->impl.get()));
                return 0;
            } else {
                return -2; // Error: invalid number of Stokes parameters
            }
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
                                 OutputC* output) {
    return engine->calculate_radiance(atmosphere, output);
}
}
