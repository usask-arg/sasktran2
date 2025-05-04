#include "brdf.h"
#include "internal_types.h"
#include "sasktran2.h"

extern "C" {

BRDF* sk_brdf_create_lambertian(int nstokes) {
    BRDF* brdf = new BRDF();

    if (nstokes == 1) {
        brdf->impl =
            std::make_shared<sasktran2::atmosphere::brdf::Lambertian<1>>();
    } else if (nstokes == 3) {
        brdf->impl =
            std::make_shared<sasktran2::atmosphere::brdf::Lambertian<3>>();
    } else {
        delete brdf;
        return nullptr;
    }

    return brdf;
}

BRDF* sk_brdf_create_kokhanovsky(int nstokes) {
    BRDF* brdf = new BRDF();

    if (nstokes == 1) {
        brdf->impl =
            std::make_shared<sasktran2::atmosphere::brdf::SnowKokhanovsky<1>>();
    } else if (nstokes == 3) {
        brdf->impl =
            std::make_shared<sasktran2::atmosphere::brdf::SnowKokhanovsky<3>>();
    } else {
        delete brdf;
        return nullptr;
    }

    return brdf;
}

BRDF* sk_brdf_create_modis(int nstokes) {
    BRDF* brdf = new BRDF();

    if (nstokes == 1) {
        brdf->impl = std::make_shared<sasktran2::atmosphere::brdf::MODIS<1>>();
    } else if (nstokes == 3) {
        brdf->impl = std::make_shared<sasktran2::atmosphere::brdf::MODIS<3>>();
    } else {
        delete brdf;
        return nullptr;
    }

    return brdf;
}

void sk_brdf_destroy(BRDF* brdf) { delete brdf; }

int sk_brdf_get_num_deriv(BRDF* config, int* num_deriv) {
    if (config == nullptr || num_deriv == nullptr) {
        return -1; // Error: invalid arguments
    }

    auto* impl1 =
        dynamic_cast<sasktran2::atmosphere::brdf::BRDF<1>*>(config->impl.get());
    auto* impl3 =
        dynamic_cast<sasktran2::atmosphere::brdf::BRDF<3>*>(config->impl.get());

    if (impl1) {
        *num_deriv = impl1->num_deriv();
        return 0;
    }

    if (impl3) {
        *num_deriv = impl3->num_deriv();
        return 0;
    }

    return -2; // Error: invalid BRDF type
}

int sk_brdf_get_num_args(BRDF* config, int* num_args) {
    if (config == nullptr || num_args == nullptr) {
        return -1; // Error: invalid arguments
    }

    auto* impl1 =
        dynamic_cast<sasktran2::atmosphere::brdf::BRDF<1>*>(config->impl.get());
    auto* impl3 =
        dynamic_cast<sasktran2::atmosphere::brdf::BRDF<3>*>(config->impl.get());

    if (impl1) {
        *num_args = impl1->num_args();
        return 0;
    }

    if (impl3) {
        *num_args = impl3->num_args();
        return 0;
    }

    return -2; // Error: invalid BRDF type
}
}
