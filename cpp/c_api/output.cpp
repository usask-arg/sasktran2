#include "output.h"
#include "sasktran2/config.h"
#include "internal_types.h"
#include <cstdio>
#include <sasktran2.h>

OutputC::OutputC(double* radiance, int nrad, int nstokes) {
    Eigen::Map<Eigen::VectorXd> radiance_map(radiance, nrad);

    if (nstokes == 1) {
        impl = std::make_unique<sasktran2::OutputC<1>>(radiance_map);
    } else if (nstokes == 3) {
        impl = std::make_unique<sasktran2::OutputC<3>>(radiance_map);
    } else {
        // Handle error case
        impl = nullptr;
    }
}

int OutputC::assign_derivative_memory(const char* name,
                                      double* derivative_mapping, int nrad,
                                      int nstokes, int nderiv) {
    if (impl == nullptr) {
        return -1; // Error: Output not initialized
    }

    // Memory structure is (nrad * nstokes, nderiv)
    Eigen::Map<Eigen::MatrixXd> derivative_map(derivative_mapping,
                                               nrad * nstokes, nderiv);

    auto* impl1 = dynamic_cast<sasktran2::OutputC<1>*>(impl.get());
    auto* impl3 = dynamic_cast<sasktran2::OutputC<3>*>(impl.get());

    if (impl1) {
        impl1->set_derivative_mapping_memory(name, derivative_map);
        return 0;
    } else if (impl3) {
        impl3->set_derivative_mapping_memory(name, derivative_map);
        return 0;
    } else {
        // Handle error case
        return -1;
    }
}

int OutputC::assign_surface_derivative_memory(const char* name,
                                              double* derivative_mapping,
                                              int nrad, int nstokes) {
    if (impl == nullptr) {
        return -1; // Error: Output not initialized
    }

    // Memory structure is (nrad * nstokes, 1)
    Eigen::Map<Eigen::MatrixXd> derivative_map(derivative_mapping,
                                               nrad * nstokes, 1);

    auto* impl1 = dynamic_cast<sasktran2::OutputC<1>*>(impl.get());
    auto* impl3 = dynamic_cast<sasktran2::OutputC<3>*>(impl.get());

    if (impl1) {
        impl1->set_surface_derivative_mapping_memory(name, derivative_map);
        return 0;
    } else if (impl3) {
        impl3->set_surface_derivative_mapping_memory(name, derivative_map);
        return 0;
    } else {
        // Handle error case
        return -1;
    }
}

extern "C" {
OutputC* sk_output_create(double* radiance, int nrad, int nstokes) {
    return new OutputC(radiance, nrad, nstokes);
}

void sk_output_destroy(OutputC* output) { delete output; }

int sk_output_assign_derivative_memory(OutputC* output, const char* name,
                                       double* derivative_mapping, int nrad,
                                       int nstokes, int nderiv) {
    if (output->impl == nullptr) {
        return -1; // Error: Output not initialized
    }

    return output->assign_derivative_memory(name, derivative_mapping, nrad,
                                            nstokes, nderiv);
}

int sk_output_assign_surface_derivative_memory(OutputC* output,
                                               const char* name,
                                               double* derivative_mapping,
                                               int nrad, int nstokes) {
    if (output->impl == nullptr) {
        return -1; // Error: Output not initialized
    }

    return output->assign_surface_derivative_memory(name, derivative_mapping,
                                                    nrad, nstokes);
}

int sk_output_get_los_optical_depth(OutputC* output, double** od) {
    if (output->impl == nullptr) {
        return -1; // Error: Output not initialized
    }

    auto* impl1 = dynamic_cast<sasktran2::OutputC<1>*>(output->impl.get());
    auto* impl3 = dynamic_cast<sasktran2::OutputC<3>*>(output->impl.get());

    if (impl1) {
        *od = impl1->los_optical_depth().data();
        return 0;
    } else if (impl3) {
        *od = impl3->los_optical_depth().data();
        return 0;
    } else {
        // Handle error case
        return -1;
    }
}
}
