#include "output.h"
#include "sasktran2/config.h"
#include "internal_types.h"
#include <cstdio>
#include <sasktran2.h>

OutputC::OutputC(double* radiance, int nrad, int nstokes, double* flux,
                 int nflux) {
    Eigen::Map<Eigen::VectorXd> radiance_map(radiance, nrad);
    Eigen::Map<Eigen::VectorXd> flux_map(flux, nflux);

    if (nstokes == 1) {
        impl = std::make_unique<sasktran2::OutputC<1>>(radiance_map, flux_map);
    } else if (nstokes == 3) {
        impl = std::make_unique<sasktran2::OutputC<3>>(radiance_map, flux_map);
    } else {
        // Handle error case
        impl = nullptr;
    }
}

OutputJVP::OutputJVP(double* radiance, double* jvp, int nrad, int nstokes) {
    Eigen::Map<Eigen::VectorXd> radiance_map(radiance, nrad * nstokes);
    Eigen::Map<Eigen::VectorXd> jvp_map(jvp, nrad * nstokes);
    if (nstokes == 1) {
        impl = std::make_unique<sasktran2::OutputJVP<1>>(radiance_map, jvp_map);
    } else if (nstokes == 3) {
        impl = std::make_unique<sasktran2::OutputJVP<3>>(radiance_map, jvp_map);
    }
}

int OutputJVP::assign_derivative_tangent(const char* name,
                                         const double* tangent, int nparam) {
    if (impl == nullptr || tangent == nullptr || nparam < 0) {
        return -1;
    }
    Eigen::Map<const Eigen::VectorXd> tangent_map(tangent, nparam);
    if (auto* output = dynamic_cast<sasktran2::OutputJVP<1>*>(impl.get())) {
        output->set_derivative_tangent(name, tangent_map);
        return 0;
    }
    if (auto* output = dynamic_cast<sasktran2::OutputJVP<3>*>(impl.get())) {
        output->set_derivative_tangent(name, tangent_map);
        return 0;
    }
    return -2;
}

int OutputJVP::assign_surface_tangent(const char* name, const double* tangent,
                                      int nparam) {
    if (impl == nullptr || tangent == nullptr || nparam < 0) {
        return -1;
    }
    Eigen::Map<const Eigen::VectorXd> tangent_map(tangent, nparam);
    if (auto* output = dynamic_cast<sasktran2::OutputJVP<1>*>(impl.get())) {
        output->set_surface_tangent(name, tangent_map);
        return 0;
    }
    if (auto* output = dynamic_cast<sasktran2::OutputJVP<3>*>(impl.get())) {
        output->set_surface_tangent(name, tangent_map);
        return 0;
    }
    return -2;
}

OutputVJP::OutputVJP(double* radiance, const double* cotangent, int nrad,
                     int nstokes) {
    Eigen::Map<Eigen::VectorXd> radiance_map(radiance, nrad * nstokes);
    Eigen::Map<const Eigen::VectorXd> cotangent_map(cotangent, nrad * nstokes);
    if (nstokes == 1) {
        impl = std::make_unique<sasktran2::OutputVJP<1>>(radiance_map,
                                                         cotangent_map);
    } else if (nstokes == 3) {
        impl = std::make_unique<sasktran2::OutputVJP<3>>(radiance_map,
                                                         cotangent_map);
    }
}

int OutputVJP::assign_derivative_gradient(const char* name, double* gradient,
                                          int nparam) {
    if (impl == nullptr || gradient == nullptr || nparam < 0) {
        return -1;
    }
    Eigen::Map<Eigen::VectorXd> gradient_map(gradient, nparam);
    if (auto* output = dynamic_cast<sasktran2::OutputVJP<1>*>(impl.get())) {
        output->set_derivative_gradient_memory(name, gradient_map);
        return 0;
    }
    if (auto* output = dynamic_cast<sasktran2::OutputVJP<3>*>(impl.get())) {
        output->set_derivative_gradient_memory(name, gradient_map);
        return 0;
    }
    return -2;
}

int OutputVJP::assign_surface_gradient(const char* name, double* gradient,
                                       int nparam) {
    if (impl == nullptr || gradient == nullptr || nparam < 0) {
        return -1;
    }
    Eigen::Map<Eigen::VectorXd> gradient_map(gradient, nparam);
    if (auto* output = dynamic_cast<sasktran2::OutputVJP<1>*>(impl.get())) {
        output->set_surface_gradient_memory(name, gradient_map);
        return 0;
    }
    if (auto* output = dynamic_cast<sasktran2::OutputVJP<3>*>(impl.get())) {
        output->set_surface_gradient_memory(name, gradient_map);
        return 0;
    }
    return -2;
}

int OutputVJP::finalize() {
    if (impl == nullptr) {
        return -1;
    }
    if (auto* output = dynamic_cast<sasktran2::OutputVJP<1>*>(impl.get())) {
        output->finalize();
        return 0;
    }
    if (auto* output = dynamic_cast<sasktran2::OutputVJP<3>*>(impl.get())) {
        output->finalize();
        return 0;
    }
    return -2;
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

int OutputC::assign_flux_derivative_memory(const char* name,
                                           double* derivative_mapping, int nrad,
                                           int nderiv) {
    if (impl == nullptr) {
        return -1; // Error: Output not initialized
    }

    // Memory structure is (nrad , nderiv)
    Eigen::Map<Eigen::MatrixXd> derivative_map(derivative_mapping, nrad,
                                               nderiv);

    auto* impl1 = dynamic_cast<sasktran2::OutputC<1>*>(impl.get());
    auto* impl3 = dynamic_cast<sasktran2::OutputC<3>*>(impl.get());

    if (impl1) {
        impl1->set_flux_derivative_mapping_memory(name, derivative_map);
        return 0;
    } else if (impl3) {
        impl3->set_flux_derivative_mapping_memory(name, derivative_map);
        return 0;
    } else {
        // Handle error case
        return -1;
    }
}

int OutputC::assign_surface_flux_derivative_memory(const char* name,
                                                   double* derivative_mapping,
                                                   int nrad) {
    if (impl == nullptr) {
        return -1; // Error: Output not initialized
    }

    // Memory structure is (nrad * nstokes, 1)
    Eigen::Map<Eigen::MatrixXd> derivative_map(derivative_mapping, nrad, 1);

    auto* impl1 = dynamic_cast<sasktran2::OutputC<1>*>(impl.get());
    auto* impl3 = dynamic_cast<sasktran2::OutputC<3>*>(impl.get());

    if (impl1) {
        impl1->set_flux_surface_derivative_mapping_memory(name, derivative_map);
        return 0;
    } else if (impl3) {
        impl3->set_flux_surface_derivative_mapping_memory(name, derivative_map);
        return 0;
    } else {
        // Handle error case
        return -1;
    }
}

extern "C" {
OutputC* sk_output_create(double* radiance, int nrad, int nstokes, double* flux,
                          int nflux) {
    return new OutputC(radiance, nrad, nstokes, flux, nflux);
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

int sk_output_assign_flux_derivative_memory(OutputC* output, const char* name,
                                            double* derivative_mapping,
                                            int nrad, int nderiv) {
    if (output->impl == nullptr) {
        return -1; // Error: Output not initialized
    }

    return output->assign_flux_derivative_memory(name, derivative_mapping, nrad,
                                                 nderiv);
}

int sk_output_assign_surface_flux_derivative_memory(OutputC* output,
                                                    const char* name,
                                                    double* derivative_mapping,
                                                    int nrad) {
    if (output->impl == nullptr) {
        return -1; // Error: Output not initialized
    }

    return output->assign_surface_flux_derivative_memory(
        name, derivative_mapping, nrad);
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

OutputJVP* sk_output_jvp_create(double* radiance, double* jvp, int nrad,
                                int nstokes) {
    return new OutputJVP(radiance, jvp, nrad, nstokes);
}

void sk_output_jvp_destroy(OutputJVP* output) { delete output; }

int sk_output_jvp_assign_derivative_tangent(OutputJVP* output, const char* name,
                                            const double* tangent, int nparam) {
    return output == nullptr
               ? -1
               : output->assign_derivative_tangent(name, tangent, nparam);
}

int sk_output_jvp_assign_surface_tangent(OutputJVP* output, const char* name,
                                         const double* tangent, int nparam) {
    return output == nullptr
               ? -1
               : output->assign_surface_tangent(name, tangent, nparam);
}

OutputVJP* sk_output_vjp_create(double* radiance, const double* cotangent,
                                int nrad, int nstokes) {
    return new OutputVJP(radiance, cotangent, nrad, nstokes);
}

void sk_output_vjp_destroy(OutputVJP* output) { delete output; }

int sk_output_vjp_assign_derivative_gradient(OutputVJP* output,
                                             const char* name, double* gradient,
                                             int nparam) {
    return output == nullptr
               ? -1
               : output->assign_derivative_gradient(name, gradient, nparam);
}

int sk_output_vjp_assign_surface_gradient(OutputVJP* output, const char* name,
                                          double* gradient, int nparam) {
    return output == nullptr
               ? -1
               : output->assign_surface_gradient(name, gradient, nparam);
}

int sk_output_vjp_finalize(OutputVJP* output) {
    return output == nullptr ? -1 : output->finalize();
}
}
