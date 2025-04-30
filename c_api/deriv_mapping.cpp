#include "deriv_mapping.h"
#include "internal_types.h"
#include "sasktran2/derivative_mapping.h"

extern "C" {

int sk_deriv_mapping_destroy(DerivativeMapping* mapping) {
    if (mapping == nullptr) {
        return -1; // Error: mapping is null
    }
    // Don't need to delete impl because it is owned by the atmosphere storage
    delete mapping;
    return 0;
}

int sk_deriv_mapping_set_zero(DerivativeMapping* mapping) {
    if (mapping == nullptr) {
        return -1; // Error: mapping is null
    }
    mapping->impl->set_zero();
    return 0;
}

int sk_deriv_mapping_get_d_ssa(DerivativeMapping* mapping, double** ssa) {
    if (mapping == nullptr) {
        return -1; // Error: mapping is null
    }

    if (!mapping->impl->native_mapping().d_ssa.has_value()) {
        mapping->impl->allocate_ssa_derivatives();
    }
    *ssa = mapping->impl->native_mapping().d_ssa.value().data();
    return 0;
}

int sk_deriv_mapping_get_d_extinction(DerivativeMapping* mapping,
                                      double** extinction) {
    if (mapping == nullptr) {
        return -1; // Error: mapping is null
    }

    if (!mapping->impl->native_mapping().d_extinction.has_value()) {
        mapping->impl->allocate_extinction_derivatives();
    }
    *extinction = mapping->impl->native_mapping().d_extinction.value().data();
    return 0;
}

int sk_deriv_mapping_get_scat_factor(DerivativeMapping* mapping,
                                     double** scat_factor) {
    if (mapping == nullptr) {
        return -1; // Error: mapping is null
    }

    if (!mapping->impl->native_mapping().scat_factor.has_value()) {
        mapping->impl->allocate_legendre_derivatives();
    }
    *scat_factor = mapping->impl->native_mapping().scat_factor.value().data();
    return 0;
}

int sk_deriv_mapping_get_d_legendre(DerivativeMapping* mapping,
                                    double** d_legendre) {
    if (mapping == nullptr) {
        return -1; // Error: mapping is null
    }

    if (!mapping->impl->native_mapping().d_legendre.has_value()) {
        mapping->impl->allocate_legendre_derivatives();
    }
    *d_legendre = mapping->impl->native_mapping().d_legendre.value().data();
    return 0;
}

int sk_deriv_mapping_get_d_emission(DerivativeMapping* mapping,
                                    double** d_emission) {
    if (mapping == nullptr) {
        return -1; // Error: mapping is null
    }

    if (!mapping->impl->native_mapping().d_emission.has_value()) {
        mapping->impl->allocate_emission_derivatives();
    }
    *d_emission = mapping->impl->native_mapping().d_emission.value().data();
    return 0;
}

int sk_deriv_mapping_get_scat_deriv_index(DerivativeMapping* mapping,
                                          int* scat_deriv_index) {
    if (mapping == nullptr) {
        return -1; // Error: mapping is null
    }

    *scat_deriv_index = mapping->impl->get_scattering_index();
    return 0;
}

int sk_deriv_mapping_set_scat_deriv_index(DerivativeMapping* mapping,
                                          int scat_deriv_index) {
    if (mapping == nullptr) {
        return -1; // Error: mapping is null
    }

    mapping->impl->set_scattering_index(scat_deriv_index);
    return 0;
}

int sk_deriv_mapping_set_interp_dim(DerivativeMapping* mapping,
                                    const char* name) {
    if (mapping == nullptr) {
        return -1; // Error: mapping is null
    }

    mapping->impl->set_interp_dim(name);
    return 0;
}

int sk_deriv_mapping_set_assign_name(DerivativeMapping* mapping,
                                     const char* name) {
    if (mapping == nullptr) {
        return -1; // Error: mapping is null
    }

    mapping->impl->set_assign_name(name);
    return 0;
}

int sk_deriv_mapping_set_log_radiance_space(DerivativeMapping* mapping,
                                            int log_radiance_space) {
    if (mapping == nullptr) {
        return -1; // Error: mapping is null
    }

    mapping->impl->set_log_radiance_space(log_radiance_space != 0);
    return 0;
}

int sk_deriv_mapping_is_scattering_derivative(DerivativeMapping* mapping,
                                              int* is_scattering_derivative) {
    if (mapping == nullptr) {
        return -1; // Error: mapping is null
    }

    *is_scattering_derivative =
        mapping->impl->is_scattering_derivative() ? 1 : 0;
    return 0;
}

int sk_deriv_mapping_get_num_legendre(DerivativeMapping* mapping,
                                      int* num_legendre) {
    if (mapping == nullptr) {
        return -1; // Error: mapping is null
    }

    *num_legendre = mapping->impl->num_legendre();
    return 0;
}

int sk_deriv_mapping_get_num_location(DerivativeMapping* mapping,
                                      int* num_location) {
    if (mapping == nullptr) {
        return -1; // Error: mapping is null
    }

    *num_location = mapping->impl->num_location();
    return 0;
}

int sk_deriv_mapping_get_num_wavel(DerivativeMapping* mapping, int* num_wavel) {
    if (mapping == nullptr) {
        return -1; // Error: mapping is null
    }

    *num_wavel = mapping->impl->num_wavel();
    return 0;
}

int sk_deriv_mapping_get_num_output(DerivativeMapping* mapping,
                                    int* num_output) {
    if (mapping == nullptr) {
        return -1; // Error: mapping is null
    }

    *num_output = mapping->impl->num_output();
    return 0;
}

int sk_deriv_mapping_get_assign_name(DerivativeMapping* mapping,
                                     const char** name) {
    if (mapping == nullptr) {
        return -1; // Error: mapping is null
    }

    *name = mapping->impl->get_assign_name().c_str();
    return 0;
}

int sk_deriv_mapping_get_interp_dim(DerivativeMapping* mapping,
                                    const char** name) {
    if (mapping == nullptr) {
        return -1; // Error: mapping is null
    }

    *name = mapping->impl->get_interp_dim().c_str();
    return 0;
}

int sk_deriv_mapping_set_interpolator(DerivativeMapping* mapping,
                                      double* interpolator, int dim1,
                                      int dim2) {
    if (mapping == nullptr) {
        return -1; // Error: mapping is null
    }

    Eigen::MatrixXd interpolator_matrix =
        Eigen::Map<Eigen::MatrixXd>(interpolator, dim1, dim2);

    mapping->impl->set_interpolator(interpolator_matrix);
    return 0;
}

int sk_deriv_mapping_get_interpolator(DerivativeMapping* mapping,
                                      double** interpolator, int* dim1,
                                      int* dim2) {
    if (mapping == nullptr) {
        return -1; // Error: mapping is null
    }
    if (!mapping->impl->get_interpolator().has_value()) {
        return -2; // Error: interpolator is not set
    }

    Eigen::MatrixXd& interpolator_matrix =
        mapping->impl->get_interpolator().value();
    *dim1 = interpolator_matrix.rows();
    *dim2 = interpolator_matrix.cols();
    *interpolator = interpolator_matrix.data();

    return 0;
}

int sk_surface_deriv_mapping_get_d_emission(SurfaceDerivativeMapping* mapping,
                                            double** emission) {
    if (mapping == nullptr) {
        return -1; // Error: mapping is null
    }

    if (!mapping->impl->native_surface_mapping().d_emission.has_value()) {
        mapping->impl->allocate_emission_derivatives();
    }
    *emission =
        mapping->impl->native_surface_mapping().d_emission.value().data();
    return 0;
}

int sk_surface_deriv_mapping_get_d_brdf(SurfaceDerivativeMapping* mapping,
                                        double** brdf) {
    if (mapping == nullptr) {
        return -1; // Error: mapping is null
    }

    if (!mapping->impl->native_surface_mapping().d_brdf.has_value()) {
        mapping->impl->allocate_brdf_derivatives();
    }
    *brdf = mapping->impl->native_surface_mapping().d_brdf.value().data();
    return 0;
}

int sk_surface_deriv_mapping_set_interpolator(SurfaceDerivativeMapping* mapping,
                                              double* interpolator, int dim1,
                                              int dim2) {
    if (mapping == nullptr) {
        return -1; // Error: mapping is null
    }

    Eigen::MatrixXd interpolator_matrix =
        Eigen::Map<Eigen::MatrixXd>(interpolator, dim1, dim2);
    mapping->impl->set_interpolator(interpolator_matrix);
    return 0;
}

int sk_surface_deriv_mapping_get_interpolator(SurfaceDerivativeMapping* mapping,
                                              double** interpolator, int* dim1,
                                              int* dim2) {
    if (mapping == nullptr) {
        return -1; // Error: mapping is null
    }
    if (!mapping->impl->get_interpolator().has_value()) {
        return -2; // Error: interpolator is not set
    }

    Eigen::MatrixXd& interpolator_matrix =
        mapping->impl->get_interpolator().value();
    *dim1 = interpolator_matrix.rows();
    *dim2 = interpolator_matrix.cols();
    *interpolator = interpolator_matrix.data();

    return 0;
}

int sk_surface_deriv_mapping_get_interp_dim(SurfaceDerivativeMapping* mapping,
                                            const char** name) {
    if (mapping == nullptr) {
        return -1; // Error: mapping is null
    }

    *name = mapping->impl->get_interp_dim().c_str();
    return 0;
}

int sk_surface_deriv_mapping_set_interp_dim(SurfaceDerivativeMapping* mapping,
                                            const char* name) {
    if (mapping == nullptr) {
        return -1; // Error: mapping is null
    }

    mapping->impl->set_interp_dim(name);
    return 0;
}

int sk_surface_deriv_mapping_set_zero(SurfaceDerivativeMapping* mapping) {
    if (mapping == nullptr) {
        return -1; // Error: mapping is null
    }

    mapping->impl->set_zero();
    return 0;
}

int sk_surface_deriv_mapping_get_num_wavel(SurfaceDerivativeMapping* mapping,
                                           int* num_wavel) {
    if (mapping == nullptr) {
        return -1; // Error: mapping is null
    }

    *num_wavel = mapping->impl->num_wavel();
    return 0;
}

int sk_surface_deriv_mapping_get_num_brdf_args(
    SurfaceDerivativeMapping* mapping, int* num_brdf_args) {
    if (mapping == nullptr) {
        return -1; // Error: mapping is null
    }

    *num_brdf_args = mapping->impl->num_brdf_args();
    return 0;
}

int sk_surface_deriv_mapping_destroy(SurfaceDerivativeMapping* mapping) {
    if (mapping == nullptr) {
        return -1; // Error: mapping is null
    }
    // Don't need to delete impl because it is owned by the atmosphere storage
    delete mapping;
    return 0;
}
}
