#include "sasktran2/atmosphere/grid_storage.h"
#include "sasktran2/config.h"
#include <sasktran2/atmosphere/atmosphere.h>

namespace sasktran2::atmosphere {

    template <int NSTOKES>
    Atmosphere<NSTOKES>::Atmosphere(
        AtmosphereGridStorageFull<NSTOKES>&& storage,
        Surface<NSTOKES>&& surface, bool calculate_derivatives)
        : m_storage_holder(std::make_shared<AtmosphereGridStorageFull<NSTOKES>>(
              std::move(storage))),
          m_surface_holder(
              std::make_shared<Surface<NSTOKES>>(std::move(surface))),
          m_storage(*m_storage_holder), m_surface(*m_surface_holder),
          m_calculate_derivatives(calculate_derivatives),
          m_include_emission_derivatives(false) {}

    template <int NSTOKES>
    Atmosphere<NSTOKES>::Atmosphere(AtmosphereGridStorageFull<NSTOKES>& stor,
                                    Surface<NSTOKES>& surf,
                                    bool calculate_derivatives,
                                    bool include_emission_derivatives)
        : m_storage(stor), m_surface(surf),
          m_calculate_derivatives(calculate_derivatives),
          m_include_emission_derivatives(include_emission_derivatives) {}

    template <int NSTOKES>
    Atmosphere<NSTOKES>::Atmosphere(int nwavel,
                                    const sasktran2::Geometry1D& geometry,
                                    const sasktran2::Config& config,
                                    bool calculate_derivatives)
        : m_storage_holder(std::make_shared<AtmosphereGridStorageFull<NSTOKES>>(
              nwavel, geometry.size(), config.num_singlescatter_moments())),
          m_calculate_derivatives(calculate_derivatives),
          m_surface_holder(std::make_shared<Surface<NSTOKES>>(nwavel)),
          m_storage(*m_storage_holder), m_surface(*m_surface_holder),
          m_include_emission_derivatives(
              (config.emission_source() !=
               sasktran2::Config::EmissionSource::none) &&
              calculate_derivatives) {}

    template <int NSTOKES> int Atmosphere<NSTOKES>::num_deriv() const {
        if (m_calculate_derivatives) {
            return (int)(m_storage.total_extinction.rows() *
                             (2 + num_source_deriv_groups()) +
                         m_surface.num_deriv() +
                         int(m_include_emission_derivatives));
        } else {
            return 0;
        }
    }

    template <int NSTOKES> int Atmosphere<NSTOKES>::num_output_deriv() const {
        int n = 0;

        if (m_calculate_derivatives) {
            for (auto& [name, mapping] :
                 m_storage.derivative_mappings_const()) {
                n += mapping.num_output();
            }
        }

        return n;
    }

    template <int NSTOKES>
    void Atmosphere<NSTOKES>::apply_delta_m_scaling(int order) {
        if (order >= m_storage.max_stored_legendre()) {
            spdlog::warn("Trying to delta scale without the correct number of "
                         "legendre orders");
            return;
        }

        // Copy the unscaled extinction/ssa
        m_storage.unscaled_total_extinction = m_storage.total_extinction;
        m_storage.unscaled_ssa = m_storage.ssa;

        int stokes_stack;
        if (NSTOKES == 1) {
            stokes_stack = 1;
        } else if (NSTOKES == 3) {
            stokes_stack = 4;
        }

        m_storage.applied_f_location = order * stokes_stack;
        m_storage.applied_f_order = order;

        // First copy the requested order into the f storage
        // TODO: Figure out how to do this with eigen operations?
        for (int i = 0; i < m_storage.f.rows(); ++i) {
            for (int j = 0; j < m_storage.f.cols(); ++j) {
                m_storage.f(i, j) =
                    m_storage.leg_coeff(m_storage.applied_f_location, i, j) /
                    (2 * order + 1);
            }
        }

        // Now scale k* = (1-wf) k
        m_storage.total_extinction.array() =
            m_storage.total_extinction.array() *
            (1 - m_storage.ssa.array() *
                     m_storage.f.template cast<double>().array());

        // And then w* = (1 - f) / (1-wf) * w
        m_storage.ssa.array() =
            (1 - m_storage.f.template cast<double>().array()) /
            (1 - m_storage.ssa.array() *
                     m_storage.f.template cast<double>().array()) *
            m_storage.ssa.array();

        // Start by assignind d_f
        for (int w = 0; w < m_storage.f.cols(); ++w) {     // wavel loop
            for (int i = 0; i < m_storage.f.rows(); ++i) { // location loop
                for (int d = 0; d < m_storage.numscatderiv; ++d) { // deriv loop
                    m_storage.d_f(i, w, d) =
                        m_storage.d_leg_coeff(m_storage.applied_f_location, i,
                                              w, d) /
                        (2 * order + 1);
                }
            }
        }

        // Lastly we need to scale the legendre coefficients, b = b / (1-f)
        // Note that this is the scaling for single scatter TMS correction, we
        // still have to subtract f/1-f for multiple scatter
        for (int w = 0; w < m_storage.f.cols(); ++w) {     // wavel loop
            for (int i = 0; i < m_storage.f.rows(); ++i) { // location loop
                for (int j = 0; j < m_storage.max_stored_legendre();
                     ++j) { // legendre loop
                    m_storage.leg_coeff(j, i, w) /= (1 - m_storage.f(i, w));

                    // Also have to scale the derivatives
                    for (int d = 0; d < m_storage.numscatderiv;
                         ++d) { // deriv loop

                        // db* = db / 1-f + (b*) * df / (1-f)
                        m_storage.d_leg_coeff(j, i, w, d) +=
                            m_storage.leg_coeff(j, i, w) *
                            m_storage.d_f(i, w, d);
                        m_storage.d_leg_coeff(j, i, w, d) /=
                            1 - m_storage.f(i, w);
                    }
                }
            }
        }
        // Now we can scale the derivative mappings based on our new values

        for (auto& [name, mapping] : m_storage.derivative_mappings()) {
            // Start with d_k
            if (mapping.native_mapping().d_extinction.has_value()) {
                // Step 1, scale the extinction
                mapping.native_mapping().d_extinction.value().array() *=
                    (1 - m_storage.unscaled_ssa.array() *
                             m_storage.f.template cast<double>().array());

                // Step 2, subtract how the term due to d_ssa
                mapping.native_mapping().d_extinction.value().array() -=
                    m_storage.unscaled_total_extinction.array() *
                    m_storage.f.template cast<double>().array() *
                    mapping.native_mapping().d_ssa.value().array();

                // and d_ssa
                // Only have to scale d_ssa
                mapping.native_mapping().d_ssa.value().array() *=
                    (1 - m_storage.f.template cast<double>().array() *
                             (1 - m_storage.ssa.array())) /
                    (1 - m_storage.unscaled_ssa.array() *
                             m_storage.f.template cast<double>().array());
            }

            // If we are a scattering derivative, we also have to calculate the
            // changes in d_k and d_ssa due to d_f
            if (mapping.is_scattering_derivative()) {
                // Now we have to scale the d_k and d_ssa
                if (mapping.native_mapping().d_extinction.has_value()) {
                    // TODO: Rewrite this with Eigen ops
                    for (int w = 0; w < m_storage.f.cols(); ++w) { // wavel loop
                        for (int i = 0; i < m_storage.f.rows();
                             ++i) { // location loop
                            double d_f =
                                m_storage.d_f(i, w,
                                              mapping.get_scattering_index()) *
                                mapping.native_mapping().scat_factor.value()(i,
                                                                             w);

                            mapping.native_mapping().d_extinction.value()(i,
                                                                          w) -=
                                m_storage.unscaled_ssa(i, w) *
                                m_storage.unscaled_total_extinction(i, w) * d_f;

                            mapping.native_mapping().d_ssa.value()(i, w) +=
                                d_f * m_storage.unscaled_ssa(i, w) /
                                (1 - m_storage.unscaled_ssa(i, w) *
                                         m_storage.f(i, w)) *
                                (m_storage.ssa(i, w) - 1);
                        }
                    }
                }
            }
        }
    }

    template class Atmosphere<1>;
    template class Atmosphere<3>;

} // namespace sasktran2::atmosphere
