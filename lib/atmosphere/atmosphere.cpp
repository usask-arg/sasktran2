#include <sasktran2/atmosphere/atmosphere.h>

namespace sasktran2::atmosphere {

    template <int NSTOKES>
    Atmosphere<NSTOKES>::Atmosphere(
        const std::vector<Constituent>& constituents,
        Surface<NSTOKES>&& surface, const Eigen::VectorXd& wavelengths,
        const sasktran2::Geometry1D& geometry, const sasktran2::Config& config,
        bool calculate_derivatives)
        : m_calculate_derivatives(calculate_derivatives),
          m_surface(std::move(surface)),
          m_storage(wavelengths.size(), geometry.size(),
                    config.num_singlescatter_moments()) {}

    template <int NSTOKES>
    Atmosphere<NSTOKES>::Atmosphere(
        AtmosphereGridStorageFull<NSTOKES>&& storage,
        Surface<NSTOKES>&& surface, bool calculate_derivatives)
        : m_storage(storage), m_surface(std::move(surface)),
          m_calculate_derivatives(calculate_derivatives) {}

    template <int NSTOKES>
    Atmosphere<NSTOKES>::Atmosphere(int nwavel,
                                    const sasktran2::Geometry1D& geometry,
                                    const sasktran2::Config& config,
                                    bool calculate_derivatives)
        : m_storage(nwavel, geometry.size(),
                    config.num_singlescatter_moments()),
          m_calculate_derivatives(calculate_derivatives), m_surface(nwavel) {}

    template <int NSTOKES> int Atmosphere<NSTOKES>::num_deriv() const {
        if (m_calculate_derivatives) {
            return (int)(m_storage.total_extinction.rows() *
                             (2 + m_storage.numscatderiv) +
                         m_surface.num_deriv());
        } else {
            return 0;
        }
    }

    template <int NSTOKES>
    void Atmosphere<NSTOKES>::apply_delta_m_scaling(int order) {
        if (order >= m_storage.max_stored_legendre()) {
            spdlog::warn("Trying to delta scale without the correct number of "
                         "legendre orders");
            return;
        }

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
    }

    template class Atmosphere<1>;
    template class Atmosphere<3>;

} // namespace sasktran2::atmosphere
