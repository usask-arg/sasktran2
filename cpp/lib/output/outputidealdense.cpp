#include <sasktran2/output.h>

namespace sasktran2 {
    template <int NSTOKES> void OutputIdealDense<NSTOKES>::resize() {
        m_radiance.resize(NSTOKES * this->m_nwavel * this->m_nlos,
                          this->m_nderiv, false);
    }

    template <int NSTOKES>
    void OutputIdealDense<NSTOKES>::assign(
        const sasktran2::Dual<double, sasktran2::dualstorage::dense, NSTOKES>&
            radiance,
        int losidx, int wavelidx, int threadidx) {

        int linear_index = NSTOKES * this->m_nlos * wavelidx + NSTOKES * losidx;

        // Quantities necessary for derivative propagation

        if constexpr (NSTOKES >= 1) {
            m_radiance.value(linear_index) = radiance.value(0);
            m_radiance.deriv(linear_index, Eigen::all) =
                radiance.deriv(0, Eigen::all);
        }

        if constexpr (NSTOKES >= 3) {
            // Q/U components have a rotation
            m_radiance.value(linear_index + 1) =
                this->m_stokes_C[losidx] * radiance.value(1) -
                this->m_stokes_S[losidx] * radiance.value(2);
            m_radiance.deriv(linear_index + 1, Eigen::all) =
                this->m_stokes_C[losidx] * radiance.deriv(1, Eigen::all) -
                this->m_stokes_S[losidx] * radiance.deriv(2, Eigen::all);

            m_radiance.value(linear_index + 2) =
                this->m_stokes_S[losidx] * radiance.value(1) +
                this->m_stokes_C[losidx] * radiance.value(2);
            m_radiance.deriv(linear_index + 2, Eigen::all) =
                this->m_stokes_S[losidx] * radiance.deriv(1, Eigen::all) +
                this->m_stokes_C[losidx] * radiance.deriv(2, Eigen::all);
        }

        if constexpr (NSTOKES == 4) {
            // V component is a strict copy
            m_radiance.value(linear_index + 3) = radiance.value(3);
            m_radiance.deriv(linear_index + 3, Eigen::all) =
                radiance.deriv(3, Eigen::all);
        }
    }

    template class OutputIdealDense<1>;
    template class OutputIdealDense<3>;

} // namespace sasktran2
