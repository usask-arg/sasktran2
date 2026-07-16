#include <sasktran2/output.h>

namespace sasktran2 {
    template <int NSTOKES> void OutputIdealDense<NSTOKES>::resize() {
        m_radiance.resize(NSTOKES * this->m_nwavel * this->m_nlos,
                          this->m_nderiv, false);
    }

    template <int NSTOKES>
    template <typename Radiance>
    void OutputIdealDense<NSTOKES>::assign_lane(const Radiance& radiance,
                                                int losidx, int wavelidx) {

        int linear_index = NSTOKES * this->m_nlos * wavelidx + NSTOKES * losidx;

        // Quantities necessary for derivative propagation

        if constexpr (NSTOKES >= 1) {
            m_radiance.value(linear_index) = radiance.value(0);
            m_radiance.deriv(linear_index, Eigen::placeholders::all) =
                radiance.deriv(0, Eigen::placeholders::all);
        }

        if constexpr (NSTOKES >= 3) {
            // Q/U components have a rotation
            m_radiance.value(linear_index + 1) =
                this->m_stokes_C[losidx] * radiance.value(1) -
                this->m_stokes_S[losidx] * radiance.value(2);
            m_radiance.deriv(linear_index + 1, Eigen::placeholders::all) =
                this->m_stokes_C[losidx] *
                    radiance.deriv(1, Eigen::placeholders::all) -
                this->m_stokes_S[losidx] *
                    radiance.deriv(2, Eigen::placeholders::all);

            m_radiance.value(linear_index + 2) =
                this->m_stokes_S[losidx] * radiance.value(1) +
                this->m_stokes_C[losidx] * radiance.value(2);
            m_radiance.deriv(linear_index + 2, Eigen::placeholders::all) =
                this->m_stokes_S[losidx] *
                    radiance.deriv(1, Eigen::placeholders::all) +
                this->m_stokes_C[losidx] *
                    radiance.deriv(2, Eigen::placeholders::all);
        }

        if constexpr (NSTOKES == 4) {
            // V component is a strict copy
            m_radiance.value(linear_index + 3) = radiance.value(3);
            m_radiance.deriv(linear_index + 3, Eigen::placeholders::all) =
                radiance.deriv(3, Eigen::placeholders::all);
        }
    }

    template <int NSTOKES>
    void OutputIdealDense<NSTOKES>::assign(
        const sasktran2::WavelengthBlock<>& block,
        const sasktran2::WavelengthBlockDual<NSTOKES>& radiance, int losidx,
        int) {
        if (block.count == 1) {
            const sasktran2::WavelengthBlockConstLaneDualView<NSTOKES>
                radiance_lane(radiance, 0);
            assign_lane(radiance_lane, losidx, block.start);
            return;
        }

        const int stride = NSTOKES * this->m_nlos;
        const int first = stride * block.start + NSTOKES * losidx;
        const auto output_indices = [&](int stokes) {
            return Eigen::seqN(first + stokes, block.count, stride);
        };

        m_radiance.value(output_indices(0)) =
            radiance.value.row(0).head(block.count).transpose();
        if (this->m_nderiv > 0) {
            m_radiance.deriv(output_indices(0), Eigen::placeholders::all) =
                radiance.derivative_stokes(0, this->m_nderiv, 0, block.count)
                    .transpose();
        }

        if constexpr (NSTOKES >= 3) {
            const double stokes_c = this->m_stokes_C[losidx];
            const double stokes_s = this->m_stokes_S[losidx];
            m_radiance.value(output_indices(1)) =
                (stokes_c * radiance.value.row(1).head(block.count) -
                 stokes_s * radiance.value.row(2).head(block.count))
                    .transpose();
            m_radiance.value(output_indices(2)) =
                (stokes_s * radiance.value.row(1).head(block.count) +
                 stokes_c * radiance.value.row(2).head(block.count))
                    .transpose();
            if (this->m_nderiv > 0) {
                const auto q = radiance.derivative_stokes(0, this->m_nderiv, 1,
                                                          block.count);
                const auto u = radiance.derivative_stokes(0, this->m_nderiv, 2,
                                                          block.count);
                m_radiance.deriv(output_indices(1), Eigen::placeholders::all) =
                    (stokes_c * q - stokes_s * u).transpose();
                m_radiance.deriv(output_indices(2), Eigen::placeholders::all) =
                    (stokes_s * q + stokes_c * u).transpose();
            }
        }

        if constexpr (NSTOKES == 4) {
            m_radiance.value(output_indices(3)) =
                radiance.value.row(3).head(block.count).transpose();
            if (this->m_nderiv > 0) {
                m_radiance.deriv(output_indices(3), Eigen::placeholders::all) =
                    radiance
                        .derivative_stokes(0, this->m_nderiv, 3, block.count)
                        .transpose();
            }
        }
    }

    template class OutputIdealDense<1>;
    template class OutputIdealDense<3>;

} // namespace sasktran2
