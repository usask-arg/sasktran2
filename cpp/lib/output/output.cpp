#include "sasktran2/atmosphere/atmosphere.h"
#include "sasktran2/config.h"
#include "sasktran2/geometry.h"
#include <sasktran2/output.h>
#include <sasktran2/math/scattering.h>

namespace sasktran2 {

    template <int NSTOKES>
    void Output<NSTOKES>::initialize(
        const sasktran2::Config& config, const sasktran2::Geometry1D& geometry,
        const std::vector<sasktran2::raytracing::TracedRay>& rays,
        const sasktran2::atmosphere::Atmosphere<NSTOKES>& atmosphere) {
        m_nlos = rays.size();
        m_nwavel = atmosphere.num_wavel();
        m_nderiv = atmosphere.num_deriv();
        m_ngeometry = atmosphere.storage().total_extinction.rows();

        m_atmosphere = &atmosphere;
        m_config = &config;

        this->resize();

        if constexpr (NSTOKES > 1) {
            m_stokes_C.resize(rays.size());
            m_stokes_S.resize(rays.size());
            m_stokes_C.setOnes();
            m_stokes_S.setZero();

            if (config.stokes_basis() ==
                sasktran2::Config::StokesBasis::solar) {
                for (int i = 0; i < rays.size(); ++i) {
                    auto CS = geometry.coordinates().stokes_standard_to_solar(
                        rays[i].observer_and_look.look_away);

                    m_stokes_C[i] = CS.first;
                    m_stokes_S[i] = CS.second;
                }
            } else if (config.stokes_basis() ==
                       sasktran2::Config::StokesBasis::observer) {
                for (int i = 0; i < rays.size(); ++i) {
                    auto CS =
                        geometry.coordinates().stokes_standard_to_observer(
                            rays[i].observer_and_look.look_away,
                            rays[i].observer_and_look.observer.position);

                    m_stokes_C[i] = CS.first;
                    m_stokes_S[i] = CS.second;
                }
            }
        }

        if (m_config->output_los_optical_depth()) {
            m_los_optical_depth.resize(m_nwavel, m_nlos);
            m_los_optical_depth.setZero();
        }
    }

    template class Output<1>;
    template class Output<3>;

} // namespace sasktran2
