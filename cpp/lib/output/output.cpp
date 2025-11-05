#include "sasktran2/atmosphere/atmosphere.h"
#include "sasktran2/config.h"
#include "sasktran2/geometry.h"
#include "sasktran2/viewinggeometry_internal.h"
#include <sasktran2/output.h>
#include <sasktran2/math/scattering.h>

namespace sasktran2 {

    template <int NSTOKES>
    void Output<NSTOKES>::initialize(
        const sasktran2::Config& config, const sasktran2::Geometry1D& geometry,
        const sasktran2::viewinggeometry::InternalViewingGeometry&
            internal_viewing,
        const sasktran2::atmosphere::Atmosphere<NSTOKES>& atmosphere) {
        m_nlos = internal_viewing.traced_rays.size();
        m_nfluxpos = internal_viewing.flux_observers.size();
        m_nfluxtype = (int)config.get_flux_types().size();
        m_nwavel = atmosphere.num_wavel();
        m_nderiv = atmosphere.num_deriv();
        m_ngeometry = atmosphere.storage().total_extinction.rows();

        m_atmosphere = &atmosphere;
        m_config = &config;

        this->resize();

        if constexpr (NSTOKES > 1) {
            m_stokes_C.resize(m_nlos);
            m_stokes_S.resize(m_nlos);
            m_stokes_C.setOnes();
            m_stokes_S.setZero();

            if (config.stokes_basis() ==
                sasktran2::Config::StokesBasis::solar) {
                for (int i = 0; i < m_nlos; ++i) {
                    auto CS = geometry.coordinates().stokes_standard_to_solar(
                        internal_viewing.traced_rays[i]
                            .observer_and_look.look_away);

                    m_stokes_C[i] = CS.first;
                    m_stokes_S[i] = CS.second;
                }
            } else if (config.stokes_basis() ==
                       sasktran2::Config::StokesBasis::observer) {
                for (int i = 0; i < m_nlos; ++i) {
                    auto CS =
                        geometry.coordinates().stokes_standard_to_observer(
                            internal_viewing.traced_rays[i]
                                .observer_and_look.look_away,
                            internal_viewing.traced_rays[i]
                                .observer_and_look.observer.position);

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
