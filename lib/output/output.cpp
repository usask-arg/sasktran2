#include "sasktran2/config.h"
#include "sasktran2/geometry.h"
#include <sasktran2/output.h>
#include <sasktran2/math/scattering.h>

namespace sasktran2 {

    template <int NSTOKES>
    void Output<NSTOKES>::initialize(
        const sasktran2::Config& config, const sasktran2::Geometry1D& geometry,
        const std::vector<sasktran2::raytracing::TracedRay>& rays) {
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
    }

    template class Output<1>;
    template class Output<3>;

} // namespace sasktran2
