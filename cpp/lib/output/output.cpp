#include "sasktran2/atmosphere/atmosphere.h"
#include "sasktran2/config.h"
#include "sasktran2/geometry.h"
#include "sasktran2/viewinggeometry_internal.h"
#include <sasktran2/output.h>
#include <sasktran2/math/scattering.h>

namespace sasktran2 {

    template <int NSTOKES>
    void Output<NSTOKES>::initialize(
        const sasktran2::Config& config, const sasktran2::Geometry& geometry,
        const sasktran2::viewinggeometry::InternalViewingGeometry&
            internal_viewing,
        const sasktran2::atmosphere::Atmosphere<NSTOKES>& atmosphere) {
        m_nlos = static_cast<int>(internal_viewing.num_rays());
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

            // The sources are calculated with the Stokes reference in the
            // plane of the observer position and the look vector, so only the
            // standard and solar bases need a rotation
            if (config.stokes_basis() !=
                sasktran2::Config::StokesBasis::observer) {
                const auto& coords = geometry.coordinates();
                const Eigen::Vector3d& to_reference =
                    config.stokes_basis() ==
                            sasktran2::Config::StokesBasis::solar
                        ? coords.sun_unit()
                        : coords.reference_z();

                for (int i = 0; i < m_nlos; ++i) {
                    const auto& ray = internal_viewing.viewing_ray(i);
                    auto CS = coords.stokes_rotation(
                        ray.look_away, ray.observer.position, to_reference);

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
