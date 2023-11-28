#include "sasktran2/geometry.h"
#include <sasktran2/output.h>
#include <sasktran2/math/scattering.h>

namespace sasktran2 {

    std::pair<double, double>
    solar_rotation_factors(double ref_cos_sza,
                           const Eigen::Vector3d& look_vector) {
        // Project the sun and the *unrotated* sun into perpindicular componets
        // to the look vector
        Eigen::Vector3d sun(0, 0, 1);
        Eigen::Vector3d true_sun(sqrt(1 - ref_cos_sza * ref_cos_sza), 0,
                                 ref_cos_sza);

        if ((abs(sun.dot(look_vector)) >= 1) ||
            (abs(true_sun.dot(look_vector)) >= 1)) {
            // Parallel sun, not sure what to do...
            // TODO: CHeck this
            return std::make_pair(1.0, 0.0);
        }

        auto perp_sun = (sun - sun.dot(look_vector) * look_vector).normalized();
        auto perp_true_sun =
            (true_sun - true_sun.dot(look_vector) * look_vector).normalized();

        // Find the angle between them and use that as the Stokes rotation angle
        double rot_rangle = acos(perp_sun.dot(perp_true_sun));

        std::pair<double, double> result;

        result.first = cos(2 * rot_rangle);
        result.second = sin(2 * rot_rangle);

        return result;
    }

    std::pair<double, double>
    observer_rotation_factors(const Eigen::Vector3d& observer,
                              const Eigen::Vector3d& look_vector) {
        Eigen::Vector3d z(0, 0, 1);

        if ((abs(z.dot(look_vector)) >= 1) ||
            (abs(observer.dot(look_vector)) >= 1)) {
            // Parallel sun, not sure what to do...
            // TODO: CHeck this
            return std::make_pair(1.0, 0.0);
        }

        auto perp_z = (z - z.dot(look_vector) * look_vector).normalized();
        auto perp_obs =
            (observer - observer.dot(look_vector) * look_vector).normalized();

        // Find the angle between them and use that as the Stokes rotation angle
        double rot_rangle = acos(perp_z.dot(perp_obs));

        std::pair<double, double> result;

        result.first = cos(2 * rot_rangle);
        result.second = sin(2 * rot_rangle);

        return result;
    }

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
                sasktran2::Config::StokesBasis::standard) {
                const auto& sun = geometry.coordinates().sun_unit();

                if (geometry.coordinates().sun_forced_z()) {
                    for (int i = 0; i < rays.size(); ++i) {
                        auto CS = solar_rotation_factors(
                            geometry.coordinates().cos_sza_at_reference(),
                            rays[i].observer_and_look.look_away);

                        m_stokes_C[i] = CS.first;
                        m_stokes_S[i] = CS.second;
                    }
                }
            } else if (config.stokes_basis() ==
                       sasktran2::Config::StokesBasis::observer) {
                for (int i = 0; i < rays.size(); ++i) {
                    auto CS = observer_rotation_factors(
                        rays[i].observer_and_look.observer.position,
                        rays[i].observer_and_look.look_away);

                    m_stokes_C[i] = CS.first;
                    m_stokes_S[i] = CS.second;
                }
            }
        }
    }

    template class Output<1>;
    template class Output<3>;

} // namespace sasktran2
