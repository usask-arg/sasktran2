#include <sasktran2/output.h>

namespace sasktran2 {
    template <int NSTOKES> void OutputDerivMapped<NSTOKES>::resize() {
        m_radiance.resize(NSTOKES * this->m_nwavel * this->m_nlos, 0, false);

        int i = 0;
        for (auto& [name, deriv] :
             this->m_atmosphere->storage().derivative_mappings_const()) {
            m_derivatives[name].resize(NSTOKES * this->m_nwavel * this->m_nlos,
                                       deriv.num_output());
        }

        for (auto& [name, deriv] :
             this->m_atmosphere->surface().derivative_mappings()) {
            m_surface_derivatives[name].resize(
                NSTOKES * this->m_nwavel * this->m_nlos, 1);
        }

        m_native_thread_storage.resize(this->m_config->num_threads());
        for (auto& storage : m_native_thread_storage) {
            storage.resize(NSTOKES, this->m_ngeometry);
        }
    }

    template <int NSTOKES>
    void OutputDerivMapped<NSTOKES>::assign(
        const sasktran2::Dual<double, sasktran2::dualstorage::dense, NSTOKES>&
            radiance,
        int losidx, int wavelidx, int threadidx) {

        int linear_index = NSTOKES * this->m_nlos * wavelidx + NSTOKES * losidx;

        auto& deriv_storage = m_native_thread_storage[threadidx];

        // Quantities necessary for derivative propagation

        if constexpr (NSTOKES >= 1) {
            m_radiance.value(linear_index) = radiance.value(0);
        }

        Eigen::Ref<const Eigen::Matrix<double, NSTOKES, -1>> d_rad_by_d_ssa =
            radiance.d_ssa(this->m_ngeometry); // Col vector [stokes X N]
        Eigen::Ref<const Eigen::Matrix<double, NSTOKES, -1>> d_rad_by_d_k =
            radiance.d_extinction(this->m_ngeometry); // col vector [stokes X N]

        // Do the atmosphere mappings
        for (auto& [name, deriv] : m_derivatives) {
            const auto& mapping =
                this->m_atmosphere->storage().derivative_mappings_const().at(
                    name);

            const auto& d_ssa =
                mapping.native_mapping().d_ssa.value(); // [N x wavelength]
            const auto& d_extinction =
                mapping.native_mapping()
                    .d_extinction.value(); // [N x wavelength]

            for (int i = 0; i < NSTOKES; ++i) {
                deriv_storage(i, Eigen::all).array() =
                    d_ssa.col(wavelidx).transpose().array() *
                        d_rad_by_d_ssa(i, Eigen::all).array() +
                    d_extinction.col(wavelidx).transpose().array() *
                        d_rad_by_d_k(i, Eigen::all).array();
            }

            if (mapping.is_scattering_derivative()) {
                // Have to include the scattering terms
                Eigen::Ref<const Eigen::Matrix<double, NSTOKES, -1>>
                    d_rad_by_d_scat = radiance.d_scatterer(
                        this->m_ngeometry, mapping.get_scattering_index());
                const auto& scat_factor =
                    mapping.native_mapping()
                        .scat_factor.value(); // [N x wavelength]

                for (int i = 0; i < NSTOKES; ++i) {
                    deriv_storage(i, Eigen::all).array() +=
                        scat_factor.col(wavelidx).transpose().array() *
                        d_rad_by_d_scat(i, Eigen::all).array();
                }
            }

            if (mapping.native_mapping().d_emission.has_value()) {
                // Include the emission terms
                Eigen::Ref<const Eigen::Matrix<double, NSTOKES, -1>>
                    d_rad_by_d_emission = radiance.d_emission(
                        this->m_ngeometry,
                        this->m_atmosphere->num_scattering_deriv_groups());

                const auto& d_emission =
                    mapping.native_mapping().d_emission.value();

                // Emission terms only ever affect the stokes = 0 component
                deriv_storage(0, Eigen::all).array() +=
                    d_emission.col(wavelidx).transpose().array() *
                    d_rad_by_d_emission(0, Eigen::all).array();
            }

            // If we are interpolating in geometry, then apply the interpolator
            // before assigning
            if (mapping.get_interpolator_const().has_value()) {
                const auto& interpolator =
                    mapping.get_interpolator_const().value(); // [N x M]

                // Always have to assign the stokes = 0 component
                deriv(linear_index, Eigen::all).array() =
                    deriv_storage(0, Eigen::all) * interpolator;

                if constexpr (NSTOKES >= 3) {
                    // Have to assign rotated Q/U components
                    deriv(linear_index + 1, Eigen::all).array() =
                        (this->m_stokes_C[losidx] *
                             deriv_storage(1, Eigen::all) -
                         this->m_stokes_S[losidx] *
                             deriv_storage(2, Eigen::all)) *
                        interpolator;

                    deriv(linear_index + 2, Eigen::all).array() =
                        (this->m_stokes_S[losidx] *
                             deriv_storage(1, Eigen::all) +
                         this->m_stokes_C[losidx] *
                             deriv_storage(2, Eigen::all)) *
                        interpolator;
                }
            } else {
                // Similarly always assign stokes = 0
                deriv(linear_index, Eigen::all).array() =
                    deriv_storage(0, Eigen::all).array();

                // And do rotate Q/U if necessary
                if constexpr (NSTOKES >= 3) {
                    deriv(linear_index + 1, Eigen::all).array() =
                        this->m_stokes_C[losidx] *
                            deriv_storage(1, Eigen::all).array() -
                        this->m_stokes_S[losidx] *
                            deriv_storage(2, Eigen::all).array();
                    deriv(linear_index + 2, Eigen::all).array() =
                        this->m_stokes_S[losidx] *
                            deriv_storage(1, Eigen::all).array() +
                        this->m_stokes_C[losidx] *
                            deriv_storage(2, Eigen::all).array();
                }
            }
            if (mapping.log_radiance_space()) {
                deriv(linear_index, Eigen::all).array() /= radiance.value(0);
                if constexpr (NSTOKES >= 3) {
                    deriv(linear_index + 1, Eigen::all).array() /=
                        radiance.value(0);
                    deriv(linear_index + 2, Eigen::all).array() /=
                        radiance.value(0);
                }
            }
        }

        Eigen::Ref<const Eigen::Matrix<double, NSTOKES, -1>>
            d_rad_by_d_surface = radiance.d_brdf(
                this->m_ngeometry,
                this->m_atmosphere->num_scattering_deriv_groups(),
                this->m_atmosphere->surface()
                    .num_deriv()); // [stokes X num_brdf_deriv]
        // Then do the surface mappings
        for (auto& [name, deriv] : m_surface_derivatives) {
            const auto& mapping =
                this->m_atmosphere->surface().derivative_mappings().at(name);

            // Do the BRDF derivatives
            if (mapping.native_surface_mapping().d_brdf.has_value()) {
                // Just assign the stokes = 0 component
                // TODO: Revisit when we have polarized brdfs?
                deriv(linear_index, 0) =
                    d_rad_by_d_surface(0, Eigen::all)
                        .dot(
                            mapping.native_surface_mapping().d_brdf.value().row(
                                wavelidx));
            }

            // Then do the emission derivatives
            if (mapping.native_surface_mapping().d_emission.has_value()) {
                double d_rad_by_d_emission = radiance.deriv(
                    0,
                    this->m_atmosphere->surface_emission_deriv_start_index());
                // Just assign the stokes = 0 component
                deriv(linear_index, 0) =
                    d_rad_by_d_emission *
                    mapping.native_surface_mapping().d_emission.value()(
                        wavelidx, 0);
            }

            if constexpr (NSTOKES >= 3) {
                // Temporaries for dQ, dU
                double dQ =
                    d_rad_by_d_surface(1, Eigen::all)
                        .dot(
                            mapping.native_surface_mapping().d_brdf.value().row(
                                wavelidx));
                double dU =
                    d_rad_by_d_surface(2, Eigen::all)
                        .dot(
                            mapping.native_surface_mapping().d_brdf.value().row(
                                wavelidx));

                // And assign with rotation
                deriv(linear_index + 1, 0) = this->m_stokes_C[losidx] * dQ -
                                             this->m_stokes_S[losidx] * dU;
                deriv(linear_index + 2, 0) = this->m_stokes_S[losidx] * dQ +
                                             this->m_stokes_C[losidx] * dU;
            }
        }

        if constexpr (NSTOKES >= 3) {
            // Q/U components have a rotation
            m_radiance.value(linear_index + 1) =
                this->m_stokes_C[losidx] * radiance.value(1) -
                this->m_stokes_S[losidx] * radiance.value(2);

            m_radiance.value(linear_index + 2) =
                this->m_stokes_S[losidx] * radiance.value(1) +
                this->m_stokes_C[losidx] * radiance.value(2);
        }

        if constexpr (NSTOKES == 4) {
            // V component is a strict copy
            m_radiance.value(linear_index + 3) = radiance.value(3);
        }
    }

    template class OutputDerivMapped<1>;
    template class OutputDerivMapped<3>;

} // namespace sasktran2
