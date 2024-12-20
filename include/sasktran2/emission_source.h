#pragma once

#include <sasktran2/source_interface.h>
#include <sasktran2/raytracing.h>
#include <sasktran2/atmosphere/atmosphere.h>
#include <sasktran2/config.h>
#include <sasktran2/dual.h>

namespace sasktran2::emission {

    template <int NSTOKES>
    class EmissionSource : public SourceTermInterface<NSTOKES> {
      private:
        const sasktran2::atmosphere::Atmosphere<NSTOKES>* m_atmosphere;
        const std::vector<sasktran2::raytracing::TracedRay>* m_los_rays;

        void integrated_source_constant(
            int wavelidx, int losidx, int layeridx, int wavel_threadidx,
            int threadidx, const sasktran2::raytracing::SphericalLayer& layer,
            const sasktran2::SparseODDualView& shell_od,
            sasktran2::Dual<double, sasktran2::dualstorage::dense, NSTOKES>&
                source) const;

      public:
        /** Here the emission source term saves the los_rays, which are
         * needed to detect ground hits and whether to include
         * surface emissions at the end of the ray.
         *
         *  @param los_rays The traced line of sight rays
         */
        void initialize_geometry(
            const std::vector<sasktran2::raytracing::TracedRay>& los_rays)
            override;

        /**
         *
         */
        void initialize_atmosphere(
            const sasktran2::atmosphere::Atmosphere<NSTOKES>& atmosphere)
            override;

        /** Calculates the integrated source term for a given layer.
         *
         * @param losidx Raw index pointing to the ray that was previously
         * passed in initialize_geometry
         * @param layeridx Raw index pointing to the layer that was previosuly
         * passed in initialize_geometry
         * @param layer The layer that we are integrating over
         * @param source The returned source term
         */
        void integrated_source(
            int wavelidx, int losidx, int layeridx, int wavel_threadidx,
            int threadidx, const sasktran2::raytracing::SphericalLayer& layer,
            const sasktran2::SparseODDualView& shell_od,
            sasktran2::Dual<double, sasktran2::dualstorage::dense, NSTOKES>&
                source) const override;

        /** Calculates the source term at the end of the ray.  Common examples
         * of this are ground scattering, ground emission, or the solar radiance
         * if looking directly at the sun.
         *
         * @param wavelidx Raw index for the wavelength we are calculating
         * @param losidx Raw index pointing to the ray that was previously
         * passed in initialize_geometry
         * @param source The returned source term
         */
        void end_of_ray_source(
            int wavelidx, int losidx, int wavel_threadidx, int threadidx,
            sasktran2::Dual<double, sasktran2::dualstorage::dense, NSTOKES>&
                source) const override;

        /** Calculates the radiance at the start of the ray, i.e., the source
         * term has done the equivalent of the integration along the ray.  This
         * is useful if the source term has it's own way of performing
         * integration that is different than the standard method used in the
         * model.  It can also be used if the source term calculates quantities
         * that can be used to get the radiance directly instead of doing
         * integration.
         *
         *  Typically source terms will only either implement
         * start_of_ray_source, or integrated_source + end_of_ray_source and not
         * both
         *
         * @brief Not used for the Emission source.
         *
         * @param wavelidx
         * @param losidx
         * @param wavel_threadidx
         * @param threadidx
         * @param source
         */
        virtual void start_of_ray_source(
            int wavelidx, int losidx, int wavel_threadidx, int threadidx,
            sasktran2::Dual<double, sasktran2::dualstorage::dense, NSTOKES>&
                source) const override{};
    };

} // namespace sasktran2::emission
