#pragma once

#include <sasktran2/source_interface.h>
#include <sasktran2/raytracing.h>
#include <sasktran2/atmosphere/atmosphere.h>
#include <sasktran2/config.h>
#include <sasktran2/dual.h>

namespace sasktran2::emission {

    template <int NSTOKES, Config::EmissionSource EMISSION_SOURCE_TYPE =
                               Config::EmissionSource::standard>
    class EmissionSource : public SourceTermInterface<NSTOKES> {
      private:
        const sasktran2::atmosphere::Atmosphere<NSTOKES>* m_atmosphere =
            nullptr;
        std::vector<bool> m_ground_is_hit;

        void integrated_source_constant(
            int wavelidx, int losidx, int layeridx, int wavel_threadidx,
            int threadidx, const sasktran2::raytracing::LayerGeometry& layer,
            const sasktran2::raytracing::GridWeightStencilView&
                entrance_weights,
            const sasktran2::raytracing::GridWeightStencilView& exit_weights,
            const sasktran2::SparseODDualView& shell_od,
            sasktran2::Dual<double, sasktran2::dualstorage::dense, NSTOKES>&
                source) const;

      public:
        bool supports_wavelength_blocks() const override { return true; }

        void calculate(const sasktran2::WavelengthBlock&, int) override {}
        /** Here the emission source term saves the los_rays, which are
         * needed to detect ground hits and whether to include
         * surface emissions at the end of the ray.
         *
         *  @param internal_viewing Information on the internal viewing
         * geometry, los_rays and flux observers
         */
        void initialize_geometry(
            const sasktran2::viewinggeometry::InternalViewingGeometry&
                internal_viewing) override;

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
      private:
        void integrated_source(
            int wavelidx, int losidx, int layeridx, int wavel_threadidx,
            int threadidx, const sasktran2::raytracing::TracedLayer& layer,
            const sasktran2::raytracing::GridWeightStencilView&
                entrance_weights,
            const sasktran2::raytracing::GridWeightStencilView& exit_weights,
            const sasktran2::SparseODDualView& shell_od,
            sasktran2::Dual<double, sasktran2::dualstorage::dense, NSTOKES>&
                source,
            typename SourceTermInterface<
                NSTOKES>::IntegrationDirection direction =
                SourceTermInterface<NSTOKES>::IntegrationDirection::none) const;

        void integrated_source_block(
            const sasktran2::WavelengthBlock& batch, int losidx, int layeridx,
            int wavel_threadidx, int threadidx,
            const sasktran2::raytracing::TracedLayer& layer,
            const sasktran2::raytracing::GridWeightStencilView&
                entrance_weights,
            const sasktran2::raytracing::GridWeightStencilView& exit_weights,
            const sasktran2::WavelengthBlockODView& shell_od,
            sasktran2::WavelengthBlockDual<NSTOKES>& source,
            typename SourceTermInterface<
                NSTOKES>::IntegrationDirection direction =
                SourceTermInterface<NSTOKES>::IntegrationDirection::none) const;

      public:
        void integrated_source(
            const sasktran2::WavelengthBlock& block, int losidx, int layeridx,
            int wavel_threadidx, int threadidx,
            const sasktran2::raytracing::TracedLayer& layer,
            const sasktran2::raytracing::GridWeightStencilView&
                entrance_weights,
            const sasktran2::raytracing::GridWeightStencilView& exit_weights,
            const sasktran2::WavelengthBlockODView& shell_od,
            sasktran2::WavelengthBlockDualView<NSTOKES>& source,
            typename SourceTermInterface<NSTOKES>::IntegrationDirection
                direction =
                    SourceTermInterface<NSTOKES>::IntegrationDirection::none)
            const override {
            if (source.is_scalar()) {
                integrated_source(
                    block.start, losidx, layeridx, wavel_threadidx, threadidx,
                    layer, entrance_weights, exit_weights, shell_od.scalar(),
                    source.scalar(), direction);
            } else {
                integrated_source_block(block, losidx, layeridx,
                                        wavel_threadidx, threadidx, layer,
                                        entrance_weights, exit_weights,
                                        shell_od, source.block(), direction);
            }
        }

        bool supports_geometry_dimension(int dimension) const override {
            return dimension >= 1;
        }

        /** Calculates the source term at the end of the ray.  Common examples
         * of this are ground scattering, ground emission, or the solar radiance
         * if looking directly at the sun.
         *
         * @param wavelidx Raw index for the wavelength we are calculating
         * @param losidx Raw index pointing to the ray that was previously
         * passed in initialize_geometry
         * @param source The returned source term
         */
      private:
        void end_of_ray_source(
            int wavelidx, int losidx, int wavel_threadidx, int threadidx,
            sasktran2::Dual<double, sasktran2::dualstorage::dense, NSTOKES>&
                source) const;

        void end_of_ray_source_block(
            const sasktran2::WavelengthBlock& batch, int losidx,
            int wavel_threadidx, int threadidx,
            sasktran2::WavelengthBlockDual<NSTOKES>& source) const;

      public:
        void end_of_ray_source(const sasktran2::WavelengthBlock& block,
                               int losidx, int wavel_threadidx, int threadidx,
                               sasktran2::WavelengthBlockDualView<NSTOKES>&
                                   source) const override {
            if (source.is_scalar()) {
                end_of_ray_source(block.start, losidx, wavel_threadidx,
                                  threadidx, source.scalar());
            } else {
                end_of_ray_source_block(block, losidx, wavel_threadidx,
                                        threadidx, source.block());
            }
        }

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
        void start_of_ray_source(
            const sasktran2::WavelengthBlock&, int, int, int,
            sasktran2::WavelengthBlockDualView<NSTOKES>&) const override {}
    };

} // namespace sasktran2::emission
