#pragma once

#include "sasktran2/atmosphere/atmosphere.h"
#include "sasktran2/config.h"
#include "sasktran2/geometry.h"
#include "sasktran2/raytracing.h"
#include <sasktran2/viewinggeometry_internal.h>
#include <Eigen/src/Core/Matrix.h>
#include <sasktran2/internal_common.h>
#include <sasktran2/dual.h>
#include <sasktran2/wavelength_block.h>

namespace sasktran2 {

    /** Essentially a pure virtual void class to interface with SWIG, removing
     * the NSTOKES template.
     *
     */
    class OutputInterface {
      public:
        virtual ~OutputInterface() {}
    };

    /** Base class for the output container.  Provides storage for the output
     * values, as well as defines what quantities are actually output.  The
     * Sasktran2 engine both solves the radiative transfer equation to get the
     * full radiance field, \f$I(location, solid angle, wavelength)\f$ and
     * integrates this along the specified LOS vectors to get
     *  \f$(I_{los_i}(wavelength)\f$.  As well as corresponding derivatives.
     *
     *  The user output is going to be a function of these two parameters.
     *
     * @tparam NSTOKES
     */
    template <int NSTOKES> class Output : public OutputInterface {
      private:
        virtual void resize() = 0;

      protected:
        int m_nlos;
        int m_nfluxpos;
        int m_nfluxtype;
        int m_nwavel;
        int m_nderiv;
        int m_ngeometry;
        int m_wavelength_block_capacity = 1;

        Eigen::VectorXd m_stokes_C;
        Eigen::VectorXd m_stokes_S;

        // Diagnostic qunatities
        Eigen::MatrixXd m_los_optical_depth;

        const sasktran2::atmosphere::Atmosphere<NSTOKES>* m_atmosphere;
        const sasktran2::Config* m_config;

      public:
        Output(){};
        virtual ~Output() {}

        /**
         *  Initializes the output container with the necessary information.
         * This is called by the engine during calculate_radiance
         */
        virtual void initialize(
            const sasktran2::Config& config,
            const sasktran2::Geometry& geometry,
            const sasktran2::viewinggeometry::InternalViewingGeometry&
                internal_viewing,
            const sasktran2::atmosphere::Atmosphere<NSTOKES>& atmosphere);

        /** Sets the wavelength block capacity negotiated by the engine. */
        void set_wavelength_block_capacity(int block_capacity) {
            if (block_capacity < 1) {
                throw std::invalid_argument(
                    "Output wavelength block capacity must be positive");
            }
            m_wavelength_block_capacity = block_capacity;
        }

        /** Method the Sasktran2 engine calls for each integrated line of sight
         * and contiguous wavelength block.
         *
         * @param block Contiguous wavelengths represented by radiance
         * @param radiance The final calculated radiance and corresponding
         * derivatives
         * @param losidx The index of this line of sight
         * @param threadidx The source-integration thread index
         */
        virtual void
        assign(const sasktran2::WavelengthBlock<>& block,
               const sasktran2::WavelengthBlockDual<NSTOKES>& radiance,
               int losidx, int threadidx) = 0;

        virtual void assign_flux(
            const sasktran2::Dual<double, sasktran2::dualstorage::dense, 1>&
                flux,
            int fluxidx, int wavelidx, int threadidx, int flux_type_idx) {
            spdlog::error(
                "Flux assignment not implemented for this output type");
        };

        /**
         *
         * @return The engine number of lines of sight
         */
        int num_los() const { return m_nlos; }

        /**
         *
         * @return The engine number of wavelengths
         */
        int num_wavel() const { return m_nwavel; }

        /**
         *
         * @return The engine number of derivatives
         */
        int num_deriv() const { return m_nderiv; }

        const Eigen::MatrixXd& los_optical_depth() const {
            return m_los_optical_depth;
        }

        Eigen::MatrixXd& los_optical_depth() { return m_los_optical_depth; }
    };

    /** An idealized output container where only the line of sight radiances are
     * stored, and are stored on the native model calculation grid.
     *
     * @tparam NSTOKES
     */
    template <int NSTOKES> class OutputIdealDense : public Output<NSTOKES> {
      private:
        sasktran2::Dual<double, sasktran2::dualstorage::dense>
            m_radiance; /**< Internal storage */

        void resize();

        template <typename Radiance>
        void assign_lane(const Radiance& radiance, int losidx, int wavelidx);

      public:
        OutputIdealDense(){};

        void assign(const sasktran2::WavelengthBlock<>& block,
                    const sasktran2::WavelengthBlockDual<NSTOKES>& radiance,
                    int losidx, int threadidx) override;

        /**
         *
         * @return The stored radiance container
         */
        sasktran2::Dual<double, sasktran2::dualstorage::dense>& radiance() {
            return m_radiance;
        }
    };

    /** An idealized output container where only the line of sight radiances are
     * stored, and the derivatives are mapped based on the atmosphere derivative
     * mappings
     *
     * @tparam NSTOKES
     */
    template <int NSTOKES> class OutputDerivMapped : public Output<NSTOKES> {
      private:
        sasktran2::Dual<double, sasktran2::dualstorage::dense>
            m_radiance; /**< Internal storage, only radiances are stored no
                           derivatives */

        std::map<std::string, Eigen::MatrixXd> m_derivatives;
        std::map<std::string, Eigen::MatrixXd> m_surface_derivatives;
        std::vector<Eigen::MatrixXd> m_native_thread_storage;
        std::vector<Eigen::MatrixXd> m_mapped_thread_storage;

        void resize();

        template <typename Radiance>
        void assign_lane(const Radiance& radiance, int losidx, int wavelidx,
                         int threadidx);

      public:
        OutputDerivMapped(){};

        void assign(const sasktran2::WavelengthBlock<>& block,
                    const sasktran2::WavelengthBlockDual<NSTOKES>& radiance,
                    int losidx, int threadidx) override;

        /**
         *
         * @return The stored radiance container
         */
        sasktran2::Dual<double, sasktran2::dualstorage::dense>& radiance() {
            return m_radiance;
        }

        const std::map<std::string, Eigen::MatrixXd>& derivatives() const {
            return m_derivatives;
        }

        const std::map<std::string, Eigen::MatrixXd>&
        surface_derivatives() const {
            return m_surface_derivatives;
        }
    };

    /** An idealized output container where only the line of sight radiances are
     * stored, and the derivatives are mapped based on the atmosphere derivative
     * mappings.  This is the output class that is used by the C api, it differs
     * from the OutputDerivMapped class in that memory is passed in through the
     * API for the derivatives and radiances.  Therefore we don't provide any
     * accessors
     *
     * @tparam NSTOKES
     */
    template <int NSTOKES> class OutputC : public Output<NSTOKES> {
      private:
        Eigen::Map<Eigen::VectorXd> m_radiance;
        Eigen::Map<Eigen::VectorXd> m_flux;

        std::map<std::string, Eigen::Map<Eigen::MatrixXd>> m_derivatives;
        std::map<std::string, Eigen::Map<Eigen::MatrixXd>>
            m_surface_derivatives;

        std::map<std::string, Eigen::Map<Eigen::MatrixXd>> m_flux_derivatives;
        std::map<std::string, Eigen::Map<Eigen::MatrixXd>>
            m_flux_surface_derivatives;
        std::vector<Eigen::MatrixXd> m_native_thread_storage;
        std::vector<Eigen::MatrixXd> m_mapped_thread_storage;

        void resize();

        template <typename Radiance>
        void assign_lane(const Radiance& radiance, int losidx, int wavelidx,
                         int threadidx);

      public:
        OutputC(Eigen::Map<Eigen::VectorXd> radiance,
                Eigen::Map<Eigen::VectorXd> flux)
            : m_radiance(radiance), m_flux(flux){};

        void set_derivative_mapping_memory(
            const std::string& name,
            Eigen::Map<Eigen::MatrixXd> derivative_mapping) {
            m_derivatives.insert(
                {name, Eigen::Map<Eigen::MatrixXd>(nullptr, 0,
                                                   0)}); // create a null map
            // then placement new into the map
            Eigen::Map<Eigen::MatrixXd>* ref = &m_derivatives.at(name);

            new (ref) Eigen::Map<Eigen::MatrixXd>(derivative_mapping.data(),
                                                  derivative_mapping.rows(),
                                                  derivative_mapping.cols());
        }

        void set_surface_derivative_mapping_memory(
            const std::string& name,
            Eigen::Map<Eigen::MatrixXd> derivative_mapping) {
            m_surface_derivatives.insert(
                {name, Eigen::Map<Eigen::MatrixXd>(nullptr, 0,
                                                   0)}); // create a null map
            // then placement new into the map
            Eigen::Map<Eigen::MatrixXd>* ref = &m_surface_derivatives.at(name);

            new (ref) Eigen::Map<Eigen::MatrixXd>(derivative_mapping.data(),
                                                  derivative_mapping.rows(),
                                                  derivative_mapping.cols());
        }

        void set_flux_derivative_mapping_memory(
            const std::string& name,
            Eigen::Map<Eigen::MatrixXd> derivative_mapping) {
            m_flux_derivatives.insert(
                {name, Eigen::Map<Eigen::MatrixXd>(nullptr, 0,
                                                   0)}); // create a null map
            // then placement new into the map
            Eigen::Map<Eigen::MatrixXd>* ref = &m_flux_derivatives.at(name);

            new (ref) Eigen::Map<Eigen::MatrixXd>(derivative_mapping.data(),
                                                  derivative_mapping.rows(),
                                                  derivative_mapping.cols());
        }

        void set_flux_surface_derivative_mapping_memory(
            const std::string& name,
            Eigen::Map<Eigen::MatrixXd> derivative_mapping) {
            m_flux_surface_derivatives.insert(
                {name, Eigen::Map<Eigen::MatrixXd>(nullptr, 0,
                                                   0)}); // create a null map
            // then placement new into the map
            Eigen::Map<Eigen::MatrixXd>* ref =
                &m_flux_surface_derivatives.at(name);

            new (ref) Eigen::Map<Eigen::MatrixXd>(derivative_mapping.data(),
                                                  derivative_mapping.rows(),
                                                  derivative_mapping.cols());
        }

        void assign(const sasktran2::WavelengthBlock<>& block,
                    const sasktran2::WavelengthBlockDual<NSTOKES>& radiance,
                    int losidx, int threadidx) override;

        void
        assign_flux(const sasktran2::Dual<double, sasktran2::dualstorage::dense,
                                          1>& flux,
                    int fluxidx, int wavelidx, int threadidx,
                    int flux_type_idx) override;
    };

    /** Output adapter that contracts each native radiance derivative row with
     * one labeled parameter direction without storing a Jacobian. */
    template <int NSTOKES> class OutputJVP : public Output<NSTOKES> {
      private:
        Eigen::Map<Eigen::VectorXd> m_radiance;
        Eigen::Map<Eigen::VectorXd> m_jvp;
        std::map<std::string, Eigen::VectorXd> m_derivative_tangents;
        std::map<std::string, Eigen::VectorXd> m_surface_tangents;
        std::vector<Eigen::MatrixXd> m_native_thread_storage;
        std::vector<Eigen::MatrixXd> m_mapped_thread_storage;

        void resize();

        template <typename Radiance>
        void assign_lane(const Radiance& radiance, int losidx, int wavelidx,
                         int threadidx);

      public:
        OutputJVP(Eigen::Map<Eigen::VectorXd> radiance,
                  Eigen::Map<Eigen::VectorXd> jvp)
            : m_radiance(radiance), m_jvp(jvp) {}

        void set_derivative_tangent(const std::string& name,
                                    Eigen::Ref<const Eigen::VectorXd> tangent) {
            m_derivative_tangents[name] = tangent;
        }

        void set_surface_tangent(const std::string& name,
                                 Eigen::Ref<const Eigen::VectorXd> tangent) {
            m_surface_tangents[name] = tangent;
        }

        void assign(const sasktran2::WavelengthBlock<>& block,
                    const sasktran2::WavelengthBlockDual<NSTOKES>& radiance,
                    int losidx, int threadidx) override;

        void native_tangent(int wavelidx, Eigen::VectorXd& tangent) const;

        void assign_native(int losidx, int wavelidx,
                           const Eigen::Vector<double, NSTOKES>& value,
                           const Eigen::Vector<double, NSTOKES>& jvp);

        void assign_flux(
            const sasktran2::Dual<double, sasktran2::dualstorage::dense, 1>&,
            int, int, int, int) override {}
    };

    /** Output adapter that contracts radiance cotangents with each native
     * derivative row and accumulates labeled parameter gradients. */
    template <int NSTOKES> class OutputVJP : public Output<NSTOKES> {
      private:
        Eigen::Map<Eigen::VectorXd> m_radiance;
        Eigen::Map<const Eigen::VectorXd> m_cotangent;
        std::map<std::string, Eigen::Map<Eigen::VectorXd>>
            m_derivative_gradients;
        std::map<std::string, Eigen::Map<Eigen::VectorXd>> m_surface_gradients;
        std::vector<std::map<std::string, Eigen::VectorXd>>
            m_thread_derivative_gradients;
        std::vector<std::map<std::string, Eigen::VectorXd>>
            m_thread_surface_gradients;
        std::vector<Eigen::MatrixXd> m_native_thread_storage;

        void resize();

        template <typename Radiance>
        void assign_lane(const Radiance& radiance, int losidx, int wavelidx,
                         int threadidx);

      public:
        OutputVJP(Eigen::Map<Eigen::VectorXd> radiance,
                  Eigen::Map<const Eigen::VectorXd> cotangent)
            : m_radiance(radiance), m_cotangent(cotangent) {}

        void set_derivative_gradient_memory(
            const std::string& name,
            Eigen::Map<Eigen::VectorXd> derivative_gradient);

        void set_surface_gradient_memory(
            const std::string& name,
            Eigen::Map<Eigen::VectorXd> derivative_gradient);

        void assign(const sasktran2::WavelengthBlock<>& block,
                    const sasktran2::WavelengthBlockDual<NSTOKES>& radiance,
                    int losidx, int threadidx) override;

        Eigen::Vector<double, NSTOKES> native_cotangent(int losidx,
                                                        int wavelidx) const;

        void assign_native_value(int losidx, int wavelidx,
                                 const Eigen::Vector<double, NSTOKES>& value);

        void accumulate_native_gradient(
            int wavelidx, int threadidx,
            Eigen::Ref<const Eigen::VectorXd> native_gradient);

        void assign_flux(
            const sasktran2::Dual<double, sasktran2::dualstorage::dense, 1>&,
            int, int, int, int) override {}

        /** Reduces thread-local gradient accumulation into caller memory. */
        void finalize();
    };
} // namespace sasktran2
