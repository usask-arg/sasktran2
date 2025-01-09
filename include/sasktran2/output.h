#pragma once

#include "sasktran2/atmosphere/atmosphere.h"
#include "sasktran2/config.h"
#include "sasktran2/sensor.h"
#include "sasktran2/geometry.h"
#include "sasktran2/raytracing.h"
#include <Eigen/src/Core/Matrix.h>
#include <sasktran2/internal_common.h>
#include <sasktran2/dual.h>

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
        int m_nwavel;
        int m_nderiv;
        int m_ngeometry;

        Eigen::VectorXd m_stokes_C;
        Eigen::VectorXd m_stokes_S;

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
            const sasktran2::Geometry1D& geometry,
            const std::vector<sasktran2::raytracing::TracedRay>& rays,
            const sasktran2::atmosphere::Atmosphere<NSTOKES>& atmosphere);

        /** Method the Sasktran2 engine calls for each integrated line of
         * sight/wavelength
         *
         * @param radiance The final calculated radiance and corresponding
         * derivatives
         * @param losidx The index of this line of sight
         * @param wavelidx The index of this wavelength
         */
        virtual void
        assign(const sasktran2::Dual<double, sasktran2::dualstorage::dense,
                                     NSTOKES>& radiance,
               int losidx, int wavelidx, int threadidx) = 0;

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

      public:
        OutputIdealDense(){};

        void assign(const sasktran2::Dual<double, sasktran2::dualstorage::dense,
                                          NSTOKES>& radiance,
                    int losidx, int wavelidx, int threadidx);

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

        void resize();

      public:
        OutputDerivMapped(){};

        void assign(const sasktran2::Dual<double, sasktran2::dualstorage::dense,
                                          NSTOKES>& radiance,
                    int losidx, int wavelidx, int threadidx);

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

    /**
     */
    template <int NSTOKES> class OutputSpectralSensor : public Output<NSTOKES> {
      private:
        const sasktran2::sensor::Sensor& m_sensor;

        sasktran2::Dual<double, sasktran2::dualstorage::dense>
            m_radiance; /**< Internal storage */

        void resize();

      public:
        OutputSpectralSensor(const sasktran2::sensor::Sensor& sensor)
            : m_sensor(sensor){};

        void assign(const sasktran2::Dual<double, sasktran2::dualstorage::dense,
                                          NSTOKES>& radiance,
                    int losidx, int wavelidx);

        /**
         *
         * @return The stored radiance container
         */
        sasktran2::Dual<double, sasktran2::dualstorage::dense>& radiance() {
            return m_radiance;
        }
    };
} // namespace sasktran2
