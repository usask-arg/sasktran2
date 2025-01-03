#pragma once

#include <sasktran2/internal_common.h>

namespace sasktran2::sensor {
    /**
     *  @brief  A class to represent a spectral sensor.
     *
     *  This class is a base class for spectral sensors. It provides a mapping
     * from the high-resolution calculation space to a reduced observation
     * space.  I.e. if I is the radiance at a high-resolution grid of shape (n,)
     * then a spectral sensor is a matrix (n, m) such that I' = I @ M where I'
     * is the radiance at the sensor's wavelengths.
     */
    class Sensor {
      private:
      public:
        Sensor();
        virtual ~Sensor();

        virtual const Eigen::MatrixXd& spectral_mapping(int ray_idx) const;
    };

    class ManualSensor : public Sensor {
      private:
        Eigen::MatrixXd m_spectral_mapping;

      public:
        ManualSensor(int num_samples, int num_wavelengths);
        virtual ~ManualSensor();

        const Eigen::MatrixXd& spectral_mapping(int ray_idx) const override;
    };
} // namespace sasktran2::sensor
