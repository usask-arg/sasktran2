#pragma once
#include <sasktran2/internal_common.h>


namespace sasktran2::mie::distribution {
    struct Distribution {
        virtual double pdf(double x) = 0;

        virtual double left_bound() {
            return 0.0;
        }

        virtual double right_bound() {
            return std::numeric_limits<double>::infinity();
        }
    };

    struct LogNormal : public Distribution {
        private:
            double m_mean;
            double m_stddev;
        public:

        LogNormal(double mean, double stddev) : m_mean(mean), m_stddev(stddev) {}

        double pdf(double x) {
            if (x <= 0.0) {
                return 0.0;
            }
            double coeff = 1.0 / (x * m_stddev * sqrt(2.0 * M_PI));
            double exponent = -pow(log(x) - m_mean, 2) / (2.0 * m_stddev * m_stddev);
            return coeff * exp(exponent);
        }
    };
}