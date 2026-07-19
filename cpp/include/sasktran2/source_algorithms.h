#pragma once

#include <sasktran2/internal_common.h>
#include <sasktran2/dual.h>

#include <cmath>

namespace sasktran2::sourcealgo {

    /** Attenuated endpoint weights for a source interpolated across a layer. */
    struct LinearSourceWeights {
        double start;
        double end;
        double d_start_d_od;
        double d_end_d_od;
    };

    inline LinearSourceWeights
    linear_source_weights(double od, double exp_minus_od,
                          double start_fraction = 0.5,
                          double end_fraction = 0.5) {
        double phi;
        double d_phi;
        double linear_start;
        double linear_end;
        double d_linear_start;
        double d_linear_end;

        if (std::abs(od) < 1e-2) {
            const double od2 = od * od;
            const double od3 = od2 * od;
            const double od4 = od3 * od;
            const double od5 = od4 * od;
            const double od6 = od5 * od;
            linear_start = od / 2.0 - od2 / 6.0 + od3 / 24.0 - od4 / 120.0 +
                           od5 / 720.0 - od6 / 5040.0;
            linear_end = od / 2.0 - od2 / 3.0 + od3 / 8.0 - od4 / 30.0 +
                         od5 / 144.0 - od6 / 840.0;
            d_linear_start = 0.5 - od / 3.0 + od2 / 8.0 - od3 / 30.0 +
                             od4 / 144.0 - od5 / 840.0;
            d_linear_end = 0.5 - 2.0 * od / 3.0 + 3.0 * od2 / 8.0 -
                           2.0 * od3 / 15.0 + 5.0 * od4 / 144.0 - od5 / 140.0;
        } else {
            phi = (1.0 - exp_minus_od) / od;
            d_phi = 1.0 / od - phi * (1.0 + 1.0 / od);
            linear_start = 1.0 - phi;
            linear_end = phi - exp_minus_od;
            d_linear_start = -d_phi;
            d_linear_end = d_phi + exp_minus_od;
        }

        const double start_correction = start_fraction - 0.5;
        const double end_correction = end_fraction - 0.5;
        const double attenuation_integral = std::abs(od) < 1e-2
                                                ? linear_start + linear_end
                                                : 1.0 - exp_minus_od;
        return {
            linear_start + start_correction * attenuation_integral,
            linear_end + end_correction * attenuation_integral,
            d_linear_start + start_correction * exp_minus_od,
            d_linear_end + end_correction * exp_minus_od,
        };
    }

    /** Linear-source weights normalized to a unit-distance coordinate. */
    inline LinearSourceWeights
    normalized_linear_source_weights(double od, double exp_minus_od,
                                     double start_fraction = 0.5,
                                     double end_fraction = 0.5) {
        if (std::abs(od) >= 1e-2) {
            const auto weights = linear_source_weights(
                od, exp_minus_od, start_fraction, end_fraction);
            const double inverse_od = 1.0 / od;
            return {
                weights.start * inverse_od,
                weights.end * inverse_od,
                (weights.d_start_d_od * od - weights.start) * inverse_od *
                    inverse_od,
                (weights.d_end_d_od * od - weights.end) * inverse_od *
                    inverse_od,
            };
        }

        const double od2 = od * od;
        const double od3 = od2 * od;
        const double od4 = od3 * od;
        const double od5 = od4 * od;
        const double od6 = od5 * od;
        const double phi = 1.0 - od / 2.0 + od2 / 6.0 - od3 / 24.0 +
                           od4 / 120.0 - od5 / 720.0 + od6 / 5040.0;
        const double d_phi = -0.5 + od / 3.0 - od2 / 8.0 + od3 / 30.0 -
                             od4 / 144.0 + od5 / 840.0;
        const double start_correction = start_fraction - 0.5;
        const double end_correction = end_fraction - 0.5;

        return {
            0.5 - od / 6.0 + od2 / 24.0 - od3 / 120.0 + od4 / 720.0 -
                od5 / 5040.0 + start_correction * phi,
            0.5 - od / 3.0 + od2 / 8.0 - od3 / 30.0 + od4 / 144.0 -
                od5 / 840.0 + end_correction * phi,
            -1.0 / 6.0 + od / 12.0 - od2 / 40.0 + od3 / 180.0 - od4 / 1008.0 +
                start_correction * d_phi,
            -1.0 / 3.0 + od / 4.0 - od2 / 10.0 + od3 / 36.0 - od4 / 168.0 +
                end_correction * d_phi,
        };
    }

    /** Endpoint weights for exact single-scatter layer integration.
     *
     * Solar transmission is represented as log-linear between the two layer
     * boundaries while the solar-free scattering source is represented as
     * linear. The returned weights include the observer-facing solar
     * transmission. Their derivatives are with respect to the viewing optical
     * depth and the solar optical depths at the observer-facing (`start`) and
     * far (`end`) boundaries.
     */
    struct SingleScatterSourceWeights {
        double start;
        double end;
        double d_start_d_view_od;
        double d_end_d_view_od;
        double d_start_d_solar_start_od;
        double d_end_d_solar_start_od;
        double d_start_d_solar_end_od;
        double d_end_d_solar_end_od;
    };

    inline SingleScatterSourceWeights single_scatter_source_weights(
        double view_od, double exp_minus_view_od,
        double solar_transmission_start, double solar_transmission_end,
        double start_fraction = 0.5, double end_fraction = 0.5) {
        const double transmission_ratio =
            solar_transmission_start / solar_transmission_end;
        const double solar_od =
            transmission_ratio > 0.0 && std::isfinite(transmission_ratio)
                ? std::log(transmission_ratio)
                : std::log(solar_transmission_start) -
                      std::log(solar_transmission_end);
        const double combined_od = view_od + solar_od;
        double start;
        double end;
        double d_start_d_combined_od;
        double d_end_d_combined_od;

        if (std::abs(combined_od) < 1e-2) {
            const double far_attenuated_transmission =
                solar_transmission_end * exp_minus_view_od;
            const auto weights = normalized_linear_source_weights(
                combined_od,
                far_attenuated_transmission / solar_transmission_start,
                start_fraction, end_fraction);
            start = solar_transmission_start * weights.start;
            end = solar_transmission_start * weights.end;
            d_start_d_combined_od =
                solar_transmission_start * weights.d_start_d_od;
            d_end_d_combined_od = solar_transmission_start * weights.d_end_d_od;
        } else {
            // Scale the exponential terms by the endpoint transmissions before
            // forming the weights. This remains finite when exp(-combined_od)
            // would overflow but its product with the start transmission would
            // not.
            const double far_attenuated_transmission =
                solar_transmission_end * exp_minus_view_od;
            const double scaled_phi =
                (solar_transmission_start - far_attenuated_transmission) /
                combined_od;
            const double scaled_d_phi = solar_transmission_start / combined_od -
                                        scaled_phi * (1.0 + 1.0 / combined_od);
            const double start_standard =
                (solar_transmission_start - scaled_phi) / combined_od;
            const double end_standard =
                (scaled_phi - far_attenuated_transmission) / combined_od;
            const double start_correction = start_fraction - 0.5;
            const double end_correction = end_fraction - 0.5;

            start = start_standard + start_correction * scaled_phi;
            end = end_standard + end_correction * scaled_phi;
            d_start_d_combined_od =
                (-scaled_d_phi - start_standard) / combined_od +
                start_correction * scaled_d_phi;
            d_end_d_combined_od =
                (scaled_d_phi + far_attenuated_transmission - end_standard) /
                    combined_od +
                end_correction * scaled_d_phi;
        }

        return {
            start,
            end,
            d_start_d_combined_od,
            d_end_d_combined_od,
            -start - d_start_d_combined_od,
            -end - d_end_d_combined_od,
            d_start_d_combined_od,
            d_end_d_combined_od,
        };
    }

    template <int NSTOKES, typename SourceType>
    void add_integrated_constant_source(
        const sasktran2::SparseODDualView& shell_od,
        const SourceType& start_source, const SourceType& end_source,
        double start_weight, double end_weight,
        sasktran2::Dual<double, sasktran2::dualstorage::dense, NSTOKES>&
            result_source);

    template <int NSTOKES, typename SourceType>
    void add_integrated_exponential_source(
        const sasktran2::SparseODDualView& shell_od,
        const SourceType& start_source, const SourceType& end_source,
        sasktran2::Dual<double, sasktran2::dualstorage::dense, NSTOKES>&
            result_source);

} // namespace sasktran2::sourcealgo
