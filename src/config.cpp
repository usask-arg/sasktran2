#include "sasktran2/config.h"
#include <pybind11/pybind11.h>
#include <sasktran2.h>

namespace py = pybind11;

void init_config(py::module_& m) {
    py::enum_<sasktran2::Config::MultipleScatterSource>(m,
                                                        "MultipleScatterSource")
        .value("DiscreteOrdinates",
               sasktran2::Config::MultipleScatterSource::discrete_ordinates)
        .value("TwoStream", sasktran2::Config::MultipleScatterSource::twostream)
        .value("SuccessiveOrders", sasktran2::Config::MultipleScatterSource::hr)
        .value("NoSource", sasktran2::Config::MultipleScatterSource::none)
        .export_values();

    py::enum_<sasktran2::Config::SingleScatterSource>(m, "SingleScatterSource")
        .value("NoSource", sasktran2::Config::SingleScatterSource::none)
        .value("Exact", sasktran2::Config::SingleScatterSource::exact)
        .value("Table", sasktran2::Config::SingleScatterSource::solartable)
        .value("DiscreteOrdinates",
               sasktran2::Config::SingleScatterSource::discrete_ordinates)
        .export_values();

    py::enum_<sasktran2::Config::OccultationSource>(m, "OccultationSource")
        .value("NoSource", sasktran2::Config::OccultationSource::none)
        .value("Standard", sasktran2::Config::OccultationSource::standard)
        .export_values();

    py::enum_<sasktran2::Config::EmissionSource>(m, "EmissionSource")
        .value("NoSource", sasktran2::Config::EmissionSource::none)
        .value("Standard", sasktran2::Config::EmissionSource::standard)
        .export_values();

    py::enum_<sasktran2::Config::StokesBasis>(m, "StokesBasis")
        .value("Standard", sasktran2::Config::StokesBasis::standard)
        .value("Solar", sasktran2::Config::StokesBasis::solar)
        .value("Observer", sasktran2::Config::StokesBasis::observer)
        .export_values();

    py::enum_<sasktran2::Config::ThreadingModel>(m, "ThreadingModel")
        .value("Wavelength", sasktran2::Config::ThreadingModel::wavelength)
        .value("Source", sasktran2::Config::ThreadingModel::source)
        .export_values();

    py::enum_<sasktran2::Config::InputValidationMode>(m, "InputValidationMode")
        .value("Strict", sasktran2::Config::InputValidationMode::strict)
        .value("Standard", sasktran2::Config::InputValidationMode::standard)
        .value("Disabled", sasktran2::Config::InputValidationMode::disabled)
        .export_values();

    py::class_<sasktran2::Config>(m, "Config")
        .def(py::init<>(),
             R"(
                Object which stores all of the configuration settings for the radiative transfer calculation.
            )")
        .def_property("num_threads", &sasktran2::Config::num_threads,
                      &sasktran2::Config::set_num_threads,
                      R"(
                Controls the number of threads used in the calculation.  For maximum performance it is
                recommended to set this to the number of physical cores on your machine.  Defaults to
                1
            )")
        .def_property("threading_model", &sasktran2::Config::threading_model,
                      &sasktran2::Config::set_threading_model,
                      R"(
                Sets the multi-threading mode to use in the calculation.

                `sasktran2.ThreadingModel.Wavelength` (Default)
                    Calculation is multi-threaded over the wavelength (batch) dimension only. This method
                    works very well when this dimension is large, but may increase memory usage. It also
                    is not very effective for a small number of wavelengths.

                `sasktran2.ThreadingModel.Source`
                    Calculation is multi-threaded individually by each source function for each wavelength.
                    This method is recommended when memory is a concern, or when the number of wavelengths
                    is small.
            )")
        .def_property("input_validation_mode",
                      &sasktran2::Config::input_validation_mode,
                      &sasktran2::Config::set_input_validation_mode,
                      R"(
                Sets the input validation mode to use in the calculation.

                `sasktran2.InputValidationMode.Strict` (Default)
                    All input validation checks are performed. This is the recommended mode for most users.

                `sasktran2.InputValidationMode.Standard`
                    Only the most important input validation checks are performed. This mode is recommended
                    for advanced users who are confident in the input data.

                `sasktran2.InputValidationMode.Disabled`
                    No input validation checks are performed. This mode is recommended for advanced users who
                    are confident in the input data and want to maximize performance.
            )")
        .def_property("num_stokes", &sasktran2::Config::num_stokes,
                      &sasktran2::Config::set_num_stokes,
                      R"(
                Sets the number of Stokes parameters used in the calculation. 1 is equivalent to the scalar approximation.
                Currently the only supported values are 1, and 3.
            )")
        .def_property("single_scatter_source",
                      &sasktran2::Config::single_scatter_source,
                      &sasktran2::Config::set_single_scatter_source,
                      R"(
                Sets which (if any) single scatter source is to be used inside the calculation.

                `sasktran2.SingleScatterSource.Exact` (Default)
                    A single scatter source where exact ray tracing is performed at each quadrature
                    point along the observer lines of sight towards the sun

                `sasktran2.SingleScatterSource.Table`
                    A single scatter source where a pre-computed table is used to calculate solar
                    transmission to quadrature points along the line of sight.

                `sasktran2.SingleScatterSource.DiscreteOrdinates`
                    Lets the discrete ordinates source function calculate the single scatter source. Only
                    has an effect if the geometry mode is set to PlaneParallel or PseudoSpherical, and
                    if the DiscreteOrdinates source function is also used for multiple scatter.

                `sasktran2.SingleScatterSource.NoSource`
                    Disables the single scatter source
            )")
        .def_property("occultation_source",
                      &sasktran2::Config::occultation_source,
                      &sasktran2::Config::set_occultation_source,
                      R"(
                Sets which (if any) occultation source is to be used inside the calculation.

                `sasktran2.OccultationSource.NoSource` (Default)
                    No occultation source included

                `sasktran2.OccultationSource.Standard`
                    A constant source of 1 is placed at the end of every individual line of sight.
            )")
        .def_property("emission_source", &sasktran2::Config::emission_source,
                      &sasktran2::Config::set_emission_source,
                      R"(
                Sets which (if any) emission source is to be used inside the calculation.

                `sasktran2.EmissionSource.NoSource` (Default)
                    No emission source included

                `sasktran2.EmissionSource.Standard`
                    An emission source defined on the atmosphere grid is enabled.
            )")
        .def_property("multiple_scatter_source",
                      &sasktran2::Config::multiple_scatter_source,
                      &sasktran2::Config::set_multiple_scatter_source,
                      R"(
                Sets which (if any) multiple scatter source is to be used inside the calculation.

                `sasktran2.MultipleScatterSource.NoSource` (Default)
                    Multiple scattering is disabled

                `sasktran2.MultipleScatterSource.DiscreteOrdinates`
                    The discrete ordinates technique is used to estimate the multiple scatter signal

                `sasktran2.MultipleScatterSource.SuccessiveOrders`
                    The successive orders of scattering method is used to estimate the multiple scatter
                    signal

            )")
        .def_property("stokes_basis", &sasktran2::Config::stokes_basis,
                      &sasktran2::Config::set_stokes_basis,
                      R"(

                      )")
        .def_property("delta_m_scaling",
                      &sasktran2::Config::apply_delta_scaling,
                      &sasktran2::Config::set_apply_delta_scaling,
                      R"(
                Controls whether the delta-M scaling is applied to the calculation.  Defaults to False.
            )")
        .def_property("los_refraction", &sasktran2::Config::los_refraction,
                      &sasktran2::Config::set_los_refraction,
                      R"(
                Controls whether or not refraction is enabled for the observer line of sight rays. Requires
                the refractive index to be set in the Geometry object for refraction to work.  Defaults to False.
            )")
        .def_property("solar_refraction", &sasktran2::Config::solar_refraction,
                      &sasktran2::Config::set_solar_refraction,
                      R"(
                      Controls whether or not refraction is enabled for the solar line of sight rays. Requires
                        the refractive index to be set in the Geometry object for refraction to work.  Only has an effect
                        when the single scatter source term is set to Table.  Defaults to False.
                      )")
        .def_property("multiple_scatter_refraction",
                      &sasktran2::Config::multiple_scatter_refraction,
                      &sasktran2::Config::set_multiple_scatter_refraction,
                      R"(
                    Controls whether or not refraction is enabled for the multiple scatter source. Requires
                    the refractive index to be set in the Geometry object for refraction to work.
                    Only has an effect when the SuccessiveOrders multiple scatter source term is being used.
                    Defaults to False.
                )")
        .def_property("num_sza", &sasktran2::Config::num_do_sza,
                      &sasktran2::Config::set_num_do_sza,
                      R"(
                The number of solar zenith angle discretizations to use when calculating the multiple scatter source.
                For the discrete ordinates source, this determines the number of independent discrete ordinates calculations to perform.
                In the successive orders of scattering source, this is directly the number of discretizations.
                Defaults to 1, indicating that the multiple scatter source is estimated only at the reference point.
            )")
        .def_property("num_successive_orders_iterations",
                      &sasktran2::Config::num_hr_spherical_iterations,
                      &sasktran2::Config::set_num_hr_spherical_iterations,
                      R"(
                The number of iterations to perform when using the successive orders of scattering multiple scatter
                source.
            )")
        .def_property("init_successive_orders_with_discrete_ordinates",
                      &sasktran2::Config::initialize_hr_with_do,
                      &sasktran2::Config::set_initialize_hr_with_do,
                      R"(
                If set to true, when using the successive orders source, it will be initialized with a source calculated with
                the discrete ordinates source instead of a single scattering source.  This greatly reduces the number
                of iterations required for the method, as well as provides a better weighting function approximation.
            )")
        .def_property("num_streams", &sasktran2::Config::num_do_streams,
                      &sasktran2::Config::set_num_do_streams,
                      R"(
                The number of streams to use in the discrete ordinates method. This is the number of streams
                in full space, i.e. each hemisphere has num_streams / 2 angular discretizations.  Must
                be an even number. Default to 16.
            )")
        .def_property("num_forced_azimuth",
                      &sasktran2::Config::num_do_forced_azimuth,
                      &sasktran2::Config::set_num_do_forced_azimuth,
                      R"(
                             If set to a value greater than 0, the discrete ordinates method will use this number of azimuth terms independent of convergence.
                             Defaults to -1, which means to use the number of azimuth terms required for convergence.
                              )")
        .def_property("do_backprop", &sasktran2::Config::do_backprop,
                      &sasktran2::Config::set_do_backprop,
                      R"(
                            Enables backpropagation for the weighting functions when using the DO source in plane parallel or pseudo-spherical geometry.
                            Can greatly improve the computation speed of the calculation when the number of lines of sight is small. Default to True
                      )")
        .def_property("num_successive_orders_points",
                      &sasktran2::Config::num_hr_full_incoming_points,
                      &sasktran2::Config::set_num_hr_full_incoming_points,
                      R"(
                The number of incoming points to use in the successive orders calculation for each solar
                zenith angle.  Must be equal to or less than the number of atmosphere altitude grid points.
                Default is -1 which means to use every altitude grid point.
            )")
        .def_property("num_singlescatter_moments",
                      &sasktran2::Config::num_singlescatter_moments,
                      &sasktran2::Config::set_num_singlescatter_moments,
                      R"(
                The number of Legendre expansion moments to use in the single scatter calculation.
                Must be greater or equal to num_streams. Default to 16.
            )")
        .def_property("num_successive_orders_incoming",
                      &sasktran2::Config::num_hr_incoming,
                      &sasktran2::Config::set_num_hr_incoming,
                      R"(
                The number of integration nodes to use in the successive orders algorithm when calculating the incoming
                radiance at each grid point.  Must be one of [6, 14, 26, 38, 50, 74, 86, 110, 146, 170, 194, 230, 266, 302, 350,
                434, 590, 770, 974, 1202, 1454, 1730, 2030, 2354, 2702, 3074, 3470, 3890]. Default is 110.
            )")
        .def_property("num_successive_orders_outgoing",
                      &sasktran2::Config::num_hr_outgoing,
                      &sasktran2::Config::set_num_hr_outgoing,
                      R"(
                The number of sample points to use in the successive orders algorithm to calculate the outgoing source function on
                radiance at each grid point.  Must be one of [6, 14, 26, 38, 50, 74, 86, 110, 146, 170, 194, 230, 266, 302, 350,
                434, 590, 770, 974, 1202, 1454, 1730, 2030, 2354, 2702, 3074, 3470, 3890]. Default is 110.
            )");
}
