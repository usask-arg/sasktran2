#include "sasktran2/config.h"
#include <pybind11/pybind11.h>
#include <sasktran2.h>

namespace py = pybind11;

void init_config(py::module_& m) {
    py::enum_<sasktran2::Config::MultipleScatterSource>(m,
                                                        "MultipleScatterSource")
        .value("DiscreteOrdinates",
               sasktran2::Config::MultipleScatterSource::discrete_ordinates)
        .value("SuccessiveOrders", sasktran2::Config::MultipleScatterSource::hr)
        .value("NoSource", sasktran2::Config::MultipleScatterSource::none)
        .export_values();

    py::enum_<sasktran2::Config::SingleScatterSource>(m, "SingleScatterSource")
        .value("NoSource", sasktran2::Config::SingleScatterSource::none)
        .value("Exact", sasktran2::Config::SingleScatterSource::exact)
        .value("Table", sasktran2::Config::SingleScatterSource::solartable)
        .export_values();

    py::enum_<sasktran2::Config::OccultationSource>(m, "OccultationSource")
        .value("NoSource", sasktran2::Config::OccultationSource::none)
        .value("Standard", sasktran2::Config::OccultationSource::standard)
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
            )");
}
