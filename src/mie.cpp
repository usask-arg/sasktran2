#include <pybind11/pybind11.h>
#include <pybind11/stl.h>
#include <pybind11/eigen.h>
#include <sasktran2.h>

namespace py = pybind11;
using namespace pybind11::literals;

void init_mie(py::module_& m) {
    // MieData
    py::class_<sasktran2::mie::MieData>(m, "MieData")
        .def_readonly("Qext", &sasktran2::mie::MieData::Qext, R"(
                Calculated Extinction Efficiency factor [unitless] for given size parameters and refractive index. Shape (size).
            )")
        .def_readonly("Qsca", &sasktran2::mie::MieData::Qsca, R"(
                Calculated Scattering Efficiency factor [unitless] for given size parameters and refractive index. Shape (size).
            )")
        .def_readonly("S1", &sasktran2::mie::MieData::S1, R"(
                Calculated Complex Scattering Amplitude [unitless] in first direction of incident polarization for given size parameters, cos(scattering angles) and refractive index. Shape (size, angle).
            )")
        .def_readonly("S2", &sasktran2::mie::MieData::S2, R"(
                Calculated Complex Scattering Amplitude [unitless] in second direction of incident polarization for given size parameters, cos(scattering angles) and refractive index. Shape (size, angle).
            )");

    // MieOutput
    py::class_<sasktran2::mie::MieOutput>(m, "MieOutput")
        .def_readonly("size_parameter", &sasktran2::mie::MieOutput::size_param,
                      R"(
                Array containing size parameters of spheres (2pi*radius/wavelength). Shape (size).
            )")
        .def_readonly("cos_angles", &sasktran2::mie::MieOutput::cos_angles, R"(
                Array containing the cosine of the scattering angles. Shape (angle).
            )")
        .def_readonly("refractive_index",
                      &sasktran2::mie::MieOutput::refractive_index, R"(
                Complex refractive index of spheres.
            )")
        .def_readonly("values", &sasktran2::mie::MieOutput::values, R"(
                MieData structure containing Extinction Efficiency, Scattering Efficiency and Scattering Amplitudes.
            )");

    py::class_<sasktran2::mie::LinearizedMie>(m, "LinearizedMie")
        .def(py::init<int>(), R"(
            A Mie object created with no input parameters.

            Standard usage is to create a Mie object, and then calculate mie parameters using
            `calculate` method.

            Parameters
            ----------
            num_threads : int
                Number of threads to use for the Mie calculation. Default is 1.

        )",
             "num_threads"_a = 1)
        .def("calculate", &sasktran2::mie::LinearizedMie::calculate,
             R"(
                Performs the Mie computation for an array of size parameters, a single refractive index, and an array that is the cosine of the scattering angles.

                Parameters
                ----------
                size_param : np.ndarray
                    Array of Mie size parameters. Shape (size).
                refractive_index : complex
                    Complex Mie refractive index
                cos_angles : np.ndarray
                    Array of cosine of angles to calculate the scattering amplitude at. Shape (angle).
                calculate_derivative : bool, optional
                    Optional parameter, initiates calculations of derivatives for size parameter and refractive index (not implemented at the moment), by default False

                Returns
                -------
                MieOutput
                    MieOutput that contains the original size parameters, cosine of angles, and refractive index, as well as the calculated mie parameters.

                Examples
                --------

                >>> import sasktran2 as sk
                >>> import numpy as np
                >>> mie = sk.mie.LinearizedMie()
                >>> size_param = np.array([3.0, 4.0, 5.0])
                >>> cos_angles = np.linspace(-1, 1, 100)
                >>> refractive_index = 1.5 + 0.0j
                >>> output = mie.calculate(size_param, refractive_index, cos_angles, True)

                >>> print(output.values.Qext)
                [3.41805617 4.05245221 3.92782673]
                >>> print(output.values.Qsca)
                [3.41805617 4.05245221 3.92782673]

             )",
             "size_param"_a, "refractive_index"_a, "cos_angles"_a,
             "calculate_derivative"_a);
}
