#include "sasktran2/math/wigner.h"
#include <pybind11/pybind11.h>
#include <pybind11/stl.h>
#include <pybind11/eigen.h>
#include <sasktran2.h>

namespace py = pybind11;
using namespace pybind11::literals;

void init_math(py::module_& m) {
    py::class_<sasktran2::math::WignerDCalculator>(m, "WignerD")
        .def(py::init<int, int>(), R"(
            Performs calculations of the Wigner (small) d function, :math:`d^l_{m, n}(\theta)`.

            First, this class is constructed for a given `m` and `n`, and then :py:meth:`d` is called
            for a given `l`.

            The Wigner functions are closely related to the associated Legendre polynomials,
            .. math::

                d^l_{m, 0}(\theta) = \sqrt{\frac{(l-m)!}{(l+m)!}} P^m_l(\cos \theta)

            and the regular Legendre polynomials,
            .. math::

                d^l_{0, 0}(\theta) = P_l(\cos \theta)

            Parameters
            ----------
            m: int
                The parameter `m` in :math:`d^l_{m, n}`

            n: int
                The parameter `n` in :math:`d^l_{m, n}`

    )",
             "m"_a, "n"_a)
        .def("d", py::vectorize(&sasktran2::math::WignerDCalculator::d), R"(
            Calculates :math:`d^l_{m, n}(\theta)` for a given `l`, and `m`, `n` provided in the constructor.
            Note that only one of `theta`, `l` can be array-like, one must be scalar.

            Parameters
            ----------
            theta: numpy.ndarray[numpy.float64]
                Angles (in radians) to calculate the function at

            l: numpy.ndarray[numpy.int32]
                The parameter `n` in :math:`d^l_{m, n}`

            Returns
            -------
            np.array
                The calculated Wigner function, either scalar or the same size as `theta` or `l`, whichever is array-like.

        )",
             "theta"_a, "l"_a)

        ;

    m.def("voigt_broaden", &sasktran2::math::spectroscopy::voigt_broaden, R"(
            Calculates the absorption coefficient spectrum for a given set of lines using the Voigt
            line shape.   Most of these parameters are taken directly from the HITRAN database.

            Parameters
            ----------
            line_center: numpy.ndarray[numpy.float64]
                The center of the lines in wavenumbers (cm^-1)

            line_intensity: numpy.ndarray[numpy.float64]
                The intensity of the lines

            lower_energy: numpy.ndarray[numpy.float64]
                The lower energy level of the lines

            gamma_air: numpy.ndarray[numpy.float64]
                The Lorentz broadening due to air

            gamma_self: numpy.ndarray[numpy.float64]
                The Lorentz broadening due to self

            delta_air: numpy.ndarray[numpy.float64]
                The pressure shift due to air

            n_air: numpy.ndarray[numpy.float64]

            iso_id: numpy.ndarray[numpy.int32]
                The identifier for the isotopalog

            partitions: numpy.ndarray[numpy.float64]
                The partition function ratios at the specified temperatures, dimension [ngeo, num_isotop]

            molecular_mass: numpy.ndarray[numpy.float64]
                The molecular mass of each isotopalog, [num_isotop]

            pressure: numpy.ndarray[numpy.float64]
                Partial pressure (1 at 101.325 kPa) at each geometry [ngeo]

            pself: numpy.ndarray[numpy.float64]
                Self partial pressure at each geometry [ngeo], only required for self broadening, but can inform
                the determination of which lines have relevant contributions

            temperature: numpy.ndarray[numpy.float64]
                Temperature in K at each geometry [ngeo]

            wavenumber_grid: numpy.ndarray[numpy.float64]
                The wavenumber grid to produce the result on [wavenumber]

            result: numpy.ndarray[numpy.float64]
                The result matrix, [wavenumber, geometry]

            line_contribution_width: float
                The width of the line of each line to consider, default 25.0

            cull_factor: float
                Lines with less than this total optical depth are culled, default 0.0 indicating no culling
                Should be experimented with to see what is appropriate for the application

            num_threads: int
                The number of threads to use if OMP is enabled, default 1
        )",
          "line_center"_a, "line_intensity"_a, "lower_energy"_a, "gamma_air"_a,
          "gamma_self"_a, "delta_air"_a, "n_air"_a, "iso_id"_a, "partitions"_a,
          "molecular_mass"_a, "pressure"_a, "pself"_a, "temperature"_a,
          "wavenumber_grid"_a, "result"_a, "line_contribution_width"_a = 25.0,
          "cull_factor"_a = 0.0, "num_threads"_a = 1,
          "interpolation_delta"_a = 0.0);
}
