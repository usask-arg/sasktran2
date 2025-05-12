from __future__ import annotations

from sasktran2._core_rust import (
    EmissionSource,
    InputValidationMode,
    MultipleScatterSource,
    OccultationSource,
    PyConfig,
    SingleScatterSource,
    StokesBasis,
    ThreadingLib,
    ThreadingModel,
)


class Config:
    _config: PyConfig

    def __init__(self):
        """
        Object which stores all of the configuration settings for the radiative transfer calculation.
        """
        self._config = PyConfig()

    def _into_rust_object(self):
        """
        Returns the rust object which contains all of the configuration settings for the radiative transfer calculation.
        """
        return self._config

    @property
    def num_threads(self) -> int:
        """
        Controls the number of threads used in the calculation.  For maximum performance it is
        recommended to set this to the number of physical cores on your machine.  Defaults to
        1
        """
        return self._config.num_threads

    @num_threads.setter
    def num_threads(self, value: int):
        """
        Controls the number of threads used in the calculation.  For maximum performance it is
        recommended to set this to the number of physical cores on your machine.  Defaults to
        1
        """
        self._config.num_threads = value

    @property
    def threading_model(self) -> ThreadingModel:
        """
        Sets the multi-threading mode to use in the calculation.

        `sasktran2.ThreadingModel.Wavelength` (Default)
            Calculation is multi-threaded over the wavelength (batch) dimension only. This method
            works very well when this dimension is large, but may increase memory usage. It also
            is not very effective for a small number of wavelengths.

        `sasktran2.ThreadingModel.Source`
            Calculation is multi-threaded individually by each source function for each wavelength.
            This method is recommended when memory is a concern, or when the number of wavelengths
            is small.
        """
        return self._config.threading_model

    @threading_model.setter
    def threading_model(self, value: ThreadingModel):
        """
        Sets the multi-threading mode to use in the calculation.

        `sasktran2.ThreadingModel.Wavelength` (Default)
            Calculation is multi-threaded over the wavelength (batch) dimension only. This method
            works very well when this dimension is large, but may increase memory usage. It also
            is not very effective for a small number of wavelengths.

        `sasktran2.ThreadingModel.Source`
            Calculation is multi-threaded individually by each source function for each wavelength.
            This method is recommended when memory is a concern, or when the number of wavelengths
            is small.
        """
        self._config.threading_model = value

    @property
    def input_validation_mode(self) -> InputValidationMode:
        """
        Sets the input validation mode to use in the calculation.

        `sasktran2.InputValidationMode.Strict` (Default)
            All input validation checks are performed. This is the recommended mode for most users.

        `sasktran2.InputValidationMode.Standard`
            Only the most important input validation checks are performed. This mode is recommended
            for advanced users who are confident in the input data.

        `sasktran2.InputValidationMode.Disabled`
            No input validation checks are performed. This mode is recommended for advanced users who
            are confident in the input data and want to maximize performance.
        """
        return self._config.input_validation_mode

    @input_validation_mode.setter
    def input_validation_mode(self, value: InputValidationMode):
        """
        Sets the input validation mode to use in the calculation.

        `sasktran2.InputValidationMode.Strict` (Default)
            All input validation checks are performed. This is the recommended mode for most users.

        `sasktran2.InputValidationMode.Standard`
            Only the most important input validation checks are performed. This mode is recommended
            for advanced users who are confident in the input data.

        `sasktran2.InputValidationMode.Disabled`
            No input validation checks are performed. This mode is recommended for advanced users who
            are confident in the input data and want to maximize performance.
        """
        self._config.input_validation_mode = value

    @property
    def num_stokes(self) -> int:
        """
        Sets the number of Stokes parameters used in the calculation. 1 is equivalent to the scalar approximation.
        Currently the only supported values are 1, and 3.
        """
        return self._config.num_stokes

    @num_stokes.setter
    def num_stokes(self, value: int):
        """
        Sets the number of Stokes parameters used in the calculation. 1 is equivalent to the scalar approximation.
        Currently the only supported values are 1, and 3.
        """
        self._config.num_stokes = value

    @property
    def single_scatter_source(self) -> SingleScatterSource:
        """
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
        """
        return self._config.single_scatter_source

    @single_scatter_source.setter
    def single_scatter_source(self, value: SingleScatterSource):
        """
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
        """
        self._config.single_scatter_source = value

    @property
    def occultation_source(self) -> OccultationSource:
        """
        Sets which (if any) occultation source is to be used inside the calculation.

        `sasktran2.OccultationSource.NoSource` (Default)
            No occultation source included

        `sasktran2.OccultationSource.Standard`
            A constant source of 1 is placed at the end of every individual line of sight.
        """
        return self._config.occultation_source

    @occultation_source.setter
    def occultation_source(self, value: OccultationSource):
        """
        Sets which (if any) occultation source is to be used inside the calculation.

        `sasktran2.OccultationSource.NoSource` (Default)
            No occultation source included

        `sasktran2.OccultationSource.Standard`
            A constant source of 1 is placed at the end of every individual line of sight.
        """
        self._config.occultation_source = value

    @property
    def emission_source(self) -> EmissionSource:
        """
        Sets which (if any) emission source is to be used inside the calculation.

        `sasktran2.EmissionSource.NoSource` (Default)
            No emission source included

        `sasktran2.EmissionSource.Standard`
            An emission source defined on the atmosphere grid is enabled.
        """
        return self._config.emission_source

    @emission_source.setter
    def emission_source(self, value: EmissionSource):
        """
        Sets which (if any) emission source is to be used inside the calculation.

        `sasktran2.EmissionSource.NoSource` (Default)
            No emission source included

        `sasktran2.EmissionSource.Standard`
            An emission source defined on the atmosphere grid is enabled.
        """
        self._config.emission_source = value

    @property
    def multiple_scatter_source(self) -> MultipleScatterSource:
        """
        Sets which (if any) multiple scatter source is to be used inside the calculation.

        `sasktran2.MultipleScatterSource.NoSource` (Default)
            Multiple scattering is disabled

        `sasktran2.MultipleScatterSource.DiscreteOrdinates`
            The discrete ordinates technique is used to estimate the multiple scatter signal

        `sasktran2.MultipleScatterSource.SuccessiveOrders`
            The successive orders of scattering method is used to estimate the multiple scatter
            signal

        """
        return self._config.multiple_scatter_source

    @multiple_scatter_source.setter
    def multiple_scatter_source(self, value: MultipleScatterSource):
        """
        Sets which (if any) multiple scatter source is to be used inside the calculation.

        `sasktran2.MultipleScatterSource.NoSource` (Default)
            Multiple scattering is disabled

        `sasktran2.MultipleScatterSource.DiscreteOrdinates`
            The discrete ordinates technique is used to estimate the multiple scatter signal

        `sasktran2.MultipleScatterSource.SuccessiveOrders`
            The successive orders of scattering method is used to estimate the multiple scatter
            signal

        """
        self._config.multiple_scatter_source = value

    @property
    def stokes_basis(self) -> StokesBasis:
        """ """
        return self._config.stokes_basis

    @stokes_basis.setter
    def stokes_basis(self, value: StokesBasis):
        """ """
        self._config.stokes_basis = value

    @property
    def delta_m_scaling(self) -> bool:
        """
        Controls whether the delta-M scaling is applied to the calculation.  Defaults to False.
        """
        return self._config.delta_m_scaling

    @delta_m_scaling.setter
    def delta_m_scaling(self, value: bool):
        """
        Controls whether the delta-M scaling is applied to the calculation.  Defaults to False.
        """
        self._config.delta_m_scaling = value

    @property
    def los_refraction(self) -> bool:
        """
        Controls whether or not refraction is enabled for the observer line of sight rays. Requires
        the refractive index to be set in the Geometry object for refraction to work.  Defaults to False.
        """
        return self._config.los_refraction

    @los_refraction.setter
    def los_refraction(self, value: bool):
        """
        Controls whether or not refraction is enabled for the observer line of sight rays. Requires
        the refractive index to be set in the Geometry object for refraction to work.  Defaults to False.
        """
        self._config.los_refraction = value

    @property
    def solar_refraction(self) -> bool:
        """
        Controls whether or not refraction is enabled for the solar rays. Requires
        the refractive index to be set in the Geometry object for refraction to work.  Defaults to False.
        """
        return self._config.solar_refraction

    @solar_refraction.setter
    def solar_refraction(self, value: bool):
        """
        Controls whether or not refraction is enabled for the solar rays. Requires
        the refractive index to be set in the Geometry object for refraction to work.  Defaults to False.
        """
        self._config.solar_refraction = value

    @property
    def output_los_optical_depth(self) -> bool:
        """
        If True, then the output from a calculate_radiance call will also contain an "optical_depth" variable
        with dimensions ["los", "wavelength"] that is the optical depth for each line of sight.
        """
        return self._config.output_los_optical_depth

    @output_los_optical_depth.setter
    def output_los_optical_depth(self, value: bool):
        """
        If True, then the output from a calculate_radiance call will also contain an "optical_depth" variable
        with dimensions ["los", "wavelength"] that is the optical depth for each line of sight.
        """
        self._config.output_los_optical_depth = value

    @property
    def multiple_scatter_refraction(self) -> bool:
        """
        Controls whether or not refraction is enabled for the solar line of sight rays. Requires
        the refractive index to be set in the Geometry object for refraction to work.  Only has an effect
        when the single scatter source term is set to Table.  Defaults to False.
        """
        return self._config.multiple_scatter_refraction

    @multiple_scatter_refraction.setter
    def multiple_scatter_refraction(self, value: bool):
        """
        Controls whether or not refraction is enabled for the solar line of sight rays. Requires
        the refractive index to be set in the Geometry object for refraction to work.  Only has an effect
        when the single scatter source term is set to Table.  Defaults to False.
        """
        self._config.multiple_scatter_refraction = value

    @property
    def num_sza(self) -> int:
        """
        The number of solar zenith angle discretizations to use when calculating the multiple scatter source.
        For the discrete ordinates source, this determines the number of independent discrete ordinates calculations to perform.
        In the successive orders of scattering source, this is directly the number of discretizations.
        Defaults to 1, indicating that the multiple scatter source is estimated only at the reference point.
        """
        return self._config.num_sza

    @num_sza.setter
    def num_sza(self, value: int):
        """
        The number of solar zenith angle discretizations to use when calculating the multiple scatter source.
        For the discrete ordinates source, this determines the number of independent discrete ordinates calculations to perform.
        In the successive orders of scattering source, this is directly the number of discretizations.
        Defaults to 1, indicating that the multiple scatter source is estimated only at the reference point.
        """
        self._config.num_sza = value

    @property
    def num_successive_orders_iterations(self) -> int:
        """
        The number of iterations to perform when using the successive orders of scattering multiple scatter
        source.
        """
        return self._config.num_successive_orders_iterations

    @num_successive_orders_iterations.setter
    def num_successive_orders_iterations(self, value: int):
        """
        The number of iterations to perform when using the successive orders of scattering multiple scatter
        source.
        """
        self._config.num_successive_orders_iterations = value

    @property
    def init_successive_orders_with_discrete_ordinates(self) -> bool:
        """
        If set to true, when using the successive orders source, it will be initialized with a source calculated with
        the discrete ordinates source instead of a single scattering source.  This greatly reduces the number
        of iterations required for the method, as well as provides a better weighting function approximation.
        """
        return self._config.init_successive_orders_with_discrete_ordinates

    @init_successive_orders_with_discrete_ordinates.setter
    def init_successive_orders_with_discrete_ordinates(self, value: bool):
        """
        If set to true, when using the successive orders source, it will be initialized with a source calculated with
        the discrete ordinates source instead of a single scattering source.  This greatly reduces the number
        of iterations required for the method, as well as provides a better weighting function approximation.
        """
        self._config.init_successive_orders_with_discrete_ordinates = value

    @property
    def num_streams(self) -> int:
        """
        The number of streams to use in the discrete ordinates method. This is the number of streams
        in full space, i.e. each hemisphere has num_streams / 2 angular discretizations.  Must
        be an even number. Default to 16.
        """
        return self._config.num_streams

    @num_streams.setter
    def num_streams(self, value: int):
        """
        The number of streams to use in the discrete ordinates method. This is the number of streams
        in full space, i.e. each hemisphere has num_streams / 2 angular discretizations.  Must
        be an even number. Default to 16.
        """
        self._config.num_streams = value

    @property
    def num_forced_azimuth(self) -> int:
        """
        If set to a value greater than 0, the discrete ordinates method will use this number of azimuth terms independent of convergence.
        Defaults to -1, which means to use the number of azimuth terms required for convergence.
        """
        return self._config.num_forced_azimuth

    @num_forced_azimuth.setter
    def num_forced_azimuth(self, value: int):
        """
        If set to a value greater than 0, the discrete ordinates method will use this number of azimuth terms independent of convergence.
        Defaults to -1, which means to use the number of azimuth terms required for convergence.
        """
        self._config.num_forced_azimuth = value

    @property
    def do_backprop(self) -> bool:
        """
        Enables backpropagation for the weighting functions when using the DO source in plane parallel or pseudo-spherical geometry.
        Can greatly improve the computation speed of the calculation when the number of lines of sight is small. Default to True
        """
        return self._config.do_backprop

    @do_backprop.setter
    def do_backprop(self, value: bool):
        """
        Enables backpropagation for the weighting functions when using the DO source in plane parallel or pseudo-spherical geometry.
        Can greatly improve the computation speed of the calculation when the number of lines of sight is small. Default to True
        """
        self._config.do_backprop = value

    @property
    def num_successive_order_points(self) -> int:
        """
        The number of incoming points to use in the successive orders calculation for each solar
        zenith angle.  Must be equal to or less than the number of atmosphere altitude grid points.
        Default is -1 which means to use every altitude grid point.
        """
        return self._config.num_successive_order_points

    @num_successive_order_points.setter
    def num_successive_order_points(self, value: int):
        """
        The number of incoming points to use in the successive orders calculation for each solar
        zenith angle.  Must be equal to or less than the number of atmosphere altitude grid points.
        Default is -1 which means to use every altitude grid point.
        """
        self._config.num_successive_order_points = value

    @property
    def num_singlescatter_moments(self) -> int:
        """
        The number of Legendre expansion moments to use in the single scatter calculation.
        Must be greater or equal to num_streams. Default to 16.
        """
        return self._config.num_singlescatter_moments

    @num_singlescatter_moments.setter
    def num_singlescatter_moments(self, value: int):
        """
        The number of Legendre expansion moments to use in the single scatter calculation.
        Must be greater or equal to num_streams. Default to 16.
        """
        self._config.num_singlescatter_moments = value

    @property
    def num_successive_orders_incoming(self) -> int:
        """
        The number of integration nodes to use in the successive orders algorithm when calculating the incoming
        radiance at each grid point.  Must be one of [6, 14, 26, 38, 50, 74, 86, 110, 146, 170, 194, 230, 266, 302, 350,
        434, 590, 770, 974, 1202, 1454, 1730, 2030, 2354, 2702, 3074, 3470, 3890]. Default is 110.
        """
        return self._config.num_successive_orders_incoming

    @num_successive_orders_incoming.setter
    def num_successive_orders_incoming(self, value: int):
        """
        The number of integration nodes to use in the successive orders algorithm when calculating the incoming
        radiance at each grid point.  Must be one of [6, 14, 26, 38, 50, 74, 86, 110, 146, 170, 194, 230, 266, 302, 350,
        434, 590, 770, 974, 1202, 1454, 1730, 2030, 2354, 2702, 3074, 3470, 3890]. Default is 110.
        """
        self._config.num_successive_orders_incoming = value

    @property
    def num_successive_orders_outgoing(self) -> int:
        """
        The number of sample points to use in the successive orders algorithm to calculate the outgoing source function on
        radiance at each grid point.  Must be one of [6, 14, 26, 38, 50, 74, 86, 110, 146, 170, 194, 230, 266, 302, 350,
        434, 590, 770, 974, 1202, 1454, 1730, 2030, 2354, 2702, 3074, 3470, 3890]. Default is 110.
        """
        return self._config.num_successive_orders_outgoing

    @num_successive_orders_outgoing.setter
    def num_successive_orders_outgoing(self, value: int):
        """
        The number of sample points to use in the successive orders algorithm to calculate the outgoing source function on
        radiance at each grid point.  Must be one of [6, 14, 26, 38, 50, 74, 86, 110, 146, 170, 194, 230, 266, 302, 350,
        434, 590, 770, 974, 1202, 1454, 1730, 2030, 2354, 2702, 3074, 3470, 3890]. Default is 110.
        """
        self._config.num_successive_orders_outgoing = value

    @property
    def threading_lib(self) -> ThreadingLib:
        """
        Sets the threading library to use in the calculation.  Default is ThreadingLib.Rayon.
        """
        return self._config.threading_lib

    @threading_lib.setter
    def threading_lib(self, value: ThreadingLib):
        """
        Sets the threading library to use in the calculation.  Default is ThreadingLib.Rayon.
        """
        self._config.threading_lib = value
