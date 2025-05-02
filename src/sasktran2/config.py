from __future__ import annotations

from sasktran2._core_rust import (
    InputValidationMode,
    PyConfig,
    SingleScatterSource,
    ThreadingModel,
)


class Config:
    _config: PyConfig

    def __init__(self):
        """
        Object which stores all of the configuration settings for the radiative transfer calculation.
        """
        self._config = PyConfig()

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
