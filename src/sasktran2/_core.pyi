from __future__ import annotations
import numpy
import typing

__all__ = [
    "AltitudeGrid",
    "AtmosphereStokes_1",
    "AtmosphereStokes_3",
    "AtmosphereStorageStokes_1",
    "AtmosphereStorageStokes_3",
    "BRDFStokes_1",
    "BRDFStokes_3",
    "Config",
    "ConstantSpacing",
    "Coordinates",
    "DerivativeMapping",
    "Disabled",
    "DiscreteOrdinates",
    "Ellipsoidal",
    "EmissionSource",
    "EngineStokes_1",
    "EngineStokes_3",
    "Exact",
    "Geodetic",
    "Geometry1D",
    "GeometryType",
    "GridSpacing",
    "GroundViewingSolar",
    "InputValidationMode",
    "InterpolationMethod",
    "LambertianStokes_1",
    "LambertianStokes_3",
    "LinearInterpolation",
    "LinearSpacing",
    "LinearizedMie",
    "LowerInterpolation",
    "MODISStokes_1",
    "MODISStokes_3",
    "MieData",
    "MieOutput",
    "MultipleScatterSource",
    "NoSource",
    "Observer",
    "OccultationSource",
    "OutOfBoundsExtend",
    "OutOfBoundsPolicy",
    "OutOfBoundsSetZero",
    "OutputDerivMappedStokes_1",
    "OutputDerivMappedStokes_3",
    "OutputIdealStokes_1",
    "OutputIdealStokes_3",
    "OutputStokes_1",
    "OutputStokes_3",
    "PlaneParallel",
    "PseudoSpherical",
    "ShellInterpolation",
    "SingleScatterSource",
    "SnowKokhanovskyStokes_1",
    "SnowKokhanovskyStokes_3",
    "Solar",
    "Source",
    "Spherical",
    "Standard",
    "StokesBasis",
    "Strict",
    "SuccessiveOrders",
    "SurfaceDerivativeMapping",
    "SurfaceStokes_1",
    "SurfaceStokes_3",
    "Table",
    "TangentAltitudeSolar",
    "ThreadingModel",
    "TwoStream",
    "ViewingGeometry",
    "ViewingGeometryBase",
    "Wavelength",
    "WignerD",
    "voigt_broaden",
]

class AltitudeGrid:
    @staticmethod
    def _pybind11_conduit_v1_(*args, **kwargs): ...
    def __init__(
        self,
        arg0: numpy.ndarray[numpy.float64[m, 1]],
        arg1: GridSpacing,
        arg2: OutOfBoundsPolicy,
        arg3: InterpolationMethod,
    ) -> None: ...

class AtmosphereStokes_1:
    storage: ...
    @staticmethod
    def _pybind11_conduit_v1_(*args, **kwargs): ...
    def __init__(
        self, arg0: int, arg1: Geometry1D, arg2: Config, arg3: bool
    ) -> None: ...
    def apply_delta_m_scaling(self, arg0: int) -> None: ...
    @property
    def surface(self) -> ...: ...

class AtmosphereStokes_3:
    storage: ...
    @staticmethod
    def _pybind11_conduit_v1_(*args, **kwargs): ...
    def __init__(
        self, arg0: int, arg1: Geometry1D, arg2: Config, arg3: bool
    ) -> None: ...
    def apply_delta_m_scaling(self, arg0: int) -> None: ...
    @property
    def surface(self) -> ...: ...

class AtmosphereStorageStokes_1:
    d_leg_coeff: numpy.ndarray[numpy.float64[..., ..., ..., ...]]
    emission_source: numpy.ndarray[numpy.float64[m, n]]
    leg_coeff: numpy.ndarray[numpy.float64[..., ..., ...]]
    solar_irradiance: numpy.ndarray[numpy.float64[m, 1]]
    ssa: numpy.ndarray[numpy.float64[m, n]]
    total_extinction: numpy.ndarray[numpy.float64[m, n]]
    @staticmethod
    def _pybind11_conduit_v1_(*args, **kwargs): ...
    def __init__(self, arg0: int, arg1: int, arg2: int) -> None: ...
    def get_derivative_mapping(self, name: str) -> DerivativeMapping: ...
    def resize_derivatives(self, num_deriv: int) -> None: ...
    @property
    def d_f(self) -> numpy.ndarray[numpy.float64[..., ..., ...]]: ...
    @property
    def derivative_mappings(self) -> dict[str, DerivativeMapping]: ...
    @property
    def f(self) -> numpy.ndarray[numpy.float64[m, n]]: ...

class AtmosphereStorageStokes_3:
    d_leg_coeff: numpy.ndarray[numpy.float64[..., ..., ..., ...]]
    emission_source: numpy.ndarray[numpy.float64[m, n]]
    leg_coeff: numpy.ndarray[numpy.float64[..., ..., ...]]
    solar_irradiance: numpy.ndarray[numpy.float64[m, 1]]
    ssa: numpy.ndarray[numpy.float64[m, n]]
    total_extinction: numpy.ndarray[numpy.float64[m, n]]
    @staticmethod
    def _pybind11_conduit_v1_(*args, **kwargs): ...
    def __init__(self, arg0: int, arg1: int, arg2: int) -> None: ...
    def get_derivative_mapping(self, name: str) -> DerivativeMapping: ...
    def resize_derivatives(self, num_deriv: int) -> None: ...
    @property
    def d_f(self) -> numpy.ndarray[numpy.float64[..., ..., ...]]: ...
    @property
    def derivative_mappings(self) -> dict[str, DerivativeMapping]: ...
    @property
    def f(self) -> numpy.ndarray[numpy.float64[m, n]]: ...

class BRDFStokes_1:
    @staticmethod
    def _pybind11_conduit_v1_(*args, **kwargs): ...

class BRDFStokes_3:
    @staticmethod
    def _pybind11_conduit_v1_(*args, **kwargs): ...

class Config:
    stokes_basis: StokesBasis
    @staticmethod
    def _pybind11_conduit_v1_(*args, **kwargs): ...
    def __init__(self) -> None:
        """
        Object which stores all of the configuration settings for the radiative transfer calculation.
        """

    @property
    def delta_m_scaling(self) -> bool:
        """
        Controls whether the delta-M scaling is applied to the calculation.  Defaults to False.
        """

    @delta_m_scaling.setter
    def delta_m_scaling(self, arg1: bool) -> None: ...
    @property
    def do_backprop(self) -> bool:
        """
        Enables backpropagation for the weighting functions when using the DO source in plane parallel or pseudo-spherical geometry.
        Can greatly improve the computation speed of the calculation when the number of lines of sight is small. Default to True
        """

    @do_backprop.setter
    def do_backprop(self, arg1: bool) -> None: ...
    @property
    def emission_source(self) -> EmissionSource:
        """
        Sets which (if any) emission source is to be used inside the calculation.

        `sasktran2.EmissionSource.NoSource` (Default)
            No emission source included

        `sasktran2.EmissionSource.Standard`
            An emission source defined on the atmosphere grid is enabled.
        """

    @emission_source.setter
    def emission_source(self, arg1: EmissionSource) -> None: ...
    @property
    def init_successive_orders_with_discrete_ordinates(self) -> bool:
        """
        If set to true, when using the successive orders source, it will be initialized with a source calculated with
        the discrete ordinates source instead of a single scattering source.  This greatly reduces the number
        of iterations required for the method, as well as provides a better weighting function approximation.
        """

    @init_successive_orders_with_discrete_ordinates.setter
    def init_successive_orders_with_discrete_ordinates(self, arg1: bool) -> None: ...
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

    @input_validation_mode.setter
    def input_validation_mode(self, arg1: InputValidationMode) -> None: ...
    @property
    def los_refraction(self) -> bool:
        """
        Controls whether or not refraction is enabled for the observer line of sight rays. Requires
        the refractive index to be set in the Geometry object for refraction to work.  Defaults to False.
        """

    @los_refraction.setter
    def los_refraction(self, arg1: bool) -> None: ...
    @property
    def multiple_scatter_refraction(self) -> bool:
        """
        Controls whether or not refraction is enabled for the multiple scatter source. Requires
        the refractive index to be set in the Geometry object for refraction to work.
        Only has an effect when the SuccessiveOrders multiple scatter source term is being used.
        Defaults to False.
        """

    @multiple_scatter_refraction.setter
    def multiple_scatter_refraction(self, arg1: bool) -> None: ...
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

    @multiple_scatter_source.setter
    def multiple_scatter_source(self, arg1: MultipleScatterSource) -> None: ...
    @property
    def num_forced_azimuth(self) -> int:
        """
        If set to a value greater than 0, the discrete ordinates method will use this number of azimuth terms independent of convergence.
        Defaults to -1, which means to use the number of azimuth terms required for convergence.
        """

    @num_forced_azimuth.setter
    def num_forced_azimuth(self, arg1: int) -> None: ...
    @property
    def num_singlescatter_moments(self) -> int:
        """
        The number of Legendre expansion moments to use in the single scatter calculation.
        Must be greater or equal to num_streams. Default to 16.
        """

    @num_singlescatter_moments.setter
    def num_singlescatter_moments(self, arg1: int) -> None: ...
    @property
    def num_stokes(self) -> int:
        """
        Sets the number of Stokes parameters used in the calculation. 1 is equivalent to the scalar approximation.
        Currently the only supported values are 1, and 3.
        """

    @num_stokes.setter
    def num_stokes(self, arg1: int) -> None: ...
    @property
    def num_streams(self) -> int:
        """
        The number of streams to use in the discrete ordinates method. This is the number of streams
        in full space, i.e. each hemisphere has num_streams / 2 angular discretizations.  Must
        be an even number. Default to 16.
        """

    @num_streams.setter
    def num_streams(self, arg1: int) -> None: ...
    @property
    def num_successive_orders_incoming(self) -> int:
        """
        The number of integration nodes to use in the successive orders algorithm when calculating the incoming
        radiance at each grid point.  Must be one of [6, 14, 26, 38, 50, 74, 86, 110, 146, 170, 194, 230, 266, 302, 350,
        434, 590, 770, 974, 1202, 1454, 1730, 2030, 2354, 2702, 3074, 3470, 3890]. Default is 110.
        """

    @num_successive_orders_incoming.setter
    def num_successive_orders_incoming(self, arg1: int) -> None: ...
    @property
    def num_successive_orders_iterations(self) -> int:
        """
        The number of iterations to perform when using the successive orders of scattering multiple scatter
        source.
        """

    @num_successive_orders_iterations.setter
    def num_successive_orders_iterations(self, arg1: int) -> None: ...
    @property
    def num_successive_orders_outgoing(self) -> int:
        """
        The number of sample points to use in the successive orders algorithm to calculate the outgoing source function on
        radiance at each grid point.  Must be one of [6, 14, 26, 38, 50, 74, 86, 110, 146, 170, 194, 230, 266, 302, 350,
        434, 590, 770, 974, 1202, 1454, 1730, 2030, 2354, 2702, 3074, 3470, 3890]. Default is 110.
        """

    @num_successive_orders_outgoing.setter
    def num_successive_orders_outgoing(self, arg1: int) -> None: ...
    @property
    def num_successive_orders_points(self) -> int:
        """
        The number of incoming points to use in the successive orders calculation for each solar
        zenith angle.  Must be equal to or less than the number of atmosphere altitude grid points.
        Default is -1 which means to use every altitude grid point.
        """

    @num_successive_orders_points.setter
    def num_successive_orders_points(self, arg1: int) -> None: ...
    @property
    def num_sza(self) -> int:
        """
        The number of solar zenith angle discretizations to use when calculating the multiple scatter source.
        For the discrete ordinates source, this determines the number of independent discrete ordinates calculations to perform.
        In the successive orders of scattering source, this is directly the number of discretizations.
        Defaults to 1, indicating that the multiple scatter source is estimated only at the reference point.
        """

    @num_sza.setter
    def num_sza(self, arg1: int) -> None: ...
    @property
    def num_threads(self) -> int:
        """
        Controls the number of threads used in the calculation.  For maximum performance it is
        recommended to set this to the number of physical cores on your machine.  Defaults to
        1
        """

    @num_threads.setter
    def num_threads(self, arg1: int) -> None: ...
    @property
    def occultation_source(self) -> OccultationSource:
        """
        Sets which (if any) occultation source is to be used inside the calculation.

        `sasktran2.OccultationSource.NoSource` (Default)
            No occultation source included

        `sasktran2.OccultationSource.Standard`
            A constant source of 1 is placed at the end of every individual line of sight.
        """

    @occultation_source.setter
    def occultation_source(self, arg1: OccultationSource) -> None: ...
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

    @single_scatter_source.setter
    def single_scatter_source(self, arg1: SingleScatterSource) -> None: ...
    @property
    def solar_refraction(self) -> bool:
        """
        Controls whether or not refraction is enabled for the solar line of sight rays. Requires
          the refractive index to be set in the Geometry object for refraction to work.  Only has an effect
          when the single scatter source term is set to Table.  Defaults to False.
        """

    @solar_refraction.setter
    def solar_refraction(self, arg1: bool) -> None: ...
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

    @threading_model.setter
    def threading_model(self, arg1: ThreadingModel) -> None: ...

class Coordinates:
    @staticmethod
    def _pybind11_conduit_v1_(*args, **kwargs): ...
    def __init__(
        self, arg0: float, arg1: float, arg2: float, arg3: GeometryType, arg4: bool
    ) -> None: ...

class DerivativeMapping:
    assign_name: str
    interp_dim: str
    log_radiance_space: bool
    scat_deriv_index: int
    @staticmethod
    def _pybind11_conduit_v1_(*args, **kwargs): ...
    def __init__(self, arg0: int, arg1: int, arg2: int) -> None: ...
    def set_zero(self) -> None: ...
    @property
    def d_emission(self) -> numpy.ndarray[numpy.float64[m, n]]: ...
    @property
    def d_extinction(self) -> numpy.ndarray[numpy.float64[m, n]]: ...
    @property
    def d_leg_coeff(self) -> numpy.ndarray[numpy.float64[..., ..., ...]]: ...
    @property
    def d_ssa(self) -> numpy.ndarray[numpy.float64[m, n]]: ...
    @property
    def interpolator(self) -> numpy.ndarray[numpy.float64[m, n]] | None: ...
    @interpolator.setter
    def interpolator(self, arg1: numpy.ndarray[numpy.float64[m, n]]) -> None: ...
    @property
    def is_scattering_derivative(self) -> bool: ...
    @property
    def num_output(self) -> int: ...
    @property
    def scat_factor(self) -> numpy.ndarray[numpy.float64[m, n]]: ...

class EmissionSource:
    """
    Members:

      NoSource

      Standard
    """

    NoSource: typing.ClassVar[EmissionSource]  # value = <EmissionSource.NoSource: 1>
    Standard: typing.ClassVar[EmissionSource]  # value = <EmissionSource.Standard: 0>
    __members__: typing.ClassVar[
        dict[str, EmissionSource]
    ]  # value = {'NoSource': <EmissionSource.NoSource: 1>, 'Standard': <EmissionSource.Standard: 0>}
    @staticmethod
    def _pybind11_conduit_v1_(*args, **kwargs): ...
    def __eq__(self, other: typing.Any) -> bool: ...
    def __getstate__(self) -> int: ...
    def __hash__(self) -> int: ...
    def __index__(self) -> int: ...
    def __init__(self, value: int) -> None: ...
    def __int__(self) -> int: ...
    def __ne__(self, other: typing.Any) -> bool: ...
    def __repr__(self) -> str: ...
    def __setstate__(self, state: int) -> None: ...
    def __str__(self) -> str: ...
    @property
    def name(self) -> str: ...
    @property
    def value(self) -> int: ...

class EngineStokes_1:
    @staticmethod
    def _pybind11_conduit_v1_(*args, **kwargs): ...
    def __init__(
        self,
        config: Config,
        model_geometry: Geometry1D,
        viewing_geometry: ViewingGeometry,
    ) -> None:
        """
        Internal SASKTRAN2 object which handles the radiative transfer calculation.

        Parameters
        ----------
        config: sasktran2.Config
            Configuration settings
        model_geometry: sasktran2.Geometry1D
            The model geometry
        viewing_geometry: sasktan2.ViewingGeometry
            The viewing geometry
        """

    def calculate_radiance(
        self, atmosphere: AtmosphereStokes_1, output: OutputStokes_1
    ) -> None:
        """
        Performs the radiative transfer calculation for the given atmosphere, placing the result in output

        Parameters
        ----------
        atmosphere: sasktran2.Atmosphere
            The atmosphere object

        output: sasktran2.Output
            The result to place the output inside
        """

class EngineStokes_3:
    @staticmethod
    def _pybind11_conduit_v1_(*args, **kwargs): ...
    def __init__(
        self,
        config: Config,
        model_geometry: Geometry1D,
        viewing_geometry: ViewingGeometry,
    ) -> None:
        """
        Internal SASKTRAN2 object which handles the radiative transfer calculation.

        Parameters
        ----------
        config: sasktran2.Config
            Configuration settings
        model_geometry: sasktran2.Geometry1D
            The model geometry
        viewing_geometry: sasktan2.ViewingGeometry
            The viewing geometry
        """

    def calculate_radiance(
        self, atmosphere: AtmosphereStokes_3, output: OutputStokes_3
    ) -> None:
        """
        Performs the radiative transfer calculation for the given atmosphere, placing the result in output

        Parameters
        ----------
        atmosphere: sasktran2.Atmosphere
            The atmosphere object

        output: sasktran2.Output
            The result to place the output inside
        """

class Geodetic:
    @staticmethod
    def _pybind11_conduit_v1_(*args, **kwargs): ...
    def __init__(self, equatorial_radius: float, flattening_factor: float) -> None:
        """
        A geodetic object based on a given equatorial (semi-major) radius and flattening factor.

        Standard usage is to create a geodetic object, and then initialize it through one of the
        `from_*` methods.

        Parameters
        ----------
        equatorial_radius: float
            Radius at the equator (semi-major) of the ellipsoid.
        flattening_factor: float
            Flattening factor of the ellipsoid. This is defined as (a-b)/a, where a is the semi-major axis
            and b is the semi-minor radius
        """

    def altitude_intercepts(
        self,
        altitude: float,
        observer: numpy.ndarray[numpy.float64[3, 1]],
        look_vector: numpy.ndarray[numpy.float64[3, 1]],
    ) -> tuple[numpy.ndarray[numpy.float64[3, 1]], numpy.ndarray[numpy.float64[3, 1]]]:
        """
                        Calculate the two intersections of a line of sight and an altitude.

                        Parameters
                        ----------
                        altitude : float
                            Altitude in meters.
                        observer : np.ndarray
                            Three element array containing the obsever position in geocentric coordinates.
                        look_vector : np.ndarray
                            Three element array containing a normalized look vector.

                        Returns
                        -------
                        np.ndarray
                            Three element array containing the first (entering) intercept in geocentric coordinates.
                        np.ndarray
                            Three element array containing the second (exiting) intercept in geocentric coordinates.

                        Examples
                        --------
                        >>> import sasktran2 as sk
                        >>> import numpy as np
                        >>> geodetic = sk.WGS84()
                        >>> look = geodetic.from_tangent_altitude(15322, [3.676013154788849600e+005, 1.009976313640051500e+006, \
                                                                    -6.871601202127538600e+006], [0, 0, 1])
                        >>> obs = geodetic.location
                        >>> intercept1, intercept2 = geodetic.altitude_intercepts(16000, obs, look)
                        >>> print(np.array_str(intercept1, precision=3))
                        [ 1147302.059  3152186.5   -5425360.027]
                        >>> print(np.array_str(intercept2, precision=3))
                        [ 1201098.489  3299990.978 -5325574.803]
        """

    def from_lat_lon_alt(
        self, latitude: float, longitude: float, altitude: float
    ) -> None:
        """
        Initializes the Geodetic based on a specifiec latitude, longitude, and altitude.

        Parameters
        ----------
        latitude : float
            Latitude in degrees (-90 to 90)
        longitude : float
            Longitude in degrees (0 to 360 or -180 to 180)
        altitude : float
            Altitude above the geoid in metres

        Examples
        --------
        >>> import sasktran2 as sk
        >>> geodetic = sk.WGS84()
        >>> geodetic.from_lat_lon_alt(latitude=-15, longitude=-20, altitude=7342)
        >>> print(geodetic)
        WGS84 Location:
        Latitude: -15.0, Longitude: 340.0, Altitude: 7342.0
        """

    def from_tangent_altitude(
        self,
        altitude: float,
        observer: numpy.ndarray[numpy.float64[3, 1]],
        boresight: numpy.ndarray[numpy.float64[3, 1]],
    ) -> numpy.ndarray[numpy.float64[3, 1]]:
        """
                        Initialized the Geodetic from a specified tangent altitude, obsever location, and bore sight plane.

                        Parameters
                        ----------
                        altitude : float
                            Tangent altitude in meters
                        observer : np.ndarray
                            Three element array containing the obsever position in geocentric coordinates
                        boresight : np.ndarray
                            Three element array containing a normalized look vector that is within the bore sight plane.

                        Returns
                        -------
                        np.ndarray
                            Three element array containing the normalized look vector to the tangent point.

                        Examples
                        --------
                        >>> import sasktran2 as sk
                        >>> geodetic = sk.WGS84()
                        >>> look = geodetic.from_tangent_altitude(15322, [ 3.676013154788849600e+005, 1.009976313640051500e+006,\
                                                                          -6.871601202127538600e+006], [0, 0, 1])
                        >>> print(look)
                        [0.28880556 0.79348676 0.53569591]
                        >>> print(geodetic)
                        WGS84 Location:
                        Latitude: -57.60888188776806, Longitude: 70.00000000000001, Altitude: 15321.971935882739
        """

    def from_tangent_point(
        self,
        observer: numpy.ndarray[numpy.float64[3, 1]],
        look_vector: numpy.ndarray[numpy.float64[3, 1]],
    ) -> None:
        """
                        Initializes  the Geodetic by calculating the tangent point from an observer position and look vector

                        Parameters
                        ----------
                        observer : np.ndarray
                            Three element array containing the observer position in geocentric coordinates
                        look_vector : np.ndarray
                            Three element array containing a normalized look vector

                        Examples
                        --------
                        >>> import sasktran2 as sk
                        >>> geodetic = sk.WGS84()
                        >>> geodetic.from_tangent_point([ 3.676013154788849600e+005, 1.009976313640051500e+006,\
                                                         -6.871601202127538600e+006], [ 2.884568631765662100e-001,\
                                                          7.925287180643269000e-001,  5.372996083468238900e-001])
                        >>> print(geodetic)
                        WGS84 Location:
                        Latitude: -57.500000192733594, Longitude: 70.0, Altitude: 10002.99586173162
        """

    def from_xyz(self, location: numpy.ndarray[numpy.float64[3, 1]]) -> None:
        """
        Initializes the Geodetic from a geocentric location

        Parameters
        ----------
        location : np.ndarray
            Three element vector containing a location in geocentric coordinates

        Examples
        --------
        >>> import sasktran2 as sk
        >>> geodetic = sk.WGS84()
        >>> geodetic.from_xyz([ 5797230.47518212, -2110019.3341472, -1642001.16317228])
        >>> print(geodetic)
        WGS84 Location:
        Latitude: -14.999999973747736, Longitude: 340.00000000000006, Altitude: 7344.999610390202
        """

    @property
    def altitude(self) -> float:
        """
        Altitude in [m] above the surface of the ellipsoid.
        """

    @property
    def latitude(self) -> float:
        """
        Geodetic latitude in degrees
        """

    @property
    def local_south(self) -> numpy.ndarray[numpy.float64[3, 1]]:
        """
        A unit vector pointing in the local south direction
        """

    @property
    def local_up(self) -> numpy.ndarray[numpy.float64[3, 1]]:
        """
        A unit vector pointing up (perpindicular to the ellipsoidal surface)
        """

    @property
    def local_west(self) -> numpy.ndarray[numpy.float64[3, 1]]:
        """
        A unit vector pointing in the local west direction
        """

    @property
    def location(self) -> numpy.ndarray[numpy.float64[3, 1]]:
        """
        Geocentric location in cartesian coordinates
        """

    @property
    def longitude(self) -> float:
        """
        Geodetic longitude in degrees
        """

    @property
    def valid(self) -> bool:
        """
        True if the geodetic object has been initialized, False otherwise.
        """

class Geometry1D:
    @staticmethod
    def _pybind11_conduit_v1_(*args, **kwargs): ...
    def __init__(
        self,
        cos_sza: float,
        solar_azimuth: float,
        earth_radius_m: float,
        altitude_grid_m: numpy.ndarray[numpy.float64[m, 1]],
        interpolation_method: InterpolationMethod,
        geometry_type: GeometryType,
    ) -> None:
        """
        Initializes a geometry where the atmosphere varies only in 1 dimension (altitude).  The reference point
        is defined by solar angles at the reference point.

        Parameters
        ----------
        cos_sza: float
            Cosine of solar zenith angle at the reference point
        solar_azimuth: float
            Solar azimuth angle at the reference point.
        earth_radius_m: float
            Radius of the earth.  Only has an effect if geometry_type is not set to PlaneParallel
        altitude_grid_m: np.array
            One dimensional altitude grid
        interpolation_method: sasktran2.InterpolationMethod
            The interpolation method to use in-between geometry grid points

            `sasktran2.InterpolationMethod.LinearInterpolation`
                In-between grid points, linear interpolation is assumed.  This means that Atmospheric quantities
                such as extinction, single scatter albedo, should be thought of as sampled on the geometry grid points.

            `sasktran2.InterpolationMethod.ShellInterpolation`
                Atmospheric quantities such as extinction, single scatter albedo, are assumed to be constant in-between
                geometry grid points.

        geometry_type: sasktran2.GeometryType
            The global geometry type used inside the radiative transfer calculation.

            `sasktran2.GeometryType.Spherical`
                All aspects of the calculation are done using spherical geometry.

            `sasktran2.GeometryType.PlaneParallel`
                All aspects of the calculation are done using plane-parallel geometry.

            `sasktran2.GeometryType.PseudoSpherical`
                Line of sight integration and the multiple scatter calculation is done using
                plane parallel geometry, however the initial solar source function is calculated
                using a spherical geometry.
        """

    def altitudes(self) -> numpy.ndarray[numpy.float64[m, 1]]: ...
    @property
    def refractive_index(self) -> numpy.ndarray[numpy.float64[m, 1]]:
        """
        The refractive index of the atmosphere.  This is used to calculate refraction in the radiative transfer calculation.
        Defaults to 1.0 which indicates no refractive effects.  Only has an effect if the refraction configuration options are
        enabled in the `sasktran2.Config` object.
        """

    @refractive_index.setter
    def refractive_index(self, arg1: numpy.ndarray[numpy.float64[m, 1]]) -> None: ...

class GeometryType:
    """
    Members:

      PlaneParallel

      Spherical

      PseudoSpherical

      Ellipsoidal
    """

    Ellipsoidal: typing.ClassVar[GeometryType]  # value = <GeometryType.Ellipsoidal: 3>
    PlaneParallel: typing.ClassVar[
        GeometryType
    ]  # value = <GeometryType.PlaneParallel: 0>
    PseudoSpherical: typing.ClassVar[
        GeometryType
    ]  # value = <GeometryType.PseudoSpherical: 1>
    Spherical: typing.ClassVar[GeometryType]  # value = <GeometryType.Spherical: 2>
    __members__: typing.ClassVar[
        dict[str, GeometryType]
    ]  # value = {'PlaneParallel': <GeometryType.PlaneParallel: 0>, 'Spherical': <GeometryType.Spherical: 2>, 'PseudoSpherical': <GeometryType.PseudoSpherical: 1>, 'Ellipsoidal': <GeometryType.Ellipsoidal: 3>}
    @staticmethod
    def _pybind11_conduit_v1_(*args, **kwargs): ...
    def __eq__(self, other: typing.Any) -> bool: ...
    def __getstate__(self) -> int: ...
    def __hash__(self) -> int: ...
    def __index__(self) -> int: ...
    def __init__(self, value: int) -> None: ...
    def __int__(self) -> int: ...
    def __ne__(self, other: typing.Any) -> bool: ...
    def __repr__(self) -> str: ...
    def __setstate__(self, state: int) -> None: ...
    def __str__(self) -> str: ...
    @property
    def name(self) -> str: ...
    @property
    def value(self) -> int: ...

class GridSpacing:
    """
    Members:

      ConstantSpacing

      LinearSpacing
    """

    ConstantSpacing: typing.ClassVar[
        GridSpacing
    ]  # value = <GridSpacing.ConstantSpacing: 0>
    LinearSpacing: typing.ClassVar[
        GridSpacing
    ]  # value = <GridSpacing.LinearSpacing: 1>
    __members__: typing.ClassVar[
        dict[str, GridSpacing]
    ]  # value = {'ConstantSpacing': <GridSpacing.ConstantSpacing: 0>, 'LinearSpacing': <GridSpacing.LinearSpacing: 1>}
    @staticmethod
    def _pybind11_conduit_v1_(*args, **kwargs): ...
    def __eq__(self, other: typing.Any) -> bool: ...
    def __getstate__(self) -> int: ...
    def __hash__(self) -> int: ...
    def __index__(self) -> int: ...
    def __init__(self, value: int) -> None: ...
    def __int__(self) -> int: ...
    def __ne__(self, other: typing.Any) -> bool: ...
    def __repr__(self) -> str: ...
    def __setstate__(self, state: int) -> None: ...
    def __str__(self) -> str: ...
    @property
    def name(self) -> str: ...
    @property
    def value(self) -> int: ...

class GroundViewingSolar(ViewingGeometryBase):
    @staticmethod
    def _pybind11_conduit_v1_(*args, **kwargs): ...
    def __init__(
        self,
        cos_sza: float,
        relative_azimuth: float,
        cos_viewing_zenith: float,
        observer_altitude_m: float,
    ) -> None:
        """
        Defines a viewing ray that is looking at the ground from angles defined at the ground location. Note that
        all of these parameters assumes straight line paths (i.e. no atmospheric refraction)

        Parameters
        ----------
        cos_sza: float
            Cosine of solar zenith angle at the ground point [unitless]
        relative_azimuth: float
            Relative azimuth angle to the sun [rad] at the ground point. An angle of 0 degrees corresponds to the forward scattering plane.
        observer_altitude_m: float
            Observer altitude relative to the earth [m]
        cos_viewing_zenith: float
            Cosine of the viewing zenith angle at the ground point [unitless]
        """

class InputValidationMode:
    """
    Members:

      Strict

      Standard

      Disabled
    """

    Disabled: typing.ClassVar[
        InputValidationMode
    ]  # value = <InputValidationMode.Disabled: 2>
    Standard: typing.ClassVar[
        InputValidationMode
    ]  # value = <InputValidationMode.Standard: 1>
    Strict: typing.ClassVar[
        InputValidationMode
    ]  # value = <InputValidationMode.Strict: 0>
    __members__: typing.ClassVar[
        dict[str, InputValidationMode]
    ]  # value = {'Strict': <InputValidationMode.Strict: 0>, 'Standard': <InputValidationMode.Standard: 1>, 'Disabled': <InputValidationMode.Disabled: 2>}
    @staticmethod
    def _pybind11_conduit_v1_(*args, **kwargs): ...
    def __eq__(self, other: typing.Any) -> bool: ...
    def __getstate__(self) -> int: ...
    def __hash__(self) -> int: ...
    def __index__(self) -> int: ...
    def __init__(self, value: int) -> None: ...
    def __int__(self) -> int: ...
    def __ne__(self, other: typing.Any) -> bool: ...
    def __repr__(self) -> str: ...
    def __setstate__(self, state: int) -> None: ...
    def __str__(self) -> str: ...
    @property
    def name(self) -> str: ...
    @property
    def value(self) -> int: ...

class InterpolationMethod:
    """
    Members:

      ShellInterpolation

      LinearInterpolation

      LowerInterpolation
    """

    LinearInterpolation: typing.ClassVar[
        InterpolationMethod
    ]  # value = <InterpolationMethod.LinearInterpolation: 1>
    LowerInterpolation: typing.ClassVar[
        InterpolationMethod
    ]  # value = <InterpolationMethod.LowerInterpolation: 2>
    ShellInterpolation: typing.ClassVar[
        InterpolationMethod
    ]  # value = <InterpolationMethod.ShellInterpolation: 0>
    __members__: typing.ClassVar[
        dict[str, InterpolationMethod]
    ]  # value = {'ShellInterpolation': <InterpolationMethod.ShellInterpolation: 0>, 'LinearInterpolation': <InterpolationMethod.LinearInterpolation: 1>, 'LowerInterpolation': <InterpolationMethod.LowerInterpolation: 2>}
    @staticmethod
    def _pybind11_conduit_v1_(*args, **kwargs): ...
    def __eq__(self, other: typing.Any) -> bool: ...
    def __getstate__(self) -> int: ...
    def __hash__(self) -> int: ...
    def __index__(self) -> int: ...
    def __init__(self, value: int) -> None: ...
    def __int__(self) -> int: ...
    def __ne__(self, other: typing.Any) -> bool: ...
    def __repr__(self) -> str: ...
    def __setstate__(self, state: int) -> None: ...
    def __str__(self) -> str: ...
    @property
    def name(self) -> str: ...
    @property
    def value(self) -> int: ...

class LambertianStokes_1(BRDFStokes_1):
    @staticmethod
    def _pybind11_conduit_v1_(*args, **kwargs): ...
    def __init__(self) -> None: ...
    @property
    def num_deriv(self) -> int:
        """
        Number of derivatives this BRDF will calculate.
        """

class LambertianStokes_3(BRDFStokes_3):
    @staticmethod
    def _pybind11_conduit_v1_(*args, **kwargs): ...
    def __init__(self) -> None: ...
    @property
    def num_deriv(self) -> int:
        """
        Number of derivatives this BRDF will calculate.
        """

class LinearizedMie:
    @staticmethod
    def _pybind11_conduit_v1_(*args, **kwargs): ...
    def __init__(self, num_threads: int = 1) -> None:
        """
        A Mie object created with no input parameters.

        Standard usage is to create a Mie object, and then calculate mie parameters using
        `calculate` method.

        Parameters
        ----------
        num_threads : int
            Number of threads to use for the Mie calculation. Default is 1.
        """

    def calculate(
        self,
        size_param: numpy.ndarray[numpy.float64[m, 1]],
        refractive_index: complex,
        cos_angles: numpy.ndarray[numpy.float64[m, 1]],
        calculate_derivative: bool,
    ) -> MieOutput:
        """
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
        """

class MODISStokes_1(BRDFStokes_1):
    @staticmethod
    def _pybind11_conduit_v1_(*args, **kwargs): ...
    def __init__(self) -> None: ...
    @property
    def num_deriv(self) -> int:
        """
        Number of derivatives this BRDF will calculate.
        """

class MODISStokes_3(BRDFStokes_3):
    @staticmethod
    def _pybind11_conduit_v1_(*args, **kwargs): ...
    def __init__(self) -> None: ...
    @property
    def num_deriv(self) -> int:
        """
        Number of derivatives this BRDF will calculate.
        """

class MieData:
    @staticmethod
    def _pybind11_conduit_v1_(*args, **kwargs): ...
    @property
    def Qext(self) -> numpy.ndarray[numpy.float64[m, 1]]:
        """
        Calculated Extinction Efficiency factor [unitless] for given size parameters and refractive index. Shape (size).
        """

    @property
    def Qsca(self) -> numpy.ndarray[numpy.float64[m, 1]]:
        """
        Calculated Scattering Efficiency factor [unitless] for given size parameters and refractive index. Shape (size).
        """

    @property
    def S1(self) -> numpy.ndarray[numpy.complex128[m, n]]:
        """
        Calculated Complex Scattering Amplitude [unitless] in first direction of incident polarization for given size parameters, cos(scattering angles) and refractive index. Shape (size, angle).
        """

    @property
    def S2(self) -> numpy.ndarray[numpy.complex128[m, n]]:
        """
        Calculated Complex Scattering Amplitude [unitless] in second direction of incident polarization for given size parameters, cos(scattering angles) and refractive index. Shape (size, angle).
        """

class MieOutput:
    @staticmethod
    def _pybind11_conduit_v1_(*args, **kwargs): ...
    @property
    def cos_angles(self) -> numpy.ndarray[numpy.float64[m, 1]]:
        """
        Array containing the cosine of the scattering angles. Shape (angle).
        """

    @property
    def refractive_index(self) -> complex:
        """
        Complex refractive index of spheres.
        """

    @property
    def size_parameter(self) -> numpy.ndarray[numpy.float64[m, 1]]:
        """
        Array containing size parameters of spheres (2pi*radius/wavelength). Shape (size).
        """

    @property
    def values(self) -> MieData:
        """
        MieData structure containing Extinction Efficiency, Scattering Efficiency and Scattering Amplitudes.
        """

class MultipleScatterSource:
    """
    Members:

      DiscreteOrdinates

      TwoStream

      SuccessiveOrders

      NoSource
    """

    DiscreteOrdinates: typing.ClassVar[
        MultipleScatterSource
    ]  # value = <MultipleScatterSource.DiscreteOrdinates: 0>
    NoSource: typing.ClassVar[
        MultipleScatterSource
    ]  # value = <MultipleScatterSource.NoSource: 3>
    SuccessiveOrders: typing.ClassVar[
        MultipleScatterSource
    ]  # value = <MultipleScatterSource.SuccessiveOrders: 1>
    TwoStream: typing.ClassVar[
        MultipleScatterSource
    ]  # value = <MultipleScatterSource.TwoStream: 2>
    __members__: typing.ClassVar[
        dict[str, MultipleScatterSource]
    ]  # value = {'DiscreteOrdinates': <MultipleScatterSource.DiscreteOrdinates: 0>, 'TwoStream': <MultipleScatterSource.TwoStream: 2>, 'SuccessiveOrders': <MultipleScatterSource.SuccessiveOrders: 1>, 'NoSource': <MultipleScatterSource.NoSource: 3>}
    @staticmethod
    def _pybind11_conduit_v1_(*args, **kwargs): ...
    def __eq__(self, other: typing.Any) -> bool: ...
    def __getstate__(self) -> int: ...
    def __hash__(self) -> int: ...
    def __index__(self) -> int: ...
    def __init__(self, value: int) -> None: ...
    def __int__(self) -> int: ...
    def __ne__(self, other: typing.Any) -> bool: ...
    def __repr__(self) -> str: ...
    def __setstate__(self, state: int) -> None: ...
    def __str__(self) -> str: ...
    @property
    def name(self) -> str: ...
    @property
    def value(self) -> int: ...

class OccultationSource:
    """
    Members:

      NoSource

      Standard
    """

    NoSource: typing.ClassVar[
        OccultationSource
    ]  # value = <OccultationSource.NoSource: 1>
    Standard: typing.ClassVar[
        OccultationSource
    ]  # value = <OccultationSource.Standard: 0>
    __members__: typing.ClassVar[
        dict[str, OccultationSource]
    ]  # value = {'NoSource': <OccultationSource.NoSource: 1>, 'Standard': <OccultationSource.Standard: 0>}
    @staticmethod
    def _pybind11_conduit_v1_(*args, **kwargs): ...
    def __eq__(self, other: typing.Any) -> bool: ...
    def __getstate__(self) -> int: ...
    def __hash__(self) -> int: ...
    def __index__(self) -> int: ...
    def __init__(self, value: int) -> None: ...
    def __int__(self) -> int: ...
    def __ne__(self, other: typing.Any) -> bool: ...
    def __repr__(self) -> str: ...
    def __setstate__(self, state: int) -> None: ...
    def __str__(self) -> str: ...
    @property
    def name(self) -> str: ...
    @property
    def value(self) -> int: ...

class OutOfBoundsPolicy:
    """
    Members:

      OutOfBoundsExtend

      OutOfBoundsSetZero
    """

    OutOfBoundsExtend: typing.ClassVar[
        OutOfBoundsPolicy
    ]  # value = <OutOfBoundsPolicy.OutOfBoundsExtend: 0>
    OutOfBoundsSetZero: typing.ClassVar[
        OutOfBoundsPolicy
    ]  # value = <OutOfBoundsPolicy.OutOfBoundsSetZero: 1>
    __members__: typing.ClassVar[
        dict[str, OutOfBoundsPolicy]
    ]  # value = {'OutOfBoundsExtend': <OutOfBoundsPolicy.OutOfBoundsExtend: 0>, 'OutOfBoundsSetZero': <OutOfBoundsPolicy.OutOfBoundsSetZero: 1>}
    @staticmethod
    def _pybind11_conduit_v1_(*args, **kwargs): ...
    def __eq__(self, other: typing.Any) -> bool: ...
    def __getstate__(self) -> int: ...
    def __hash__(self) -> int: ...
    def __index__(self) -> int: ...
    def __init__(self, value: int) -> None: ...
    def __int__(self) -> int: ...
    def __ne__(self, other: typing.Any) -> bool: ...
    def __repr__(self) -> str: ...
    def __setstate__(self, state: int) -> None: ...
    def __str__(self) -> str: ...
    @property
    def name(self) -> str: ...
    @property
    def value(self) -> int: ...

class OutputDerivMappedStokes_1(OutputStokes_1):
    @staticmethod
    def _pybind11_conduit_v1_(*args, **kwargs): ...
    def __init__(self) -> None: ...
    @property
    def deriv_map(self) -> dict[str, numpy.ndarray[numpy.float64[m, n]]]: ...
    @property
    def radiance(self) -> numpy.ndarray[numpy.float64[m, 1]]: ...
    @property
    def surface_deriv_map(self) -> dict[str, numpy.ndarray[numpy.float64[m, n]]]: ...

class OutputDerivMappedStokes_3(OutputStokes_3):
    @staticmethod
    def _pybind11_conduit_v1_(*args, **kwargs): ...
    def __init__(self) -> None: ...
    @property
    def deriv_map(self) -> dict[str, numpy.ndarray[numpy.float64[m, n]]]: ...
    @property
    def radiance(self) -> numpy.ndarray[numpy.float64[m, 1]]: ...
    @property
    def surface_deriv_map(self) -> dict[str, numpy.ndarray[numpy.float64[m, n]]]: ...

class OutputIdealStokes_1(OutputStokes_1):
    @staticmethod
    def _pybind11_conduit_v1_(*args, **kwargs): ...
    def __init__(self) -> None: ...
    @property
    def d_radiance(self) -> numpy.ndarray[numpy.float64[m, n]]: ...
    @property
    def radiance(self) -> numpy.ndarray[numpy.float64[m, 1]]: ...

class OutputIdealStokes_3(OutputStokes_3):
    @staticmethod
    def _pybind11_conduit_v1_(*args, **kwargs): ...
    def __init__(self) -> None: ...
    @property
    def d_radiance(self) -> numpy.ndarray[numpy.float64[m, n]]: ...
    @property
    def radiance(self) -> numpy.ndarray[numpy.float64[m, 1]]: ...

class OutputStokes_1:
    @staticmethod
    def _pybind11_conduit_v1_(*args, **kwargs): ...

class OutputStokes_3:
    @staticmethod
    def _pybind11_conduit_v1_(*args, **kwargs): ...

class SingleScatterSource:
    """
    Members:

      NoSource

      Exact

      Table

      DiscreteOrdinates
    """

    DiscreteOrdinates: typing.ClassVar[
        SingleScatterSource
    ]  # value = <SingleScatterSource.DiscreteOrdinates: 2>
    Exact: typing.ClassVar[
        SingleScatterSource
    ]  # value = <SingleScatterSource.Exact: 0>
    NoSource: typing.ClassVar[
        SingleScatterSource
    ]  # value = <SingleScatterSource.NoSource: 3>
    Table: typing.ClassVar[
        SingleScatterSource
    ]  # value = <SingleScatterSource.Table: 1>
    __members__: typing.ClassVar[
        dict[str, SingleScatterSource]
    ]  # value = {'NoSource': <SingleScatterSource.NoSource: 3>, 'Exact': <SingleScatterSource.Exact: 0>, 'Table': <SingleScatterSource.Table: 1>, 'DiscreteOrdinates': <SingleScatterSource.DiscreteOrdinates: 2>}
    @staticmethod
    def _pybind11_conduit_v1_(*args, **kwargs): ...
    def __eq__(self, other: typing.Any) -> bool: ...
    def __getstate__(self) -> int: ...
    def __hash__(self) -> int: ...
    def __index__(self) -> int: ...
    def __init__(self, value: int) -> None: ...
    def __int__(self) -> int: ...
    def __ne__(self, other: typing.Any) -> bool: ...
    def __repr__(self) -> str: ...
    def __setstate__(self, state: int) -> None: ...
    def __str__(self) -> str: ...
    @property
    def name(self) -> str: ...
    @property
    def value(self) -> int: ...

class SnowKokhanovskyStokes_1(BRDFStokes_1):
    @staticmethod
    def _pybind11_conduit_v1_(*args, **kwargs): ...
    def __init__(self) -> None: ...
    @property
    def num_deriv(self) -> int:
        """
        Number of derivatives this BRDF will calculate.
        """

class SnowKokhanovskyStokes_3(BRDFStokes_3):
    @staticmethod
    def _pybind11_conduit_v1_(*args, **kwargs): ...
    def __init__(self) -> None: ...
    @property
    def num_deriv(self) -> int:
        """
        Number of derivatives this BRDF will calculate.
        """

class StokesBasis:
    """
    Members:

      Standard

      Solar

      Observer
    """

    Observer: typing.ClassVar[StokesBasis]  # value = <StokesBasis.Observer: 2>
    Solar: typing.ClassVar[StokesBasis]  # value = <StokesBasis.Solar: 1>
    Standard: typing.ClassVar[StokesBasis]  # value = <StokesBasis.Standard: 0>
    __members__: typing.ClassVar[
        dict[str, StokesBasis]
    ]  # value = {'Standard': <StokesBasis.Standard: 0>, 'Solar': <StokesBasis.Solar: 1>, 'Observer': <StokesBasis.Observer: 2>}
    @staticmethod
    def _pybind11_conduit_v1_(*args, **kwargs): ...
    def __eq__(self, other: typing.Any) -> bool: ...
    def __getstate__(self) -> int: ...
    def __hash__(self) -> int: ...
    def __index__(self) -> int: ...
    def __init__(self, value: int) -> None: ...
    def __int__(self) -> int: ...
    def __ne__(self, other: typing.Any) -> bool: ...
    def __repr__(self) -> str: ...
    def __setstate__(self, state: int) -> None: ...
    def __str__(self) -> str: ...
    @property
    def name(self) -> str: ...
    @property
    def value(self) -> int: ...

class SurfaceDerivativeMapping:
    interp_dim: str
    @staticmethod
    def _pybind11_conduit_v1_(*args, **kwargs): ...
    def __init__(self, arg0: int, arg1: int) -> None: ...
    def set_zero(self) -> None: ...
    @property
    def d_brdf(self) -> numpy.ndarray[numpy.float64[m, n]]: ...
    @property
    def d_emission(self) -> numpy.ndarray[numpy.float64[m, n]]: ...
    @property
    def interpolator(self) -> numpy.ndarray[numpy.float64[m, n]] | None: ...
    @interpolator.setter
    def interpolator(self, arg1: numpy.ndarray[numpy.float64[m, n]]) -> None: ...

class SurfaceStokes_1:
    brdf: BRDFStokes_1
    @staticmethod
    def _pybind11_conduit_v1_(*args, **kwargs): ...
    def get_derivative_mapping(self, name: str) -> SurfaceDerivativeMapping: ...
    @property
    def albedo(self) -> numpy.ndarray[numpy.float64[m, n]]: ...
    @property
    def brdf_args(self) -> numpy.ndarray[numpy.float64[m, n]]: ...
    @property
    def d_brdf_args(self) -> list[numpy.ndarray[numpy.float64[m, n]]]: ...
    @property
    def derivative_mappings(self) -> dict[str, SurfaceDerivativeMapping]: ...
    @property
    def emission(self) -> numpy.ndarray[numpy.float64[m, 1]]: ...
    @property
    def max_azimuthal_order(self) -> int: ...

class SurfaceStokes_3:
    brdf: BRDFStokes_3
    @staticmethod
    def _pybind11_conduit_v1_(*args, **kwargs): ...
    def get_derivative_mapping(self, name: str) -> SurfaceDerivativeMapping: ...
    @property
    def albedo(self) -> numpy.ndarray[numpy.float64[m, n]]: ...
    @property
    def brdf_args(self) -> numpy.ndarray[numpy.float64[m, n]]: ...
    @property
    def d_brdf_args(self) -> list[numpy.ndarray[numpy.float64[m, n]]]: ...
    @property
    def derivative_mappings(self) -> dict[str, SurfaceDerivativeMapping]: ...
    @property
    def emission(self) -> numpy.ndarray[numpy.float64[m, 1]]: ...
    @property
    def max_azimuthal_order(self) -> int: ...

class TangentAltitudeSolar(ViewingGeometryBase):
    @staticmethod
    def _pybind11_conduit_v1_(*args, **kwargs): ...
    def __init__(
        self,
        tangent_altitude_m: float,
        relative_azimuth: float,
        observer_altitude_m: float,
        cos_sza: float,
    ) -> None:
        """
        Defines a viewing ray from the observer altitude, and tangent point parameters. Note that all of
        these parameters assume straight line paths (i.e. no atmospheric refraction)

        Parameters
        ----------
        tangent_altitude_m: float
            Tangent altitude in [m]
        relative_azimuth: float
            Relative azimuth angle to the sun. An angle of 0 degrees corresponds to the forward scattering plane. [rad]
        observer_altitude_m: float
            Observer altitude relative to the earth [m]
        cos_sza: float
            Cosine of the solar zenith angle at the tangent point [unitless]
        """

class ThreadingModel:
    """
    Members:

      Wavelength

      Source
    """

    Source: typing.ClassVar[ThreadingModel]  # value = <ThreadingModel.Source: 1>
    Wavelength: typing.ClassVar[
        ThreadingModel
    ]  # value = <ThreadingModel.Wavelength: 0>
    __members__: typing.ClassVar[
        dict[str, ThreadingModel]
    ]  # value = {'Wavelength': <ThreadingModel.Wavelength: 0>, 'Source': <ThreadingModel.Source: 1>}
    @staticmethod
    def _pybind11_conduit_v1_(*args, **kwargs): ...
    def __eq__(self, other: typing.Any) -> bool: ...
    def __getstate__(self) -> int: ...
    def __hash__(self) -> int: ...
    def __index__(self) -> int: ...
    def __init__(self, value: int) -> None: ...
    def __int__(self) -> int: ...
    def __ne__(self, other: typing.Any) -> bool: ...
    def __repr__(self) -> str: ...
    def __setstate__(self, state: int) -> None: ...
    def __str__(self) -> str: ...
    @property
    def name(self) -> str: ...
    @property
    def value(self) -> int: ...

class ViewingGeometry:
    @staticmethod
    def _pybind11_conduit_v1_(*args, **kwargs): ...
    def __init__(self) -> None: ...
    def add_ray(self, arg0: ViewingGeometryBase) -> None: ...
    @property
    def observer_rays(self) -> list: ...

class ViewingGeometryBase:
    @staticmethod
    def _pybind11_conduit_v1_(*args, **kwargs): ...
    def __repr__(self) -> str: ...

class WignerD:
    @staticmethod
    def _pybind11_conduit_v1_(*args, **kwargs): ...
    def __init__(self, m: int, n: int) -> None:
        """
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
        """

    def d(
        self, theta: numpy.ndarray[numpy.float64], l: numpy.ndarray[numpy.int32]
    ) -> typing.Any:
        """
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
        """

def voigt_broaden(
    line_center: numpy.ndarray[numpy.float64[m, 1]],
    line_intensity: numpy.ndarray[numpy.float64[m, 1]],
    lower_energy: numpy.ndarray[numpy.float64[m, 1]],
    gamma_air: numpy.ndarray[numpy.float64[m, 1]],
    gamma_self: numpy.ndarray[numpy.float64[m, 1]],
    delta_air: numpy.ndarray[numpy.float64[m, 1]],
    n_air: numpy.ndarray[numpy.float64[m, 1]],
    iso_id: numpy.ndarray[numpy.int32[m, 1]],
    partitions: numpy.ndarray[numpy.float64[m, n], numpy.ndarray.flags.f_contiguous],
    molecular_mass: numpy.ndarray[numpy.float64[m, 1]],
    pressure: numpy.ndarray[numpy.float64[m, 1]],
    pself: numpy.ndarray[numpy.float64[m, 1]],
    temperature: numpy.ndarray[numpy.float64[m, 1]],
    wavenumber_grid: numpy.ndarray[numpy.float64[m, 1]],
    result: numpy.ndarray[
        numpy.float64[m, n],
        numpy.ndarray.flags.writeable,
        numpy.ndarray.flags.f_contiguous,
    ],
    line_contribution_width: float = 10.0,
    cull_factor: float = 0.0,
    num_threads: int = 1,
) -> None: ...

ConstantSpacing: GridSpacing  # value = <GridSpacing.ConstantSpacing: 0>
Disabled: InputValidationMode  # value = <InputValidationMode.Disabled: 2>
DiscreteOrdinates: (
    SingleScatterSource  # value = <SingleScatterSource.DiscreteOrdinates: 2>
)
Ellipsoidal: GeometryType  # value = <GeometryType.Ellipsoidal: 3>
Exact: SingleScatterSource  # value = <SingleScatterSource.Exact: 0>
LinearInterpolation: (
    InterpolationMethod  # value = <InterpolationMethod.LinearInterpolation: 1>
)
LinearSpacing: GridSpacing  # value = <GridSpacing.LinearSpacing: 1>
LowerInterpolation: (
    InterpolationMethod  # value = <InterpolationMethod.LowerInterpolation: 2>
)
NoSource: EmissionSource  # value = <EmissionSource.NoSource: 1>
Observer: StokesBasis  # value = <StokesBasis.Observer: 2>
OutOfBoundsExtend: OutOfBoundsPolicy  # value = <OutOfBoundsPolicy.OutOfBoundsExtend: 0>
OutOfBoundsSetZero: (
    OutOfBoundsPolicy  # value = <OutOfBoundsPolicy.OutOfBoundsSetZero: 1>
)
PlaneParallel: GeometryType  # value = <GeometryType.PlaneParallel: 0>
PseudoSpherical: GeometryType  # value = <GeometryType.PseudoSpherical: 1>
ShellInterpolation: (
    InterpolationMethod  # value = <InterpolationMethod.ShellInterpolation: 0>
)
Solar: StokesBasis  # value = <StokesBasis.Solar: 1>
Source: ThreadingModel  # value = <ThreadingModel.Source: 1>
Spherical: GeometryType  # value = <GeometryType.Spherical: 2>
Standard: InputValidationMode  # value = <InputValidationMode.Standard: 1>
Strict: InputValidationMode  # value = <InputValidationMode.Strict: 0>
SuccessiveOrders: (
    MultipleScatterSource  # value = <MultipleScatterSource.SuccessiveOrders: 1>
)
Table: SingleScatterSource  # value = <SingleScatterSource.Table: 1>
TwoStream: MultipleScatterSource  # value = <MultipleScatterSource.TwoStream: 2>
Wavelength: ThreadingModel  # value = <ThreadingModel.Wavelength: 0>
