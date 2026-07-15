from __future__ import annotations

import numpy as np

import sasktran2 as sk
from sasktran2.constants import K_BOLTZMANN
from sasktran2.database.o3 import O3BirkWagnerDatabase, O3SerdyuchenkoDatabase
from sasktran2.database.web import StandardDatabase

from . import database, refraction  # noqa: F401
from .aerline import AERLineAbsorber  # noqa: F401
from .baum import BaumIceCrystal  # noqa: F401
from .henyey import HenyeyGreenstein  # noqa: F401
from .hitran import (
    HITRANAbsorber,  # noqa: F401
)
from .mie import Mie  # noqa: F401
from .xsec_absorber import (  # noqa: F401
    IOBremen,
    IOGeisa,
    O2LymanAlpha,
    O2SchumannRunge,
    OClOGeisa,
    XsecAbsorber,
)


class O3DBM(database.OpticalDatabaseGenericAbsorber):
    def __init__(self) -> None:
        """
        Tabulated high resolution cross-sections of O3 measured by Daumont, Brion and Malicet in the early 1990's [1].
        The wavelength range slightly varies with temperature but covers the entire UV to NIR region, from 194.50 nm to
        830.00 nm at 0.01 to 0.02 nm resolution. The cross-section data were collected at 0.01-0.02 nm resolution and
        each wavelength/cross-section table varies in size from 22,052 to 63,501 entries. The data consists of 5 tables
        of wavelength versus cross-section for 5 temperatures.

        Notes
        -----
        Temerature Range
            Measurements are provided at 5 temperatures covering typical stratospheric and
            tropospheric conditions:

                | 218 K
                | 228 K
                | 243 K
                | 273 K
                | 295 K

        Wavelength Range
            The wavelength range of each temperature table is slightly different and is given below.
            Note that most of the temperature variation occurs in the Huggins band
            between 315 and 360 nm:

                | 218K -> 194.50nm to 650.01nm
                | 228K -> 194.50nm to 520.01nm
                | 243K -> 194.50nm to 519.01nm
                | 273K -> 299.50nm to 520.01nm
                | 295K -> 195.00nm to 830.00nm

            We looked into temperature interpolation and while DBM suggest that a quadratic interpolation scheme [3] they do
            not indicate an explicit technique. We tested several quadratic fitting routines and found that a truncated linear
            fit in temperature was visually more appealing than any of the quadratic fits and had none of the undesirable
            artifacts (excessive curvature etc.) that naturally arise with quadratic curve fitting. Consequently this object
            uses a truncated linear fit in temperature.

        Data Source
            These data are an exact replication of the data files:

                | O3_CRS_BDM_218K.dat
                | O3_CRS_BDM_228K.dat
                | O3_CRS_BDM_243K.dat
                | O3_CRS_BDM_273K.dat
                | O3_CRS_BDM_295K.dat

            Data is from the IGACO site, http://igaco-o3.fmi.fi/ACSO/cross_sections.html. The files were copied on
            July 16-July 25 2012.

        References
        ----------
        .. [1] Daumont, D., et al. "Ozone UV spectroscopy I: Absorption cross-sections at room temperature."
                Journal of Atmospheric Chemistry 15.2 (1992): 145-155.
        .. [2] Brion, J., et al. "High-resolution laboratory absorption cross section of O3. Temperature effect."
                Chemical physics letters 213.5 - 6(1993): 610-612.
        .. [3] Malicet, J., et al. "Ozone UV spectroscopy. II. Absorption cross-sections and temperature dependence."
                Journal of Atmospheric Chemistry 21.3 (1995): 263-273.
        .. [4] Brion, J., et al. "Absorption spectra measurements for the ozone molecule in the 350-830 nm region."
                Journal of Atmospheric Chemistry 30.2 (1998): 291-299.

        Raises
        ------
        OSError
            If the file could not be found
        """
        dbm_file = StandardDatabase().path("cross_sections/o3/dbm.nc")

        if dbm_file.exists():
            super().__init__(dbm_file)
        else:
            msg = "Could not find DBM file"
            raise OSError(msg)


class O3BirkWagner(database.OpticalDatabaseGenericAbsorber):
    def __init__(self) -> None:
        """
        Tabulated high resolution UV cross-sections of O3 measured by Birk and Wagner.

        This database is useful for UV ozone absorption calculations in the
        Hartley/Huggins region. The raw Zenodo archive contains measured
        absorption cross sections and uncertainties for six temperatures:

            | 193 K
            | 213 K
            | 233 K
            | 253 K
            | 273 K
            | 293 K

        The wavelength grid is common to all six tables and spans approximately
        243.9 nm to 346.0 nm. The wavelength grid is used directly as a vacuum
        wavelength grid; an overlap comparison with the independently
        vacuum-labelled Serdyuchenko dataset shows the Huggins-band structure is
        aligned as-is and shifted if an air-to-vacuum conversion is applied.
        Cross sections are provided by the source files in cm^2/molecule and
        converted internally to m^2/molecule.

        Data Source
            The data are from the ESA SEOM-IAS ozone UV absorption cross-section
            database created by M. Birk and G. Wagner at DLR for improved
            atmospheric spectroscopy in support of Sentinel-5P/TROPOMI. The
            source archive also contains uncertainty columns and polynomial
            coefficients; this class uses the measured tabulated cross sections.

        The source files are downloaded from Zenodo record 1485588 when the local
        NetCDF cache is missing. Downloading requires the optional ``zenodo-get``
        package, available through the ``sasktran2[zenodo]`` extra.

        References
        ----------
        Birk, M. and Wagner, G. "ESA SEOM-IAS - Measurement and ACS database
        O3 UV region." Zenodo, 2018. doi:10.5281/zenodo.1485588.

        Raises
        ------
        ImportError
            If the data must be downloaded and zenodo-get is not installed
        OSError
            If the data could not be downloaded or converted
        """
        super().__init__(O3BirkWagnerDatabase().path())


class O3Serdyuchenko(database.OpticalDatabaseGenericAbsorber):
    def __init__(self) -> None:
        """
        Tabulated UV/VIS/NIR cross-sections of O3 measured by Serdyuchenko et al.

        This database covers a broader spectral range than O3DBM or
        O3BirkWagner, extending from 213.33 nm to 1100.00 nm on a 0.01 nm
        wavelength grid. Cross sections are available at 11 temperatures:

            | 193 K
            | 203 K
            | 213 K
            | 223 K
            | 233 K
            | 243 K
            | 253 K
            | 263 K
            | 273 K
            | 283 K
            | 293 K

        The measurements combine a Bruker HR 120 Fourier transform spectrometer
        and an ESA 400 Echelle spectrometer. Cross sections are provided by the
        source file in cm^2/molecule and converted internally to m^2/molecule.

        Data Source
            The source file is the Serdyuchenko-Gorshelev UV/VIS/NIR ozone
            absorption cross-section dataset. The original file stores
            temperatures in descending order from 293 K to 193 K; SASKTRAN2
            sorts them into ascending order before interpolation.

        The source file is downloaded from Zenodo record 5793207 when the local
        NetCDF cache is missing. Downloading requires the optional ``zenodo-get``
        package, available through the ``sasktran2[zenodo]`` extra.

        References
        ----------
        Gorshelev, V., Serdyuchenko, A., Weber, M., Chehade, W., and Burrows,
        J. P. "High spectral resolution ozone absorption cross-sections -
        Part 1: Measurements, data analysis and comparison with previous
        measurements around 293 K." Atmospheric Measurement Techniques 7,
        609-624, 2014. doi:10.5194/amt-7-609-2014.

        Serdyuchenko, A., Gorshelev, V., Weber, M., Chehade, W., and Burrows,
        J. P. "High spectral resolution ozone absorption cross-sections -
        Part 2: Temperature dependence." Atmospheric Measurement Techniques 7,
        625-636, 2014. doi:10.5194/amt-7-625-2014.

        Serdyuchenko, A., Gorshelev, V., Weber, M., and Burrows, J. P.
        "Serdyuchenko-Gorshelev UV/VIS/NIR ozone absorption cross-section."
        Zenodo, 2021. doi:10.5281/zenodo.5793207.

        Raises
        ------
        ImportError
            If the data must be downloaded and zenodo-get is not installed
        OSError
            If the data could not be downloaded or converted
        """
        super().__init__(O3SerdyuchenkoDatabase().path())


class NO2Vandaele(database.OpticalDatabaseGenericAbsorber):
    def __init__(self) -> None:
        """
        Calculates the absorption cross section of NO2 molecules from 230 nm to 1000 nm at 220 K to 294 K following [1]

        .. [1] Vandaele, Ann Carine, et al. "Measurements of the NO2 absorption cross-section from 42 000 cm-1 to
            10 000 cm-1 (238-1000 nm) at 220 K and 294 K." Journal of Quantitative Spectroscopy and Radiative Transfer
            59.3-5 (1998): 171-184.

        Raises
        ------
        OSError
            If the Vandaele file cannot be found
        """
        v_file = StandardDatabase().path("cross_sections/no2/vandaele.nc")

        if v_file.exists():
            super().__init__(v_file)
        else:
            msg = "Could not find Vandaele file"
            raise OSError(msg)


class HITRANUV(database.OpticalDatabaseGenericAbsorber):
    def __init__(self, name: str, version: str = "2022") -> None:
        """
        HITRAN UV cross-sections for a given constituent.  The constituent must be one of the following:
        ["O3", "NO2", "BrO", "SO2"]

        These are the files that you obtain from the HITRAN website when selecting UV cross sections.
        Note that usually these are not very good, and you should use a specific database for each constituent instead
        such as O3DBM.  The HITRAN UV cross sections are only provided for convenience.

        Specific Issues
            O3
                We have noticed the Ozone cross sections are unreasonably large in the region < 310 nm, and does not match the source data.
                We do not recommend using this.
            NO2
                Even though the HITRAN database specifically says all cross sections are given in vacuum wavenumbers, this one
                seems to be specified in air wavenumbers.  We have converted it to vacuum wavenumbers manually. After doing so
                it matches the Vandaele source files relatively well.

        Parameters
        ----------
        name : str
            Constituent
        version : str, optional
            HITRAN version number, currently only "2022" is supported, by default "2022"

        Raises
        ------
        OSError
            If the databases are not installed
        """
        data_file = StandardDatabase().path(f"cross_sections/{name}/hitran{version}.nc")

        if data_file.exists():
            super().__init__(data_file)
        else:
            msg = f"Could not find HITRAN UV database for {name} at {data_file}"
            raise OSError(msg)


class HITRANTabulated(database.OpticalDatabaseGenericAbsorber):
    def __init__(self, name: str, res="01nm") -> None:
        """
        Loads in a database tabulated from HITRAN line entries as a function of pressure/temperature that have been
        reduced to a given resolution.

        This requires the extended databases to be downloaded, e.g. by running sk.appconfig.download_extended_databases()

        Currently supported species and resolutions are:

            | H2O : 01nm
            | O2 : 01nm

        Parameters
        ----------
        name : str
            Species name
        res : str, optional
            Resolution of the database to load in, by default "01nm"

        Raises
        ------
        OSError
            If the specified file cannot be found
        """
        data_file = sk.appconfig.database_root().joinpath(
            f"cross_sections/{name.lower()}/hitran_{res}_res.nc"
        )

        if data_file.exists():
            super().__init__(data_file)
        else:
            msg = f"Could not find HITRAN database for {name} at {data_file}"
            raise OSError(msg)


class HITRANCollision(database.OpticalDatabaseGenericAbsorber):
    def __init__(self, name: str) -> None:
        """
        Loads collision induced absorption (CIA) cross sections compiled from data found at https://hitran.org/cia/.

        This requires the extended databases to be downloaded, e.g. by running sk.appconfig.download_extended_databases()

        Currently supported species are:

            | O2O2

        Parameters
        ----------
        name : str
            Species name
        """

        data_file = StandardDatabase().path(
            f"cross_sections/{name.lower()}/hitran_cia.nc"
        )

        if data_file.exists():
            super().__init__(data_file)
        else:
            msg = f"Could not find HITRAN CIA database for {name} at {data_file}"
            raise OSError(msg)


def pressure_temperature_to_numberdensity(
    pressure_pa: np.array, temperature_k: np.array, include_derivatives=False
) -> np.array:
    """
    Converts pressure and temperature to number density using the ideal gas law

    Parameters
    ----------
    pressure_pa : np.array
        Pressure in [Pa]
    temperature_k : np.array
        Temperature in [K]
    include_derivatives: bool, optional
        If set to true, the derivative of the number density with respect to pressure and temperature will be returned
        alongside the numberdensity.  The return signature is then np.array, np.array, np.array

    Returns
    -------
    np.array, if include_derivatives is False
        Number density in [molecules / m^3]

    np.array, np.array, np.array if include_derivatives is True
        Number density, dN/dP, dN/dT in [molecules / m^3], [molecules / m^3 / Pa], [molecules / m^3 / K]
    """
    # Ideal gas law is PV = N k T
    # Or (N/V) = P / (k T)

    N = pressure_pa / (K_BOLTZMANN * temperature_k)

    if not include_derivatives:
        return N

    dN_dP = 1 / (K_BOLTZMANN * temperature_k)
    dN_dT = -pressure_pa / (K_BOLTZMANN * temperature_k**2)

    return N, dN_dP, dN_dT


def air_wavelength_to_vacuum_wavelength(wavelength_nm: np.array) -> np.array:
    """
    Converts wavelength specified in Air at STP to wavelengths at vacuum

    Parameters
    ----------
    wavelength_nm : np.array
        Wavelength in air at STP

    Returns
    -------
    np.array
        Vacuum wavelengths
    """
    s = 1e4 / (wavelength_nm * 10)  # Convert to angstroms

    refac_index = (
        1
        + 0.00008336624212083
        + 0.02408926869968 / (130.1065924522 - s**2)
        + 0.0001599740894897 / (38.92568793293 - s**2)
    )

    return wavelength_nm * refac_index


def vacuum_wavelength_to_air_wavelength(wavelength_nm: np.array) -> np.array:
    """
    Converts wavelength specified in vacuum to wavelengths at air at STP

    Parameters
    ----------
    wavelength_nm : np.array
        Wavelengths in vacuum

    Returns
    -------
    np.array
        Wavelengths in air at STP
    """
    s = 1e4 / (wavelength_nm * 10)  # Convert to angstroms

    n = 1 + 0.0000834254 + 0.02406147 / (130 - s**2) + 0.00015998 / (38.9 - s**2)

    return wavelength_nm / n
