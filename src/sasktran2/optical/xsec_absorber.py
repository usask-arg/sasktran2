from __future__ import annotations

from pathlib import Path

import numpy as np
import xarray as xr

from sasktran2._core_rust import PyXsecAbsorber
from sasktran2.database.web import StandardDatabase, WebDatabase
from sasktran2.optical import database
from sasktran2.optical.base import OpticalProperty

_CM2_TO_M2 = 1e-4
_IO_BREMEN_CACHE_DIR = Path("cross_sections/io/bremen")
_IO_BREMEN_HIGH_RESOLUTION_URL = (
    "https://www.iup.uni-bremen.de/gruppen/molspec/IO/"
    "RefSpec_0.07nmFWHM_IO_IUPBremen.txt"
)
_IO_BREMEN_TEMPERATURE_URL = (
    "https://www.iup.uni-bremen.de/gruppen/molspec/downloads/"
    "sigmaio1.3nmfwhm233k298k.txt"
)
_IO_BREMEN_TEMPERATURES_K = np.array(
    [233.0, 243.0, 263.0, 273.0, 283.0, 298.0],
    dtype=np.float64,
)
_IO_BREMEN_TEMPERATURE_SIGMA_COLUMNS = np.array([3, 5, 7, 9, 11, 13])
_IO_BREMEN_REFERENCE_TEMPERATURE_K = 298.0
_IO_BREMEN_REFERENCE_MIN_FRACTION = 0.01


class XsecAbsorber(OpticalProperty):
    """
    Cross section absorber loaded from HITRAN fixed-width format files.

    This absorber loads cross section data from .xsc files and performs
    interpolation in temperature (linear) and uses the nearest pressure level.
    Wavenumber interpolation uses zero padding outside the available range.

    Parameters
    ----------
    source : str or list of str or Path
        Can be:
        - Path to a single .xsc file
        - Path to a folder containing .xsc files
        - List of paths to .xsc files

    Examples
    --------
    Load from a folder:

    >>> absorber = XsecAbsorber("/path/to/hitran/xsec/CCl4")

    Load from a single file:

    >>> absorber = XsecAbsorber("/path/to/hitran/xsec/CCl4_v1.xsc")

    Load from multiple files:

    >>> absorber = XsecAbsorber(["/path/to/file1.xsc", "/path/to/file2.xsc"])
    """

    _internal: PyXsecAbsorber

    def __init__(self, source: str | Path | list[str] | list[Path]):
        if isinstance(source, str | Path):
            source_path = Path(source)
            if source_path.is_dir():
                # Load from folder
                self._internal = PyXsecAbsorber.from_folder(source_path.as_posix())
            else:
                # Load from single file
                self._internal = PyXsecAbsorber.from_file(source_path.as_posix())
        elif isinstance(source, list):
            # Load from list of files
            file_paths = [Path(f).as_posix() for f in source]
            self._internal = PyXsecAbsorber.from_files(file_paths)
        else:
            msg = "source must be a path to a file, folder, or list of file paths"
            raise TypeError(msg)

    @classmethod
    def from_lbl_database(cls, species_name: str) -> XsecAbsorber:
        """
        Load cross sections from the LBLRTM database for a given species.

        This method reads the FSCDXS index file to determine which cross section
        files to load for the requested species, then loads them.

        Parameters
        ----------
        species_name : str
            Name of the species. Case insensitive.

            Supported species:
            ACET, ACETICACI, BRO, C2CL2F4 (F114), C2CL3F3 (F113), C2CLF5 (F115),
            C2F6, C2HCl2F3CF2, CCL2FCH3, CCL4, CCLF3 (F13), CF3CH2CF3, CF3CH3,
            CF4 (F14), CFH2CF3, CH2F2, CH3CCLF2, CH3CHF2, CH3CN, CHCl2C2F5,
            CHCl2CF3, CHCL2F, CHClF2, CHClFCF3, CHF2CF3, CHF2CH2CF3, CHF3,
            CLONO2, F11, F12, FURAN, GLYCOLALD, HCHO, HNO3, HNO4, ISOP, N2O5,
            NF3, NO2, PAN, PROPENE, SF6, SO2

        Returns
        -------
        XsecAbsorber
            XsecAbsorber instance with the loaded cross section data

        Raises
        ------
        ValueError
            If the species name is not found in the LBLRTM database

        Examples
        --------
        >>> absorber = XsecAbsorber.from_lbl_database('CCL4')
        >>> absorber = XsecAbsorber.from_lbl_database('HNO3')
        >>> absorber = XsecAbsorber.from_lbl_database('f11')  # Case insensitive
        """
        db = StandardDatabase()
        fscdxs_path = db.path("cross_sections/lblrtm/FSCDXS")

        # Parse the FSCDXS file to find entries for this species
        file_list = []
        species_upper = species_name.upper()

        with fscdxs_path.open() as f:
            # Skip header line
            next(f)

            for line in f:
                # Skip empty lines
                if not line.strip():
                    continue

                # Parse the molecule name (first field, whitespace separated)
                parts = line.split()
                if not parts:
                    continue

                molecule = parts[0]

                # Check if this line matches our species
                if molecule.upper() == species_upper:
                    # Extract file names from the line
                    # Files can be concatenated like: xs/HNO3AT6xs/HNO3AT5xs/HNO3AT4...
                    # Split the entire line on 'xs/' to extract all filenames
                    line_parts = line.split("xs/")
                    for part in line_parts[1:]:  # Skip first part (before first 'xs/')
                        # Extract just the filename (everything before next space or end)
                        filename = part.split()[0] if part.strip() else ""
                        if filename:
                            file_path = db.path(f"cross_sections/lblrtm/xs/{filename}")
                            file_list.append(file_path)

        if not file_list:
            msg = f"No cross section files found for species '{species_name}' in LBLRTM database"
            raise ValueError(msg)

        # Create and return the XsecAbsorber from the file list
        return cls(file_list)

    def atmosphere_quantities(self, atmo, **kwargs):
        """
        Calculate optical quantities for a given atmosphere.

        Parameters
        ----------
        atmo : Atmosphere
            Atmosphere object containing temperature and pressure profiles
        **kwargs
            Must include wavenumbers_cminv as an array

        Returns
        -------
        OpticalQuantities
            Object containing cross sections and other optical properties
        """
        return self._internal.atmosphere_quantities(atmo, **kwargs)

    def optical_derivatives(self, atmo, **kwargs):
        """
        Calculate optical derivatives for a given atmosphere.

        Parameters
        ----------
        atmo : Atmosphere
            Atmosphere object containing temperature and pressure profiles
        **kwargs
            Must include wavenumbers_cminv as an array

        Returns
        -------
        dict
            Dictionary with "temperature_k" key containing temperature derivatives
        """
        return self._internal.optical_derivatives(atmo, **kwargs)

    def _into_rust_object(self):
        """Return the internal Rust object."""
        return self._internal


class O2SchumannRunge(XsecAbsorber):
    def __init__(self) -> None:
        """
        Schumann runge absorption taken from https://lweb.cfa.harvard.edu/amp/ampdata/cfamols.html.  Includes both the
        Continuum from 130 to 175 nm and the Schumann-Runge bands from 175 to 203 nm.  The data is reduced down to 0.1 nm resolution and is provided in vacuum wavelengths.

        Raises
        ------
        OSError
            If the specified file cannot be found
        """
        data_file = StandardDatabase().path("cross_sections/o2/O2SCHRUNG")

        if data_file.exists():
            super().__init__(data_file)
        else:
            msg = f"Could not find O2 Schumann-Runge database at {data_file}"
            raise OSError(msg)


class O2LymanAlpha(database.OpticalDatabaseGenericAbsorber):
    """
    O2 Lyman-alpha effective absorption cross section for attenuation at 121.567 nm.

    This is represented as a narrow pseudo-line optical database so it can be used by
    the existing cross-section interpolation machinery. The default cross section is
    derived from the Yankovsky O2 Lyman-alpha TOA photolysis rate,
    ``3.40e-9 s^-1``, and a nominal line-integrated solar Lyman-alpha flux,
    ``3.2e15 photons m^-2 s^-1``.
    """

    WAVELENGTH_NM = 121.567
    TOA_RATE_S = 3.40e-9
    TOA_FLUX_PHOTONS_M2_S = 3.2e15
    EFFECTIVE_CROSS_SECTION_M2 = TOA_RATE_S / TOA_FLUX_PHOTONS_M2_S

    def __init__(
        self,
        cross_section_m2: float = EFFECTIVE_CROSS_SECTION_M2,
        center_wavelength_nm: float = WAVELENGTH_NM,
        half_width_nm: float = 5.0e-4,
    ) -> None:
        if cross_section_m2 < 0:
            msg = "cross_section_m2 must be non-negative"
            raise ValueError(msg)
        if half_width_nm <= 0:
            msg = "half_width_nm must be positive"
            raise ValueError(msg)

        wavelength_nm = np.array(
            [
                center_wavelength_nm - half_width_nm,
                center_wavelength_nm,
                center_wavelength_nm + half_width_nm,
            ],
            dtype=np.float64,
        )
        xs = np.array([0.0, cross_section_m2, 0.0], dtype=np.float64)
        ds = xr.Dataset(
            data_vars={"xs": (["wavelength_nm"], xs)},
            coords={"wavelength_nm": wavelength_nm},
        )

        database.OpticalDatabase.__init__(self, db=ds)


class OClOGeisa(XsecAbsorber):
    def __init__(self) -> None:
        """
        OClO absorption cross sections from the GEISA database, reduced to 0.1 nm
        resolution from the 20 cm-1 source files.

        Raises
        ------
        OSError
            If the specified file cannot be found
        """
        data_file = StandardDatabase().path("cross_sections/oclo/OCLO20")

        if data_file.exists():
            super().__init__(data_file)
        else:
            msg = f"Could not find OClO GEISA database at {data_file}"
            raise OSError(msg)


class IOGeisa(XsecAbsorber):
    """
    IO (Iodine Oxide) cross sections from the GEISA database.

    These are fixed-width format cross sections (0.1 nm averaged)
    at 298 K temperature.

    Raises
    ------
    OSError
        If the IO file is not found in the database
    """

    def __init__(self):
        """Load IO cross sections from the standard database."""
        db = StandardDatabase()
        file_path = db.path("cross_sections/io/IO")
        super().__init__(file_path)


def _download_bremen_io_file(url: str) -> Path:
    db = WebDatabase(url, rel_path=_IO_BREMEN_CACHE_DIR)
    db.load()

    output_file = db.output_file()
    if not output_file.exists():
        msg = f"Could not download Bremen IO cross section data from {url}"
        raise OSError(msg)

    return output_file


def _read_bremen_io_table(path: Path, min_columns: int) -> np.ndarray:
    rows = []

    with path.open(encoding="latin-1") as f:
        for line in f:
            try:
                values = [float(value) for value in line.split()]
            except ValueError:
                continue

            if len(values) >= min_columns:
                rows.append(values[:min_columns])

    if not rows:
        msg = f"No numeric Bremen IO data rows found in {path}"
        raise ValueError(msg)

    return np.array(rows, dtype=np.float64)


def _bremen_io_high_resolution_cross_section(
    path: Path,
) -> tuple[np.ndarray, np.ndarray]:
    data = _read_bremen_io_table(path, min_columns=4)

    wavelength_nm = data[:, 1]
    xs_m2 = data[:, 3] * _CM2_TO_M2

    return wavelength_nm, xs_m2


def _bremen_io_temperature_scale_factors(
    high_resolution_wavelength_nm: np.ndarray,
    temperature_path: Path,
) -> np.ndarray:
    data = _read_bremen_io_table(temperature_path, min_columns=15)

    low_resolution_wavelength_nm = data[:, 1]
    low_resolution_sigma = data[:, _IO_BREMEN_TEMPERATURE_SIGMA_COLUMNS]
    reference_index = int(
        np.where(_IO_BREMEN_TEMPERATURES_K == _IO_BREMEN_REFERENCE_TEMPERATURE_K)[0][0]
    )
    reference_sigma = low_resolution_sigma[:, reference_index]

    scaling = np.ones(
        (_IO_BREMEN_TEMPERATURES_K.size, high_resolution_wavelength_nm.size),
        dtype=np.float64,
    )
    max_reference_sigma = np.nanmax(reference_sigma)
    if not np.isfinite(max_reference_sigma) or max_reference_sigma <= 0:
        return scaling

    reference_threshold = max_reference_sigma * _IO_BREMEN_REFERENCE_MIN_FRACTION

    for index, sigma in enumerate(low_resolution_sigma.T):
        if index == reference_index:
            continue

        valid = (
            np.isfinite(sigma)
            & np.isfinite(reference_sigma)
            & (sigma > 0)
            & (reference_sigma > reference_threshold)
        )

        if np.count_nonzero(valid) < 2:
            continue

        ratio = sigma[valid] / reference_sigma[valid]
        scaling[index] = np.interp(
            high_resolution_wavelength_nm,
            low_resolution_wavelength_nm[valid],
            ratio,
            left=ratio[0],
            right=ratio[-1],
        )

    return scaling


def _bremen_io_dataset(
    high_resolution_path: Path,
    temperature_path: Path | None = None,
) -> xr.Dataset:
    wavelength_nm, high_resolution_xs_m2 = _bremen_io_high_resolution_cross_section(
        high_resolution_path
    )

    if temperature_path is None:
        xs = high_resolution_xs_m2
        coords = {"wavelength_nm": wavelength_nm}
        dims = ["wavelength_nm"]
    else:
        scaling = _bremen_io_temperature_scale_factors(wavelength_nm, temperature_path)
        xs = scaling * high_resolution_xs_m2[np.newaxis, :]
        coords = {
            "temperature_k": _IO_BREMEN_TEMPERATURES_K,
            "wavelength_nm": wavelength_nm,
        }
        dims = ["temperature_k", "wavelength_nm"]

    return xr.Dataset(
        data_vars={"xs": (dims, xs)},
        coords=coords,
        attrs={
            "source": "IUP Bremen IO reference spectra",
            "reference_temperature_k": _IO_BREMEN_REFERENCE_TEMPERATURE_K,
            "high_resolution_fwhm_nm": 0.07,
        },
    )


class IOBremen(database.OpticalDatabaseGenericAbsorber):
    """
    IO (Iodine Oxide) absorption cross sections from the IUP Bremen database.

    By default this uses the 298 K IO spectrum measured at 0.07 nm FWHM. If
    ``include_temperature_dependence`` is true, the high-resolution spectrum is
    scaled at each temperature using ratios from the 1.3 nm FWHM measurements at
    233 K, 243 K, 263 K, 273 K, 283 K, and 298 K.

    Parameters
    ----------
    include_temperature_dependence : bool, optional
        Include the Bremen 1.3 nm temperature dependence as wavelength-dependent
        scaling factors applied to the 0.07 nm 298 K spectrum, by default False.
    """

    def __init__(self, include_temperature_dependence: bool = False) -> None:
        high_resolution_path = _download_bremen_io_file(_IO_BREMEN_HIGH_RESOLUTION_URL)

        if include_temperature_dependence:
            temperature_path = _download_bremen_io_file(_IO_BREMEN_TEMPERATURE_URL)
        else:
            temperature_path = None

        db = _bremen_io_dataset(high_resolution_path, temperature_path)
        database.OpticalDatabase.__init__(self, db=db)
