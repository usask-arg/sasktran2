from __future__ import annotations

from pathlib import Path

from sasktran2._core_rust import PyXsecAbsorber
from sasktran2.database.web import StandardDatabase
from sasktran2.optical.base import OpticalProperty


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
