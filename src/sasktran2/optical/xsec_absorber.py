from __future__ import annotations

from pathlib import Path

from sasktran2._core_rust import PyXsecAbsorber
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
        if isinstance(source, (str, Path)):
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
