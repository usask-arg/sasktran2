use crate::constituent::atmo_storage::AtmosphereStorage;
use crate::optical::optical_quantities::PyOpticalQuantities;
use pyo3::exceptions::PyValueError;
use pyo3::types::PyDict;
use pyo3::{IntoPyObjectExt, prelude::*};
use sasktran2_rs::optical::traits::{OpticalProperty, OpticalPropertyExt};
use sasktran2_rs::optical::types::xsec_absorber::{self, XsecDatabase};

use super::xsec_dbase::{HasDb, PyDictWrapper};

#[pyclass(unsendable)]
/// A cross section absorber database loaded from HITRAN fixed-width format files.
///
/// This database stores cross sections as a function of wavenumber at various
/// temperature and pressure conditions. It supports interpolation in temperature
/// (linear between two nearest points) and uses the nearest pressure level.
///
/// Examples
/// --------
/// Create from a single file:
///
/// >>> absorber = XsecAbsorber.from_file("path/to/file.xsc")
///
/// Create from multiple files:
///
/// >>> absorber = XsecAbsorber.from_files(["file1.xsc", "file2.xsc"])
///
/// Create from a folder (reads all .xsc files):
///
/// >>> absorber = XsecAbsorber.from_folder("path/to/folder")
pub struct PyXsecAbsorber {
    pub xsec_absorber: xsec_absorber::XsecDatabase,
}

#[pymethods]
impl PyXsecAbsorber {
    #[staticmethod]
    /// Create a cross section database from a single file.
    ///
    /// Parameters
    /// ----------
    /// filepath : str
    ///     Path to a HITRAN fixed-width format cross section file (.xsc)
    ///
    /// Returns
    /// -------
    /// XsecAbsorber
    ///     The loaded cross section database
    pub fn from_file(filepath: &str) -> PyResult<Self> {
        let path = std::path::PathBuf::from(filepath);

        let db = xsec_absorber::read_fwf_xsec(path)
            .ok_or_else(|| PyValueError::new_err("Failed to read cross section file"))?;

        Ok(PyXsecAbsorber { xsec_absorber: db })
    }

    #[staticmethod]
    /// Create a cross section database from multiple files.
    ///
    /// Parameters
    /// ----------
    /// filepaths : list of str
    ///     List of paths to HITRAN fixed-width format cross section files (.xsc)
    ///
    /// Returns
    /// -------
    /// XsecAbsorber
    ///     The combined cross section database
    pub fn from_files(filepaths: Vec<String>) -> PyResult<Self> {
        let mut combined_db = XsecDatabase(vec![]);

        for filepath in filepaths {
            let path = std::path::PathBuf::from(filepath);
            let db = xsec_absorber::read_fwf_xsec(path.clone()).ok_or_else(|| {
                PyValueError::new_err(format!("Failed to read file: {}", path.display()))
            })?;

            combined_db += db;
        }

        if combined_db.0.is_empty() {
            return Err(PyValueError::new_err("No cross section data was loaded"));
        }

        Ok(PyXsecAbsorber {
            xsec_absorber: combined_db,
        })
    }

    #[staticmethod]
    /// Create a cross section database from all .xsc files in a folder.
    ///
    /// Parameters
    /// ----------
    /// folder : str
    ///     Path to a folder containing HITRAN fixed-width format cross section files
    ///
    /// Returns
    /// -------
    /// XsecAbsorber
    ///     The combined cross section database from all .xsc files in the folder
    pub fn from_folder(folder: &str) -> PyResult<Self> {
        let path = std::path::PathBuf::from(folder);

        let db = xsec_absorber::read_fwf_folder(path)
            .ok_or_else(|| PyValueError::new_err("Failed to read cross section folder"))?;

        if db.0.is_empty() {
            return Err(PyValueError::new_err(
                "No cross section data was loaded from folder",
            ));
        }

        Ok(PyXsecAbsorber { xsec_absorber: db })
    }

    #[pyo3(signature = (atmo, **kwargs))]
    /// Calculate optical quantities for a given atmosphere.
    ///
    /// Parameters
    /// ----------
    /// atmo : Atmosphere
    ///     Atmosphere object containing temperature_k and pressure_pa profiles
    /// **kwargs
    ///     Additional keyword arguments. Must include:
    ///     - wavenumbers_cminv : ndarray
    ///         Wavenumbers in cm^-1 at which to calculate cross sections
    ///
    /// Returns
    /// -------
    /// OpticalQuantities
    ///     Object containing cross_section and other optical properties
    fn atmosphere_quantities<'py>(
        &self,
        atmo: Bound<'py, PyAny>,
        kwargs: Option<&Bound<'py, PyDict>>,
    ) -> PyResult<Bound<'py, PyAny>> {
        let rust_atmo = AtmosphereStorage::new(&atmo)?;
        let aux_inputs = PyDictWrapper(kwargs);

        let oq = self
            .xsec_absorber
            .optical_quantities(&rust_atmo.inputs, &aux_inputs)
            .map_err(|e| PyValueError::new_err(format!("Failed to get optical quantities: {e}")))?;

        PyOpticalQuantities::new(oq).into_bound_py_any(atmo.py())
    }

    #[pyo3(signature = (atmo, **kwargs))]
    /// Calculate optical derivatives for a given atmosphere.
    ///
    /// Parameters
    /// ----------
    /// atmo : Atmosphere
    ///     Atmosphere object containing temperature_k and pressure_pa profiles
    /// **kwargs
    ///     Additional keyword arguments. Must include:
    ///     - wavenumbers_cminv : ndarray
    ///         Wavenumbers in cm^-1 at which to calculate cross sections
    ///
    /// Returns
    /// -------
    /// dict
    ///     Dictionary with key "temperature_k" containing OpticalQuantities
    ///     with derivatives of cross section with respect to temperature
    fn optical_derivatives<'py>(
        &self,
        atmo: Bound<'py, PyAny>,
        kwargs: Option<&Bound<'py, PyDict>>,
    ) -> PyResult<Bound<'py, PyDict>> {
        let rust_atmo = AtmosphereStorage::new(&atmo)?;
        let aux_inputs = PyDictWrapper(kwargs);

        let oq = self
            .xsec_absorber
            .optical_derivatives(&rust_atmo.inputs, &aux_inputs)
            .map_err(|e| {
                PyValueError::new_err(format!("Failed to get optical derivatives: {e}"))
            })?;

        let py = atmo.py();
        let py_dict = PyDict::new(py);

        for (key, oq) in oq {
            let py_oq = PyOpticalQuantities::new(oq).into_bound_py_any(py)?;
            py_dict.set_item(key.as_str(), py_oq)?;
        }

        Ok(py_dict)
    }
}

impl HasDb for PyXsecAbsorber {
    fn db(&self) -> &dyn OpticalProperty {
        &self.xsec_absorber
    }
}
