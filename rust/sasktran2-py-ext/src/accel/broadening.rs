use numpy::*;
use pyo3::prelude::*;
use sasktran2_rs::optical::broaden;

use crate::prelude::IntoPyResult;

#[pyfunction]
#[allow(clippy::too_many_arguments)]
pub fn voigt_broaden_uniform<'py>(
    line_center: PyReadonlyArray1<'py, f64>,
    line_intensity: PyReadonlyArray1<'py, f64>,
    lower_energy: PyReadonlyArray1<'py, f64>,
    gamma_air: PyReadonlyArray1<'py, f64>,
    gamma_self: PyReadonlyArray1<'py, f64>,
    delta_air: PyReadonlyArray1<'py, f64>,
    n_air: PyReadonlyArray1<'py, f64>,
    iso_id: PyReadonlyArray1<'py, i32>,
    partitions: PyReadonlyArray2<'py, f64>,
    mol_mass: PyReadonlyArray1<'py, f64>,
    pressure: PyReadonlyArray1<'py, f64>,
    pself: PyReadonlyArray1<'py, f64>,
    temperature: PyReadonlyArray1<'py, f64>,
    first_wavenumber: f64,
    wavenumber_spacing: f64,
    mut result: PyReadwriteArray2<'py, f64>, // [geometry, wavenumber]
    line_contribution_width: f64,
    cull_factor: f64,
    num_threads: usize,
    subtract_pedastal: bool,
) -> PyResult<()> {
    broaden::voigt_broaden_uniform(
        line_center.as_array().view(),
        line_intensity.as_array().view(),
        lower_energy.as_array().view(),
        gamma_air.as_array().view(),
        gamma_self.as_array().view(),
        delta_air.as_array().view(),
        n_air.as_array().view(),
        iso_id.as_array().view(),
        partitions.as_array().view(),
        mol_mass.as_array().view(),
        pressure.as_array().view(),
        pself.as_array().view(),
        temperature.as_array().view(),
        first_wavenumber,
        wavenumber_spacing,
        result.as_array_mut().view_mut(), // [geometry, wavenumber]
        line_contribution_width,
        cull_factor,
        num_threads,
        subtract_pedastal,
    )
    .into_pyresult()?;

    Ok(())
}

#[pyfunction]
#[allow(clippy::too_many_arguments)]
pub fn voigt_broaden<'py>(
    line_center: PyReadonlyArray1<'py, f64>,
    line_intensity: PyReadonlyArray1<'py, f64>,
    lower_energy: PyReadonlyArray1<'py, f64>,
    gamma_air: PyReadonlyArray1<'py, f64>,
    gamma_self: PyReadonlyArray1<'py, f64>,
    delta_air: PyReadonlyArray1<'py, f64>,
    n_air: PyReadonlyArray1<'py, f64>,
    iso_id: PyReadonlyArray1<'py, i32>,
    partitions: PyReadonlyArray2<'py, f64>,
    mol_mass: PyReadonlyArray1<'py, f64>,
    pressure: PyReadonlyArray1<'py, f64>,
    pself: PyReadonlyArray1<'py, f64>,
    temperature: PyReadonlyArray1<'py, f64>,
    wavenumber_grid: PyReadonlyArray1<'py, f64>,
    mut result: PyReadwriteArray2<'py, f64>, // [geometry, wavenumber]
    line_contribution_width: f64,
    cull_factor: f64,
    num_threads: usize,
    subtract_pedastal: bool,
) -> PyResult<()> {
    broaden::voigt_broaden(
        line_center.as_array().view(),
        line_intensity.as_array().view(),
        lower_energy.as_array().view(),
        gamma_air.as_array().view(),
        gamma_self.as_array().view(),
        delta_air.as_array().view(),
        n_air.as_array().view(),
        iso_id.as_array().view(),
        partitions.as_array().view(),
        mol_mass.as_array().view(),
        pressure.as_array().view(),
        pself.as_array().view(),
        temperature.as_array().view(),
        wavenumber_grid.as_array().view(),
        result.as_array_mut().view_mut(), // [geometry, wavenumber]
        line_contribution_width,
        cull_factor,
        num_threads,
        subtract_pedastal,
    )
    .into_pyresult()?;

    Ok(())
}

#[pyfunction]
#[allow(clippy::too_many_arguments)]
pub fn voigt_broaden_with_line_coupling<'py>(
    line_center: PyReadonlyArray1<'py, f64>,
    line_intensity: PyReadonlyArray1<'py, f64>,
    lower_energy: PyReadonlyArray1<'py, f64>,
    gamma_air: PyReadonlyArray1<'py, f64>,
    gamma_self: PyReadonlyArray1<'py, f64>,
    delta_air: PyReadonlyArray1<'py, f64>,
    n_air: PyReadonlyArray1<'py, f64>,
    iso_id: PyReadonlyArray1<'py, i32>,
    partitions: PyReadonlyArray2<'py, f64>,
    y_coupling: PyReadonlyArray2<'py, f64>,
    g_coupling: PyReadonlyArray2<'py, f64>,
    mol_mass: PyReadonlyArray1<'py, f64>,
    pressure: PyReadonlyArray1<'py, f64>,
    pself: PyReadonlyArray1<'py, f64>,
    temperature: PyReadonlyArray1<'py, f64>,
    wavenumber_grid: PyReadonlyArray1<'py, f64>,
    mut result: PyReadwriteArray2<'py, f64>, // [geometry, wavenumber]
    line_contribution_width: f64,
    cull_factor: f64,
    num_threads: usize,
    subtract_pedastal: bool,
) -> PyResult<()> {
    broaden::voigt_broaden_with_line_coupling(
        line_center.as_array().view(),
        line_intensity.as_array().view(),
        lower_energy.as_array().view(),
        gamma_air.as_array().view(),
        gamma_self.as_array().view(),
        delta_air.as_array().view(),
        n_air.as_array().view(),
        iso_id.as_array().view(),
        partitions.as_array().view(),
        y_coupling.as_array().view(),
        g_coupling.as_array().view(),
        mol_mass.as_array().view(),
        pressure.as_array().view(),
        pself.as_array().view(),
        temperature.as_array().view(),
        wavenumber_grid.as_array().view(),
        result.as_array_mut().view_mut(), // [geometry, wavenumber]
        line_contribution_width,
        cull_factor,
        num_threads,
        subtract_pedastal,
    )
    .into_pyresult()?;

    Ok(())
}
