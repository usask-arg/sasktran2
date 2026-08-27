use ndarray::Array2;
use numpy::{IntoPyArray, PyArray2, PyReadonlyArray1};
use pyo3::exceptions::PyValueError;
use pyo3::prelude::*;
use rayon::prelude::*;
use sasktran2_rs::continuum::mt_ckd::{
    AtmosphericState, LinearizedSpectrum, MtCkd, OUTPUT_POINTS, calculate, calculate_absorption,
    calculate_absorption_linearized, calculate_linearized,
};
use std::sync::OnceLock;

fn model(data_file: &str) -> PyResult<&'static MtCkd> {
    static MODEL: OnceLock<MtCkd> = OnceLock::new();

    if MODEL.get().is_none() {
        let loaded = MtCkd::from_file(data_file).map_err(|error| {
            PyValueError::new_err(format!("failed to load MT_CKD coefficient data: {error}"))
        })?;
        let _ = MODEL.set(loaded);
    }

    Ok(MODEL.get().expect("MT_CKD model must be initialized"))
}

fn states(
    pressure_pa: PyReadonlyArray1<'_, f64>,
    temperature_k: PyReadonlyArray1<'_, f64>,
    vmr_h2o: PyReadonlyArray1<'_, f64>,
    vmr_co2: PyReadonlyArray1<'_, f64>,
    vmr_o3: PyReadonlyArray1<'_, f64>,
) -> PyResult<Vec<AtmosphericState>> {
    let pressure_pa = pressure_pa.as_array().to_vec();
    let temperature_k = temperature_k.as_array().to_vec();
    let vmr_h2o = vmr_h2o.as_array().to_vec();
    let vmr_co2 = vmr_co2.as_array().to_vec();
    let vmr_o3 = vmr_o3.as_array().to_vec();
    let length = pressure_pa.len();
    if [
        temperature_k.len(),
        vmr_h2o.len(),
        vmr_co2.len(),
        vmr_o3.len(),
    ]
    .into_iter()
    .any(|candidate| candidate != length)
    {
        return Err(PyValueError::new_err(
            "pressure_pa, temperature_k, vmr_h2o, vmr_co2, and vmr_o3 must have the same length",
        ));
    }

    Ok((0..length)
        .map(|index| AtmosphericState {
            pressure_pa: pressure_pa[index],
            temperature_k: temperature_k[index],
            h2o_vmr: vmr_h2o[index],
            co2_vmr: vmr_co2[index],
            o3_vmr: vmr_o3[index],
        })
        .collect())
}

fn array_from_values(values: Vec<f64>, row_count: usize) -> Array2<f64> {
    Array2::from_shape_vec((row_count, OUTPUT_POINTS), values).unwrap()
}

fn array_from_rows(rows: impl IntoIterator<Item = Vec<f64>>, row_count: usize) -> Array2<f64> {
    array_from_values(rows.into_iter().flatten().collect(), row_count)
}

/// MT_CKD 4.3 extinction coefficient on the fixed 1,991-point grid.
#[pyfunction(name = "_mt_ckd")]
#[allow(clippy::too_many_arguments)]
pub fn mt_ckd_py<'py>(
    py: Python<'py>,
    pressure_pa: PyReadonlyArray1<'py, f64>,
    temperature_k: PyReadonlyArray1<'py, f64>,
    vmr_h2o: PyReadonlyArray1<'py, f64>,
    vmr_co2: PyReadonlyArray1<'py, f64>,
    vmr_o3: PyReadonlyArray1<'py, f64>,
    data_file: &str,
    include_rayleigh: bool,
) -> PyResult<Bound<'py, PyArray2<f64>>> {
    let model = model(data_file)?;
    let states = states(pressure_pa, temperature_k, vmr_h2o, vmr_co2, vmr_o3)?;
    let row_count = states.len();
    let spectra: Vec<_> = py.detach(|| {
        states
            .into_par_iter()
            .map(|state| {
                if include_rayleigh {
                    calculate(model, state)
                } else {
                    calculate_absorption(model, state)
                }
            })
            .collect()
    });
    Ok(array_from_rows(spectra, row_count).into_pyarray(py))
}

type LinearizedArrays<'py> = (
    Bound<'py, PyArray2<f64>>,
    Bound<'py, PyArray2<f64>>,
    Bound<'py, PyArray2<f64>>,
    Bound<'py, PyArray2<f64>>,
    Bound<'py, PyArray2<f64>>,
    Bound<'py, PyArray2<f64>>,
);

/// Extinction coefficient followed by analytic Jacobians for P, T, H2O, CO2, and O3.
#[pyfunction(name = "_mt_ckd_linearized")]
#[allow(clippy::too_many_arguments)]
pub fn mt_ckd_linearized_py<'py>(
    py: Python<'py>,
    pressure_pa: PyReadonlyArray1<'py, f64>,
    temperature_k: PyReadonlyArray1<'py, f64>,
    vmr_h2o: PyReadonlyArray1<'py, f64>,
    vmr_co2: PyReadonlyArray1<'py, f64>,
    vmr_o3: PyReadonlyArray1<'py, f64>,
    data_file: &str,
    include_rayleigh: bool,
) -> PyResult<LinearizedArrays<'py>> {
    let model = model(data_file)?;
    let states = states(pressure_pa, temperature_k, vmr_h2o, vmr_co2, vmr_o3)?;
    let row_count = states.len();
    let spectra: Vec<LinearizedSpectrum> = py.detach(|| {
        states
            .into_par_iter()
            .map(|state| {
                if include_rayleigh {
                    calculate_linearized(model, state)
                } else {
                    calculate_absorption_linearized(model, state)
                }
            })
            .collect()
    });

    let mut fields: [Vec<f64>; 6] =
        std::array::from_fn(|_| Vec::with_capacity(row_count * OUTPUT_POINTS));
    for spectrum in spectra {
        fields[0].extend(spectrum.extinction);
        fields[1].extend(spectrum.d_pressure_pa);
        fields[2].extend(spectrum.d_temperature_k);
        fields[3].extend(spectrum.d_h2o_vmr);
        fields[4].extend(spectrum.d_co2_vmr);
        fields[5].extend(spectrum.d_o3_vmr);
    }
    let [extinction, d_pressure, d_temperature, d_h2o, d_co2, d_o3] = fields;
    Ok((
        array_from_values(extinction, row_count).into_pyarray(py),
        array_from_values(d_pressure, row_count).into_pyarray(py),
        array_from_values(d_temperature, row_count).into_pyarray(py),
        array_from_values(d_h2o, row_count).into_pyarray(py),
        array_from_values(d_co2, row_count).into_pyarray(py),
        array_from_values(d_o3, row_count).into_pyarray(py),
    ))
}
