use numpy::*;
use pyo3::prelude::*;
use sasktran2_rs::optical::storage::OpticalQuantities;
use sasktran2_rs::optical::util::legendre::LegendreAccess;

#[pyclass]
pub struct PyOpticalQuantities {
    oq: OpticalQuantities,
}

impl PyOpticalQuantities {
    pub fn new(oq: OpticalQuantities) -> Self {
        Self { oq }
    }
}

#[pymethods]
impl PyOpticalQuantities {
    #[new]
    #[pyo3(signature = (num_geometry, num_wavelengths, fortran_ordering = false))]
    fn py_new(num_geometry: usize, num_wavelengths: usize, fortran_ordering: bool) -> Self {
        let oq = OpticalQuantities::new(num_geometry, num_wavelengths, fortran_ordering);
        Self { oq }
    }

    #[getter]
    fn get_cross_section<'py>(this: Bound<'py, Self>) -> Bound<'py, PyArray2<f64>> {
        let array = &this.borrow().oq.cross_section;

        unsafe { PyArray2::borrow_from_array(array, this.into_any()) }
    }

    #[getter]
    // This is just for interfacing
    fn get_extinction<'py>(this: Bound<'py, Self>) -> Bound<'py, PyArray2<f64>> {
        let array = &this.borrow().oq.cross_section;

        unsafe { PyArray2::borrow_from_array(array, this.into_any()) }
    }

    #[getter]
    fn get_ssa<'py>(this: Bound<'py, Self>) -> Bound<'py, PyArray2<f64>> {
        let array = &this.borrow().oq.ssa;

        unsafe { PyArray2::borrow_from_array(array, this.into_any()) }
    }

    // Native optical derivatives use this same storage type. These aliases
    // implement the NativeGridDerivative interface used by Python wrappers.
    #[getter]
    fn get_d_extinction<'py>(this: Bound<'py, Self>) -> Bound<'py, PyArray2<f64>> {
        Self::get_extinction(this)
    }

    #[getter]
    fn get_d_ssa<'py>(this: Bound<'py, Self>) -> Bound<'py, PyArray2<f64>> {
        Self::get_ssa(this)
    }

    #[getter]
    fn get_d_leg_coeff<'py>(this: Bound<'py, Self>) -> Option<Bound<'py, PyArray3<f64>>> {
        let binding = this.borrow();
        let array = binding.oq.legendre.as_ref()?;
        // NativeGridDerivative uses (coefficient, location, wavelength), while
        // the native optical storage keeps the coefficient dimension last.
        let view = array.view().permuted_axes([2, 0, 1]);
        unsafe { Some(PyArray3::borrow_from_array(&view, this.into_any())) }
    }

    #[getter]
    fn get_a1<'py>(this: Bound<'py, Self>) -> Option<Bound<'py, PyArray3<f64>>> {
        let binding = this.borrow();
        let view = binding.oq.legendre.as_ref();

        if let Some(view) = view {
            let a1 = view.view();
            let a1 = LegendreAccess::<Ix3, 4>::new(a1).a1;

            unsafe { Some(PyArray3::borrow_from_array(&a1, this.into_any())) }
        } else {
            None
        }
    }

    #[getter]
    fn get_leg_coeff<'py>(this: Bound<'py, Self>) -> Option<Bound<'py, PyArray3<f64>>> {
        let array = &this.borrow().oq.legendre;

        if let Some(array) = array {
            unsafe { Some(PyArray3::borrow_from_array(array, this.into_any())) }
        } else {
            None
        }
    }
}
