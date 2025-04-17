#![allow(non_snake_case)]

use crate::mie::Mie;

use num::complex::Complex64;
use std::os::raw::c_double;

#[unsafe(no_mangle)]
pub extern "C" fn mie_new() -> *mut Mie {
    let mie = Mie::new();
    Box::into_raw(Box::new(mie))
}

#[unsafe(no_mangle)]
pub extern "C" fn mie_with_cos_angles(mie_ptr: *mut Mie, angles: *const c_double, len: usize) {
    if mie_ptr.is_null() || angles.is_null() {
        return;
    }

    let angles_slice = unsafe { std::slice::from_raw_parts(angles, len) };
    let mie = unsafe { &mut *mie_ptr };
    *mie = std::mem::replace(mie, Mie::new()).with_cos_angles(angles_slice.to_vec());
}

#[unsafe(no_mangle)]
pub extern "C" fn mie_calculate(
    mie_ptr: *mut Mie,
    size_param: c_double,
    refractive_real: c_double,
    refractive_imag: c_double,
    Qext: *mut c_double,
    Qsca: *mut c_double,
) {
    if mie_ptr.is_null() || Qext.is_null() || Qsca.is_null() {
        return;
    }

    let mie = unsafe { &mut *mie_ptr };
    let refractive_index = Complex64::new(refractive_real, refractive_imag);
    let (ext, sca) = mie.calculate(size_param, refractive_index);

    unsafe {
        *Qext = ext;
        *Qsca = sca;
    }
}

#[unsafe(no_mangle)]
pub extern "C" fn mie_free(mie_ptr: *mut Mie) {
    if !mie_ptr.is_null() {
        unsafe {
            drop(Box::from_raw(mie_ptr));
        }
    }
}
