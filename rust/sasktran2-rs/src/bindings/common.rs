use sasktran2_sys::ffi;

pub fn openmp_support_enabled() -> bool {
    let openmp = unsafe { ffi::sk_openmp_support_enabled() };
    openmp != 0
}
