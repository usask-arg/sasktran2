use sasktran2_sys::ffi;

/// Returns true if the openmp support is enabled in the sasktran2 library
/// Controlled with a compile time flag when building the library
pub fn openmp_support_enabled() -> bool {
    let openmp = unsafe { ffi::sk_openmp_support_enabled() };
    openmp != 0
}
