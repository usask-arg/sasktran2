pub mod ffi;

#[cfg(force_link_scipy_blas)]
#[link(name = "scipy_openblas64_")]
unsafe extern "C" {
    fn scipy_dgbsv_64_();
}

#[cfg(force_link_scipy_blas)]
#[used]
static _FORCE_LINK_BLAS: unsafe extern "C" fn() = scipy_dgbsv_64_;
