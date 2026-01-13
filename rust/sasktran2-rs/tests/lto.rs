use sasktran2_sys::ffi::run_lto_tests;

#[test]
fn test_lto() {
    unsafe {
        run_lto_tests();
    }
}
