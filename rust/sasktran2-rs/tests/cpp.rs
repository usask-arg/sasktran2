use sasktran2_sys::ffi::{run_catch2_tests};

#[test]
fn run_all_cpp_tests() {
    // Catch2 expects argc>=1 and argv[0] set to a program name
    let prog = std::ffi::CString::new("sasktran2_tests").unwrap();
    let mut argv = vec![prog.as_ptr() as *mut i8, std::ptr::null_mut()];
    let result = unsafe { run_catch2_tests(1, argv.as_mut_ptr()) };
    assert_eq!(result, 0, "C++ tests failed");
}