use sasktran2_sys::cxxffi;

#[test]
fn test_integration() {
    let result = cxxffi::new_wigner_d_calculator(1, 0);
}
