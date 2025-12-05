pub mod ffi;

#[cxx::bridge(namespace="sasktran2::math")]
pub mod cxxffi {

    extern "Rust" {
        fn test() -> f64;
    }

    unsafe extern "C++" {
        include!("sasktran2/math/wigner.h");

        type WignerDCalculator;

        fn new_wigner_d_calculator(m: i32, n: i32) -> UniquePtr<WignerDCalculator>;
    }
}


fn test() -> f64 {
    42.0
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_wigner_d_basic() {
        let calc = cxxffi::new_wigner_d_calculator(1, 1);

    }
}