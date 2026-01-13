pub mod math;

#[cxx::bridge(namespace = "sasktran2::rust::testing")]
pub mod ffi {
    extern "Rust" {
        fn rust_echo(input: i32) -> i32;
    }
}

#[inline(always)]
pub fn rust_echo(input: i32) -> i32 {
    input
}
