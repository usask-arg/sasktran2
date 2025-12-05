// src/ffi.rs
#[allow(non_camel_case_types, non_snake_case, non_upper_case_globals)]
mod bindings {
    include!("bindings.rs");

    unsafe extern "C" {
        pub unsafe fn run_catch2_tests(argc: i32, argv: *mut *mut i8) -> i32;
    }
}

pub use bindings::*; // optional, or only re-export some items
