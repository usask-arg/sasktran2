// src/ffi.rs
#[allow(non_camel_case_types, non_snake_case, non_upper_case_globals)]
mod bindings {
    include!("bindings.rs");
}

pub use bindings::*; // optional, or only re-export some items
