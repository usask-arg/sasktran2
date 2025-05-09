pub mod prelude;
pub mod traits;
pub mod types;
pub mod raytracers;

#[cxx::bridge]
mod ffi {
    extern "Rust" {
        fn trace_ray(x: f64, y: f64, z: f64) -> f64;
    }
}

pub fn trace_ray(x: f64, y: f64, z: f64) -> f64 {
    0.0
}