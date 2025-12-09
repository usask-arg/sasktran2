// build.rs

fn main() {
    let _ = cxx_build::bridge("src/math/wigner.rs");
}
