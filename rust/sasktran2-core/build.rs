fn main() {
    let _ = cxx_build::bridge("src/math/wigner.rs");
    let _ = cxx_build::bridge("src/math/gaussquad.rs");
    //.file("src/math/gaussquad.rs");
}
