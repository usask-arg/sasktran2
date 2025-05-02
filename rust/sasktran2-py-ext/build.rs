fn main() {
    if let Ok(openblas_path) = std::env::var("DEP_SASKTRAN2_OPENBLAS_PATH") {
        println!("cargo:rustc-link-arg={}", openblas_path);
    }
}
