use std::env;

fn main() {
    let out_dir = env::var("OUT_DIR").unwrap();
    let install_prefix = std::path::Path::new(&out_dir).join("install");

    let use_omp = env::var("USE_OMP").unwrap_or_else(|_| "OFF".to_string());
    let sktran_blas_vendor = env::var("SKTRAN_BLAS_VENDOR").unwrap_or_else(|_| "Apple".to_string());
    let do_stream_templates = env::var("DO_STREAM_TEMPLATES").unwrap_or_else(|_| "OFF".to_string());

    let dst = cmake::Config::new("../../")
        .define("BUILD_SHARED_LIBS", "OFF")
        .define("CMAKE_INSTALL_PREFIX", &install_prefix)
        .define("DO_STREAM_TEMPLATES", do_stream_templates)
        .define("USE_OMP", use_omp)
        .define("SKTRAN_BLAS_VENDOR", sktran_blas_vendor)
        .build();

    println!("cargo:root={}", install_prefix.display());

    println!("cargo:rustc-link-search=native={}/build/lib", dst.display());
    println!(
        "cargo:rustc-link-search=native={}/build/c_api",
        dst.display()
    );
    println!("cargo:rustc-link-lib=static=csasktran2");
    println!("cargo:rustc-link-lib=static=sasktran2");

    println!("cargo:rerun-if-changed=../../include");
    println!("cargo:rerun-if-changed=../../lib");
    println!("cargo:rerun-if-changed=../../CMakeLists.txt");

    let bindings = bindgen::Builder::default()
        .header("../../include/c_api/sasktran2.h")
        .generate_inline_functions(true)
        .generate()
        .expect("Unable to generate bindings");

    bindings
        .write_to_file("src/bindings.rs")
        .expect("Couldn't write bindings!");

    // Also link C++ stdlib if needed
    println!("cargo:rustc-link-lib=stdc++");
    // println!("cargo:rustc-link-lib=dylib=stdc++"); // Or libc++ on macOS
}
