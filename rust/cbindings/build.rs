use std::env;

fn main() {
    let out_dir = env::var("OUT_DIR").unwrap();
    let install_prefix = std::path::Path::new(&out_dir).join("install");

    let dst = cmake::Config::new("../../")
        .define("BUILD_SHARED_LIBS", "OFF")
        .define("CMAKE_INSTALL_PREFIX", &install_prefix)
        .define("USE_OMP", "OFF")
        .define("SKTRAN_BLAS_VENDOR", "Apple")
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

    let bindings = bindgen::Builder::default()
        .header("../../include/c_api/sasktran2.h")
        .generate_inline_functions(true)
        .generate()
        .expect("Unable to generate bindings");

    bindings
        .write_to_file("src/bindings.rs")
        .expect("Couldn't write bindings!");

    // Also link C++ stdlib if needed
    // println!("cargo:rustc-link-lib=dylib=stdc++"); // Or libc++ on macOS
}
