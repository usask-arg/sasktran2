use std::env;
use std::fs;
use std::path::{Path, PathBuf};


fn main() {
    // print all environment variables
    for (key, value) in env::vars() {
        println!("{}: {}", key, value);
    }

    let out_dir = env::var("OUT_DIR").unwrap();
    let install_prefix = std::path::Path::new(&out_dir).join("install");

    let use_omp = env::var("USE_OMP").unwrap_or_else(|_| "OFF".to_string());
    let sktran_blas_vendor = env::var("SKTRAN_BLAS_VENDOR").unwrap_or_else(|_| "OpenBLAS".to_string());
    let do_stream_templates = env::var("DO_STREAM_TEMPLATES").unwrap_or_else(|_| "OFF".to_string());

    let dst = cmake::Config::new("../../")
        .define("BUILD_SHARED_LIBS", "OFF")
        .define("CMAKE_INSTALL_PREFIX", &install_prefix)
        .define("DO_STREAM_TEMPLATES", do_stream_templates)
        .define("USE_OMP", use_omp)
        .define("SKTRAN_BLAS_VENDOR", sktran_blas_vendor);

    if cfg!(target_os = "windows") {
        let build_type = if cfg!(debug_assertions) {
            "Debug"
        } else {
            dst.define("CMAKE_MSVC_RUNTIME_LIBRARY", "MultiThreadedDebugDLL"); // /MDd for Debu
            "Release"
        };
        dst.define("CMAKE_BUILD_TYPE", build_type);
    }

    let dst = dst.build();

    let lib_file_path = Path::new(&out_dir).join("build").join("libs_to_link.txt");

    // Read the file content
    let lib_contents =
        fs::read_to_string(lib_file_path).expect("Failed to read library paths file");

    // Iterate over each path in the file (assuming semicolon-separated)
    for lib_path in lib_contents.trim().split(';') {
        if !lib_path.is_empty() {
            let path = PathBuf::from(lib_path);
            // Extract the directory where the library is located
            if let Some(parent) = path.parent() {
                println!("cargo:rustc-link-search=native={}", parent.display());
            }
            // Extract the library name
            if let Some(lib_name) = path.file_stem() {
                // Assumes the library name starts with 'lib' as in 'libopenblas'
                if let Some(name) = lib_name.to_str() {

                    if let Some(extension) = path.extension() {
                        if let Some(extension) = extension.to_str() {
                            if extension == "framework" {
                                println!("cargo:rustc-link-lib=framework={}", name);
                            } else {
                                if name.starts_with("lib") {
                                    // if we are on windows, we keep the 'lib' prefix
                                    // otherwise we remove it
                                    if cfg!(target_os = "windows") {
                                        println!("cargo:rustc-link-lib=dylib={}", name);
                                    } else {
                                        println!("cargo:rustc-link-lib=dylib={}", &name[3..]);
                                    }
                                } else {
                                    println!("cargo:rustc-link-lib=dylib={}", name);
                                }
                            }
                        }
                    }
                }
            }
        }
    }

    println!("cargo:root={}", install_prefix.display());

    println!("cargo:rustc-link-search=native={}/install/lib", dst.display());
    println!(
        "cargo:rustc-link-search=native={}/install/c_api",
        dst.display()
    );
    println!("cargo:rustc-link-lib=static=csasktran2");
    println!("cargo:rustc-link-lib=static=sasktran2");

    println!("cargo:rerun-if-changed=../../include");
    println!("cargo:rerun-if-changed=../../lib");
    println!("cargo:rerun-if-changed=../../CMakeLists.txt");
    #[cfg(feature = "build-bindings")]
    {
        let bindings = bindgen::Builder::default()
            .header("../../include/c_api/sasktran2.h")
            .generate_inline_functions(true)
            .generate()
            .expect("Unable to generate bindings");

        bindings
            .write_to_file("src/bindings.rs")
            .expect("Couldn't write bindings!");
    }

    // Also link C++ stdlib if needed
    if cfg!(target_os = "linux") {
        println!("cargo:rustc-link-lib=stdc++");
    }

    println!("Build artifacts at: {}", dst.display());
    for entry in std::fs::read_dir(dst.join("build/lib")).unwrap() {
        println!("Found lib: {:?}", entry.unwrap().path());
    }

    // println!("cargo:rustc-link-lib=dylib=stdc++"); // Or libc++ on macOS
}
