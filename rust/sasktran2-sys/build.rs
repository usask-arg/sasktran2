use std::env;
use std::fs;
use std::path::{Path, PathBuf};

fn cpp_src_path() -> (PathBuf, bool) {
    // if the vendor folder exists, we are in a sdist mode, use that
    // otherwise we are in a dev mode, use the cpp folder
    let vendor_path = PathBuf::from("vendor/");
    if vendor_path.exists() {
        return (vendor_path, true);
    }

    (PathBuf::from("../../cpp"), false)
}

fn main() {
    let (cpp_src, vendored) = cpp_src_path();

    let vendored = match vendored {
        true => "ON",
        false => "OFF",
    };

    let out_dir = env::var("OUT_DIR").unwrap();
    let install_prefix = std::path::Path::new(&out_dir).join("install");

    let use_omp = env::var("USE_OMP").unwrap_or_else(|_| "OFF".to_string());

    let default_blas = if cfg!(target_os = "macos") {
        "Apple"
    } else {
        "OpenBLAS"
    };

    let sktran_blas_vendor =
        env::var("SKTRAN_BLAS_VENDOR").unwrap_or_else(|_| default_blas.to_string());
    let do_stream_templates = env::var("DO_STREAM_TEMPLATES").unwrap_or_else(|_| "OFF".to_string());

    let mut binding = cmake::Config::new(&cpp_src);

    let target = std::env::var("TARGET").unwrap();

    binding
        .define("BUILD_SHARED_LIBS", "OFF")
        .define("CMAKE_INSTALL_PREFIX", &install_prefix)
        .define("DO_STREAM_TEMPLATES", do_stream_templates)
        .define("USE_OMP", use_omp)
        .define("Rust_CARGO_TARGET", target)
        .define("VENDORED", vendored)
        .define("SKTRAN_BLAS_VENDOR", sktran_blas_vendor);

    if cfg!(target_os = "windows") {
        binding
            .define("CMAKE_CONFIGURATION_TYPES", "Release")
            .profile("Release");
        // Use release-mode MSVC runtime
        binding.define("CMAKE_MSVC_RUNTIME_LIBRARY", "MultiThreadedDLL");
    }

    let dst = binding.build();

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
                if let Some(name) = lib_name.to_str()
                    && let Some(extension) = path.extension()
                    && let Some(extension) = extension.to_str()
                {
                    if extension == "framework" {
                        println!("cargo:rustc-link-lib=framework={name}");
                    } else if name.starts_with("lib") {
                        // if we are on windows, we keep the 'lib' prefix
                        // otherwise we remove it
                        if cfg!(target_os = "windows") {
                            println!("cargo:rustc-link-lib=dylib={name}");
                        } else {
                            println!(
                                "cargo:rustc-link-lib=dylib={}",
                                name.strip_prefix("lib").unwrap_or(name)
                            );
                        }
                    } else {
                        println!("cargo:rustc-link-lib=dylib={name}");
                    }
                }
            }
        }
    }

    println!("cargo:root={}", install_prefix.display());

    println!(
        "cargo:rustc-link-search=native={}/install/lib",
        dst.display()
    );
    println!(
        "cargo:rustc-link-search=native={}/install/c_api",
        dst.display()
    );
    println!("cargo:rustc-link-lib=static=csasktran2");
    println!("cargo:rustc-link-lib=static=sasktran2");

    println!("cargo:rerun-if-changed={}/include", cpp_src.display());
    println!("cargo:rerun-if-changed={}/lib", cpp_src.display());
    println!("cargo:rerun-if-changed={}/c_api", cpp_src.display());
    println!(
        "cargo:rerun-if-changed={}/CMakeLists.txt",
        cpp_src.display()
    );
    #[cfg(feature = "build-bindings")]
    {
        let bindings = bindgen::Builder::default()
            .header("../../cpp/include/c_api/sasktran2.h")
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
