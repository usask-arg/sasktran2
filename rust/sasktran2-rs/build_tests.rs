// This is a hack to link the test library for integration tests
// Integration tests in tests/ are separate binaries that don't automatically
// get the link directives from sasktran2-sys

fn main() {
    // Get the path from sasktran2-sys's build
    let out_dir = std::env::var("OUT_DIR").unwrap();
    let sys_build_dir = std::path::Path::new(&out_dir)
        .parent().unwrap()
        .parent().unwrap()
        .parent().unwrap()
        .join("build/sasktran2-sys-*/out/install/lib");
    
    // Find the actual build directory (glob pattern)
    if let Ok(entries) = std::fs::read_dir(sys_build_dir.parent().unwrap().parent().unwrap()) {
        for entry in entries.filter_map(Result::ok) {
            let path = entry.path();
            if path.file_name().unwrap().to_string_lossy().starts_with("sasktran2-sys-") {
                let lib_dir = path.join("out/install/lib");
                if lib_dir.exists() {
                    println!("cargo:rustc-link-search=native={}", lib_dir.display());
                    println!("cargo:rustc-link-lib=static=sasktran2_tests");
                    break;
                }
            }
        }
    }
}
