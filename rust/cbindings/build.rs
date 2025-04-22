fn main() {
    let bindings = bindgen::Builder::default()
        .header("/Users/djz828/dev/sasktran2/include/c_api/sasktran2.h")
        .generate_inline_functions(true)
        .generate()
        .expect("Unable to generate bindings");

    bindings
        .write_to_file("src/bindings.rs")
        .expect("Couldn't write bindings!");

    println!("cargo:rustc-link-search=native=/Users/djz828/dev/sasktran2/build/c_api");
    println!("cargo:rustc-link-search=native=/Users/djz828/dev/sasktran2/build/lib");
    println!("cargo:rustc-link-lib=static=csasktran2");
    println!("cargo:rustc-link-lib=static=sasktran2");

    // Also link C++ stdlib if needed
    // println!("cargo:rustc-link-lib=dylib=stdc++"); // Or libc++ on macOS
}
