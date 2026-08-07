#include <sasktran2/testing/lto.h>
#include <chrono>
#include <iostream>
#ifdef SKTRAN_RUST_SUPPORT
#include <sasktran2-core/src/lib.rs.h>
#endif

int cpp_echo(int val) { return val; }

int test_fun() {
    int sum = 0;
    for (int i = 0; i < 1000000; ++i) {
        sum += cpp_echo(i);
    }
    return sum;
}

#ifdef SKTRAN_RUST_SUPPORT
int test_rust_fun() {
    int sum = 0;
    for (int i = 0; i < 1000000; ++i) {
        sum += sasktran2::rust::testing::rust_echo(i);
    }
    return sum;
}
#endif

void sasktran2::testing::run_lto_tests() {
    auto t1 = std::chrono::high_resolution_clock::now();
    auto sum = test_fun();
    auto t2 = std::chrono::high_resolution_clock::now();

    auto duration =
        std::chrono::duration_cast<std::chrono::nanoseconds>(t2 - t1).count();

    std::cout << "C++ sum: " << sum << " Time (ns): " << duration << std::endl;

#ifdef SKTRAN_RUST_SUPPORT
    auto t1_r = std::chrono::high_resolution_clock::now();
    auto sum_r = test_rust_fun();
    auto t2_r = std::chrono::high_resolution_clock::now();

    auto duration_r =
        std::chrono::duration_cast<std::chrono::nanoseconds>(t2_r - t1_r)
            .count();

    std::cout << "Rust sum: " << sum_r << " Time (ns): " << duration_r
              << std::endl;
#endif
}
