#include <sasktran2/test_helper.h>

#include <sasktran2.h>

#ifdef SASKTRAN_RUST_SUPPORT

TEST_CASE("Test Rust integration", "[sasktran2][rust]") {

    double result = sasktran2::rust::test();

    auto test = sasktran2::rust::new (2, 3);

    REQUIRE(result == 42.0);
}

#endif
