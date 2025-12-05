
#ifdef SKTRAN_CATCH2_VERSION3
#define CATCH_CONFIG_RUNNER
#else
#define CATCH_CONFIG_MAIN
#endif

#include <sasktran2/test_helper.h>
#include <catch2/catch_session.hpp>

TEST_CASE("Main") {}


extern "C" int run_catch2_tests(int argc, char** argv) {
    int result = Catch::Session().run( argc, argv );
    return result;
}
