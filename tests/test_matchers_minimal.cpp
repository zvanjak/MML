// Minimal test to isolate the compilation issue
#include <catch2/catch_test_macros.hpp>

TEST_CASE("Minimal test", "[minimal]") {
    REQUIRE(1 == 1);
}
