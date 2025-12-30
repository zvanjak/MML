// Test Catch2 with just matchers included
#include <catch2/catch_test_macros.hpp>
#include "TestMatchers.h"

using namespace MML::Testing::Matchers;

TEST_CASE("Simple matcher test", "[minimal]") {
    REQUIRE_THAT(REAL(1.0), RealEquals(REAL(1.0)));
}
