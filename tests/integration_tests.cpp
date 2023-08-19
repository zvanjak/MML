#include "../catch/catch.hpp"

#ifdef MML_USE_SINGLE_HEADER
#include "MML.h"
#else
#include "algorithms/Integration.h"
#endif

TEST_CASE("Test_Integration", "[simple]") {
	REQUIRE(1.0 ==  1.0);
    // row - integration method
    // column - req prec
}