#include "../catch/catch.hpp"

#ifdef MML_USE_SINGLE_HEADER
#include "MMLBasicTypes.h"
#else
#include "algorithms/Integration.h"
#endif

TEST_CASE("Test_Integration", "[simple]") {
	REQUIRE(1.0 ==  1.0);
}