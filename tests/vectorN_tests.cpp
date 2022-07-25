#include "../catch/catch.hpp"

#ifdef MML_USE_SINGLE_HEADER
#include "MML.h"
#else
#include "basic_types/VectorN.h"
#endif

TEST_CASE("Test_VectorN", "[simple]") {
	Vector2 a{1.0, 1.0};

	REQUIRE(a[0] ==  1.0);
	REQUIRE(a[1] ==  1.0);
}