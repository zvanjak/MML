#include "../catch/catch.hpp"

#ifdef MML_USE_SINGLE_HEADER
#include "MMLBasicTypes.h"
#else
#include "basic_types/Vector.h"
#endif

TEST_CASE("Test_Vector_mul_double", "[simple]") {
    MML::Vector a({1.0, 1.0});

	REQUIRE(1.0 ==  a[0]);
	REQUIRE(1.0 ==  a[1]);
}

// zbrajanje, oduzimanje, op. sa skalaraom

// kako reagira na index violation?
