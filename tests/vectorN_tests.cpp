#include "../catch/catch.hpp"

#ifdef MML_USE_SINGLE_HEADER
#include "MML.h"
#else
#include "basic_types/VectorN.h"
#endif

TEST_CASE("VectorN_default_ctor_init_to_zero", "[simple]") {
    MML::VectorN<3> a;

	REQUIRE(0.0 ==  a[0]);
	REQUIRE(0.0 ==  a[1]);
	REQUIRE(0.0 ==  a[2]);
}

TEST_CASE("VectorN_initializer_list_ctor", "[simple]") {
    MML::VectorN<3> a{1.0, 2.0, 3.0};

	REQUIRE(1.0 ==  a[0]);
	REQUIRE(2.0 ==  a[1]);
	REQUIRE(3.0 ==  a[2]);
}

TEST_CASE("Test_VectorN", "[simple]") {
	Vector2 a{1.0, 1.0};

	REQUIRE(a[0] ==  1.0);
	REQUIRE(a[1] ==  1.0);
}