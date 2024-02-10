#include "../catch/catch.hpp"

#ifdef MML_USE_SINGLE_HEADER
#include "MML.h"
#else
#include "base/VectorN.h"
#endif

using namespace MML;

TEST_CASE("VectorN_default_ctor_init_to_zero", "[simple]") {
    VectorN<Real, 3> a;

	REQUIRE(0.0 ==  a[0]);
	REQUIRE(0.0 ==  a[1]);
	REQUIRE(0.0 ==  a[2]);
}

TEST_CASE("VectorN_initializer_list_ctor", "[simple]") {
    VectorN<Real, 3> a{1.0, 2.0, 3.0};

	REQUIRE(1.0 ==  a[0]);
	REQUIRE(2.0 ==  a[1]);
	REQUIRE(3.0 ==  a[2]);
}

TEST_CASE("Test_VectorN", "[simple]") {
	Vector2Dbl a{1.0, 1.0};

	REQUIRE(a[0] ==  1.0);
	REQUIRE(a[1] ==  1.0);
}

TEST_CASE("VectorN_init_to_value", "[simple]") {
    VectorN<Real, 3> a(1.5);

	REQUIRE(3 == a.size());

	REQUIRE(1.5 ==  a[0]);
	REQUIRE(1.5 ==  a[1]);
	REQUIRE(1.5 ==  a[2]);
}

TEST_CASE("VectorN_init_from_std::vector", "[simple]") {
	std::vector<Real> x{1.0, 2.0, 3.0};
    VectorN<Real,3> a(x);

	REQUIRE(3 == a.size());

	REQUIRE(1.0 ==  a[0]);
	REQUIRE(2.0 ==  a[1]);
	REQUIRE(3.0 ==  a[2]);
}

TEST_CASE("VectorN_init_from_double_array", "[simple]") {
	Real x[] = {1.0, 2.0, 3.0};
    VectorN<Real, 3> a(x);

	REQUIRE(3 == a.size());

	REQUIRE(1.0 ==  a[0]);
	REQUIRE(2.0 ==  a[1]);
	REQUIRE(3.0 ==  a[2]);
}

TEST_CASE("Test_VectorN_IsEqual", "[simple]") {
    Vector2Dbl a({1.0, 2.0});
    Vector2Dbl b({1.0, 2.0000001});

	REQUIRE(true == a.IsEqual(b, 1e-7));
    REQUIRE(false == a.IsEqual(b, 1e-8));
}

TEST_CASE("Test_VectorN_Op+-", "[simple]") {
    Vector2Dbl a({1.0, 2.0});
    Vector2Dbl b({1.0, 2.0});

    auto c = a + b;
    auto d = a - b;

	REQUIRE(2.0 ==  c[0]);
	REQUIRE(4.0 ==  c[1]);

	REQUIRE(0.0 ==  d[0]);
	REQUIRE(0.0 ==  d[1]);
}

TEST_CASE("Test_VectorN_mul_double", "[simple]") {
    Vector2Dbl a({1.0, 100.0});

	auto b = a * 2.0;
	auto c = 2.0 * a;

	REQUIRE(2.0 ==  b[0]);
	REQUIRE(2.0 ==  c[0]);

	REQUIRE(200.0 ==  b[1]);
	REQUIRE(200.0 ==  c[1]);	
}

TEST_CASE("Test_VectorN_div_double", "[simple]") {
    Vector2Dbl a({4.0, 400.0});

	auto b = a / 2.0;

	REQUIRE(2.0 ==  b[0]);
	REQUIRE(200.0 ==  b[1]);
}

TEST_CASE("Test_VectorN_ScalarProductCartesian", "[simple]") {
    Vector2Dbl a({1.0, 2.0});
    Vector2Dbl b({1.0, 2.0});

	REQUIRE(5.0 == a.ScalarProductCartesian(b));
}

TEST_CASE("Test_VectorN_NormL2", "[simple]") {
    Vector2Dbl a({2.0, 2.0});

	REQUIRE(sqrt(8) == a.NormL2());
}

TEST_CASE("Test_VectorN_to_string", "[simple]") {
    Vector2Dbl a({2.0, 2.0});

	REQUIRE("[    2,     2]" == a.to_string(5,3));

    Vector4Dbl b({123.0, 1.0, 10.0, -8.0});

	REQUIRE("[    123,       1,      10,      -8]" == b.to_string(7,3));

    Vector3Dbl c({123.123, 1.9876543, 10.0});

	REQUIRE("[    123.12,     1.9877,         10]" == c.to_string(10,5));
	REQUIRE("[        123.123,       1.9876543,              10]" == c.to_string(15,9));
}
