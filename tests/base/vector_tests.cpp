#include "../catch/catch.hpp"

#ifdef MML_USE_SINGLE_HEADER
#include "MML.h"
#else
#include "base/Vector.h"
#endif

using namespace MML;

namespace MML::VectorTests
{
TEST_CASE("Vector_default_ctor_init_to_zero", "[simple]") {
    Vector<Real> a(3);

    REQUIRE(3 == a.size());

	REQUIRE(0.0 ==  a[0]);
	REQUIRE(0.0 ==  a[1]);
	REQUIRE(0.0 ==  a[2]);

    Vector<Real> b(0);
    REQUIRE(0 == b.size());
}

TEST_CASE("Vector_init_to_value", "[simple]") {
    Vector<Real> a(3, 1.5);

	REQUIRE(3 == a.size());

	REQUIRE(1.5 ==  a[0]);
	REQUIRE(1.5 ==  a[1]);
	REQUIRE(1.5 ==  a[2]);
}

TEST_CASE("Vector_init_from_std::vector", "[simple]") {
	std::vector<Real> x{1.0, 2.0, 3.0};
    Vector<Real> a(x);

	REQUIRE(3 == a.size());

	REQUIRE(1.0 ==  a[0]);
	REQUIRE(2.0 ==  a[1]);
	REQUIRE(3.0 ==  a[2]);
}

TEST_CASE("Vector_init_from_double_array", "[simple]") {
	Real x[] = {1.0, 2.0, 3.0};
    Vector<Real> a(3, x);

	REQUIRE(3 == a.size());

	REQUIRE(1.0 ==  a[0]);
	REQUIRE(2.0 ==  a[1]);
	REQUIRE(3.0 ==  a[2]);
}

TEST_CASE("Vector_initializer_list_ctor", "[simple]") {
    Vector<Real> a{1.0, 2.0, 3.0};

	REQUIRE(3 == a.size());

	REQUIRE(1.0 == a[0]);
	REQUIRE(2.0 == a[1]);
	REQUIRE(3.0 == a[2]);
}

TEST_CASE("Test_Vector_IsEqual", "[simple]") {
    Vector a({1.0, 2.0});
    Vector b({1.0, 2.0000001});

	REQUIRE(true == a.IsEqual(b, 1e-7));
    REQUIRE(false == a.IsEqual(b, 1e-8));
}

TEST_CASE("Test_Vector_AreEqual", "[simple]") {
    Vector a({1.0, 2.0});
    Vector b({1.0, 2.0000001});

	REQUIRE(true == Vector<Real>::AreEqual(a, b, 1e-7));
    REQUIRE(false == Vector<Real>::AreEqual(a, b, 1e-8));
}

// zbrajanje, oduzimanje
TEST_CASE("Test_Vector_Op+-", "[simple]") {
    Vector a({1.0, 2.0});
    Vector b({1.0, 2.0});

    auto c = a + b;
    auto d = a - b;

	REQUIRE(2.0 ==  c[0]);
	REQUIRE(4.0 ==  c[1]);

	REQUIRE(0.0 ==  d[0]);
	REQUIRE(0.0 ==  d[1]);
}

// op. sa skalaraom
TEST_CASE("Test_Vector_mul_double", "[simple]") {
    Vector a({1.0, 100.0});

	auto b = a * 2.0;
	auto c = 2.0 * a;

	REQUIRE(2.0 ==  b[0]);
	REQUIRE(2.0 ==  c[0]);

	REQUIRE(200.0 ==  b[1]);
	REQUIRE(200.0 ==  c[1]);	
}

TEST_CASE("Test_Vector_div_double", "[simple]") {
    Vector a({4.0, 400.0});

	auto b = a / 2.0;

	REQUIRE(2.0 ==  b[0]);
	REQUIRE(200.0 ==  b[1]);
}

TEST_CASE("Test_Vector_ScalarProductCartesian", "[simple]") {
    Vector a({1.0, 2.0});
    Vector b({1.0, 2.0});

	REQUIRE(5.0 == a.ScalarProductCartesian(b));
}

TEST_CASE("Test_Vector_NormL2", "[simple]") {
    Vector a({2.0, 2.0});

	REQUIRE(sqrt(8) == a.NormL2());
}

TEST_CASE("Test_Vector_AngleToVector", "[simple]") {
    Vector a({0.0, 0.0, 5.0});
    Vector b({1.0, 0.0, 0.0});

	REQUIRE( Constants::PI/2 == a.AngleToVector(b) );
}

TEST_CASE("Test_Vector_AngleToVector2", "[simple]") {
    Vector a({0.0, 1.0, 1.0});
    Vector b({1.0, 0.0, 0.0});

	REQUIRE( Constants::PI/2 == a.AngleToVector(b) );
}

TEST_CASE("Test_Vector_AngleToVector3", "[simple]") {
    Vector a({1.0, 1.0, 0.0});
    Vector b({1.0, 0.0, 0.0});

    REQUIRE_THAT(Constants::PI/4, Catch::Matchers::WithinAbs(a.AngleToVector(b), 1e-16));
}

TEST_CASE("Test_Vector_to_string", "[simple]") {
    Vector a({2.0, 2.0});

	REQUIRE("[    2,     2]" == a.to_string(5,3));

    Vector b({123.0, 1.0, 10.0, -8.0});

	REQUIRE("[    123,       1,      10,      -8]" == b.to_string(7,3));

    Vector c({123.123, 1.9876543, 10.0});

	REQUIRE("[    123.12,     1.9877,         10]" == c.to_string(10,5));
	REQUIRE("[        123.123,       1.9876543,              10]" == c.to_string(15,9));
}

TEST_CASE("Test_Vector_exceptions", "[simple]") 
{
    Vector<Real> vec_dim_5(5);
    Vector<Real> vec_dim_4( {1.0, 0.0, 0.0, 1.0} );

    REQUIRE_THROWS_AS(vec_dim_4.IsEqual(vec_dim_5), VectorDimensionError); 
    REQUIRE_THROWS_AS(vec_dim_4 + vec_dim_5, VectorDimensionError); 
    REQUIRE_THROWS_AS(vec_dim_4 - vec_dim_5, VectorDimensionError); 
    REQUIRE_THROWS_AS(vec_dim_4.ScalarProductCartesian(vec_dim_5), VectorDimensionError); 
}
} // namespace MML::VectorTests