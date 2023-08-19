#include "../catch/catch.hpp"

#ifdef MML_USE_SINGLE_HEADER
#include "MML.h"
#else
#include "core/Vector.h"
#endif

TEST_CASE("Vector_default_ctor_init_to_zero", "[simple]") {
    MML::Vector<Real> a(3);

    REQUIRE(3 == a.size());

	REQUIRE(0.0 ==  a[0]);
	REQUIRE(0.0 ==  a[1]);
	REQUIRE(0.0 ==  a[2]);

    MML::Vector<Real> b(0);
    REQUIRE(0 == b.size());

    // TODO - check exception thrown when less than zero

}

TEST_CASE("Vector_init_to_value", "[simple]") {
    MML::Vector<Real> a(3, 1.5);

	REQUIRE(3 == a.size());

	REQUIRE(1.5 ==  a[0]);
	REQUIRE(1.5 ==  a[1]);
	REQUIRE(1.5 ==  a[2]);
}

TEST_CASE("Vector_init_from_std::vector", "[simple]") {
	std::vector<double> x{1.0, 2.0, 3.0};
    MML::Vector<Real> a(x);

	REQUIRE(3 == a.size());

	REQUIRE(1.0 ==  a[0]);
	REQUIRE(2.0 ==  a[1]);
	REQUIRE(3.0 ==  a[2]);
}

TEST_CASE("Vector_init_from_double_array", "[simple]") {
	double x[] = {1.0, 2.0, 3.0};
    MML::Vector<Real> a(3, x);

	REQUIRE(3 == a.size());

	REQUIRE(1.0 ==  a[0]);
	REQUIRE(2.0 ==  a[1]);
	REQUIRE(3.0 ==  a[2]);
}

TEST_CASE("Vector_initializer_list_ctor", "[simple]") {
    MML::Vector<Real> a{1.0, 2.0, 3.0};

	REQUIRE(3 == a.size());

	REQUIRE(1.0 == a[0]);
	REQUIRE(2.0 == a[1]);
	REQUIRE(3.0 == a[2]);
}

TEST_CASE("Test_Vector_IsEqual", "[simple]") {
    MML::Vector a({1.0, 2.0});
    MML::Vector b({1.0, 2.0000001});

	REQUIRE(true == a.IsEqual(b, 1e-7));
    REQUIRE(false == a.IsEqual(b, 1e-8));
}

// zbrajanje, oduzimanje
TEST_CASE("Test_Vector_Op+-", "[simple]") {
    MML::Vector a({1.0, 2.0});
    MML::Vector b({1.0, 2.0});

    auto c = a + b;
    auto d = a - b;

	REQUIRE(2.0 ==  c[0]);
	REQUIRE(4.0 ==  c[1]);

	REQUIRE(0.0 ==  d[0]);
	REQUIRE(0.0 ==  d[1]);
}

// op. sa skalaraom
TEST_CASE("Test_Vector_mul_double", "[simple]") {
    MML::Vector a({1.0, 100.0});

	auto b = a * 2.0;
	auto c = 2.0 * a;

	REQUIRE(2.0 ==  b[0]);
	REQUIRE(2.0 ==  c[0]);

	REQUIRE(200.0 ==  b[1]);
	REQUIRE(200.0 ==  c[1]);	
}

TEST_CASE("Test_Vector_div_double", "[simple]") {
    MML::Vector a({4.0, 400.0});

	auto b = a / 2.0;

	REQUIRE(2.0 ==  b[0]);
	REQUIRE(200.0 ==  b[1]);
}

TEST_CASE("Test_Vector_ScalarProductCartesian", "[simple]") {
    MML::Vector a({1.0, 2.0});
    MML::Vector b({1.0, 2.0});

	REQUIRE(5.0 == a.ScalarProductCartesian(b));
}

TEST_CASE("Test_Vector_NormL2", "[simple]") {
    MML::Vector a({2.0, 2.0});

	REQUIRE(sqrt(8) == a.NormL2());
}

TEST_CASE("Test_Vector_to_string", "[simple]") {
    MML::Vector a({2.0, 2.0});

	REQUIRE("[    2,     2]" == a.to_string(5,3));

    MML::Vector b({123.0, 1.0, 10.0, -8.0});

	REQUIRE("[    123,       1,      10,      -8]" == b.to_string(7,3));

    MML::Vector c({123.123, 1.9876543, 10.0});

	REQUIRE("[    123.12,     1.9877,         10]" == c.to_string(10,5));
	REQUIRE("[        123.123,       1.9876543,              10]" == c.to_string(15,9));
}

TEST_CASE("Test_Vector_exceptions", "[simple]") 
{
    MML::Vector<Real> vec_dim_5(5);
    MML::Vector<Real> vec_dim_4( {1.0, 0.0, 0.0, 1.0} );

    REQUIRE_THROWS_AS(vec_dim_4.IsEqual(vec_dim_5), MML::VectorDimensionError); 
    REQUIRE_THROWS_AS(vec_dim_4 + vec_dim_5, MML::VectorDimensionError); 
    REQUIRE_THROWS_AS(vec_dim_4 - vec_dim_5, MML::VectorDimensionError); 
    REQUIRE_THROWS_AS(vec_dim_4.ScalarProductCartesian(vec_dim_5), MML::VectorDimensionError); 
}