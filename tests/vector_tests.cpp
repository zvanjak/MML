#include "../catch/catch.hpp"

#ifdef MML_USE_SINGLE_HEADER
#include "MML.h"
#else
#include "basic_types/Vector.h"
#endif

TEST_CASE("Vector_default_ctor_init_to_zero", "[simple]") {
    MML::Vector<Real> a(3);

	REQUIRE(0.0 ==  a[0]);
	REQUIRE(0.0 ==  a[1]);
	REQUIRE(0.0 ==  a[2]);
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

	REQUIRE(1.0 ==  a[0]);
	REQUIRE(2.0 ==  a[1]);
	REQUIRE(3.0 ==  a[2]);
}

TEST_CASE("Test_Vector_ScalarProductCartesian", "[simple]") {
    MML::Vector a({1.0, 2.0});
    MML::Vector b({1.0, 2.0});

	REQUIRE(5.0 ==  a.ScalarProductCartesian(b));
}

TEST_CASE("Test_Vector_NormL2", "[simple]") {
    MML::Vector a({2.0, 2.0});

	REQUIRE(sqrt(8) ==  a.NormL2());
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

// kako reagira na index violation?
