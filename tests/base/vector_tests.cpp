#include "../catch/catch.hpp"

#ifdef MML_USE_SINGLE_HEADER
#include "MML.h"
#else
#include "base/Vector.h"
#endif

using namespace MML;

namespace MML::Tests::Base::VectorTests
{
	TEST_CASE("Vector::Default_ctor_has_zero_size", "[simple]")
	{
		Vector<Real> b, c(0);
		REQUIRE(0 == b.size());
		REQUIRE(0 == c.size());
	}

	TEST_CASE("Vector::Ctor_init_elems_to_zero", "[simple]") {
		Vector<Real> a(3);

		REQUIRE(3 == a.size());
		REQUIRE(0.0 == a[0]);
		REQUIRE(0.0 == a[1]);
		REQUIRE(0.0 == a[2]);

		Vector<Real> big(10000);

		REQUIRE(10000 == big.size());
		REQUIRE(0.0 == big[0]);
		REQUIRE(0.0 == big[5000]);
		REQUIRE(0.0 == big[7258]);
	}

	TEST_CASE("Vector::Ctor_init_to_value", "[simple]") {
		Vector<Real> a(3, 1.5);

		REQUIRE(3 == a.size());

		REQUIRE(1.5 == a[0]);
		REQUIRE(1.5 == a[1]);
		REQUIRE(1.5 == a[2]);
	}

	TEST_CASE("Vector::Ctor_init_from_real_array", "[simple]") {
		Real x[] = { 1.0, 2.0, 3.0 };
		Vector<Real> a(3, x);

		REQUIRE(3 == a.size());

		REQUIRE(1.0 == a[0]);
		REQUIRE(2.0 == a[1]);
		REQUIRE(3.0 == a[2]);
	}

	TEST_CASE("Vector::Ctor_negative_size_throws", "[simple]")
	{
		REQUIRE_THROWS_AS(Vector<Real>(-5), VectorInitializationError);
		REQUIRE_THROWS_AS(Vector<Real>(-10, 0.5), VectorInitializationError);

		Real x[] = { 1.0, 2.0, 3.0 };
		REQUIRE_THROWS_AS(Vector<Real>(-3, x), VectorInitializationError);
	}

	TEST_CASE("Vector::Ctor_init_from_std vector", "[simple]") {
		std::vector<Real> x{ 1.0, 2.0, 3.0 };
		Vector<Real> a(x);

		REQUIRE(3 == a.size());

		REQUIRE(1.0 == a[0]);
		REQUIRE(2.0 == a[1]);
		REQUIRE(3.0 == a[2]);
	}

	TEST_CASE("Vector::Ctor_initializer_list", "[simple]") {
		Vector<Real> a{ 1.0, 2.0, 3.0 };

		REQUIRE(3 == a.size());

		REQUIRE(1.0 == a[0]);
		REQUIRE(2.0 == a[1]);
		REQUIRE(3.0 == a[2]);
	}

	TEST_CASE("Vector::IsEqual", "[simple]") {
		Vector<Real> a({ 1.0, 2.0 });
		Vector<Real> b({ 1.0, 2.0000001 });

		REQUIRE(true == a.IsEqualTo(b, 1e-7));
		REQUIRE(false == a.IsEqualTo(b, 1e-8));
	}

	TEST_CASE("Vector::AreEqual", "[simple]") {
		Vector<Real> a({ 1.0, 2.0 });
		Vector<Real> b({ 1.0, 2.0000001 });

		REQUIRE(true == Vector<Real>::AreEqual(a, b, 1e-7));
		REQUIRE(false == Vector<Real>::AreEqual(a, b, 1e-8));
	}

	TEST_CASE("Vector::access_operators", "[simple]") {
		Vector<Real> a({ 1.0, 2.0, 3.0 });

		REQUIRE(1.0 == a[0]);
		REQUIRE(2.0 == a[1]);
		REQUIRE(3.0 == a[2]);

		REQUIRE(1.0 == a.at(0));
		REQUIRE(2.0 == a.at(1));
		REQUIRE(3.0 == a.at(2));
	}

	TEST_CASE("Vector::template_type_deduction", "[simple]") {
		Vector vec_double({ 1.0, 2.0, 3.0 });

		REQUIRE(1.0 == vec_double[0]);
		REQUIRE(2.0 == vec_double[1]);
		REQUIRE(3.0 == vec_double[2]);

		Vector vec_float({ 1.0F, 2.0F, 3.0F });

		REQUIRE(1.0F == vec_float[0]);
		REQUIRE(2.0F == vec_float[1]);
		REQUIRE(3.0F == vec_float[2]);

		Vector b_int({ 1, 2, 3 });

		REQUIRE(1 == b_int[0]);
		REQUIRE(2 == b_int[1]);
		REQUIRE(3 == b_int[2]);

		Vector vec_cmplx({ Complex(1,1), Complex(0,2) });

		REQUIRE(Complex(1, 1) == vec_cmplx[0]);
		REQUIRE(Complex(0, 2) == vec_cmplx[1]);

		int x = 5;
		Vector vec_dbl_1({ x * 1.0, sin(2.0), 3.0 });
		// CAN'T DO THIS!
		// Vector vec_dbl_2({ 1 , sin(2.0), 3.0 });
		// needs explicit template type argument
		Vector<Real> vec_dbl_2({ 1 , sin(2.0), 3.0 });

		Vector vec_cmplx_1({ 3.0 * Complex(1,1), Complex(0,2) - sin(2.5) });
		Vector<Complex> vec_cmplx_2({ 1, 3.5678, Complex(1, 1), sin(2.5) });
	}

	TEST_CASE("Vector::Op+-", "[simple]") {
		Vector a({ 1.0, 2.0 });
		Vector b({ 1.0, 2.0 });

		auto c = a + b;
		auto d = a - b;

		REQUIRE(2.0 == c[0]);
		REQUIRE(4.0 == c[1]);

		REQUIRE(0.0 == d[0]);
		REQUIRE(0.0 == d[1]);
	}

	TEST_CASE("Vector::mul_double", "[simple]") {
		Vector a({ 1.0, 100.0 });

		auto b = a * 2.0;
		auto c = 2.0 * a;

		REQUIRE(2.0 == b[0]);
		REQUIRE(2.0 == c[0]);

		REQUIRE(200.0 == b[1]);
		REQUIRE(200.0 == c[1]);
	}

	TEST_CASE("Vector::div_double", "[simple]") {
		Vector a({ 4.0, 400.0 });

		auto b = a / 2.0;

		REQUIRE(2.0 == b[0]);
		REQUIRE(200.0 == b[1]);
	}

	TEST_CASE("Vector::ScalarProductCartesian", "[simple]") {
		Vector a({ 1.0, 2.0 });
		Vector b({ 1.0, 2.0 });

		REQUIRE(5.0 == a.ScalarProductCartesian(b));
	}

	TEST_CASE("Vector::NormL2", "[simple]") {
		Vector a({ 2.0, 2.0 });

		REQUIRE(sqrt(8) == a.NormL2());
	}

	TEST_CASE("Vector::AngleToVector", "[simple]") {
		Vector a({ 0.0, 0.0, 5.0 });
		Vector b({ 1.0, 0.0, 0.0 });

		REQUIRE(Constants::PI / 2 == a.AngleToVector(b));

		Vector c({ 0.0, 1.0 });
		Vector d({ 1.0, 1.0 });

		REQUIRE_THAT(Constants::PI / 4, Catch::Matchers::WithinAbs(c.AngleToVector(d), 1e-16));
	}

	TEST_CASE("Vector::to_string", "[simple]") {
		Vector a({ 2.0, 2.0 });

		REQUIRE("[    2,     2]" == a.to_string(5, 3));

		Vector b({ 123.0, 1.0, 10.0, -8.0 });

		REQUIRE("[    123,       1,      10,      -8]" == b.to_string(7, 3));

		Vector c({ 123.123, 1.9876543, 10.0 });

		REQUIRE("[    123.12,     1.9877,         10]" == c.to_string(10, 5));
		REQUIRE("[        123.123,       1.9876543,              10]" == c.to_string(15, 9));
	}

	TEST_CASE("Vector::exceptions", "[simple]")
	{
		Vector<Real> vec_dim_5(5);
		Vector<Real> vec_dim_4({ 1.0, 0.0, 0.0, 1.0 });

		REQUIRE_THROWS_AS(Vector<Real>::GetUnitVector(-1, 5), VectorDimensionError);
		REQUIRE_THROWS_AS(Vector<Real>::GetUnitVector(5, 7), VectorDimensionError);

		REQUIRE_THROWS_AS(vec_dim_4.IsEqualTo(vec_dim_5), VectorDimensionError);

		REQUIRE_THROWS_AS(vec_dim_4.at(5), VectorDimensionError);
		REQUIRE_THROWS_AS(vec_dim_4.at(-1), VectorDimensionError);
		REQUIRE_THROWS_AS(vec_dim_4.at(4) = 5.0, VectorDimensionError);
		REQUIRE_THROWS_AS(vec_dim_4.at(-1) = 5.0, VectorDimensionError);

		REQUIRE_THROWS_AS(vec_dim_4 + vec_dim_5, VectorDimensionError);
		REQUIRE_THROWS_AS(vec_dim_4 - vec_dim_5, VectorDimensionError);
		REQUIRE_THROWS_AS(vec_dim_4 == vec_dim_5, VectorDimensionError);

		REQUIRE_THROWS_AS(vec_dim_4.ScalarProductCartesian(vec_dim_5), VectorDimensionError);
		REQUIRE_THROWS_AS(vec_dim_4.AngleToVector(vec_dim_5), VectorDimensionError);
	}
} // namespace MML::VectorTests