#include "../catch/catch.hpp"

#ifdef MML_USE_SINGLE_HEADER
#include "MML.h"
#else
#include "base/VectorN.h"
#endif

using namespace MML;

namespace MML::Tests::Base::VectorNTests
{
	TEST_CASE("VectorN::default_ctor_init_to_zero", "[VectorN]") {
		VectorN<Real, 3> a;

		REQUIRE(0.0 == a[0]);
		REQUIRE(0.0 == a[1]);
		REQUIRE(0.0 == a[2]);
	}

	TEST_CASE("VectorN::initializer_list_ctor", "[VectorN]") {
		VectorN<Real, 3> a{ 1.0, 2.0, 3.0 };

		REQUIRE(1.0 == a[0]);
		REQUIRE(2.0 == a[1]);
		REQUIRE(3.0 == a[2]);
	}

	TEST_CASE("VectorN::init_to_value", "[VectorN]") {
		VectorN<Real, 3> a(1.5);

		REQUIRE(3 == a.size());

		REQUIRE(1.5 == a[0]);
		REQUIRE(1.5 == a[1]);
		REQUIRE(1.5 == a[2]);
	}

	TEST_CASE("VectorN::init_from_std::vector", "[VectorN]") {
		std::vector<Real> x{ 1.0, 2.0, 3.0 };
		VectorN<Real, 3> a(x);

		REQUIRE(3 == a.size());

		REQUIRE(1.0 == a[0]);
		REQUIRE(2.0 == a[1]);
		REQUIRE(3.0 == a[2]);
	}

	TEST_CASE("VectorN::init_from_double_array", "[VectorN]") {
		Real x[] = { 1.0, 2.0, 3.0 };
		VectorN<Real, 3> a(x);

		REQUIRE(3 == a.size());

		REQUIRE(1.0 == a[0]);
		REQUIRE(2.0 == a[1]);
		REQUIRE(3.0 == a[2]);
	}

	TEST_CASE("VectorN::IsEqual", "[VectorN]") {
		Vec2Dbl a({ 1.0, 2.0 });
		Vec2Dbl b({ 1.0, 2.0000001 });

		REQUIRE(true == a.IsEqualTo(b, 1e-7));
		REQUIRE(false == a.IsEqualTo(b, 1e-8));
	}

	TEST_CASE("VectorN::Op+-", "[VectorN]") {
		Vec2Dbl a({ 1.0, 2.0 });
		Vec2Dbl b({ 1.0, 2.0 });

		auto c = a + b;
		auto d = a - b;

		REQUIRE(2.0 == c[0]);
		REQUIRE(4.0 == c[1]);

		REQUIRE(0.0 == d[0]);
		REQUIRE(0.0 == d[1]);
	}

	TEST_CASE("VectorN::mul_double", "[VectorN]") {
		Vec2Dbl a({ 1.0, 100.0 });

		auto b = a * 2.0;
		auto c = 2.0 * a;

		REQUIRE(2.0 == b[0]);
		REQUIRE(2.0 == c[0]);

		REQUIRE(200.0 == b[1]);
		REQUIRE(200.0 == c[1]);
	}

	TEST_CASE("VectorN::div_double", "[VectorN]") {
		Vec2Dbl a({ 4.0, 400.0 });

		auto b = a / 2.0;

		REQUIRE(2.0 == b[0]);
		REQUIRE(200.0 == b[1]);
	}

	TEST_CASE("VectorN::ScalarProductCartesian", "[VectorN]") {
		Vec2Dbl a({ 1.0, 2.0 });
		Vec2Dbl b({ 1.0, 2.0 });

		REQUIRE(5.0 == a.ScalarProductCartesian(b));
	}

	TEST_CASE("VectorN::NormL2", "[VectorN]") {
		Vec2Dbl a({ 2.0, 2.0 });

		REQUIRE(sqrt(8) == a.NormL2());
	}

	TEST_CASE("VectorN::to_string", "[VectorN]") {
		Vec2Dbl a({ 2.0, 2.0 });

		REQUIRE("[    2,     2]" == a.to_string(5, 3));

		Vec4Dbl b({ 123.0, 1.0, 10.0, -8.0 });

		REQUIRE("[    123,       1,      10,      -8]" == b.to_string(7, 3));

		Vec3Dbl c({ 123.123, 1.9876543, 10.0 });

		REQUIRE("[    123.12,     1.9877,         10]" == c.to_string(10, 5));
		REQUIRE("[        123.123,       1.9876543,              10]" == c.to_string(15, 9));
	}

	TEST_CASE("VectorN::exceptions", "[VectorN]")
	{
		Vec2Dbl vec2({ 2.0, 2.0 });
		REQUIRE_THROWS_AS(vec2.at(5), VectorDimensionError);
		REQUIRE_THROWS_AS(vec2.at(-1), VectorDimensionError);

		REQUIRE_THROWS_AS(vec2.at(4) = 5.0, VectorDimensionError);
		REQUIRE_THROWS_AS(vec2.at(-1) = 5.0, VectorDimensionError);
	}
} // namespace MML::Tests::VectorNTests