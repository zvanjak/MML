#include <catch2/catch_all.hpp>
#include "../TestPrecision.h"
#include "../TestMatchers.h"

#ifdef MML_USE_SINGLE_HEADER
#include "MML.h"
#else
#include "base/VectorN.h"
#include "base/BaseUtils.h"
#endif

using namespace MML;
using namespace MML::Testing;

namespace MML::Tests::Base::VectorNTests
{
	TEST_CASE("VectorN::default_ctor_init_to_zero", "[VectorN]") {
			TEST_PRECISION_INFO();
		VectorN<Real, 3> a;

		REQUIRE(REAL(0.0) == a[0]);
		REQUIRE(REAL(0.0) == a[1]);
		REQUIRE(REAL(0.0) == a[2]);
	}

	TEST_CASE("VectorN::initializer_list_ctor", "[VectorN]") {
			TEST_PRECISION_INFO();
		VectorN<Real, 3> a{ REAL(1.0), REAL(2.0), REAL(3.0) };

		REQUIRE(REAL(1.0) == a[0]);
		REQUIRE(REAL(2.0) == a[1]);
		REQUIRE(REAL(3.0) == a[2]);
	}

	TEST_CASE("VectorN::init_to_value", "[VectorN]") {
			TEST_PRECISION_INFO();
		VectorN<Real, 3> a(REAL(1.5));

		REQUIRE(3 == a.size());

		REQUIRE(REAL(1.5) == a[0]);
		REQUIRE(REAL(1.5) == a[1]);
		REQUIRE(REAL(1.5) == a[2]);
	}

	TEST_CASE("VectorN::init_from_ std vector", "[VectorN]") {
			TEST_PRECISION_INFO();
		std::vector<Real> x{ REAL(1.0), REAL(2.0), REAL(3.0) };
		VectorN<Real, 3> a(x);

		REQUIRE(3 == a.size());

		REQUIRE(REAL(1.0) == a[0]);
		REQUIRE(REAL(2.0) == a[1]);
		REQUIRE(REAL(3.0) == a[2]);
	}

	TEST_CASE("VectorN::init_from_double_array", "[VectorN]") {
			TEST_PRECISION_INFO();
		Real x[] = { REAL(1.0), REAL(2.0), REAL(3.0) };
		VectorN<Real, 3> a(x);

		REQUIRE(3 == a.size());

		REQUIRE(REAL(1.0) == a[0]);
		REQUIRE(REAL(2.0) == a[1]);
		REQUIRE(REAL(3.0) == a[2]);
	}

	TEST_CASE("VectorN::IsEqual", "[VectorN]") {
			TEST_PRECISION_INFO();
		Vec2Dbl a({ REAL(1.0), REAL(2.0) });
		// Difference of 1e-5 is detectable across all precisions after scaling
		Vec2Dbl b({ REAL(1.0), REAL(2.00001) });

		// 1e-4 scales to 1e-2 (float), 1e-8 scales to 1e-6 (float)
		REQUIRE(true == a.IsEqualTo(b, ScaleTolerance(REAL(1e-4))));
		REQUIRE(false == a.IsEqualTo(b, ScaleTolerance(REAL(1e-8))));
	}

	TEST_CASE("VectorN::Op+-", "[VectorN]") {
			TEST_PRECISION_INFO();
		Vec2Dbl a({ REAL(1.0), REAL(2.0) });
		Vec2Dbl b({ REAL(1.0), REAL(2.0) });

		auto c = a + b;
		auto d = a - b;

		REQUIRE(REAL(2.0) == c[0]);
		REQUIRE(REAL(4.0) == c[1]);

		REQUIRE(REAL(0.0) == d[0]);
		REQUIRE(REAL(0.0) == d[1]);
	}

	TEST_CASE("VectorN::mul_double", "[VectorN]") {
			TEST_PRECISION_INFO();
		Vec2Dbl a({ REAL(1.0), REAL(100.0) });

		auto b = a * REAL(2.0);
		auto c = REAL(2.0) * a;

		REQUIRE(REAL(2.0) == b[0]);
		REQUIRE(REAL(2.0) == c[0]);

		REQUIRE(REAL(200.0) == b[1]);
		REQUIRE(REAL(200.0) == c[1]);
	}

	TEST_CASE("VectorN::div_double", "[VectorN]") {
			TEST_PRECISION_INFO();
		Vec2Dbl a({ REAL(4.0), REAL(400.0) });

		auto b = a / REAL(2.0);

		REQUIRE(REAL(2.0) == b[0]);
		REQUIRE(REAL(200.0) == b[1]);
	}

	TEST_CASE("VectorN::NormL2", "[VectorN]") {
			TEST_PRECISION_INFO();
		Vec2Dbl a({ REAL(2.0), REAL(2.0) });

		REQUIRE(sqrt(8) == a.NormL2());
	}

	TEST_CASE("VectorN::to_string", "[VectorN]") {
			TEST_PRECISION_INFO();
		Vec2Dbl a({ REAL(2.0), REAL(2.0) });

		REQUIRE("[    2,     2]" == a.to_string(5, 3));

		Vec4Dbl b({ REAL(123.0), REAL(1.0), REAL(10.0), -REAL(8.0) });

		REQUIRE("[    123,       1,      10,      -8]" == b.to_string(7, 3));

		Vec3Dbl c({ REAL(123.123), REAL(1.9876543), REAL(10.0) });

	REQUIRE("[    123.12,     1.9877,         10]" == c.to_string(10, 5));
	REQUIRE("[        123.123,       1.9876543,              10]" == c.to_string(15, 9));
	}

	TEST_CASE("VectorN::exceptions", "[VectorN]")
	{
			TEST_PRECISION_INFO();
		Vec2Dbl vec2({ REAL(2.0), REAL(2.0) });
		REQUIRE_THROWS_AS(vec2.at(5), VectorDimensionError);
		REQUIRE_THROWS_AS(vec2.at(-1), VectorDimensionError);

		REQUIRE_THROWS_AS(vec2.at(4) = REAL(5.0), VectorDimensionError);
		REQUIRE_THROWS_AS(vec2.at(-1) = REAL(5.0), VectorDimensionError);
	}
} // namespace MML::Tests::VectorNTests