#include <catch2/catch_test_macros.hpp>
#ifdef MML_SINGLE_HEADER_INCLUDE
#include MML_SINGLE_HEADER_INCLUDE
#else
#include <mml/base/Vector/Vector.h>
#include <mml/base/Vector/VectorTypes2D.h>
#include <mml/base/Vector/VectorTypes3D.h>
#include <mml/base/Vector/VectorTypes4D.h>
#endif

#include <algorithm>
#include <array>
#include <cmath>
#include <complex>
#include <limits>
#include <type_traits>

namespace {
	using namespace MML;
	// Keep both scales finite while their naive squares overflow or underflow for every Real type.
	constexpr int ExtremeScaleExponent = std::numeric_limits<Real>::max_exponent * 2 / 3;
	const Real TinyCoordinateScale = std::scalbn(Real{1}, -ExtremeScaleExponent);
	const Real HugeCoordinateScale = std::scalbn(Real{1}, ExtremeScaleExponent);
	constexpr int DirectionScaleExponent = std::min(1000, std::numeric_limits<Real>::max_exponent - 4);

	std::array<Real, 3> Cartesian(const Vector2Polar& v) {
		return {v.R() * std::cos(v.Phi()), v.R() * std::sin(v.Phi()), 0};
	}
	std::array<Real, 3> Cartesian(const Vector3Spherical& v) {
		return {v.R() * std::sin(v.Theta()) * std::cos(v.Phi()),
			v.R() * std::sin(v.Theta()) * std::sin(v.Phi()), v.R() * std::cos(v.Theta())};
	}
	std::array<Real, 3> Cartesian(const Vector3Cylindrical& v) {
		return {v.R() * std::cos(v.Phi()), v.R() * std::sin(v.Phi()), v.Z()};
	}
	void CheckCartesian(const std::array<Real, 3>& actual, const std::array<Real, 3>& expected) {
		const Real scale = std::max({std::abs(expected[0]), std::abs(expected[1]), std::abs(expected[2]), Real{1}});
		for (size_t i = 0; i < actual.size(); ++i) {
			CHECK(std::isfinite(actual[i]));
			CHECK(std::abs(actual[i] / scale - expected[i] / scale) <= Precision::EqualityTolerance);
		}
	}
	template<class Coordinate>
	void CheckCompound(Coordinate a, Coordinate b) {
		static_assert(std::is_same_v<decltype(a += b), Coordinate&>);
		static_assert(std::is_same_v<decltype(a -= b), Coordinate&>);
		static_assert(std::is_same_v<decltype(a *= Real{2}), Coordinate&>);
		static_assert(std::is_same_v<decltype(a /= Real{2}), Coordinate&>);
		const auto first = Cartesian(a), second = Cartesian(b);
		auto compound = a;
		compound += b;
		CheckCartesian(Cartesian(compound), Cartesian(a + b));
		CheckCartesian(Cartesian(compound), {first[0] + second[0], first[1] + second[1], first[2] + second[2]});
		compound = a;
		compound -= b;
		CheckCartesian(Cartesian(compound), Cartesian(a - b));
		CheckCartesian(Cartesian(compound), {first[0] - second[0], first[1] - second[1], first[2] - second[2]});
		compound = a;
		compound *= Real{-2};
		CheckCartesian(Cartesian(compound), Cartesian(a * Real{-2}));
		compound = a;
		compound /= Real{2};
		CheckCartesian(Cartesian(compound), Cartesian(a / Real{2}));
		compound = a;
		compound += compound;
		CheckCartesian(Cartesian(compound), Cartesian(a + a));
		compound = a;
		compound -= compound;
		CheckCartesian(Cartesian(compound), {0, 0, 0});
	}
}

TEST_CASE("Coordinate scale fixtures are finite and stress squared arithmetic", "[BaseScaleContracts][CoordinateScaleContracts]")
{
	REQUIRE(std::isfinite(TinyCoordinateScale));
	REQUIRE(TinyCoordinateScale > 0);
	REQUIRE(std::isfinite(HugeCoordinateScale));
	REQUIRE(HugeCoordinateScale > 1);
	CHECK(std::isinf(HugeCoordinateScale * HugeCoordinateScale));
	CHECK(TinyCoordinateScale * TinyCoordinateScale == 0);
}

TEST_CASE("Coordinate compounds agree with geometric binary operations", "[BaseScaleContracts][CoordinateScaleContracts]")
{
	CheckCompound(Vector2Polar(2, Real{0.7}), Vector2Polar(3, Real{-0.4}));
	CheckCompound(Vector3Spherical(2, Real{1.1}, Real{0.7}), Vector3Spherical(3, Real{2.2}, Real{-0.4}));
	CheckCompound(Vector3Cylindrical(2, Real{0.7}, -4), Vector3Cylindrical(3, Real{-0.4}, 5));
	Vector3Spherical spherical(1, 1, 1);
	Vector3Cylindrical cylindrical(1, 1, 1);
	CHECK_THROWS_AS(spherical /= Real{0}, DivisionByZeroError);
	CHECK_THROWS_AS(cylindrical /= Real{0}, DivisionByZeroError);
}

TEST_CASE("Polar arithmetic resolves close angles and representable extreme radii", "[BaseScaleContracts][CoordinateScaleContracts]")
{
	for (const Real radius : {TinyCoordinateScale, Real{1}, HugeCoordinateScale}) {
		REQUIRE(std::isfinite(radius));
		REQUIRE(radius > 0);
		for (const Real firstAngle : {Real{0}, Real{0.7}, Real{10}}) {
			const Real secondAngle = std::max(firstAngle + Real{1e-10},
				std::nextafter(firstAngle, std::numeric_limits<Real>::infinity()));
			REQUIRE(secondAngle > firstAngle);
			const Vector2Polar first(radius, firstAngle), second(radius, secondAngle);
			const auto difference = first - second;
			const Real expected = radius * (Real{2} * std::abs(std::sin((secondAngle - firstAngle) / Real{2})));
			CHECK(std::isfinite(difference.R()));
			CHECK(difference.R() > 0);
			CHECK(std::abs(difference.R() / expected - Real{1}) <= Precision::EqualityTolerance);
		}
		for (const auto& angles : {std::array<Real, 2>{0, HugeCoordinateScale},
			std::array<Real, 2>{HugeCoordinateScale, -HugeCoordinateScale},
			std::array<Real, 2>{HugeCoordinateScale, HugeCoordinateScale}}) {
			const Vector2Polar first(2, angles[0]), second(3, angles[1]);
			CheckCompound(first, second);
			CHECK(std::abs((first + second).Phi()) <= Constants::PI);
			CHECK(std::abs((first - second).Phi()) <= Constants::PI);
		}
		const auto sum = Vector2Polar(radius, 0) + Vector2Polar(radius, Constants::PI / Real{2});
		CHECK(std::abs(sum.R() / radius - std::sqrt(Real{2})) <= Precision::EqualityTolerance);
	}
	CHECK_THROWS_AS((Vector2Polar(1, std::numeric_limits<Real>::quiet_NaN()) + Vector2Polar(1, 0)), NumericInputError);
}

TEST_CASE("Cartesian direction queries are scale independent and finite guarded", "[BaseScaleContracts][CoordinateScaleContracts]")
{
	for (const int exponent : {-DirectionScaleExponent, -DirectionScaleExponent * 3 / 5,
		-DirectionScaleExponent / 10, 0, DirectionScaleExponent / 10,
		DirectionScaleExponent * 3 / 5, DirectionScaleExponent}) {
		const Real scale = std::scalbn(Real{1}, exponent);
		REQUIRE(std::isfinite(scale));
		REQUIRE(scale > 0);
		Vector3Cartesian first(scale, 0, 0), diagonal(scale, scale, 0), perpendicular(0, scale, 0);
		CHECK(first.IsParallelTo(first));
		CHECK(first.IsParallelTo(-first));
		CHECK_FALSE(first.IsPerpendicularTo(first));
		CHECK(first.IsPerpendicularTo(perpendicular));
		CHECK_FALSE(first.IsParallelTo(diagonal));
		CHECK(std::abs(first.AngleToVector(diagonal) - Constants::PI / Real{4}) <= Precision::EqualityTolerance);
		Vector3Cartesian u, v;
		REQUIRE(first.GetPerpendicularVectors(u, v));
		CHECK(std::abs(u.NormL2() - Real{1}) <= Precision::EqualityTolerance);
		CHECK(std::abs(v.NormL2() - Real{1}) <= Precision::EqualityTolerance);
		CHECK(first.IsPerpendicularTo(u));
		CHECK(first.IsPerpendicularTo(v));
	}
	const Real maximum = std::numeric_limits<Real>::max();
	Vector3Cartesian huge(maximum, maximum, 0), axis(1, 0, 0), tiny(std::numeric_limits<Real>::denorm_min(), 0, 0);
	CHECK(huge.IsParallelTo(huge));
	CHECK(std::abs(huge.AngleToVector(axis) - Constants::PI / Real{4}) <= Precision::EqualityTolerance);
	CHECK(tiny.IsParallelTo(axis));
	CHECK(tiny.AngleToVector(axis) == 0);
	Vector3Cartesian zero, invalid(std::numeric_limits<Real>::quiet_NaN(), 0, 0);
	CHECK_FALSE(axis.IsParallelTo(zero));
	CHECK_FALSE(axis.IsPerpendicularTo(zero));
	CHECK_FALSE(axis.IsParallelTo(invalid));
	CHECK_FALSE(axis.IsPerpendicularTo(invalid));
	CHECK_THROWS_AS(axis.AngleToVector(zero), DivisionByZeroError);
	CHECK_THROWS_AS(axis.AngleToVector(invalid), NumericInputError);
	for (const Real eps : {Real{-1}, std::numeric_limits<Real>::infinity(), std::numeric_limits<Real>::quiet_NaN()}) {
		CHECK_FALSE(axis.IsParallelTo(axis, eps));
		CHECK_FALSE(axis.IsPerpendicularTo(axis, eps));
		CHECK_FALSE(axis.IsEqualTo(axis, eps));
	}
}

TEST_CASE("Coordinate unit directions survive enormous and tiny magnitudes", "[BaseScaleContracts][CoordinateScaleContracts]")
{
	const Vector3Spherical position(1, Real{0.4}, Real{0.9});
	const Vector3Cylindrical cylindricalPosition(1, Real{0.9}, 0);
	const auto sphericalReference = Vector3Spherical(1, Real{1.1}, Real{-0.6}).GetAsUnitVectorAtPos(position);
	const auto globalReference = Vector3Spherical(1, Real{1.1}, Real{-0.6}).GetAsUnitVectorAtPosGlobalSpher(position);
	const auto cylindricalReference = Vector3Cylindrical(1, Real{-0.6}, 2).GetAsUnitVectorAtPos(cylindricalPosition);
	for (const Real scale : {TinyCoordinateScale, HugeCoordinateScale}) {
		REQUIRE(std::isfinite(scale));
		REQUIRE(scale > 0);
		const Vector3Spherical spherical(scale, Real{1.1}, Real{-0.6});
		const Vector3Cylindrical cylindrical(scale, Real{-0.6}, Real{2} * scale);
		CHECK(spherical.GetAsUnitVectorAtPos(position).IsEqualTo(sphericalReference));
		CHECK(spherical.GetAsUnitVectorAtPosGlobalSpher(position).IsEqualTo(globalReference));
		CHECK(cylindrical.GetAsUnitVectorAtPos(cylindricalPosition).IsEqualTo(cylindricalReference));
	}
}

TEST_CASE("Integer and unsigned vector members instantiate with explicit erase semantics", "[BaseScaleContracts][CoordinateScaleContracts]")
{
	auto unit = Vector<int>::UnitVector(3, 0);
	CHECK(unit == Vector<int>{1, 0, 0});
	unit.erase(0);
	CHECK(unit == Vector<int>{0, 0});
	unit.eraseValue(0);
	CHECK(unit.empty());
	Vector<int> values{1, 2, 1, 3};
	values.eraseValue(1);
	CHECK(values == Vector<int>{2, 3});
	Vector<unsigned> dynamic{2, 4};
	CHECK(dynamic.safeDivide(2) == Vector<unsigned>{1, 2});
	CHECK_THROWS_AS(dynamic.safeDivide(0), DivisionByZeroError);
	CHECK(Vector<int>{std::numeric_limits<int>::min()}.safeDivide(std::numeric_limits<int>::min()) == Vector<int>{1});
	CHECK((VectorN<unsigned, 2>{0, 0}.IsNearZero()));
	CHECK_FALSE((VectorN<unsigned, 2>{1, 0}.IsNearZero()));
	CHECK(Vector<unsigned>::UnitVector(2, 1) == Vector<unsigned>{0, 1});
}

TEST_CASE("Vector finite predicates retain the repaired common-helper contract", "[BaseScaleContracts][CoordinateScaleContracts]")
{
	const Real nan = std::numeric_limits<Real>::quiet_NaN(), inf = std::numeric_limits<Real>::infinity();
	for (const Real value : {nan, inf, -inf}) {
		CHECK_FALSE(Vector<Real>{value}.IsEqualTo(Vector<Real>{0}));
		CHECK_FALSE(Vector<Real>{value}.IsNearZero());
		CHECK_FALSE((VectorN<Real, 1>{value}.IsEqualTo(VectorN<Real, 1>{0})));
		CHECK_FALSE((VectorN<Real, 1>{value}.IsNearZero()));
		CHECK_FALSE(Vector<Complex>{Complex(value, 0)}.IsNearZero());
	}
	for (const Real eps : {Real{-1}, nan, inf}) {
		CHECK_FALSE(Vector<Real>{0}.IsNearZero(eps));
		CHECK_FALSE((VectorN<Real, 1>{0}.IsNearZero(eps)));
		CHECK_FALSE(Vector2Cartesian(0, 0).IsEqualTo(Vector2Cartesian(0, 0), eps));
		CHECK_FALSE(Vector3Spherical(0, 0, 0).IsEqualTo(Vector3Spherical(0, 0, 0), eps));
	}
	StableNormL2Accumulator norm;
	norm.add(inf);
	norm.add(inf);
	CHECK(std::isinf(norm.value()));
	norm.add(nan);
	CHECK(std::isnan(norm.value()));
	StableNormL2Accumulator reversed;
	reversed.add(nan);
	reversed.add(inf);
	CHECK(std::isnan(reversed.value()));
}

TEST_CASE("Minkowski norms and classifications preserve scale and mixed-range cancellation", "[BaseScaleContracts][CoordinateScaleContracts]")
{
	for (const int exponent : {-DirectionScaleExponent, 0, DirectionScaleExponent}) {
		const Real scale = std::scalbn(Real{1}, exponent);
		REQUIRE(std::isfinite(scale));
		REQUIRE(scale > 0);
		for (int t = 0; t <= 4; ++t)
			for (int x = 0; x <= 4; ++x)
				for (int y = 0; y <= 4; ++y)
					for (int z = 0; z <= 2; ++z) {
						const int squared = -t * t + x * x + y * y + z * z;
						const Vector4Minkowski value{t * scale, x * scale, y * scale, z * scale};
						CHECK(value.IsTimelike() == (squared < 0));
						CHECK(value.IsSpacelike() == (squared > 0));
						CHECK(value.IsLightlike() == (squared == 0));
						const Real expected = (squared < 0 ? Real{1} : Real{-1}) * std::sqrt(Real(std::abs(squared)));
						CHECK(std::abs(value.Norm() / scale - expected) <= Precision::EqualityTolerance);
					}
	}
	for (const int exponent : {-DirectionScaleExponent, -DirectionScaleExponent * 3 / 5,
		0, DirectionScaleExponent * 3 / 5, DirectionScaleExponent}) {
		const Real scale = std::scalbn(Real{1}, exponent);
		REQUIRE(std::isfinite(scale));
		REQUIRE(scale > 0);
		const Vector4Minkowski time{Real{5} * scale, Real{3} * scale, 0, 0};
		const Vector4Minkowski space{Real{3} * scale, Real{5} * scale, 0, 0};
		const Vector4Minkowski light{Real{5} * scale, Real{3} * scale, Real{4} * scale, 0};
		CHECK(time.IsTimelike());
		CHECK_FALSE(time.IsSpacelike());
		CHECK(space.IsSpacelike());
		CHECK_FALSE(space.IsTimelike());
		CHECK(light.IsLightlike());
		CHECK(light.Norm() == 0);
		CHECK(std::abs(time.Norm() / scale - Real{4}) <= Precision::EqualityTolerance);
		CHECK(std::abs(space.Norm() / scale + Real{4}) <= Precision::EqualityTolerance);
	}
	const Real huge = HugeCoordinateScale;
	const Vector4Minkowski mixed{huge, huge, 1, 0};
	CHECK(mixed.IsSpacelike());
	CHECK(mixed.Norm() == -1);
	CHECK(ScalarProduct(mixed, mixed) == 1);
	const Vector4Minkowski smaller{huge, huge, TinyCoordinateScale, 0};
	CHECK(smaller.IsSpacelike());
	CHECK(smaller.Norm() == -TinyCoordinateScale);
	for (const Vector4Minkowski value : {
		Vector4Minkowski{huge, huge, TinyCoordinateScale, 1},
		Vector4Minkowski{huge, TinyCoordinateScale, huge, 1},
		Vector4Minkowski{huge, TinyCoordinateScale, 1, huge}}) {
		CHECK(value.IsSpacelike());
		CHECK(value.Norm() == -1);
		CHECK(ScalarProduct(value, value) == 1);
	}
	const Vector4Minkowski a{huge, huge, 3, 0}, b{huge, huge, 4, 0};
	CHECK(ScalarProduct(a, b) == 12);
	const Real epsilon = std::numeric_limits<Real>::epsilon();
	CHECK(ScalarProduct(Vector4Minkowski{1 + epsilon, 1, 0, 0},
		Vector4Minkowski{1 - epsilon, 1, 0, 0}) == epsilon * epsilon);
	const Vector4Minkowski extreme{std::numeric_limits<Real>::max(), std::numeric_limits<Real>::max(),
		std::numeric_limits<Real>::denorm_min(), 0};
	CHECK(extreme.IsSpacelike());
	CHECK(extreme.Norm() == -std::numeric_limits<Real>::denorm_min());
	const Vector4Minkowski invalid{std::numeric_limits<Real>::quiet_NaN(), 0, 0, 0};
	CHECK_FALSE(invalid.IsLightlike());
	CHECK_FALSE(invalid.IsTimelike());
	CHECK_FALSE(invalid.IsSpacelike());
	CHECK_THROWS_AS(invalid.Norm(), NumericInputError);
}
