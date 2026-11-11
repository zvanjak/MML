#include <catch2/catch_test_macros.hpp>

#ifdef MML_SINGLE_HEADER_INCLUDE
#include MML_SINGLE_HEADER_INCLUDE
#else
#include <mml/base/Quaternions.h>
#endif

#include <array>
#include <cmath>
#include <limits>

namespace MML::Tests::Base::QuaternionScaleContractTests {
namespace {
using Triple = std::array<Real, 3>;

void RequireComponents(const Quaternion& actual, const Quaternion& expected)
{
	for (int i = 0; i < 4; ++i) {
		CAPTURE(i, actual[i], expected[i]);
		REQUIRE(IsNearlyEqual(actual[i], expected[i], Precision::StrictTolerance, Precision::StrictTolerance));
	}
}

void RequireIdentity(const Quaternion& q)
{
	RequireComponents(q, Quaternion::Identity());
}

Triple RotateAxis(const Triple& v, int axis, Real angle)
{
	const Real c = std::cos(angle), s = std::sin(angle);
	if (axis == 0)
		return {v[0], c * v[1] - s * v[2], s * v[1] + c * v[2]};
	if (axis == 1)
		return {c * v[0] + s * v[2], v[1], -s * v[0] + c * v[2]};
	return {c * v[0] - s * v[1], s * v[0] + c * v[1], v[2]};
}

Triple EulerOracle(Triple v, const std::array<int, 3>& axes, const Triple& angles)
{
	// Independent elementary rotations, applied right to left for column vectors.
	for (int i = 2; i >= 0; --i)
		v = RotateAxis(v, axes[i], angles[i]);
	return v;
}

void RequireRotation(const Quaternion& q, const std::array<int, 3>& axes, const Triple& angles)
{
	const auto matrix = q.ToRotationMatrix();
	for (int column = 0; column < 3; ++column) {
		Triple basis = {0, 0, 0};
		basis[column] = 1;
		const Triple expected = EulerOracle(basis, axes, angles);
		const auto actual = q.Rotate(Vec3Cart(basis[0], basis[1], basis[2]));
		for (int row = 0; row < 3; ++row) {
			CAPTURE(column, row, angles[0], angles[1], angles[2]);
			REQUIRE(std::abs(actual[row] - expected[row]) < Precision::EqualityTolerance);
			REQUIRE(std::abs(matrix[row][column] - expected[row]) < Precision::EqualityTolerance);
		}
	}
}
}

TEST_CASE("Quaternion scaled norm normalization and inverse contracts", "[quaternion][scale-contract]")
{
	// Exercise intentional overflow at runtime rather than during constant folding.
	volatile Real runtimeMax = std::numeric_limits<Real>::max();
	const Real max = runtimeMax;
	const Real threshold = PrecisionValues<Real>::QuaternionZeroThreshold;
	const Quaternion direction(1, -2, 3, -4);
	const Quaternion expected = direction / std::sqrt(Real(30));

	for (Real scale : {threshold * Real(2), Real(1), std::sqrt(max) * Real(16), max / Real(8)}) {
		CAPTURE(scale);
		const Quaternion q = direction * scale;
		const Quaternion before = q;
		const Quaternion normalized = q.Normalized();
		REQUIRE(q == before);
		RequireComponents(normalized, expected);
		REQUIRE(normalized.IsUnit(Precision::StrictTolerance));
		Quaternion inPlace = q;
		inPlace.Normalize();
		REQUIRE(inPlace == normalized);

		const Quaternion inverse = q.Inverse();
		RequireIdentity(q * inverse);
		RequireIdentity(inverse * q);
		// Cancelling the scale checks the inverse components without underflowing the oracle.
		RequireComponents(inverse * scale, direction.Conjugate() / Real(30));
	}

	SECTION("Squared norm overflow does not imply norm overflow")
	{
		const Real scale = max / Real(4);
		const Quaternion q(scale, -scale, scale, -scale);
		REQUIRE(std::isinf(q.NormSquared()));
		REQUIRE(std::isfinite(q.Norm()));
		REQUIRE(IsNearlyEqual(q.Norm() / scale, Real(2), Precision::StrictTolerance, Precision::StrictTolerance));
	}

	SECTION("Normalization and inverse do not require a representable norm")
	{
		const Quaternion q(max, -max, max, -max);
		REQUIRE(std::isinf(q.Norm()));
		RequireComponents(q.Normalized(), Quaternion(Real(0.5), Real(-0.5), Real(0.5), Real(-0.5)));
		const Quaternion inverse = q.Inverse();
		for (int i = 0; i < 4; ++i) {
			REQUIRE(std::isfinite(inverse[i]));
			REQUIRE(inverse[i] != 0);
		}
		RequireIdentity(q * inverse);
		RequireIdentity(inverse * q);
	}

	SECTION("Original double precision report witnesses")
	{
		if constexpr (std::numeric_limits<Real>::max_exponent == 1024) {
			const Real huge = std::pow(Real(10), Real(200));
			const Quaternion q(huge, huge, 0, 0);
			REQUIRE(std::isinf(q.NormSquared())); // The old sqrt(NormSquared()) path failed here.
			const Real rootHalf = Real(1) / std::sqrt(Real(2));
			RequireComponents(q.Normalized(), Quaternion(rootHalf, rootHalf, 0, 0));
			RequireIdentity(q * q.Inverse());

			const Real small = std::pow(Real(10), Real(-7));
			REQUIRE(small > threshold);
			REQUIRE(small * small < threshold); // Old inversion compared different dimensions.
			RequireIdentity(Quaternion(small, 0, 0, 0) * Quaternion(small, 0, 0, 0).Inverse());
		}
	}
}

TEST_CASE("Quaternion norm thresholds retain their length dimension", "[quaternion][scale-contract]")
{
	const Real threshold = PrecisionValues<Real>::QuaternionZeroThreshold;
	const Real below = std::nextafter(threshold, Real(0));
	const Real above = std::nextafter(threshold, std::numeric_limits<Real>::infinity());
	for (Real norm : {Real(0), below, std::numeric_limits<Real>::denorm_min()}) {
		Quaternion q(norm, 0, 0, 0);
		const Quaternion before = q;
		REQUIRE_THROWS_AS(q.Normalized(), QuaternionError);
		REQUIRE_THROWS_AS(q.Normalize(), QuaternionError);
		REQUIRE(q == before);
		REQUIRE_THROWS_AS(q.Inverse(), QuaternionError);
	}
	for (Real norm : {threshold, above}) {
		const Quaternion q(norm, 0, 0, 0);
		CAPTURE(norm);
		RequireIdentity(q.Normalized());
		RequireIdentity(q * q.Inverse());
	}

	// Neither max-abs nor squared norm may replace the documented length comparison.
	const Real component = threshold * Real(0.75);
	const Quaternion q(component, component, 0, 0);
	const Quaternion before = q;
	const Quaternion inverse = q.Inverse();
	CAPTURE(component, inverse.w(), inverse.x());
	REQUIRE(component < threshold);
	REQUIRE(q.Norm() > threshold);
	REQUIRE(q.Normalized().IsUnit(Precision::StrictTolerance));
	REQUIRE(q == before);
	const Quaternion product = q * inverse;
	RequireIdentity(product);
	RequireIdentity(inverse * q);
	Quaternion compound = q;
	compound *= inverse;
	RequireIdentity(compound);

	const Real tiny = std::numeric_limits<Real>::denorm_min() * Real(1024);
	const Quaternion subnormal(tiny, tiny, 0, 0);
	REQUIRE(subnormal.NormSquared() == 0);
	REQUIRE(subnormal.Norm() > 0);
	REQUIRE(std::isfinite(subnormal.Norm()));
}

TEST_CASE("Quaternion direct scalar division and compound equivalence", "[quaternion][scale-contract]")
{
	const Real tiny = std::numeric_limits<Real>::denorm_min() * Real(1024);
	const Real large = std::numeric_limits<Real>::max() / Real(8);
	const Quaternion direction(1, -2, 3, -4);
	for (Real scalar : {tiny, -tiny, Real(1), Real(-2), large, -large}) {
		CAPTURE(scalar);
		const Quaternion q = direction * scalar;
		const Quaternion quotient = q / scalar;
		RequireComponents(quotient, direction);
		for (int i = 0; i < 4; ++i)
			REQUIRE(quotient[i] == q[i] / scalar);
		Quaternion compound = q;
		REQUIRE(&(compound /= scalar) == &compound);
		REQUIRE(compound == quotient);
		Quaternion product = direction;
		REQUIRE(&(product *= scalar) == &product);
		REQUIRE(product == direction * scalar);
		REQUIRE(product == scalar * direction);
	}

	// The deliberately overflowing reciprocal must be evaluated at runtime.
	volatile Real runtimeDivisor = tiny;
	REQUIRE_FALSE(std::isfinite(Real(1) / runtimeDivisor));
	RequireComponents(Quaternion(tiny, 0, 0, 0) / tiny, Quaternion::Identity());
	RequireComponents(Quaternion(0, 0, 0, 0) / tiny, Quaternion(0, 0, 0, 0));

	if constexpr (std::numeric_limits<Real>::max_exponent == 1024) {
		const Real reportDivisor = std::pow(Real(10), Real(-310));
		REQUIRE(reportDivisor > 0);
		REQUIRE_FALSE(std::isfinite(Real(1) / reportDivisor));
		RequireIdentity(Quaternion(reportDivisor, 0, 0, 0) / reportDivisor);
	}
}

TEST_CASE("Quaternion checked operations reject zero and nonfinite inputs", "[quaternion][scale-contract]")
{
	const Quaternion valid(1, -2, 3, -4);
	for (Real scalar : {Real(0), -Real(0)}) {
		Quaternion q = valid;
		REQUIRE_THROWS_AS(q / scalar, DivisionByZeroError);
		REQUIRE_THROWS_AS(q /= scalar, DivisionByZeroError);
		REQUIRE(q == valid);
	}

	for (Real invalid : {std::numeric_limits<Real>::quiet_NaN(),
		std::numeric_limits<Real>::infinity(), -std::numeric_limits<Real>::infinity()}) {
		Quaternion q = valid;
		REQUIRE_THROWS_AS(q / invalid, NumericInputError);
		REQUIRE_THROWS_AS(q /= invalid, NumericInputError);
		REQUIRE(q == valid);
		for (int component = 0; component < 4; ++component) {
			q = valid;
			q[component] = invalid;
			REQUIRE_THROWS_AS(q.Normalized(), NumericInputError);
			REQUIRE_THROWS_AS(q.Normalize(), NumericInputError);
			REQUIRE_THROWS_AS(q.Inverse(), NumericInputError);
			REQUIRE_THROWS_AS(q / Real(2), NumericInputError);
			REQUIRE_THROWS_AS(q /= Real(2), NumericInputError);
			for (int i = 0; i < 4; ++i) {
				if (i == component && std::isnan(invalid))
					REQUIRE(std::isnan(q[i]));
				else
					REQUIRE(q[i] == (i == component ? invalid : valid[i]));
			}
		}
	}
}

TEST_CASE("Quaternion Euler extraction reconstructs both conventions poles and signs", "[quaternion][euler][scale-contract]")
{
	struct Convention {
		const char* name;
		std::array<int, 3> axes;
		Quaternion (*create)(Real, Real, Real);
		Vec3Cart (Quaternion::*extract)() const;
		bool proper;
	};
	const std::array<Convention, 4> conventions = {{
		{"XYZ", {0, 1, 2}, Quaternion::FromEulerXYZ, &Quaternion::ToEulerXYZ, false},
		{"ZYX", {2, 1, 0}, Quaternion::FromEulerZYX, &Quaternion::ToEulerZYX, false},
		{"ZXZ", {2, 0, 2}, Quaternion::FromEulerZXZ, &Quaternion::ToEulerZXZ, true},
		{"ZYZ", {2, 1, 2}, Quaternion::FromEulerZYZ, &Quaternion::ToEulerZYZ, true}
	}};
	const Real pole = Precision::EqualityTolerance / Real(4);
	for (const auto& convention : conventions) {
		const Real nearPole = convention.proper ? Precision::SqrtEps : std::sqrt(Precision::EqualityTolerance);
		const std::array<Real, 7> middle = convention.proper
			? std::array<Real, 7>{0, pole, nearPole, Constants::PI / Real(3),
				Constants::PI - nearPole, Constants::PI - pole, Constants::PI}
			: std::array<Real, 7>{-Constants::PI / Real(2), -Constants::PI / Real(2) + nearPole,
				-Constants::PI / Real(3), 0, Constants::PI / Real(3),
				Constants::PI / Real(2) - nearPole, Constants::PI / Real(2)};
		for (Real first : {Real(-2), Real(0.3)})
			for (Real last : {Real(-0.7), Real(1.3)})
				for (Real second : middle) {
					const Triple angles = {first, second, last};
					const Quaternion original = convention.create(first, second, last);
					for (Real sign : {Real(1), Real(-1)}) {
						CAPTURE(convention.name, first, second, last, sign);
						const Quaternion q = original * sign;
						RequireRotation(q, convention.axes, angles);
						const Vec3Cart extracted = (q.*convention.extract)();
						for (int i = 0; i < 3; ++i)
							REQUIRE(std::isfinite(extracted[i]));
						if (convention.proper) {
							REQUIRE(extracted[1] >= 0);
							REQUIRE(extracted[1] <= Constants::PI);
						}
						const Quaternion reconstructed = convention.create(extracted[0], extracted[1], extracted[2]);
						RequireRotation(reconstructed, convention.axes, angles);
					}
				}
	}

}

TEST_CASE("Quaternion Euler extraction retains sub-epsilon middle angles",
	"[quaternion][euler][scale-contract]")
{
	for (Real sign : {Real(1), Real(-1)}) {
		const Real pitch = sign * std::numeric_limits<Real>::epsilon() / Real(8);
		const auto xyz = Quaternion::FromEulerXYZ(0, pitch, 0).ToEulerXYZ();
		const auto zyx = Quaternion::FromEulerZYX(0, pitch, 0).ToEulerZYX();
		for (Real actual : {xyz[1], zyx[1]}) {
			CHECK(actual != 0);
			CHECK(std::abs(actual / pitch - Real(1)) <= Real(8) * std::numeric_limits<Real>::epsilon());
		}
	}
}
}
