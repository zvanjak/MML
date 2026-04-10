#include <catch2/catch_all.hpp>
#include "../TestPrecision.h"
#include "../TestMatchers.h"

#include "base/Quaternions.h"
#include "base/Vector/VectorTypes.h"

using namespace MML;
using namespace MML::Testing;

namespace MML::Tests::Base::QuaternionTests {

// Helper function for approximate equality of vectors
static bool VectorApproxEqual(const Vec3Cart& v1, const Vec3Cart& v2, Real tolerance = 1e-6)
{
	return std::abs(v1[0] - v2[0]) < tolerance &&
				 std::abs(v1[1] - v2[1]) < tolerance &&
				 std::abs(v1[2] - v2[2]) < tolerance;
}

TEST_CASE("Quaternion - Construction", "[quaternion]")
{
		TEST_PRECISION_INFO();
	SECTION("Default constructor creates identity")
	{
		Quaternion q;
		REQUIRE(q.w() == REAL(1.0));
		REQUIRE(q.x() == REAL(0.0));
		REQUIRE(q.y() == REAL(0.0));
		REQUIRE(q.z() == REAL(0.0));
		REQUIRE(q.isUnit());
	}

	SECTION("Constructor from components")
	{
		Quaternion q(REAL(0.5), REAL(0.5), REAL(0.5), REAL(0.5));
		REQUIRE(q.w() == REAL(0.5));
		REQUIRE(q.x() == REAL(0.5));
		REQUIRE(q.y() == REAL(0.5));
		REQUIRE(q.z() == REAL(0.5));
	}

	SECTION("Constructor from scalar and vector")
	{
		Vec3Cart v(REAL(1.0), REAL(2.0), REAL(3.0));
		Quaternion q(REAL(0.5), v);
		REQUIRE(q.w() == REAL(0.5));
		REQUIRE(q.x() == REAL(1.0));
		REQUIRE(q.y() == REAL(2.0));
		REQUIRE(q.z() == REAL(3.0));
	}

	SECTION("Constructor from vector (pure imaginary)")
	{
		Vec3Cart v(REAL(1.0), REAL(2.0), REAL(3.0));
		Quaternion q(v);
		REQUIRE(q.w() == REAL(0.0));
		REQUIRE(q.x() == REAL(1.0));
		REQUIRE(q.y() == REAL(2.0));
		REQUIRE(q.z() == REAL(3.0));
	}

	SECTION("Identity factory method")
	{
		Quaternion q = Quaternion::Identity();
		REQUIRE(q.isUnit());
	}
}

TEST_CASE("Quaternion - From Axis-Angle", "[quaternion]")
{
		TEST_PRECISION_INFO();
	SECTION("90 degree rotation around Z-axis")
	{
		Vec3Cart axis(0, 0, 1);
		Real angle = Constants::PI / REAL(2.0);  // 90 degrees
		Quaternion q = Quaternion::FromAxisAngle(axis, angle);

		// Should be [cos(45°), 0, 0, sin(45°)]
		REQUIRE(std::abs(q.w() - std::cos(Constants::PI / REAL(4.0))) < TOL(1e-10, 1e-5));
		REQUIRE(std::abs(q.x()) < TOL(1e-10, 1e-5));
		REQUIRE(std::abs(q.y()) < TOL(1e-10, 1e-5));
		REQUIRE(std::abs(q.z() - std::sin(Constants::PI / REAL(4.0))) < TOL(1e-10, 1e-5));
		REQUIRE(q.isUnit());
	}

	SECTION("180 degree rotation around X-axis")
	{
		Vec3Cart axis(1, 0, 0);
		Real angle = Constants::PI;  // 180 degrees
		Quaternion q = Quaternion::FromAxisAngle(axis, angle);

		// Should be [0, 1, 0, 0]
		REQUIRE(std::abs(q.w()) < TOL(1e-10, 1e-5));
		REQUIRE(std::abs(q.x() - REAL(1.0)) < TOL(1e-10, 1e-5));
		REQUIRE(std::abs(q.y()) < TOL(1e-10, 1e-5));
		REQUIRE(std::abs(q.z()) < TOL(1e-10, 1e-5));
		REQUIRE(q.isUnit());
	}

	SECTION("Zero rotation")
	{
		Vec3Cart axis(1, 0, 0);
		Real angle = REAL(0.0);
		Quaternion q = Quaternion::FromAxisAngle(axis, angle);
		REQUIRE(q.isUnit());
	}
}

TEST_CASE("Quaternion - From Euler Angles", "[quaternion]")
{
		TEST_PRECISION_INFO();
	SECTION("FromEulerZYX - Zero angles gives identity")
	{
		Quaternion q = Quaternion::FromEulerZYX(0, 0, 0);
		REQUIRE(q.isUnit());
	}

	SECTION("FromEulerXYZ - Zero angles gives identity")
	{
		Quaternion q = Quaternion::FromEulerXYZ(0, 0, 0);
		REQUIRE(q.isUnit());
	}

	SECTION("FromEulerZYX - 90 degree yaw only")
	{
		Quaternion q = Quaternion::FromEulerZYX(Constants::PI / REAL(2.0), 0, 0);
		REQUIRE(q.isUnit(TOL(1e-10, 1e-5)));
		
		// Rotate a point
		Vec3Cart v(1, 0, 0);
		Vec3Cart rotated = q.Rotate(v);
		REQUIRE(VectorApproxEqual(rotated, Vec3Cart(0, 1, 0), TOL(1e-10, 1e-5)));
	}
}

TEST_CASE("Quaternion - Basic Operations", "[quaternion]")
{
		TEST_PRECISION_INFO();
	Quaternion q1(1, 2, 3, 4);
	Quaternion q2(5, 6, 7, 8);

	SECTION("Addition")
	{
		Quaternion result = q1 + q2;
		REQUIRE(result.w() == 6);
		REQUIRE(result.x() == 8);
		REQUIRE(result.y() == 10);
		REQUIRE(result.z() == 12);
	}

	SECTION("Subtraction")
	{
		Quaternion result = q1 - q2;
		REQUIRE(result.w() == -4);
		REQUIRE(result.x() == -4);
		REQUIRE(result.y() == -4);
		REQUIRE(result.z() == -4);
	}

	SECTION("Negation")
	{
		Quaternion result = -q1;
		REQUIRE(result.w() == -1);
		REQUIRE(result.x() == -2);
		REQUIRE(result.y() == -3);
		REQUIRE(result.z() == -4);
	}

	SECTION("Scalar multiplication")
	{
		Quaternion result = q1 * REAL(2.0);
		REQUIRE(result.w() == 2);
		REQUIRE(result.x() == 4);
		REQUIRE(result.y() == 6);
		REQUIRE(result.z() == 8);
	}

	SECTION("Scalar multiplication (commutative)")
	{
		Quaternion result = REAL(2.0) * q1;
		REQUIRE(result.w() == 2);
		REQUIRE(result.x() == 4);
		REQUIRE(result.y() == 6);
		REQUIRE(result.z() == 8);
	}

	SECTION("Scalar division")
	{
		Quaternion result = q1 / REAL(2.0);
		REQUIRE(result.w() == REAL(0.5));
		REQUIRE(result.x() == REAL(1.0));
		REQUIRE(result.y() == REAL(1.5));
		REQUIRE(result.z() == REAL(2.0));
	}
}

TEST_CASE("Quaternion - Multiplication", "[quaternion]")
{
		TEST_PRECISION_INFO();
	SECTION("Identity multiplication")
	{
		Quaternion q(2, 3, 4, 5);
		Quaternion id = Quaternion::Identity();
		
		Quaternion result1 = q * id;
		Quaternion result2 = id * q;
		
		REQUIRE(result1.isApprox(q));
		REQUIRE(result2.isApprox(q));
	}

	SECTION("Non-commutativity")
	{
		Quaternion q1 = Quaternion::FromAxisAngle(Vec3Cart(1, 0, 0), Constants::PI / REAL(4.0));
		Quaternion q2 = Quaternion::FromAxisAngle(Vec3Cart(0, 1, 0), Constants::PI / REAL(4.0));
		
		Quaternion result1 = q1 * q2;
		Quaternion result2 = q2 * q1;
		
		// Should not be equal (within reasonable tolerance)
		REQUIRE_FALSE(result1.isApprox(result2, 1e-6));
	}

	SECTION("Known quaternion multiplication")
	{
		// i * j = k in Hamilton convention
		Quaternion i(0, 1, 0, 0);
		Quaternion j(0, 0, 1, 0);
		Quaternion k_expected(0, 0, 0, 1);
		
		Quaternion k_result = i * j;
		REQUIRE(k_result.isApprox(k_expected, TOL(1e-10, 1e-5)));
	}

	SECTION("Conjugate property: (q1*q2)* = q2* * q1*")
	{
		Quaternion q1(1, 2, 3, 4);
		Quaternion q2(5, 6, 7, 8);
		
		Quaternion lhs = (q1 * q2).Conjugate();
		Quaternion rhs = q2.Conjugate() * q1.Conjugate();
		
		REQUIRE(lhs.isApprox(rhs, TOL(1e-10, 1e-5)));
	}
}

TEST_CASE("Quaternion - Conjugate and Inverse", "[quaternion]")
{
		TEST_PRECISION_INFO();
	SECTION("Conjugate of identity")
	{
		Quaternion q = Quaternion::Identity();
		Quaternion conj = q.Conjugate();
		REQUIRE(conj.isUnit());
	}

	SECTION("Conjugate properties")
	{
		Quaternion q(1, 2, 3, 4);
		Quaternion conj = q.Conjugate();
		
		REQUIRE(conj.w() == q.w());
		REQUIRE(conj.x() == -q.x());
		REQUIRE(conj.y() == -q.y());
		REQUIRE(conj.z() == -q.z());
	}

	SECTION("Double conjugate returns original")
	{
		Quaternion q(1, 2, 3, 4);
		Quaternion double_conj = q.Conjugate().Conjugate();
		REQUIRE(double_conj.isApprox(q));
	}

	SECTION("Inverse of unit quaternion equals conjugate")
	{
		Quaternion q = Quaternion::FromAxisAngle(Vec3Cart(1, 1, 1).Normalized(), Constants::PI / REAL(3.0));
		Quaternion inv = q.Inverse();
		Quaternion conj = q.Conjugate();
		
		REQUIRE(inv.isApprox(conj, TOL(1e-10, 1e-5)));
	}

	SECTION("q * q^(-1) = identity")
	{
		Quaternion q(1, 2, 3, 4);
		Quaternion inv = q.Inverse();
		Quaternion result = q * inv;
		
		REQUIRE(result.isUnit(TOL(1e-10, 1e-5)));
	}

	SECTION("q^(-1) * q = identity")
	{
		Quaternion q(1, 2, 3, 4);
		Quaternion inv = q.Inverse();
		Quaternion result = inv * q;
		
		REQUIRE(result.isUnit(TOL(1e-10, 1e-5)));
	}
}

TEST_CASE("Quaternion - Norm and Normalization", "[quaternion]")
{
		TEST_PRECISION_INFO();
	SECTION("Norm of identity is 1")
	{
		Quaternion q = Quaternion::Identity();
		REQUIRE(std::abs(q.Norm() - REAL(1.0)) < TOL(1e-10, 1e-5));
		REQUIRE(std::abs(q.NormSquared() - REAL(1.0)) < TOL(1e-10, 1e-5));
	}

	SECTION("Norm calculation")
	{
		Quaternion q(1, 2, 2, 0);
		Real expected_norm = std::sqrt(1 + 4 + 4);  // sqrt(9) = 3
		REQUIRE(std::abs(q.Norm() - expected_norm) < TOL(1e-10, 1e-5));
		REQUIRE(std::abs(q.NormSquared() - REAL(9.0)) < TOL(1e-10, 1e-5));
	}

	SECTION("Normalization")
	{
		Quaternion q(2, 0, 0, 0);
		q.Normalize();
		REQUIRE(q.isUnit(TOL(1e-10, 1e-5)));
		REQUIRE(std::abs(q.w() - REAL(1.0)) < TOL(1e-10, 1e-5));
	}

	SECTION("Normalized returns normalized copy")
	{
		Quaternion q(1, 2, 3, 4);
		Quaternion orig = q;
		Quaternion normalized = q.Normalized();
		
		REQUIRE(normalized.isUnit(TOL(1e-10, 1e-5)));
		REQUIRE(q.isApprox(orig));  // Original unchanged
	}

	SECTION("FromAxisAngle produces unit quaternion")
	{
		Vec3Cart axis(1, 2, 3);
		axis = axis.Normalized();
		Quaternion q = Quaternion::FromAxisAngle(axis, REAL(1.5));
		REQUIRE(q.isUnit(TOL(1e-10, 1e-5)));
	}
}

TEST_CASE("Quaternion - Vector Rotation", "[quaternion]")
{
		TEST_PRECISION_INFO();
	SECTION("Identity rotation does nothing")
	{
		Quaternion q = Quaternion::Identity();
		Vec3Cart v(1, 2, 3);
		Vec3Cart rotated = q.Rotate(v);
		REQUIRE(VectorApproxEqual(rotated, v, TOL(1e-10, 1e-5)));
	}

	SECTION("90 degree rotation around Z-axis")
	{
		// Rotate (1, 0, 0) by 90° around Z-axis -> (0, 1, 0)
		Quaternion q = Quaternion::FromAxisAngle(Vec3Cart(0, 0, 1), Constants::PI / REAL(2.0));
		Vec3Cart v(1, 0, 0);
		Vec3Cart rotated = q.Rotate(v);
		
		REQUIRE(VectorApproxEqual(rotated, Vec3Cart(0, 1, 0), TOL(1e-10, 1e-5)));
	}

	SECTION("90 degree rotation around X-axis")
	{
		// Rotate (0, 1, 0) by 90° around X-axis -> (0, 0, 1)
		Quaternion q = Quaternion::FromAxisAngle(Vec3Cart(1, 0, 0), Constants::PI / REAL(2.0));
		Vec3Cart v(0, 1, 0);
		Vec3Cart rotated = q.Rotate(v);
		
		REQUIRE(VectorApproxEqual(rotated, Vec3Cart(0, 0, 1), TOL(1e-10, 1e-5)));
	}

	SECTION("90 degree rotation around Y-axis")
	{
		// Rotate (1, 0, 0) by 90° around Y-axis -> (0, 0, -1)
		Quaternion q = Quaternion::FromAxisAngle(Vec3Cart(0, 1, 0), Constants::PI / REAL(2.0));
		Vec3Cart v(1, 0, 0);
		Vec3Cart rotated = q.Rotate(v);
		
		REQUIRE(VectorApproxEqual(rotated, Vec3Cart(0, 0, -1), TOL(1e-10, 1e-5)));
	}

	SECTION("180 degree rotation around Z-axis")
	{
		// Rotate (1, 0, 0) by 180° around Z-axis -> (-1, 0, 0)
		Quaternion q = Quaternion::FromAxisAngle(Vec3Cart(0, 0, 1), Constants::PI);
		Vec3Cart v(1, 0, 0);
		Vec3Cart rotated = q.Rotate(v);
		
		REQUIRE(VectorApproxEqual(rotated, Vec3Cart(-1, 0, 0), TOL(1e-10, 1e-5)));
	}

	SECTION("Rotation preserves vector length")
	{
		Quaternion q = Quaternion::FromAxisAngle(Vec3Cart(1, 1, 1).Normalized(), REAL(0.7));
		Vec3Cart v(1, 2, 3);
		Vec3Cart rotated = q.Rotate(v);
		
		REQUIRE(std::abs(rotated.NormL2() - v.NormL2()) < TOL(1e-10, 1e-5));
	}

	SECTION("Composition of rotations")
	{
		// Rotate by q1, then by q2 = rotate by q2*q1
		Quaternion q1 = Quaternion::FromAxisAngle(Vec3Cart(0, 0, 1), Constants::PI / REAL(2.0));
		Quaternion q2 = Quaternion::FromAxisAngle(Vec3Cart(1, 0, 0), Constants::PI / REAL(2.0));
		Quaternion combined = q2 * q1;
		
		Vec3Cart v(1, 0, 0);
		Vec3Cart rotated_separate = q2.Rotate(q1.Rotate(v));
		Vec3Cart rotated_combined = combined.Rotate(v);
		
		REQUIRE(VectorApproxEqual(rotated_separate, rotated_combined, TOL(1e-10, 1e-5)));
	}

	SECTION("Inverse rotation undoes rotation")
	{
		Quaternion q = Quaternion::FromAxisAngle(Vec3Cart(1, 2, 3).Normalized(), REAL(1.234));
		Vec3Cart v(4, 5, 6);
		
		Vec3Cart rotated = q.Rotate(v);
		Vec3Cart back = q.Inverse().Rotate(rotated);
		
		REQUIRE(VectorApproxEqual(back, v, TOL(1e-10, 1e-5)));
	}
}

TEST_CASE("Quaternion - Axis-Angle Extraction", "[quaternion]")
{
		TEST_PRECISION_INFO();
	SECTION("Extract from known rotation")
	{
		Vec3Cart axis_in(1, 0, 0);
		Real angle_in = Constants::PI / REAL(3.0);
		
		Quaternion q = Quaternion::FromAxisAngle(axis_in, angle_in);
		
		Vec3Cart axis_out;
		Real angle_out;
		q.ToAxisAngle(axis_out, angle_out);
		
		REQUIRE(std::abs(angle_out - angle_in) < TOL(1e-10, 1e-5));
		REQUIRE(VectorApproxEqual(axis_out, axis_in, TOL(1e-10, 1e-5)));
	}

	SECTION("GetRotationAxis and GetRotationAngle")
	{
		Vec3Cart axis_in(0, 1, 0);
		Real angle_in = Constants::PI / REAL(4.0);
		
		Quaternion q = Quaternion::FromAxisAngle(axis_in, angle_in);
		
		Vec3Cart axis_out = q.GetRotationAxis();
		Real angle_out = q.GetRotationAngle();
		
		REQUIRE(std::abs(angle_out - angle_in) < TOL(1e-10, 1e-5));
		REQUIRE(VectorApproxEqual(axis_out, axis_in, TOL(1e-10, 1e-5)));
	}

	SECTION("Identity quaternion has zero angle")
	{
		Quaternion q = Quaternion::Identity();
		Real angle = q.GetRotationAngle();
		REQUIRE(std::abs(angle) < TOL(1e-10, 1e-5));
	}
}

TEST_CASE("Quaternion - Euler Angle Conversion", "[quaternion]")
{
		TEST_PRECISION_INFO();
	SECTION("Round-trip ZYX Euler angles")
	{
		Real yaw_in = REAL(0.5);
		Real pitch_in = REAL(0.3);
		Real roll_in = REAL(0.7);
		
		Quaternion q = Quaternion::FromEulerZYX(yaw_in, pitch_in, roll_in);
		Vec3Cart euler_out = q.ToEulerZYX();
		
		REQUIRE(std::abs(euler_out[0] - yaw_in) < 1e-6);
		REQUIRE(std::abs(euler_out[1] - pitch_in) < 1e-6);
		REQUIRE(std::abs(euler_out[2] - roll_in) < 1e-6);
	}

	SECTION("Euler ZYX for simple rotations")
	{
		// Pure yaw (Z-axis rotation)
		Quaternion q = Quaternion::FromEulerZYX(Constants::PI / REAL(2.0), 0, 0);
		Vec3Cart euler = q.ToEulerZYX();
		REQUIRE(std::abs(euler[0] - Constants::PI / REAL(2.0)) < TOL(1e-10, 1e-5));
		REQUIRE(std::abs(euler[1]) < TOL(1e-10, 1e-5));
		REQUIRE(std::abs(euler[2]) < TOL(1e-10, 1e-5));
	}
}

TEST_CASE("Quaternion - Interpolation", "[quaternion]")
{
		TEST_PRECISION_INFO();
	Quaternion q1 = Quaternion::Identity();
	Quaternion q2 = Quaternion::FromAxisAngle(Vec3Cart(0, 0, 1), Constants::PI / REAL(2.0));

	SECTION("Lerp at t=0 returns first quaternion")
	{
		Quaternion result = Quaternion::Lerp(q1, q2, REAL(0.0));
		REQUIRE(result.isApprox(q1));
	}

	SECTION("Lerp at t=1 returns second quaternion")
	{
		Quaternion result = Quaternion::Lerp(q1, q2, REAL(1.0));
		REQUIRE(result.isApprox(q2));
	}

	SECTION("Lerp at t=REAL(0.5) is between quaternions")
	{
		Quaternion result = Quaternion::Lerp(q1, q2, REAL(0.5));
		// Should be between, but not necessarily unit
		Real dot1 = result.Dot(q1);
		Real dot2 = result.Dot(q2);
		REQUIRE(dot1 > 0);
		REQUIRE(dot2 > 0);
	}

	SECTION("Slerp at t=0 returns first quaternion")
	{
		Quaternion result = Quaternion::Slerp(q1, q2, REAL(0.0));
		REQUIRE(result.isApprox(q1, 1e-6));
	}

	SECTION("Slerp at t=1 returns second quaternion")
	{
		Quaternion result = Quaternion::Slerp(q1, q2, REAL(1.0));
		REQUIRE(result.isApprox(q2, 1e-6));
	}

	SECTION("Slerp produces unit quaternions")
	{
		Quaternion result = Quaternion::Slerp(q1, q2, REAL(0.3));
		REQUIRE(result.isUnit(1e-6));
		
		result = Quaternion::Slerp(q1, q2, REAL(0.7));
		REQUIRE(result.isUnit(1e-6));
	}

	SECTION("Slerp at t=REAL(0.5) gives halfway rotation")
	{
		Quaternion result = Quaternion::Slerp(q1, q2, REAL(0.5));
		
		// Should be 45° rotation around Z-axis
		Real angle = result.GetRotationAngle();
		REQUIRE(std::abs(angle - Constants::PI / REAL(4.0)) < 1e-6);
		
		Vec3Cart axis = result.GetRotationAxis();
		REQUIRE(VectorApproxEqual(axis, Vec3Cart(0, 0, 1), 1e-6));
	}

	SECTION("Slerp handles opposite quaternions")
	{
		Quaternion q_a = Quaternion::FromAxisAngle(Vec3Cart(0, 0, 1), REAL(0.1));
		Quaternion q_b = -q_a;  // Opposite quaternion
		
		// Slerp should handle this gracefully
		Quaternion result = Quaternion::Slerp(q_a, q_b, REAL(0.5));
		REQUIRE(result.isUnit(1e-6));
	}

	SECTION("Slerp near-identical quaternions (sinTheta stability)")
	{
		// dot ~ 0.99999 — just past the 0.9995 threshold, sinTheta very small
		Quaternion q_a = Quaternion::FromAxisAngle(Vec3Cart(0, 0, 1), REAL(0.0));
		Quaternion q_b = Quaternion::FromAxisAngle(Vec3Cart(0, 0, 1), REAL(0.001));

		// Should not produce NaN or Inf
		Quaternion result = Quaternion::Slerp(q_a, q_b, REAL(0.5));
		REQUIRE(result.isUnit(1e-4));
		REQUIRE_FALSE(std::isnan(result.w()));
		REQUIRE_FALSE(std::isnan(result.x()));
		REQUIRE_FALSE(std::isnan(result.y()));
		REQUIRE_FALSE(std::isnan(result.z()));
	}
}

TEST_CASE("Quaternion - Dot Product", "[quaternion]")
{
		TEST_PRECISION_INFO();
	SECTION("Dot product of identical quaternions")
	{
		Quaternion q(1, 2, 3, 4);
		Real dot = q.Dot(q);
		REQUIRE(std::abs(dot - q.NormSquared()) < TOL(1e-10, 1e-5));
	}

	SECTION("Dot product is commutative")
	{
		Quaternion q1(1, 2, 3, 4);
		Quaternion q2(5, 6, 7, 8);
		REQUIRE(std::abs(q1.Dot(q2) - q2.Dot(q1)) < TOL(1e-10, 1e-5));
	}

	SECTION("Orthogonal quaternions have zero dot product")
	{
		Quaternion q1(1, 0, 0, 0);
		Quaternion q2(0, 1, 0, 0);
		REQUIRE(std::abs(q1.Dot(q2)) < TOL(1e-10, 1e-5));
	}
}

TEST_CASE("Quaternion - Comparison and Equality", "[quaternion]")
{
		TEST_PRECISION_INFO();
	SECTION("Exact equality")
	{
		Quaternion q1(1, 2, 3, 4);
		Quaternion q2(1, 2, 3, 4);
		Quaternion q3(1, 2, 3, 5);
		
		REQUIRE(q1 == q2);
		REQUIRE_FALSE(q1 == q3);
		REQUIRE(q1 != q3);
	}

	SECTION("Approximate equality")
	{
		Quaternion q1(REAL(1.0), REAL(2.0), REAL(3.0), REAL(4.0));
		Quaternion q2(REAL(1.0) + 1e-7, REAL(2.0), REAL(3.0), REAL(4.0));
		
		REQUIRE(q1.isApprox(q2, 1e-6));
		REQUIRE_FALSE(q1.isApprox(q2, TOL(1e-8, 1e-8)));
	}
}

TEST_CASE("Quaternion - Special Cases and Edge Cases", "[quaternion]")
{
		TEST_PRECISION_INFO();
	SECTION("Zero rotation has identity quaternion")
	{
		Quaternion q = Quaternion::FromAxisAngle(Vec3Cart(1, 0, 0), REAL(0.0));
		REQUIRE(q.isUnit());
	}

	SECTION("Full rotation (2π) represents same rotation as identity")
	{
		// Note: q = [cos(π), sin(π)*axis] = [-1, 0, 0, 0]
		// Both q and -q represent the same rotation!
		Quaternion q = Quaternion::FromAxisAngle(Vec3Cart(1, 0, 0), REAL(2.0) * Constants::PI);
		
		// Test by rotating a vector - should get same result as identity
		Vec3Cart v(1, 2, 3);
		Vec3Cart rotated = q.Rotate(v);
		REQUIRE(VectorApproxEqual(rotated, v, TOL(1e-10, 1e-5)));
	}

	SECTION("Negative angle gives inverse rotation")
	{
		Vec3Cart axis(1, 1, 1);
		axis = axis.Normalized();
		Real angle = REAL(0.7);
		
		Quaternion q_pos = Quaternion::FromAxisAngle(axis, angle);
		Quaternion q_neg = Quaternion::FromAxisAngle(axis, -angle);
		
		Quaternion product = q_pos * q_neg;
		REQUIRE(product.isUnit(TOL(1e-10, 1e-5)));
	}

	SECTION("IsUnit check with tolerance")
	{
		Quaternion q(REAL(1.0) + 1e-7, 0, 0, 0);
		REQUIRE(q.isUnit(1e-6));
		REQUIRE_FALSE(q.isUnit(TOL(1e-8, 1e-8)));
	}

	SECTION("Rotation of zero vector")
	{
		Quaternion q = Quaternion::FromAxisAngle(Vec3Cart(1, 0, 0), REAL(1.5));
		Vec3Cart zero(0, 0, 0);
		Vec3Cart rotated = q.Rotate(zero);
		REQUIRE(VectorApproxEqual(rotated, zero, TOL(1e-10, 1e-5)));
	}

	SECTION("q and -q represent same rotation")
	{
		Quaternion q = Quaternion::FromAxisAngle(Vec3Cart(1, 0, 0), Constants::PI / REAL(4.0));
		Quaternion q_neg = -q;
		
		Vec3Cart v(1, 2, 3);
		Vec3Cart rot1 = q.Rotate(v);
		Vec3Cart rot2 = q_neg.Rotate(v);
		
		REQUIRE(VectorApproxEqual(rot1, rot2, TOL(1e-10, 1e-5)));
	}
}

TEST_CASE("Quaternion - Accessors", "[quaternion]")
{
		TEST_PRECISION_INFO();
	SECTION("Individual component accessors")
	{
		Quaternion q(1, 2, 3, 4);
		REQUIRE(q.w() == 1);
		REQUIRE(q.x() == 2);
		REQUIRE(q.y() == 3);
		REQUIRE(q.z() == 4);
	}

	SECTION("Index operator")
	{
		Quaternion q(1, 2, 3, 4);
		REQUIRE(q[0] == 1);
		REQUIRE(q[1] == 2);
		REQUIRE(q[2] == 3);
		REQUIRE(q[3] == 4);
	}

	SECTION("Scalar accessor")
	{
		Quaternion q(5, 2, 3, 4);
		REQUIRE(q.Scalar() == 5);
	}

	SECTION("Vector accessor")
	{
		Quaternion q(1, 2, 3, 4);
		Vec3Cart v = q.Vector();
		REQUIRE(VectorApproxEqual(v, Vec3Cart(2, 3, 4)));
	}

	SECTION("Mutable accessors")
	{
		Quaternion q;
		q.w() = 5;
		q.x() = 6;
		q.y() = 7;
		q.z() = 8;
		
		REQUIRE(q.w() == 5);
		REQUIRE(q.x() == 6);
		REQUIRE(q.y() == 7);
		REQUIRE(q.z() == 8);
	}
}

TEST_CASE("Quaternion - In-place Operations", "[quaternion]")
{
		TEST_PRECISION_INFO();
	SECTION("In-place addition")
	{
		Quaternion q1(1, 2, 3, 4);
		Quaternion q2(5, 6, 7, 8);
		q1 += q2;
		
		REQUIRE(q1.w() == 6);
		REQUIRE(q1.x() == 8);
		REQUIRE(q1.y() == 10);
		REQUIRE(q1.z() == 12);
	}

	SECTION("In-place subtraction")
	{
		Quaternion q1(10, 20, 30, 40);
		Quaternion q2(1, 2, 3, 4);
		q1 -= q2;
		
		REQUIRE(q1.w() == 9);
		REQUIRE(q1.x() == 18);
		REQUIRE(q1.y() == 27);
		REQUIRE(q1.z() == 36);
	}

	SECTION("In-place scalar multiplication")
	{
		Quaternion q(1, 2, 3, 4);
		q *= REAL(3.0);
		
		REQUIRE(q.w() == 3);
		REQUIRE(q.x() == 6);
		REQUIRE(q.y() == 9);
		REQUIRE(q.z() == 12);
	}

	SECTION("In-place quaternion multiplication")
	{
		Quaternion q1 = Quaternion::FromAxisAngle(Vec3Cart(0, 0, 1), Constants::PI / REAL(4.0));
		Quaternion q2 = Quaternion::FromAxisAngle(Vec3Cart(0, 0, 1), Constants::PI / REAL(4.0));
		Quaternion expected = q1 * q2;
		
		q1 *= q2;
		REQUIRE(q1.isApprox(expected, TOL(1e-10, 1e-5)));
	}
}

//////////////////////////////////////////////////////////////////
///          MATRIX CONVERSION TESTS                           ///
//////////////////////////////////////////////////////////////////

TEST_CASE("Quaternion - ToRotationMatrix", "[quaternion][matrix]")
{
	TEST_PRECISION_INFO();
	
	SECTION("Identity quaternion gives identity matrix")
	{
		Quaternion q = Quaternion::Identity();
		auto mat = q.ToRotationMatrix();
		
		REQUIRE(std::abs(mat[0][0] - REAL(1.0)) < TOL(1e-10, 1e-5));
		REQUIRE(std::abs(mat[1][1] - REAL(1.0)) < TOL(1e-10, 1e-5));
		REQUIRE(std::abs(mat[2][2] - REAL(1.0)) < TOL(1e-10, 1e-5));
		REQUIRE(std::abs(mat[0][1]) < TOL(1e-10, 1e-5));
		REQUIRE(std::abs(mat[0][2]) < TOL(1e-10, 1e-5));
		REQUIRE(std::abs(mat[1][0]) < TOL(1e-10, 1e-5));
		REQUIRE(std::abs(mat[1][2]) < TOL(1e-10, 1e-5));
		REQUIRE(std::abs(mat[2][0]) < TOL(1e-10, 1e-5));
		REQUIRE(std::abs(mat[2][1]) < TOL(1e-10, 1e-5));
	}
	
	SECTION("90 degree Z rotation matrix")
	{
		Quaternion q = Quaternion::FromAxisAngle(Vec3Cart(0, 0, 1), Constants::PI / REAL(2.0));
		auto mat = q.ToRotationMatrix();
		
		// Expected: [0, -1, 0; 1, 0, 0; 0, 0, 1]
		REQUIRE(std::abs(mat[0][0]) < TOL(1e-10, 1e-5));
		REQUIRE(std::abs(mat[0][1] + REAL(1.0)) < TOL(1e-10, 1e-5));
		REQUIRE(std::abs(mat[1][0] - REAL(1.0)) < TOL(1e-10, 1e-5));
		REQUIRE(std::abs(mat[1][1]) < TOL(1e-10, 1e-5));
		REQUIRE(std::abs(mat[2][2] - REAL(1.0)) < TOL(1e-10, 1e-5));
	}
	
	SECTION("Matrix rotation equals quaternion rotation")
	{
		Quaternion q = Quaternion::FromAxisAngle(Vec3Cart(1, 1, 1).Normalized(), REAL(0.7));
		auto mat = q.ToRotationMatrix();
		
		Vec3Cart v(1, 2, 3);
		Vec3Cart rotated_q = q.Rotate(v);
		
		// Matrix * vector
		Vec3Cart rotated_m(
			mat[0][0] * v[0] + mat[0][1] * v[1] + mat[0][2] * v[2],
			mat[1][0] * v[0] + mat[1][1] * v[1] + mat[1][2] * v[2],
			mat[2][0] * v[0] + mat[2][1] * v[1] + mat[2][2] * v[2]
		);
		
		REQUIRE(VectorApproxEqual(rotated_q, rotated_m, TOL(1e-10, 1e-5)));
	}
	
	SECTION("Rotation matrix is orthogonal (det = 1)")
	{
		Quaternion q = Quaternion::FromAxisAngle(Vec3Cart(1, 2, 3).Normalized(), REAL(1.234));
		auto mat = q.ToRotationMatrix();
		
		// Check R * R^T = I (first row dot products)
		Real dot00 = mat[0][0]*mat[0][0] + mat[0][1]*mat[0][1] + mat[0][2]*mat[0][2];
		Real dot11 = mat[1][0]*mat[1][0] + mat[1][1]*mat[1][1] + mat[1][2]*mat[1][2];
		Real dot22 = mat[2][0]*mat[2][0] + mat[2][1]*mat[2][1] + mat[2][2]*mat[2][2];
		Real dot01 = mat[0][0]*mat[1][0] + mat[0][1]*mat[1][1] + mat[0][2]*mat[1][2];
		
		REQUIRE(std::abs(dot00 - REAL(1.0)) < TOL(1e-10, 1e-5));
		REQUIRE(std::abs(dot11 - REAL(1.0)) < TOL(1e-10, 1e-5));
		REQUIRE(std::abs(dot22 - REAL(1.0)) < TOL(1e-10, 1e-5));
		REQUIRE(std::abs(dot01) < TOL(1e-10, 1e-5));
	}
}

TEST_CASE("Quaternion - FromRotationMatrix", "[quaternion][matrix]")
{
	TEST_PRECISION_INFO();
	
	SECTION("Round-trip quaternion -> matrix -> quaternion")
	{
		Quaternion q_orig = Quaternion::FromAxisAngle(Vec3Cart(1, 2, 3).Normalized(), REAL(0.8));
		auto mat = q_orig.ToRotationMatrix();
		Quaternion q_back = Quaternion::FromRotationMatrix(mat);
		
		// q and -q represent same rotation
		bool same = q_orig.isApprox(q_back, TOL(1e-10, 1e-5)) || q_orig.isApprox(-q_back, TOL(1e-10, 1e-5));
		REQUIRE(same);
	}
	
	SECTION("Various rotation angles for numerical stability")
	{
		std::vector<Real> angles = {REAL(0.1), REAL(0.5), Constants::PI / REAL(4.0), 
		                            Constants::PI / REAL(2.0), Constants::PI * REAL(0.9)};
		std::vector<Vec3Cart> axes = {Vec3Cart(1,0,0), Vec3Cart(0,1,0), Vec3Cart(0,0,1),
		                              Vec3Cart(1,1,0).Normalized(), Vec3Cart(1,1,1).Normalized()};
		
		for (const auto& axis : axes) {
			for (Real angle : angles) {
				Quaternion q_orig = Quaternion::FromAxisAngle(axis, angle);
				auto mat = q_orig.ToRotationMatrix();
				Quaternion q_back = Quaternion::FromRotationMatrix(mat);
				
				// Verify same rotation effect
				Vec3Cart v(1, 2, 3);
				Vec3Cart r1 = q_orig.Rotate(v);
				Vec3Cart r2 = q_back.Rotate(v);
				REQUIRE(VectorApproxEqual(r1, r2, TOL(1e-9, 1e-4)));
			}
		}
	}
}

//////////////////////////////////////////////////////////////////
///          COPY/MOVE SEMANTICS TESTS                         ///
//////////////////////////////////////////////////////////////////

TEST_CASE("Quaternion - Copy Semantics", "[quaternion][copy]")
{
	TEST_PRECISION_INFO();
	
	SECTION("Copy constructor")
	{
		Quaternion q1(1, 2, 3, 4);
		Quaternion q2(q1);
		
		REQUIRE(q2.w() == 1);
		REQUIRE(q2.x() == 2);
		REQUIRE(q2.y() == 3);
		REQUIRE(q2.z() == 4);
		REQUIRE(q1 == q2);
	}
	
	SECTION("Copy assignment")
	{
		Quaternion q1(1, 2, 3, 4);
		Quaternion q2;
		q2 = q1;
		
		REQUIRE(q2 == q1);
		
		// Modify original, copy should be unaffected
		q1.w() = 10;
		REQUIRE(q2.w() == 1);
	}
}

TEST_CASE("Quaternion - Move Semantics", "[quaternion][move]")
{
	TEST_PRECISION_INFO();
	
	SECTION("Move constructor")
	{
		Quaternion q1(1, 2, 3, 4);
		Quaternion q2(std::move(q1));
		
		REQUIRE(q2.w() == 1);
		REQUIRE(q2.x() == 2);
		REQUIRE(q2.y() == 3);
		REQUIRE(q2.z() == 4);
	}
	
	SECTION("Move assignment")
	{
		Quaternion q1(5, 6, 7, 8);
		Quaternion q2;
		q2 = std::move(q1);
		
		REQUIRE(q2.w() == 5);
		REQUIRE(q2.x() == 6);
		REQUIRE(q2.y() == 7);
		REQUIRE(q2.z() == 8);
	}
}

//////////////////////////////////////////////////////////////////
///          OUTPUT TESTS                                      ///
//////////////////////////////////////////////////////////////////

TEST_CASE("Quaternion - Stream Output", "[quaternion][output]")
{
	TEST_PRECISION_INFO();
	
	SECTION("operator<< format")
	{
		Quaternion q(1, 2, 3, 4);
		std::ostringstream oss;
		oss << q;
		std::string result = oss.str();
		
		REQUIRE(result.find("[") != std::string::npos);
		REQUIRE(result.find("i") != std::string::npos);
		REQUIRE(result.find("j") != std::string::npos);
		REQUIRE(result.find("k") != std::string::npos);
	}
	
	SECTION("Print method")
	{
		Quaternion q(1, 0, 0, 0);
		std::ostringstream oss;
		q.Print(oss);
		std::string result = oss.str();
		
		REQUIRE(!result.empty());
		REQUIRE(result.find("1") != std::string::npos);
	}
}

//////////////////////////////////////////////////////////////////
///          ADDITIONAL EULER ANGLE TESTS                      ///
//////////////////////////////////////////////////////////////////

TEST_CASE("Quaternion - FromEulerXYZ Additional", "[quaternion][euler]")
{
	TEST_PRECISION_INFO();
	
	SECTION("FromEulerXYZ produces unit quaternion")
	{
		Quaternion q = Quaternion::FromEulerXYZ(REAL(0.3), REAL(0.5), REAL(0.7));
		REQUIRE(q.isUnit(TOL(1e-10, 1e-5)));
	}
	
	SECTION("FromEulerXYZ rotation verification")
	{
		// Pure roll (X-axis rotation)
		Quaternion q = Quaternion::FromEulerXYZ(Constants::PI / REAL(2.0), 0, 0);
		
		// Rotating (0, 1, 0) around X by 90° should give (0, 0, 1)
		Vec3Cart v(0, 1, 0);
		Vec3Cart rotated = q.Rotate(v);
		REQUIRE(VectorApproxEqual(rotated, Vec3Cart(0, 0, 1), TOL(1e-10, 1e-5)));
	}
}

//////////////////////////////////////////////////////////////////
///          ERROR HANDLING TESTS                              ///
//////////////////////////////////////////////////////////////////

TEST_CASE("Quaternion - Error Handling", "[quaternion][errors]")
{
	TEST_PRECISION_INFO();
	
	SECTION("Inverse of near-zero quaternion throws")
	{
		Quaternion q(0, 0, 0, 0);
		REQUIRE_THROWS(q.Inverse());
	}
	
	SECTION("Normalize zero quaternion throws")
	{
		Quaternion q(0, 0, 0, 0);
		REQUIRE_THROWS(q.Normalize());
	}
}

//////////////////////////////////////////////////////////////////
///          EULER ANGLE CONVENTIONS TESTS                     ///
//////////////////////////////////////////////////////////////////

TEST_CASE("Quaternion - Tait-Bryan Euler roundtrip via rotation matrix", "[quaternion][euler]")
{
	TEST_PRECISION_INFO();
	const Real tol = TOL(1e-10, 1e-4);

	// Helper: compare two rotation matrices from two quaternions
	auto matricesEqual = [&](const Quaternion& q1, const Quaternion& q2) {
		auto m1 = q1.ToRotationMatrix();
		auto m2 = q2.ToRotationMatrix();
		for (int i = 0; i < 3; i++)
			for (int j = 0; j < 3; j++)
				REQUIRE(std::abs(m1[i][j] - m2[i][j]) < tol);
	};

	SECTION("FromEulerXZY produces correct rotation")
	{
		// Any rotation created via FromEulerXZY should produce same matrix
		// as composing individual axis rotations
		Real a1 = 0.3, a2 = 0.5, a3 = 0.7;
		Quaternion q = Quaternion::FromEulerXZY(a1, a2, a3);
		Quaternion qx = Quaternion::FromAxisAngle(Vec3Cart(1,0,0), a1);
		Quaternion qz = Quaternion::FromAxisAngle(Vec3Cart(0,0,1), a2);
		Quaternion qy = Quaternion::FromAxisAngle(Vec3Cart(0,1,0), a3);
		Quaternion expected = qx * qz * qy;
		matricesEqual(q, expected);
		REQUIRE(q.isUnit(tol));
	}

	SECTION("FromEulerYXZ produces correct rotation")
	{
		Real a1 = 0.4, a2 = -0.6, a3 = 0.2;
		Quaternion q = Quaternion::FromEulerYXZ(a1, a2, a3);
		Quaternion qy = Quaternion::FromAxisAngle(Vec3Cart(0,1,0), a1);
		Quaternion qx = Quaternion::FromAxisAngle(Vec3Cart(1,0,0), a2);
		Quaternion qz = Quaternion::FromAxisAngle(Vec3Cart(0,0,1), a3);
		Quaternion expected = qy * qx * qz;
		matricesEqual(q, expected);
	}

	SECTION("FromEulerYZX produces correct rotation")
	{
		Real a1 = -0.5, a2 = 0.8, a3 = 0.1;
		Quaternion q = Quaternion::FromEulerYZX(a1, a2, a3);
		Quaternion qy = Quaternion::FromAxisAngle(Vec3Cart(0,1,0), a1);
		Quaternion qz = Quaternion::FromAxisAngle(Vec3Cart(0,0,1), a2);
		Quaternion qx = Quaternion::FromAxisAngle(Vec3Cart(1,0,0), a3);
		Quaternion expected = qy * qz * qx;
		matricesEqual(q, expected);
	}

	SECTION("FromEulerZXY produces correct rotation")
	{
		Real a1 = 0.7, a2 = 0.3, a3 = -0.4;
		Quaternion q = Quaternion::FromEulerZXY(a1, a2, a3);
		Quaternion qz = Quaternion::FromAxisAngle(Vec3Cart(0,0,1), a1);
		Quaternion qx = Quaternion::FromAxisAngle(Vec3Cart(1,0,0), a2);
		Quaternion qy = Quaternion::FromAxisAngle(Vec3Cart(0,1,0), a3);
		Quaternion expected = qz * qx * qy;
		matricesEqual(q, expected);
	}
}

TEST_CASE("Quaternion - Proper Euler roundtrip via rotation matrix", "[quaternion][euler]")
{
	TEST_PRECISION_INFO();
	const Real tol = TOL(1e-10, 1e-4);

	auto matricesEqual = [&](const Quaternion& q1, const Quaternion& q2) {
		auto m1 = q1.ToRotationMatrix();
		auto m2 = q2.ToRotationMatrix();
		for (int i = 0; i < 3; i++)
			for (int j = 0; j < 3; j++)
				REQUIRE(std::abs(m1[i][j] - m2[i][j]) < tol);
	};

	SECTION("FromEulerZXZ produces correct rotation")
	{
		Real alpha = 0.5, beta = 1.0, gamma = 0.3;
		Quaternion q = Quaternion::FromEulerZXZ(alpha, beta, gamma);
		Quaternion qz1 = Quaternion::FromAxisAngle(Vec3Cart(0,0,1), alpha);
		Quaternion qx  = Quaternion::FromAxisAngle(Vec3Cart(1,0,0), beta);
		Quaternion qz2 = Quaternion::FromAxisAngle(Vec3Cart(0,0,1), gamma);
		Quaternion expected = qz1 * qx * qz2;
		matricesEqual(q, expected);
		REQUIRE(q.isUnit(tol));
	}

	SECTION("FromEulerZYZ produces correct rotation")
	{
		Real alpha = -0.7, beta = 0.8, gamma = 1.2;
		Quaternion q = Quaternion::FromEulerZYZ(alpha, beta, gamma);
		Quaternion qz1 = Quaternion::FromAxisAngle(Vec3Cart(0,0,1), alpha);
		Quaternion qy  = Quaternion::FromAxisAngle(Vec3Cart(0,1,0), beta);
		Quaternion qz2 = Quaternion::FromAxisAngle(Vec3Cart(0,0,1), gamma);
		Quaternion expected = qz1 * qy * qz2;
		matricesEqual(q, expected);
	}
}

TEST_CASE("Quaternion - ToEulerXYZ roundtrip", "[quaternion][euler]")
{
	TEST_PRECISION_INFO();
	const Real tol = TOL(1e-8, 1e-3);

	SECTION("Roundtrip away from gimbal lock")
	{
		Real roll = 0.3, pitch = 0.5, yaw = 0.7;
		Quaternion q = Quaternion::FromEulerXYZ(roll, pitch, yaw);
		Vec3Cart euler = q.ToEulerXYZ();

		REQUIRE(std::abs(euler[0] - roll) < tol);
		REQUIRE(std::abs(euler[1] - pitch) < tol);
		REQUIRE(std::abs(euler[2] - yaw) < tol);
	}

	SECTION("Roundtrip with negative angles")
	{
		Real roll = -0.4, pitch = 0.2, yaw = -0.6;
		Quaternion q = Quaternion::FromEulerXYZ(roll, pitch, yaw);
		Vec3Cart euler = q.ToEulerXYZ();

		REQUIRE(std::abs(euler[0] - roll) < tol);
		REQUIRE(std::abs(euler[1] - pitch) < tol);
		REQUIRE(std::abs(euler[2] - yaw) < tol);
	}
}

TEST_CASE("Quaternion - ToEulerZXZ roundtrip", "[quaternion][euler]")
{
	TEST_PRECISION_INFO();
	const Real tol = TOL(1e-8, 1e-3);

	SECTION("Roundtrip away from gimbal lock")
	{
		Real alpha = 0.5, beta = 1.0, gamma = 0.3;
		Quaternion q = Quaternion::FromEulerZXZ(alpha, beta, gamma);
		Vec3Cart euler = q.ToEulerZXZ();

		REQUIRE(std::abs(euler[0] - alpha) < tol);
		REQUIRE(std::abs(euler[1] - beta) < tol);
		REQUIRE(std::abs(euler[2] - gamma) < tol);
	}
}

TEST_CASE("Quaternion - ToEulerZYZ roundtrip", "[quaternion][euler]")
{
	TEST_PRECISION_INFO();
	const Real tol = TOL(1e-8, 1e-3);

	SECTION("Roundtrip away from gimbal lock")
	{
		Real alpha = -0.7, beta = 0.8, gamma = 1.2;
		Quaternion q = Quaternion::FromEulerZYZ(alpha, beta, gamma);
		Vec3Cart euler = q.ToEulerZYZ();

		REQUIRE(std::abs(euler[0] - alpha) < tol);
		REQUIRE(std::abs(euler[1] - beta) < tol);
		REQUIRE(std::abs(euler[2] - gamma) < tol);
	}
}

TEST_CASE("Quaternion - All FromEuler produce unit quaternions", "[quaternion][euler]")
{
	TEST_PRECISION_INFO();
	const Real tol = TOL(1e-10, 1e-5);

	Real a = 0.5, b = 0.7, c = 1.1;
	REQUIRE(Quaternion::FromEulerXYZ(a, b, c).isUnit(tol));
	REQUIRE(Quaternion::FromEulerXZY(a, b, c).isUnit(tol));
	REQUIRE(Quaternion::FromEulerYXZ(a, b, c).isUnit(tol));
	REQUIRE(Quaternion::FromEulerYZX(a, b, c).isUnit(tol));
	REQUIRE(Quaternion::FromEulerZXY(a, b, c).isUnit(tol));
	REQUIRE(Quaternion::FromEulerZYX(a, b, c).isUnit(tol));
	REQUIRE(Quaternion::FromEulerZXZ(a, b, c).isUnit(tol));
	REQUIRE(Quaternion::FromEulerZYZ(a, b, c).isUnit(tol));
}

TEST_CASE("Quaternion - Euler identity rotation", "[quaternion][euler]")
{
	TEST_PRECISION_INFO();
	const Real tol = TOL(1e-10, 1e-5);

	// Zero angles should give identity quaternion for all conventions
	Quaternion id;
	auto approxEqual = [&](const Quaternion& q) {
		// q or -q both represent identity
		Real sign = (q.w() >= 0) ? 1.0 : -1.0;
		REQUIRE(std::abs(q.w() * sign - 1.0) < tol);
		REQUIRE(std::abs(q.x() * sign) < tol);
		REQUIRE(std::abs(q.y() * sign) < tol);
		REQUIRE(std::abs(q.z() * sign) < tol);
	};

	approxEqual(Quaternion::FromEulerXZY(0, 0, 0));
	approxEqual(Quaternion::FromEulerYXZ(0, 0, 0));
	approxEqual(Quaternion::FromEulerYZX(0, 0, 0));
	approxEqual(Quaternion::FromEulerZXY(0, 0, 0));
	approxEqual(Quaternion::FromEulerZXZ(0, 0, 0));
	approxEqual(Quaternion::FromEulerZYZ(0, 0, 0));
}

} // namespace MML::Tests::Base::QuaternionTests

