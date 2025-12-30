#include <catch2/catch_all.hpp>
#include "../TestPrecision.h"
#include "../TestMatchers.h"

#include "base/Quaternions.h"
#include "base/VectorTypes.h"

using namespace MML;
using namespace MML::Testing;

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
		REQUIRE(q.IsIdentity());
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
		REQUIRE(q.IsIdentity());
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
		REQUIRE(std::abs(q.w() - std::cos(Constants::PI / REAL(4.0))) < 1e-10);
		REQUIRE(std::abs(q.x()) < 1e-10);
		REQUIRE(std::abs(q.y()) < 1e-10);
		REQUIRE(std::abs(q.z() - std::sin(Constants::PI / REAL(4.0))) < 1e-10);
		REQUIRE(q.IsUnit());
	}

	SECTION("180 degree rotation around X-axis")
	{
		Vec3Cart axis(1, 0, 0);
		Real angle = Constants::PI;  // 180 degrees
		Quaternion q = Quaternion::FromAxisAngle(axis, angle);

		// Should be [0, 1, 0, 0]
		REQUIRE(std::abs(q.w()) < 1e-10);
		REQUIRE(std::abs(q.x() - REAL(1.0)) < 1e-10);
		REQUIRE(std::abs(q.y()) < 1e-10);
		REQUIRE(std::abs(q.z()) < 1e-10);
		REQUIRE(q.IsUnit());
	}

	SECTION("Zero rotation")
	{
		Vec3Cart axis(1, 0, 0);
		Real angle = REAL(0.0);
		Quaternion q = Quaternion::FromAxisAngle(axis, angle);
		REQUIRE(q.IsIdentity());
	}
}

TEST_CASE("Quaternion - From Euler Angles", "[quaternion]")
{
		TEST_PRECISION_INFO();
	SECTION("FromEulerZYX - Zero angles gives identity")
	{
		Quaternion q = Quaternion::FromEulerZYX(0, 0, 0);
		REQUIRE(q.IsIdentity());
	}

	SECTION("FromEulerXYZ - Zero angles gives identity")
	{
		Quaternion q = Quaternion::FromEulerXYZ(0, 0, 0);
		REQUIRE(q.IsIdentity());
	}

	SECTION("FromEulerZYX - 90 degree yaw only")
	{
		Quaternion q = Quaternion::FromEulerZYX(Constants::PI / REAL(2.0), 0, 0);
		REQUIRE(q.IsUnit(1e-10));
		
		// Rotate a point
		Vec3Cart v(1, 0, 0);
		Vec3Cart rotated = q.Rotate(v);
		REQUIRE(VectorApproxEqual(rotated, Vec3Cart(0, 1, 0), 1e-10));
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
		
		REQUIRE(result1.IsApprox(q));
		REQUIRE(result2.IsApprox(q));
	}

	SECTION("Non-commutativity")
	{
		Quaternion q1 = Quaternion::FromAxisAngle(Vec3Cart(1, 0, 0), Constants::PI / REAL(4.0));
		Quaternion q2 = Quaternion::FromAxisAngle(Vec3Cart(0, 1, 0), Constants::PI / REAL(4.0));
		
		Quaternion result1 = q1 * q2;
		Quaternion result2 = q2 * q1;
		
		// Should not be equal (within reasonable tolerance)
		REQUIRE_FALSE(result1.IsApprox(result2, 1e-6));
	}

	SECTION("Known quaternion multiplication")
	{
		// i * j = k in Hamilton convention
		Quaternion i(0, 1, 0, 0);
		Quaternion j(0, 0, 1, 0);
		Quaternion k_expected(0, 0, 0, 1);
		
		Quaternion k_result = i * j;
		REQUIRE(k_result.IsApprox(k_expected, 1e-10));
	}

	SECTION("Conjugate property: (q1*q2)* = q2* * q1*")
	{
		Quaternion q1(1, 2, 3, 4);
		Quaternion q2(5, 6, 7, 8);
		
		Quaternion lhs = (q1 * q2).Conjugate();
		Quaternion rhs = q2.Conjugate() * q1.Conjugate();
		
		REQUIRE(lhs.IsApprox(rhs, 1e-10));
	}
}

TEST_CASE("Quaternion - Conjugate and Inverse", "[quaternion]")
{
		TEST_PRECISION_INFO();
	SECTION("Conjugate of identity")
	{
		Quaternion q = Quaternion::Identity();
		Quaternion conj = q.Conjugate();
		REQUIRE(conj.IsIdentity());
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
		REQUIRE(double_conj.IsApprox(q));
	}

	SECTION("Inverse of unit quaternion equals conjugate")
	{
		Quaternion q = Quaternion::FromAxisAngle(Vec3Cart(1, 1, 1).Normalized(), Constants::PI / REAL(3.0));
		Quaternion inv = q.Inverse();
		Quaternion conj = q.Conjugate();
		
		REQUIRE(inv.IsApprox(conj, 1e-10));
	}

	SECTION("q * q^(-1) = identity")
	{
		Quaternion q(1, 2, 3, 4);
		Quaternion inv = q.Inverse();
		Quaternion result = q * inv;
		
		REQUIRE(result.IsIdentity(1e-10));
	}

	SECTION("q^(-1) * q = identity")
	{
		Quaternion q(1, 2, 3, 4);
		Quaternion inv = q.Inverse();
		Quaternion result = inv * q;
		
		REQUIRE(result.IsIdentity(1e-10));
	}
}

TEST_CASE("Quaternion - Norm and Normalization", "[quaternion]")
{
		TEST_PRECISION_INFO();
	SECTION("Norm of identity is 1")
	{
		Quaternion q = Quaternion::Identity();
		REQUIRE(std::abs(q.Norm() - REAL(1.0)) < 1e-10);
		REQUIRE(std::abs(q.NormSquared() - REAL(1.0)) < 1e-10);
	}

	SECTION("Norm calculation")
	{
		Quaternion q(1, 2, 2, 0);
		Real expected_norm = std::sqrt(1 + 4 + 4);  // sqrt(9) = 3
		REQUIRE(std::abs(q.Norm() - expected_norm) < 1e-10);
		REQUIRE(std::abs(q.NormSquared() - REAL(9.0)) < 1e-10);
	}

	SECTION("Normalization")
	{
		Quaternion q(2, 0, 0, 0);
		q.Normalize();
		REQUIRE(q.IsUnit(1e-10));
		REQUIRE(std::abs(q.w() - REAL(1.0)) < 1e-10);
	}

	SECTION("Normalized returns normalized copy")
	{
		Quaternion q(1, 2, 3, 4);
		Quaternion orig = q;
		Quaternion normalized = q.Normalized();
		
		REQUIRE(normalized.IsUnit(1e-10));
		REQUIRE(q.IsApprox(orig));  // Original unchanged
	}

	SECTION("FromAxisAngle produces unit quaternion")
	{
		Vec3Cart axis(1, 2, 3);
		axis = axis.Normalized();
		Quaternion q = Quaternion::FromAxisAngle(axis, REAL(1.5));
		REQUIRE(q.IsUnit(1e-10));
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
		REQUIRE(VectorApproxEqual(rotated, v, 1e-10));
	}

	SECTION("90 degree rotation around Z-axis")
	{
		// Rotate (1, 0, 0) by 90° around Z-axis -> (0, 1, 0)
		Quaternion q = Quaternion::FromAxisAngle(Vec3Cart(0, 0, 1), Constants::PI / REAL(2.0));
		Vec3Cart v(1, 0, 0);
		Vec3Cart rotated = q.Rotate(v);
		
		REQUIRE(VectorApproxEqual(rotated, Vec3Cart(0, 1, 0), 1e-10));
	}

	SECTION("90 degree rotation around X-axis")
	{
		// Rotate (0, 1, 0) by 90° around X-axis -> (0, 0, 1)
		Quaternion q = Quaternion::FromAxisAngle(Vec3Cart(1, 0, 0), Constants::PI / REAL(2.0));
		Vec3Cart v(0, 1, 0);
		Vec3Cart rotated = q.Rotate(v);
		
		REQUIRE(VectorApproxEqual(rotated, Vec3Cart(0, 0, 1), 1e-10));
	}

	SECTION("90 degree rotation around Y-axis")
	{
		// Rotate (1, 0, 0) by 90° around Y-axis -> (0, 0, -1)
		Quaternion q = Quaternion::FromAxisAngle(Vec3Cart(0, 1, 0), Constants::PI / REAL(2.0));
		Vec3Cart v(1, 0, 0);
		Vec3Cart rotated = q.Rotate(v);
		
		REQUIRE(VectorApproxEqual(rotated, Vec3Cart(0, 0, -1), 1e-10));
	}

	SECTION("180 degree rotation around Z-axis")
	{
		// Rotate (1, 0, 0) by 180° around Z-axis -> (-1, 0, 0)
		Quaternion q = Quaternion::FromAxisAngle(Vec3Cart(0, 0, 1), Constants::PI);
		Vec3Cart v(1, 0, 0);
		Vec3Cart rotated = q.Rotate(v);
		
		REQUIRE(VectorApproxEqual(rotated, Vec3Cart(-1, 0, 0), 1e-10));
	}

	SECTION("Rotation preserves vector length")
	{
		Quaternion q = Quaternion::FromAxisAngle(Vec3Cart(1, 1, 1).Normalized(), REAL(0.7));
		Vec3Cart v(1, 2, 3);
		Vec3Cart rotated = q.Rotate(v);
		
		REQUIRE(std::abs(rotated.NormL2() - v.NormL2()) < 1e-10);
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
		
		REQUIRE(VectorApproxEqual(rotated_separate, rotated_combined, 1e-10));
	}

	SECTION("Inverse rotation undoes rotation")
	{
		Quaternion q = Quaternion::FromAxisAngle(Vec3Cart(1, 2, 3).Normalized(), REAL(1.234));
		Vec3Cart v(4, 5, 6);
		
		Vec3Cart rotated = q.Rotate(v);
		Vec3Cart back = q.Inverse().Rotate(rotated);
		
		REQUIRE(VectorApproxEqual(back, v, 1e-10));
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
		
		REQUIRE(std::abs(angle_out - angle_in) < 1e-10);
		REQUIRE(VectorApproxEqual(axis_out, axis_in, 1e-10));
	}

	SECTION("GetRotationAxis and GetRotationAngle")
	{
		Vec3Cart axis_in(0, 1, 0);
		Real angle_in = Constants::PI / REAL(4.0);
		
		Quaternion q = Quaternion::FromAxisAngle(axis_in, angle_in);
		
		Vec3Cart axis_out = q.GetRotationAxis();
		Real angle_out = q.GetRotationAngle();
		
		REQUIRE(std::abs(angle_out - angle_in) < 1e-10);
		REQUIRE(VectorApproxEqual(axis_out, axis_in, 1e-10));
	}

	SECTION("Identity quaternion has zero angle")
	{
		Quaternion q = Quaternion::Identity();
		Real angle = q.GetRotationAngle();
		REQUIRE(std::abs(angle) < 1e-10);
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
		REQUIRE(std::abs(euler[0] - Constants::PI / REAL(2.0)) < 1e-10);
		REQUIRE(std::abs(euler[1]) < 1e-10);
		REQUIRE(std::abs(euler[2]) < 1e-10);
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
		REQUIRE(result.IsApprox(q1));
	}

	SECTION("Lerp at t=1 returns second quaternion")
	{
		Quaternion result = Quaternion::Lerp(q1, q2, REAL(1.0));
		REQUIRE(result.IsApprox(q2));
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
		REQUIRE(result.IsApprox(q1, 1e-6));
	}

	SECTION("Slerp at t=1 returns second quaternion")
	{
		Quaternion result = Quaternion::Slerp(q1, q2, REAL(1.0));
		REQUIRE(result.IsApprox(q2, 1e-6));
	}

	SECTION("Slerp produces unit quaternions")
	{
		Quaternion result = Quaternion::Slerp(q1, q2, REAL(0.3));
		REQUIRE(result.IsUnit(1e-6));
		
		result = Quaternion::Slerp(q1, q2, REAL(0.7));
		REQUIRE(result.IsUnit(1e-6));
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
		REQUIRE(result.IsUnit(1e-6));
	}
}

TEST_CASE("Quaternion - Dot Product", "[quaternion]")
{
		TEST_PRECISION_INFO();
	SECTION("Dot product of identical quaternions")
	{
		Quaternion q(1, 2, 3, 4);
		Real dot = q.Dot(q);
		REQUIRE(std::abs(dot - q.NormSquared()) < 1e-10);
	}

	SECTION("Dot product is commutative")
	{
		Quaternion q1(1, 2, 3, 4);
		Quaternion q2(5, 6, 7, 8);
		REQUIRE(std::abs(q1.Dot(q2) - q2.Dot(q1)) < 1e-10);
	}

	SECTION("Orthogonal quaternions have zero dot product")
	{
		Quaternion q1(1, 0, 0, 0);
		Quaternion q2(0, 1, 0, 0);
		REQUIRE(std::abs(q1.Dot(q2)) < 1e-10);
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
		
		REQUIRE(q1.IsApprox(q2, 1e-6));
		REQUIRE_FALSE(q1.IsApprox(q2, 1e-8));
	}
}

TEST_CASE("Quaternion - Special Cases and Edge Cases", "[quaternion]")
{
		TEST_PRECISION_INFO();
	SECTION("Zero rotation has identity quaternion")
	{
		Quaternion q = Quaternion::FromAxisAngle(Vec3Cart(1, 0, 0), REAL(0.0));
		REQUIRE(q.IsIdentity());
	}

	SECTION("Full rotation (2π) represents same rotation as identity")
	{
		// Note: q = [cos(π), sin(π)*axis] = [-1, 0, 0, 0]
		// Both q and -q represent the same rotation!
		Quaternion q = Quaternion::FromAxisAngle(Vec3Cart(1, 0, 0), REAL(2.0) * Constants::PI);
		
		// Test by rotating a vector - should get same result as identity
		Vec3Cart v(1, 2, 3);
		Vec3Cart rotated = q.Rotate(v);
		REQUIRE(VectorApproxEqual(rotated, v, 1e-10));
	}

	SECTION("Negative angle gives inverse rotation")
	{
		Vec3Cart axis(1, 1, 1);
		axis = axis.Normalized();
		Real angle = REAL(0.7);
		
		Quaternion q_pos = Quaternion::FromAxisAngle(axis, angle);
		Quaternion q_neg = Quaternion::FromAxisAngle(axis, -angle);
		
		Quaternion product = q_pos * q_neg;
		REQUIRE(product.IsIdentity(1e-10));
	}

	SECTION("IsUnit check with tolerance")
	{
		Quaternion q(REAL(1.0) + 1e-7, 0, 0, 0);
		REQUIRE(q.IsUnit(1e-6));
		REQUIRE_FALSE(q.IsUnit(1e-8));
	}

	SECTION("Rotation of zero vector")
	{
		Quaternion q = Quaternion::FromAxisAngle(Vec3Cart(1, 0, 0), REAL(1.5));
		Vec3Cart zero(0, 0, 0);
		Vec3Cart rotated = q.Rotate(zero);
		REQUIRE(VectorApproxEqual(rotated, zero, 1e-10));
	}

	SECTION("q and -q represent same rotation")
	{
		Quaternion q = Quaternion::FromAxisAngle(Vec3Cart(1, 0, 0), Constants::PI / REAL(4.0));
		Quaternion q_neg = -q;
		
		Vec3Cart v(1, 2, 3);
		Vec3Cart rot1 = q.Rotate(v);
		Vec3Cart rot2 = q_neg.Rotate(v);
		
		REQUIRE(VectorApproxEqual(rot1, rot2, 1e-10));
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
		REQUIRE(q1.IsApprox(expected, 1e-10));
	}
}
