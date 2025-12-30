#include <catch2/catch_all.hpp>
#include "../TestPrecision.h"
#include "../TestMatchers.h"

#include "base/Quaternions.h"
#include "base/VectorTypes.h"
#include "core/CoordTransf/CoordTransf3D.h"

using namespace MML;
using namespace MML::Testing;

// Helper function for approximate equality of vectors
static bool VectorApproxEqual(const Vec3Cart& v1, const Vec3Cart& v2, Real tolerance = 1e-10)
{
	return std::abs(v1[0] - v2[0]) < tolerance &&
				 std::abs(v1[1] - v2[1]) < tolerance &&
				 std::abs(v1[2] - v2[2]) < tolerance;
}

TEST_CASE("Quaternion - Matrix Conversion", "[quaternion][matrix]")
{
		TEST_PRECISION_INFO();
	SECTION("Identity quaternion to identity matrix")
	{
		Quaternion q = Quaternion::Identity();
		MatrixNM<Real, 3, 3> mat = q.ToRotationMatrix();

		// Should be identity matrix
		for (int i = 0; i < 3; i++)
			for (int j = 0; j < 3; j++)
			{
				if (i == j)
					REQUIRE(std::abs(mat[i][j] - REAL(1.0)) < 1e-10);
				else
					REQUIRE(std::abs(mat[i][j]) < 1e-10);
			}
	}

	SECTION("90 degree rotation around Z-axis to matrix")
	{
		Quaternion q = Quaternion::FromAxisAngle(Vec3Cart(0, 0, 1), Constants::PI / REAL(2.0));
		MatrixNM<Real, 3, 3> mat = q.ToRotationMatrix();

		// Apply to vector (1,0,0) -> should get (0,1,0)
		Vec3Cart v(1, 0, 0);
		Vec3Cart result(mat[0][0] * v[0] + mat[0][1] * v[1] + mat[0][2] * v[2],
									mat[1][0] * v[0] + mat[1][1] * v[1] + mat[1][2] * v[2],
									mat[2][0] * v[0] + mat[2][1] * v[1] + mat[2][2] * v[2]);

		REQUIRE(VectorApproxEqual(result, Vec3Cart(0, 1, 0)));
	}

	SECTION("Matrix to quaternion and back")
	{
		// Create a known rotation quaternion
		Quaternion q_original = Quaternion::FromAxisAngle(
			Vec3Cart(1, 1, 1).Normalized(), Constants::PI / REAL(3.0));

		// Convert to matrix
		MatrixNM<Real, 3, 3> mat = q_original.ToRotationMatrix();

		// Convert back to quaternion
		Quaternion q_reconstructed = Quaternion::FromRotationMatrix(mat);

		// Quaternions q and -q represent the same rotation
		// Check if they're equal or negatives
		bool same = q_original.IsApprox(q_reconstructed, 1e-10) ||
								q_original.IsApprox(-q_reconstructed, 1e-10);
		REQUIRE(same);

		// Verify they rotate vectors the same way
		Vec3Cart v(1, 2, 3);
		Vec3Cart rot1 = q_original.Rotate(v);
		Vec3Cart rot2 = q_reconstructed.Rotate(v);
		REQUIRE(VectorApproxEqual(rot1, rot2));
	}

	SECTION("FromRotationMatrix handles all four cases")
	{
		// Case 1: trace > 0 (w largest)
		Quaternion q1 = Quaternion::FromAxisAngle(Vec3Cart(1, 1, 1).Normalized(), REAL(0.5));
		MatrixNM<Real, 3, 3> mat1 = q1.ToRotationMatrix();
		Quaternion q1_back = Quaternion::FromRotationMatrix(mat1);
		REQUIRE(VectorApproxEqual(q1.Rotate(Vec3Cart(1, 0, 0)), 
															q1_back.Rotate(Vec3Cart(1, 0, 0))));

		// Case 2: x largest (90° around X-axis)
		Quaternion q2 = Quaternion::FromAxisAngle(Vec3Cart(1, 0, 0), Constants::PI / REAL(2.0));
		MatrixNM<Real, 3, 3> mat2 = q2.ToRotationMatrix();
		Quaternion q2_back = Quaternion::FromRotationMatrix(mat2);
		REQUIRE(VectorApproxEqual(q2.Rotate(Vec3Cart(0, 1, 0)), 
															q2_back.Rotate(Vec3Cart(0, 1, 0))));

		// Case 3: y largest (90° around Y-axis)
		Quaternion q3 = Quaternion::FromAxisAngle(Vec3Cart(0, 1, 0), Constants::PI / REAL(2.0));
		MatrixNM<Real, 3, 3> mat3 = q3.ToRotationMatrix();
		Quaternion q3_back = Quaternion::FromRotationMatrix(mat3);
		REQUIRE(VectorApproxEqual(q3.Rotate(Vec3Cart(1, 0, 0)), 
															q3_back.Rotate(Vec3Cart(1, 0, 0))));

		// Case 4: z largest (90° around Z-axis)
		Quaternion q4 = Quaternion::FromAxisAngle(Vec3Cart(0, 0, 1), Constants::PI / REAL(2.0));
		MatrixNM<Real, 3, 3> mat4 = q4.ToRotationMatrix();
		Quaternion q4_back = Quaternion::FromRotationMatrix(mat4);
		REQUIRE(VectorApproxEqual(q4.Rotate(Vec3Cart(1, 0, 0)), 
															q4_back.Rotate(Vec3Cart(1, 0, 0))));
	}

	SECTION("Matrix determinant is 1 (proper rotation)")
	{
		Quaternion q = Quaternion::FromAxisAngle(Vec3Cart(1, 2, 3).Normalized(), REAL(1.234));
		MatrixNM<Real, 3, 3> mat = q.ToRotationMatrix();

		// Compute determinant (3x3)
		Real det = mat[0][0] * (mat[1][1] * mat[2][2] - mat[1][2] * mat[2][1])
						 - mat[0][1] * (mat[1][0] * mat[2][2] - mat[1][2] * mat[2][0])
						 + mat[0][2] * (mat[1][0] * mat[2][1] - mat[1][1] * mat[2][0]);

		REQUIRE(std::abs(det - REAL(1.0)) < 1e-10);
	}

	SECTION("Matrix is orthogonal (M^T * M = I)")
	{
		Quaternion q = Quaternion::FromAxisAngle(Vec3Cart(0, 1, 0), REAL(0.789));
		MatrixNM<Real, 3, 3> mat = q.ToRotationMatrix();

		// Compute M^T * M
		MatrixNM<Real, 3, 3> product;
		for (int i = 0; i < 3; i++)
			for (int j = 0; j < 3; j++)
			{
				product[i][j] = 0;
				for (int k = 0; k < 3; k++)
					product[i][j] += mat[k][i] * mat[k][j];  // M^T[i][k] * M[k][j]
			}

		// Should be identity
		for (int i = 0; i < 3; i++)
			for (int j = 0; j < 3; j++)
			{
				if (i == j)
					REQUIRE(std::abs(product[i][j] - REAL(1.0)) < 1e-10);
				else
					REQUIRE(std::abs(product[i][j]) < 1e-10);
			}
	}
}

TEST_CASE("CoordTransfCart3DRotationQuaternion - Construction", "[coordtransf][quaternion]")
{
		TEST_PRECISION_INFO();
	SECTION("Construct from quaternion")
	{
		Quaternion q = Quaternion::FromAxisAngle(Vec3Cart(0, 0, 1), Constants::PI / REAL(2.0));
		CoordTransfCart3DRotationQuaternion rot(q);

		Vec3Cart v(1, 0, 0);
		Vec3Cart result = rot.transf(v);
		REQUIRE(VectorApproxEqual(result, Vec3Cart(0, 1, 0)));
	}

	SECTION("Construct from axis-angle")
	{
		CoordTransfCart3DRotationQuaternion rot(Vec3Cart(1, 0, 0), Constants::PI / REAL(2.0));

		Vec3Cart v(0, 1, 0);
		Vec3Cart result = rot.transf(v);
		REQUIRE(VectorApproxEqual(result, Vec3Cart(0, 0, 1)));
	}

	SECTION("Construct from rotation matrix")
	{
		// Create matrix for 90° rotation around Y-axis
		MatrixNM<Real, 3, 3> mat;
		mat[0][0] = 0;  mat[0][1] = 0;  mat[0][2] = 1;
		mat[1][0] = 0;  mat[1][1] = 1;  mat[1][2] = 0;
		mat[2][0] = -1; mat[2][1] = 0;  mat[2][2] = 0;

		CoordTransfCart3DRotationQuaternion rot(mat);

		Vec3Cart v(1, 0, 0);
		Vec3Cart result = rot.transf(v);
		REQUIRE(VectorApproxEqual(result, Vec3Cart(0, 0, -1)));
	}

	SECTION("Construct from Euler angles ZYX")
	{
		CoordTransfCart3DRotationQuaternion rot = 
			CoordTransfCart3DRotationQuaternion::FromEulerZYX(
				Constants::PI / REAL(2.0), 0, 0);  // Pure yaw

		Vec3Cart v(1, 0, 0);
		Vec3Cart result = rot.transf(v);
		REQUIRE(VectorApproxEqual(result, Vec3Cart(0, 1, 0)));
	}

	SECTION("Non-unit quaternion is automatically normalized")
	{
		Quaternion q(2, 0, 0, 0);  // Not unit
		CoordTransfCart3DRotationQuaternion rot(q);

		// Should still work as identity rotation
		Vec3Cart v(1, 2, 3);
		Vec3Cart result = rot.transf(v);
		REQUIRE(VectorApproxEqual(result, v));
	}
}

TEST_CASE("CoordTransfCart3DRotationQuaternion - Transformations", "[coordtransf][quaternion]")
{
		TEST_PRECISION_INFO();
	SECTION("Identity rotation does nothing")
	{
		Quaternion q = Quaternion::Identity();
		CoordTransfCart3DRotationQuaternion rot(q);

		Vec3Cart v(1, 2, 3);
		Vec3Cart result = rot.transf(v);
		REQUIRE(VectorApproxEqual(result, v));
	}

	SECTION("90 degree rotations around principal axes")
	{
		// X-axis
		CoordTransfCart3DRotationQuaternion rotX(Vec3Cart(1, 0, 0), Constants::PI / REAL(2.0));
		REQUIRE(VectorApproxEqual(rotX.transf(Vec3Cart(0, 1, 0)), Vec3Cart(0, 0, 1)));

		// Y-axis
		CoordTransfCart3DRotationQuaternion rotY(Vec3Cart(0, 1, 0), Constants::PI / REAL(2.0));
		REQUIRE(VectorApproxEqual(rotY.transf(Vec3Cart(1, 0, 0)), Vec3Cart(0, 0, -1)));

		// Z-axis
		CoordTransfCart3DRotationQuaternion rotZ(Vec3Cart(0, 0, 1), Constants::PI / REAL(2.0));
		REQUIRE(VectorApproxEqual(rotZ.transf(Vec3Cart(1, 0, 0)), Vec3Cart(0, 1, 0)));
	}

	SECTION("Inverse transformation undoes rotation")
	{
		CoordTransfCart3DRotationQuaternion rot(
			Vec3Cart(1, 1, 1).Normalized(), REAL(1.234));

		Vec3Cart v(1, 2, 3);
		Vec3Cart rotated = rot.transf(v);
		Vec3Cart back = rot.transfInverse(rotated);

		REQUIRE(VectorApproxEqual(back, v));
	}

	SECTION("Transformation preserves vector length")
	{
		CoordTransfCart3DRotationQuaternion rot(
			Vec3Cart(1, 2, 3).Normalized(), REAL(0.789));

		Vec3Cart v(4, 5, 6);
		Vec3Cart rotated = rot.transf(v);

		REQUIRE(std::abs(rotated.NormL2() - v.NormL2()) < 1e-10);
	}
}

TEST_CASE("CoordTransfCart3DRotationQuaternion - Composition", "[coordtransf][quaternion]")
{
		TEST_PRECISION_INFO();
	SECTION("Compose two rotations")
	{
		CoordTransfCart3DRotationQuaternion rot1(Vec3Cart(0, 0, 1), Constants::PI / REAL(4.0));
		CoordTransfCart3DRotationQuaternion rot2(Vec3Cart(1, 0, 0), Constants::PI / REAL(4.0));

		CoordTransfCart3DRotationQuaternion combined = rot1.Compose(rot2);

		Vec3Cart v(1, 0, 0);
		Vec3Cart result_separate = rot2.transf(rot1.transf(v));
		Vec3Cart result_combined = combined.transf(v);

		REQUIRE(VectorApproxEqual(result_separate, result_combined));
	}

	SECTION("Composition is associative")
	{
		CoordTransfCart3DRotationQuaternion r1(Vec3Cart(1, 0, 0), REAL(0.5));
		CoordTransfCart3DRotationQuaternion r2(Vec3Cart(0, 1, 0), REAL(0.7));
		CoordTransfCart3DRotationQuaternion r3(Vec3Cart(0, 0, 1), REAL(0.3));

		CoordTransfCart3DRotationQuaternion combined1 = r1.Compose(r2).Compose(r3);
		CoordTransfCart3DRotationQuaternion combined2 = r1.Compose(r2.Compose(r3));

		Vec3Cart v(1, 2, 3);
		Vec3Cart result1 = combined1.transf(v);
		Vec3Cart result2 = combined2.transf(v);

		REQUIRE(VectorApproxEqual(result1, result2));
	}

	SECTION("Compose with identity is identity")
	{
		CoordTransfCart3DRotationQuaternion rot(Vec3Cart(1, 1, 1).Normalized(), REAL(1.0));
		CoordTransfCart3DRotationQuaternion identity(Quaternion::Identity());

		CoordTransfCart3DRotationQuaternion result1 = rot.Compose(identity);
		CoordTransfCart3DRotationQuaternion result2 = identity.Compose(rot);

		Vec3Cart v(1, 2, 3);
		REQUIRE(VectorApproxEqual(result1.transf(v), rot.transf(v)));
		REQUIRE(VectorApproxEqual(result2.transf(v), rot.transf(v)));
	}
}

TEST_CASE("CoordTransfCart3DRotationQuaternion - Interpolation", "[coordtransf][quaternion]")
{
		TEST_PRECISION_INFO();
	CoordTransfCart3DRotationQuaternion rot1(Quaternion::Identity());
	CoordTransfCart3DRotationQuaternion rot2(Vec3Cart(0, 0, 1), Constants::PI / REAL(2.0));

	SECTION("Interpolation at t=0 gives first rotation")
	{
		CoordTransfCart3DRotationQuaternion result = rot1.Interpolate(rot2, REAL(0.0));

		Vec3Cart v(1, 0, 0);
		REQUIRE(VectorApproxEqual(result.transf(v), rot1.transf(v), 1e-6));
	}

	SECTION("Interpolation at t=1 gives second rotation")
	{
		CoordTransfCart3DRotationQuaternion result = rot1.Interpolate(rot2, REAL(1.0));

		Vec3Cart v(1, 0, 0);
		REQUIRE(VectorApproxEqual(result.transf(v), rot2.transf(v), 1e-6));
	}

	SECTION("Interpolation at t=REAL(0.5) gives halfway rotation")
	{
		CoordTransfCart3DRotationQuaternion result = rot1.Interpolate(rot2, REAL(0.5));

		// Should be 45° rotation around Z-axis
		Vec3Cart v(1, 0, 0);
		Vec3Cart rotated = result.transf(v);

		// Result should be approximately (cos(45°), sin(45°), 0)
		Real sqrt2_inv = REAL(1.0) / std::sqrt(REAL(2.0));
		REQUIRE(VectorApproxEqual(rotated, Vec3Cart(sqrt2_inv, sqrt2_inv, 0), 1e-6));
	}

	SECTION("Interpolated rotations preserve vector length")
	{
		CoordTransfCart3DRotationQuaternion result = rot1.Interpolate(rot2, REAL(0.3));

		Vec3Cart v(1, 2, 3);
		Vec3Cart rotated = result.transf(v);

		REQUIRE(std::abs(rotated.NormL2() - v.NormL2()) < 1e-10);
	}
}

TEST_CASE("CoordTransfCart3DRotationQuaternion - Accessors", "[coordtransf][quaternion]")
{
		TEST_PRECISION_INFO();
	Vec3Cart axis(1, 1, 1);
	axis = axis.Normalized();
	Real angle = Constants::PI / REAL(3.0);

	CoordTransfCart3DRotationQuaternion rot(axis, angle);

	SECTION("GetQuaternion returns correct quaternion")
	{
		Quaternion q = rot.GetQuaternion();
		REQUIRE(q.IsUnit(1e-10));

		// Verify it rotates correctly
		Vec3Cart v(1, 0, 0);
		REQUIRE(VectorApproxEqual(q.Rotate(v), rot.transf(v)));
	}

	SECTION("GetRotationAxis returns correct axis")
	{
		Vec3Cart retrieved_axis = rot.GetRotationAxis();
		REQUIRE(VectorApproxEqual(retrieved_axis, axis));
	}

	SECTION("GetRotationAngle returns correct angle")
	{
		Real retrieved_angle = rot.GetRotationAngle();
		REQUIRE(std::abs(retrieved_angle - angle) < 1e-10);
	}

	SECTION("GetTransformationMatrix returns rotation matrix")
	{
		MatrixNM<Real, 3, 3> mat = rot.GetTransformationMatrix();

		// Apply matrix to vector and compare with quaternion rotation
		Vec3Cart v(1, 2, 3);
		Vec3Cart mat_result(
			mat[0][0] * v[0] + mat[0][1] * v[1] + mat[0][2] * v[2],
			mat[1][0] * v[0] + mat[1][1] * v[1] + mat[1][2] * v[2],
			mat[2][0] * v[0] + mat[2][1] * v[1] + mat[2][2] * v[2]
		);
		Vec3Cart quat_result = rot.transf(v);

		REQUIRE(VectorApproxEqual(mat_result, quat_result));
	}
}

TEST_CASE("CoordTransfCart3DRotationQuaternion - Compatibility with Matrix Rotations", "[coordtransf][quaternion]")
{
		TEST_PRECISION_INFO();
	SECTION("Matches CoordTransfCart3DRotationXAxis")
	{
		Real angle = Constants::PI / REAL(3.0);

		CoordTransfCart3DRotationXAxis matRot(angle);
		CoordTransfCart3DRotationQuaternion quatRot(Vec3Cart(1, 0, 0), angle);

		Vec3Cart v(1, 2, 3);
		Vec3Cart result_mat = matRot.transf(v);
		Vec3Cart result_quat = quatRot.transf(v);

		REQUIRE(VectorApproxEqual(result_mat, result_quat));
	}

	SECTION("Matches CoordTransfCart3DRotationYAxis")
	{
		Real angle = Constants::PI / REAL(4.0);

		CoordTransfCart3DRotationYAxis matRot(angle);
		CoordTransfCart3DRotationQuaternion quatRot(Vec3Cart(0, 1, 0), angle);

		Vec3Cart v(1, 2, 3);
		Vec3Cart result_mat = matRot.transf(v);
		Vec3Cart result_quat = quatRot.transf(v);

		REQUIRE(VectorApproxEqual(result_mat, result_quat));
	}

	SECTION("Matches CoordTransfCart3DRotationZAxis")
	{
		Real angle = Constants::PI / REAL(6.0);

		CoordTransfCart3DRotationZAxis matRot(angle);
		CoordTransfCart3DRotationQuaternion quatRot(Vec3Cart(0, 0, 1), angle);

		Vec3Cart v(1, 2, 3);
		Vec3Cart result_mat = matRot.transf(v);
		Vec3Cart result_quat = quatRot.transf(v);

		REQUIRE(VectorApproxEqual(result_mat, result_quat));
	}
}

TEST_CASE("CoordTransfCart3DRotationQuaternion - Edge Cases", "[coordtransf][quaternion]")
{
		TEST_PRECISION_INFO();
	SECTION("Zero rotation")
	{
		CoordTransfCart3DRotationQuaternion rot(Vec3Cart(1, 0, 0), REAL(0.0));

		Vec3Cart v(1, 2, 3);
		REQUIRE(VectorApproxEqual(rot.transf(v), v));
	}

	SECTION("180 degree rotation")
	{
		CoordTransfCart3DRotationQuaternion rot(Vec3Cart(0, 0, 1), Constants::PI);

		Vec3Cart v(1, 0, 0);
		Vec3Cart result = rot.transf(v);
		REQUIRE(VectorApproxEqual(result, Vec3Cart(-1, 0, 0)));
	}

	SECTION("Rotation of zero vector")
	{
		CoordTransfCart3DRotationQuaternion rot(Vec3Cart(1, 1, 1).Normalized(), REAL(1.5));

		Vec3Cart zero(0, 0, 0);
		Vec3Cart result = rot.transf(zero);
		REQUIRE(VectorApproxEqual(result, zero));
	}

	SECTION("Multiple full rotations (2π)")
	{
		CoordTransfCart3DRotationQuaternion rot(Vec3Cart(1, 0, 0), REAL(2.0) * Constants::PI);

		Vec3Cart v(1, 2, 3);
		Vec3Cart result = rot.transf(v);
		// Should be same as original (within tolerance)
		REQUIRE(VectorApproxEqual(result, v));
	}
}
