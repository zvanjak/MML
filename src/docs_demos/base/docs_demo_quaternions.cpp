///////////////////////////////////////////////////////////////////////////////////////////
///                         Docs Demo - Quaternions                                   ///
///                                                                                   ///
///  Demonstrates quaternion operations from base/Quaternions.h                       ///
///////////////////////////////////////////////////////////////////////////////////////////

#include <iostream>
#include <iomanip>
#include <cmath>

#include "MMLBase.h"
#include "base/Quaternions.h"
#include "base/VectorTypes.h"

using namespace MML;

///////////////////////////////////////////////////////////////////////////////////////////
///                         CONSTRUCTORS AND FACTORY METHODS                            ///
///////////////////////////////////////////////////////////////////////////////////////////

void Docs_Demo_Quaternions_Constructors()
{
	std::cout << "\n==========================================================================\n";
	std::cout << "Demo: Quaternion Constructors and Factory Methods\n";
	std::cout << "==========================================================================\n\n";

	// Default constructor - identity quaternion
	std::cout << "--- Constructors ---\n";
	Quaternion q1;
	std::cout << "Default (identity): " << q1 << "\n";

	// From components
	Quaternion q2(0.707, 0.707, 0, 0);
	std::cout << "From components (0.707, 0.707, 0, 0): " << q2 << "\n";

	// From scalar and vector
	Vec3Cart axis(0, 0, 1);
	Quaternion q3(0.707, Vec3Cart(0, 0, 0.707));
	std::cout << "From scalar and vector: " << q3 << "\n";

	// Pure imaginary
	Quaternion q4(Vec3Cart(1, 2, 3));
	std::cout << "Pure imaginary from (1,2,3): " << q4 << "\n";

	// Factory methods
	std::cout << "\n--- Factory Methods ---\n";
	
	// FromAxisAngle
	Vec3Cart zAxis(0, 0, 1);
	auto qRot90Z = Quaternion::FromAxisAngle(zAxis, Constants::PI / 2);
	std::cout << "90° around Z-axis: " << qRot90Z << "\n";
	std::cout << "  Expected: [cos(45°), 0, 0, sin(45°)] ≈ [0.707, 0, 0, 0.707]\n";

	// FromEulerZYX
	auto qEuler = Quaternion::FromEulerZYX(0.5, 0.2, 0.1);
	std::cout << "FromEulerZYX(yaw=0.5, pitch=0.2, roll=0.1): " << qEuler << "\n";

	// Identity
	auto qId = Quaternion::Identity();
	std::cout << "Identity(): " << qId << "\n";
}

///////////////////////////////////////////////////////////////////////////////////////////
///                         ACCESSORS                                                   ///
///////////////////////////////////////////////////////////////////////////////////////////

void Docs_Demo_Quaternions_Accessors()
{
	std::cout << "\n==========================================================================\n";
	std::cout << "Demo: Quaternion Accessors\n";
	std::cout << "==========================================================================\n\n";

	Quaternion q(0.5, 0.5, 0.5, 0.5);
	std::cout << "Quaternion: " << q << "\n\n";

	std::cout << "--- Component Access ---\n";
	std::cout << "q.w() = " << q.w() << " (scalar part)\n";
	std::cout << "q.x() = " << q.x() << " (i component)\n";
	std::cout << "q.y() = " << q.y() << " (j component)\n";
	std::cout << "q.z() = " << q.z() << " (k component)\n";

	std::cout << "\n--- Index Access ---\n";
	std::cout << "q[0] = " << q[0] << " (w)\n";
	std::cout << "q[1] = " << q[1] << " (x)\n";
	std::cout << "q[2] = " << q[2] << " (y)\n";
	std::cout << "q[3] = " << q[3] << " (z)\n";

	std::cout << "\n--- Scalar and Vector Parts ---\n";
	std::cout << "Scalar(): " << q.Scalar() << "\n";
	std::cout << "Vector(): " << q.Vector() << "\n";
}

///////////////////////////////////////////////////////////////////////////////////////////
///                         ARITHMETIC OPERATIONS                                       ///
///////////////////////////////////////////////////////////////////////////////////////////

void Docs_Demo_Quaternions_Arithmetic()
{
	std::cout << "\n==========================================================================\n";
	std::cout << "Demo: Quaternion Arithmetic\n";
	std::cout << "==========================================================================\n\n";

	Quaternion p(1, 2, 3, 4);
	Quaternion q(5, 6, 7, 8);

	std::cout << "p = " << p << "\n";
	std::cout << "q = " << q << "\n\n";

	std::cout << "--- Addition / Subtraction ---\n";
	std::cout << "p + q = " << (p + q) << "\n";
	std::cout << "p - q = " << (p - q) << "\n";
	std::cout << "-p = " << (-p) << "\n";

	std::cout << "\n--- Scalar Operations ---\n";
	std::cout << "p * 2 = " << (p * 2.0) << "\n";
	std::cout << "p / 2 = " << (p / 2.0) << "\n";
	std::cout << "2 * p = " << (2.0 * p) << "\n";

	std::cout << "\n--- Quaternion Multiplication (Hamilton Product) ---\n";
	std::cout << "p * q = " << (p * q) << "\n";
	std::cout << "q * p = " << (q * p) << "\n";
	std::cout << "NOTE: p*q ≠ q*p (non-commutative!)\n";
}

///////////////////////////////////////////////////////////////////////////////////////////
///                         QUATERNION-SPECIFIC OPERATIONS                              ///
///////////////////////////////////////////////////////////////////////////////////////////

void Docs_Demo_Quaternions_Operations()
{
	std::cout << "\n==========================================================================\n";
	std::cout << "Demo: Quaternion-Specific Operations\n";
	std::cout << "==========================================================================\n\n";

	Quaternion q(1, 2, 3, 4);
	std::cout << "q = " << q << "\n\n";

	std::cout << "--- Conjugate and Inverse ---\n";
	std::cout << "Conjugate(): " << q.Conjugate() << "\n";
	std::cout << "Inverse():   " << q.Inverse() << "\n";
	std::cout << "q * q^(-1) = " << (q * q.Inverse()) << " ≈ identity\n";

	std::cout << "\n--- Norm ---\n";
	std::cout << "NormSquared(): " << q.NormSquared() << " (w²+x²+y²+z² = 1+4+9+16 = 30)\n";
	std::cout << "Norm(): " << q.Norm() << " (√30 ≈ 5.477)\n";

	std::cout << "\n--- Normalization ---\n";
	Quaternion qNorm = q.Normalized();
	std::cout << "Normalized(): " << qNorm << "\n";
	std::cout << "||Normalized()|| = " << qNorm.Norm() << " (should be 1)\n";

	std::cout << "\n--- Unit Checks ---\n";
	std::cout << "q.IsUnit(): " << (q.IsUnit() ? "yes" : "no") << "\n";
	std::cout << "qNorm.IsUnit(): " << (qNorm.IsUnit() ? "yes" : "no") << "\n";
	std::cout << "Identity().IsIdentity(): " << (Quaternion::Identity().IsIdentity() ? "yes" : "no") << "\n";

	std::cout << "\n--- Dot Product ---\n";
	Quaternion p = Quaternion::Identity();
	auto r = Quaternion::FromAxisAngle(Vec3Cart(0, 0, 1), Constants::PI / 4);
	std::cout << "p = Identity()\n";
	std::cout << "r = 45° rotation around Z\n";
	std::cout << "p.Dot(r) = " << p.Dot(r) << " (cos of half-angle ≈ 0.924)\n";
}

///////////////////////////////////////////////////////////////////////////////////////////
///                         ROTATION OPERATIONS                                         ///
///////////////////////////////////////////////////////////////////////////////////////////

void Docs_Demo_Quaternions_Rotation()
{
	std::cout << "\n==========================================================================\n";
	std::cout << "Demo: Rotation Operations\n";
	std::cout << "==========================================================================\n\n";

	std::cout << "--- Rotating Vectors ---\n";
	// 90° rotation around Z-axis
	auto q = Quaternion::FromAxisAngle(Vec3Cart(0, 0, 1), Constants::PI / 2);
	std::cout << "Rotation: 90° around Z-axis\n";
	std::cout << "Quaternion: " << q << "\n";

	Vec3Cart vx(1, 0, 0);
	Vec3Cart vy(0, 1, 0);
	Vec3Cart vz(0, 0, 1);

	std::cout << "\nX-axis " << vx << " → " << q.Rotate(vx) << " (expect Y-axis)\n";
	std::cout << "Y-axis " << vy << " → " << q.Rotate(vy) << " (expect -X-axis)\n";
	std::cout << "Z-axis " << vz << " → " << q.Rotate(vz) << " (expect Z-axis)\n";

	std::cout << "\n--- Extracting Rotation Parameters ---\n";
	auto q2 = Quaternion::FromAxisAngle(Vec3Cart(1, 0, 0), Constants::PI / 3);
	std::cout << "Created: 60° around X-axis\n";
	std::cout << "GetRotationAxis(): " << q2.GetRotationAxis() << "\n";
	std::cout << "GetRotationAngle(): " << q2.GetRotationAngle() << " rad (≈ " 
	          << q2.GetRotationAngle() * 180.0 / Constants::PI << "°)\n";

	Vec3Cart axis;
	Real angle;
	q2.ToAxisAngle(axis, angle);
	std::cout << "ToAxisAngle(): axis=" << axis << ", angle=" << angle << " rad\n";

	std::cout << "\n--- Euler Angles Conversion ---\n";
	auto qEuler = Quaternion::FromEulerZYX(0.5, 0.2, 0.1);
	std::cout << "FromEulerZYX(0.5, 0.2, 0.1): " << qEuler << "\n";
	Vec3Cart euler = qEuler.ToEulerZYX();
	std::cout << "ToEulerZYX(): [" << euler[0] << ", " << euler[1] << ", " << euler[2] << "]\n";
	std::cout << "  (yaw=" << euler[0] << ", pitch=" << euler[1] << ", roll=" << euler[2] << ")\n";

	std::cout << "\n--- Rotation Matrix Conversion ---\n";
	auto qMat = Quaternion::FromAxisAngle(Vec3Cart(0, 1, 0), Constants::PI / 4);
	auto rotMat = qMat.ToRotationMatrix();
	std::cout << "45° rotation around Y-axis:\n";
	std::cout << rotMat << "\n";

	// Round-trip test
	auto qBack = Quaternion::FromRotationMatrix(rotMat);
	std::cout << "Round-trip: original " << qMat << "\n";
	std::cout << "           recovered " << qBack << "\n";
}

///////////////////////////////////////////////////////////////////////////////////////////
///                         INTERPOLATION                                               ///
///////////////////////////////////////////////////////////////////////////////////////////

void Docs_Demo_Quaternions_Interpolation()
{
	std::cout << "\n==========================================================================\n";
	std::cout << "Demo: Quaternion Interpolation (SLERP)\n";
	std::cout << "==========================================================================\n\n";

	// Two orientations
	auto q1 = Quaternion::Identity();
	auto q2 = Quaternion::FromAxisAngle(Vec3Cart(0, 0, 1), Constants::PI);

	std::cout << "Start: " << q1 << " (identity)\n";
	std::cout << "End:   " << q2 << " (180° around Z)\n\n";

	std::cout << "--- SLERP Interpolation ---\n";
	std::cout << "t\tQuaternion\t\t\t\tAngle (deg)\n";
	std::cout << "--------------------------------------------------------------\n";
	for (Real t = 0; t <= 1.0; t += 0.25) {
		auto q = Quaternion::Slerp(q1, q2, t);
		Real angleDeg = q.GetRotationAngle() * 180.0 / Constants::PI;
		std::cout << std::fixed << std::setprecision(2) << t << "\t" << q 
		          << "\t" << angleDeg << "°\n";
	}

	std::cout << "\n--- LERP vs SLERP ---\n";
	std::cout << "SLERP: Constant angular velocity, always unit length\n";
	std::cout << "LERP: Faster but variable speed, needs normalization\n\n";

	auto qLerp = Quaternion::Lerp(q1, q2, 0.5);
	auto qSlerp = Quaternion::Slerp(q1, q2, 0.5);
	std::cout << "At t=0.5:\n";
	std::cout << "  LERP:  " << qLerp << " ||q|| = " << qLerp.Norm() << "\n";
	std::cout << "  SLERP: " << qSlerp << " ||q|| = " << qSlerp.Norm() << "\n";
}

///////////////////////////////////////////////////////////////////////////////////////////
///                         PRACTICAL APPLICATIONS                                      ///
///////////////////////////////////////////////////////////////////////////////////////////

void Docs_Demo_Quaternions_Applications()
{
	std::cout << "\n==========================================================================\n";
	std::cout << "Demo: Practical Applications\n";
	std::cout << "==========================================================================\n\n";

	std::cout << "--- Combining Rotations ---\n";
	// Roll, then pitch, then yaw
	auto roll = Quaternion::FromAxisAngle(Vec3Cart(1, 0, 0), 0.1);
	auto pitch = Quaternion::FromAxisAngle(Vec3Cart(0, 1, 0), 0.2);
	auto yaw = Quaternion::FromAxisAngle(Vec3Cart(0, 0, 1), 0.3);

	// Combined rotation (right-to-left application)
	auto combined = yaw * pitch * roll;
	std::cout << "roll (X, 0.1 rad):  " << roll << "\n";
	std::cout << "pitch (Y, 0.2 rad): " << pitch << "\n";
	std::cout << "yaw (Z, 0.3 rad):   " << yaw << "\n";
	std::cout << "Combined (yaw*pitch*roll): " << combined << "\n";

	Vec3Cart v(1, 0, 0);
	std::cout << "Rotating X-axis: " << v << " → " << combined.Rotate(v) << "\n";

	std::cout << "\n--- Avoiding Gimbal Lock ---\n";
	std::cout << "Euler angles have gimbal lock at pitch = ±90°.\n";
	std::cout << "Quaternions handle all orientations smoothly.\n\n";

	// Near gimbal lock with Euler
	auto qNearLock = Quaternion::FromEulerZYX(0.5, Constants::PI/2 - 0.01, 0.3);
	std::cout << "Near gimbal lock (pitch ≈ 90°):\n";
	std::cout << "  Quaternion: " << qNearLock << "\n";
	std::cout << "  Can still rotate smoothly!\n";

	std::cout << "\n--- Incremental Rotations ---\n";
	Quaternion orientation = Quaternion::Identity();
	auto smallRotation = Quaternion::FromAxisAngle(Vec3Cart(0, 1, 0), 0.01);
	
	std::cout << "Applying small rotations incrementally:\n";
	for (int i = 0; i < 5; i++) {
		orientation = smallRotation * orientation;
		orientation.Normalize();  // Prevent drift
		Real angle = orientation.GetRotationAngle() * 180.0 / Constants::PI;
		std::cout << "  Step " << (i+1) << ": angle = " << std::fixed << std::setprecision(2) 
		          << angle << "°\n";
	}
}

///////////////////////////////////////////////////////////////////////////////////////////
///                         COMPARISON AND EQUALITY                                     ///
///////////////////////////////////////////////////////////////////////////////////////////

void Docs_Demo_Quaternions_Comparison()
{
	std::cout << "\n==========================================================================\n";
	std::cout << "Demo: Quaternion Comparison\n";
	std::cout << "==========================================================================\n\n";

	Quaternion p(1, 0, 0, 0);
	Quaternion q(1, 0, 0, 0);
	Quaternion r(1, 1e-12, 0, 0);

	std::cout << "p = " << p << "\n";
	std::cout << "q = " << q << "\n";
	std::cout << "r = " << r << " (tiny perturbation)\n\n";

	std::cout << "p == q: " << (p == q ? "true" : "false") << "\n";
	std::cout << "p == r: " << (p == r ? "true" : "false") << "\n";
	std::cout << "p.IsApprox(r, 1e-10): " << (p.IsApprox(r, 1e-10) ? "true" : "false") << "\n";

	std::cout << "\n--- Note: q and -q represent same rotation ---\n";
	Quaternion qPos = Quaternion::FromAxisAngle(Vec3Cart(0, 0, 1), Constants::PI / 2);
	Quaternion qNeg = -qPos;
	Vec3Cart v(1, 0, 0);
	std::cout << "q = " << qPos << "\n";
	std::cout << "-q = " << qNeg << "\n";
	std::cout << "q.Rotate(v) = " << qPos.Rotate(v) << "\n";
	std::cout << "(-q).Rotate(v) = " << qNeg.Rotate(v) << "\n";
	std::cout << "(Same result - they represent the same rotation!)\n";
}

///////////////////////////////////////////////////////////////////////////////////////////
///                         MAIN DEMO FUNCTION                                          ///
///////////////////////////////////////////////////////////////////////////////////////////

void Docs_Demo_Quaternions()
{
	std::cout << "\n##########################################################################\n";
	std::cout << "#                    QUATERNION DEMOS                                    #\n";
	std::cout << "##########################################################################\n";
	
	Docs_Demo_Quaternions_Constructors();
	Docs_Demo_Quaternions_Accessors();
	Docs_Demo_Quaternions_Arithmetic();
	Docs_Demo_Quaternions_Operations();
	Docs_Demo_Quaternions_Rotation();
	Docs_Demo_Quaternions_Interpolation();
	Docs_Demo_Quaternions_Applications();
	Docs_Demo_Quaternions_Comparison();
	
	std::cout << "\n=== All Quaternion Demos Complete ===\n";
}
