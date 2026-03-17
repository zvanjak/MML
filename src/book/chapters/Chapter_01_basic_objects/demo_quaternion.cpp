///////////////////////////////////////////////////////////////////////////////////////////
///  File:        demo_quaternion.cpp
///  Description: Comprehensive Quaternion demonstrations for Chapter 01
///               Construction, operations, rotations, and practical applications
///////////////////////////////////////////////////////////////////////////////////////////

#ifdef MML_USE_SINGLE_HEADER
#include "MML.h"
#else
#include "MMLBase.h"

#include "mml/base/Quaternions.h"
#include "mml/base/Vector/VectorTypes.h"
#endif

using namespace MML;

///////////////////////////////////////////////////////////////////////////////////////////
///                        Quaternion Basics                                            ///
///////////////////////////////////////////////////////////////////////////////////////////

void Demo_Quaternion_Basics()
{
    std::cout << "\n=== Quaternion Basics ===\n\n";

    std::cout << "Quaternions: q = w + xi + yj + zk\n";
    std::cout << "  - w is the scalar (real) part\n";
    std::cout << "  - (x, y, z) is the vector (imaginary) part\n";
    std::cout << "  - i² = j² = k² = ijk = -1\n\n";

    // Construction methods
    std::cout << "--- Constructors ---\n";
    
    Quaternion q1;                              // Identity (1, 0, 0, 0)
    std::cout << "Default (identity): " << q1 << "\n";

    Quaternion q2(0.707, 0.707, 0, 0);          // From components
    std::cout << "From components (0.707, 0.707, 0, 0): " << q2 << "\n";

    Quaternion q3(Vec3Cart(1, 2, 3));           // Pure imaginary (w=0)
    std::cout << "Pure imaginary from (1,2,3): " << q3 << "\n";

    // Factory methods
    std::cout << "\n--- Factory Methods ---\n";
    
    auto qId = Quaternion::Identity();
    std::cout << "Identity(): " << qId << "\n";

    Vec3Cart zAxis(0, 0, 1);
    auto qRot90Z = Quaternion::FromAxisAngle(zAxis, Constants::PI / 2);
    std::cout << "FromAxisAngle(Z, 90°): " << qRot90Z << "\n";
    std::cout << "  Expected: [cos(45°), 0, 0, sin(45°)] ≈ [0.707, 0, 0, 0.707]\n";

    auto qEuler = Quaternion::FromEulerZYX(0.5, 0.2, 0.1);
    std::cout << "FromEulerZYX(yaw=0.5, pitch=0.2, roll=0.1): " << qEuler << "\n";
}

///////////////////////////////////////////////////////////////////////////////////////////
///                        Quaternion Components                                        ///
///////////////////////////////////////////////////////////////////////////////////////////

void Demo_Quaternion_Components()
{
    std::cout << "\n=== Quaternion Components ===\n\n";

    Quaternion q(0.5, 0.5, 0.5, 0.5);
    std::cout << "q = " << q << "\n\n";

    // Named accessors
    std::cout << "--- Named Accessors ---\n";
    std::cout << "q.w() = " << q.w() << " (scalar part)\n";
    std::cout << "q.x() = " << q.x() << " (i component)\n";
    std::cout << "q.y() = " << q.y() << " (j component)\n";
    std::cout << "q.z() = " << q.z() << " (k component)\n";

    // Index access
    std::cout << "\n--- Index Access ---\n";
    std::cout << "q[0] = " << q[0] << " (w)\n";
    std::cout << "q[1] = " << q[1] << " (x)\n";
    std::cout << "q[2] = " << q[2] << " (y)\n";
    std::cout << "q[3] = " << q[3] << " (z)\n";

    // Scalar and vector parts
    std::cout << "\n--- Scalar and Vector Parts ---\n";
    std::cout << "Scalar(): " << q.Scalar() << "\n";
    std::cout << "Vector(): " << q.Vector() << "\n";
}

///////////////////////////////////////////////////////////////////////////////////////////
///                        Quaternion Arithmetic                                        ///
///////////////////////////////////////////////////////////////////////////////////////////

void Demo_Quaternion_Arithmetic()
{
    std::cout << "\n=== Quaternion Arithmetic ===\n\n";

    Quaternion p(1, 2, 3, 4);
    Quaternion q(5, 6, 7, 8);

    std::cout << "p = " << p << "\n";
    std::cout << "q = " << q << "\n\n";

    // Addition and subtraction
    std::cout << "--- Addition / Subtraction ---\n";
    std::cout << "p + q = " << (p + q) << "\n";
    std::cout << "p - q = " << (p - q) << "\n";
    std::cout << "-p    = " << (-p) << "\n";

    // Scalar operations
    std::cout << "\n--- Scalar Operations ---\n";
    std::cout << "p * 2 = " << (p * 2.0) << "\n";
    std::cout << "p / 2 = " << (p / 2.0) << "\n";
    std::cout << "2 * p = " << (2.0 * p) << "\n";

    // Quaternion multiplication (Hamilton product)
    std::cout << "\n--- Hamilton Product (Non-commutative!) ---\n";
    std::cout << "p * q = " << (p * q) << "\n";
    std::cout << "q * p = " << (q * p) << "\n";
    std::cout << "NOTE: p*q ≠ q*p (quaternion multiplication is NOT commutative)\n";
}

///////////////////////////////////////////////////////////////////////////////////////////
///                        Quaternion Operations                                        ///
///////////////////////////////////////////////////////////////////////////////////////////

void Demo_Quaternion_Operations()
{
    std::cout << "\n=== Quaternion Operations ===\n\n";

    Quaternion q(1, 2, 3, 4);
    std::cout << "q = " << q << "\n\n";

    // Conjugate and inverse
    std::cout << "--- Conjugate and Inverse ---\n";
    std::cout << "q.Conjugate(): " << q.Conjugate() << " (negates vector part)\n";
    std::cout << "q.Inverse():   " << q.Inverse() << "\n";
    std::cout << "q * q^(-1) = " << (q * q.Inverse()) << " ≈ identity\n";

    // Norm
    std::cout << "\n--- Norm ---\n";
    std::cout << "NormSquared(): " << q.NormSquared() << " (1+4+9+16 = 30)\n";
    std::cout << "Norm(): " << q.Norm() << " (√30 ≈ 5.477)\n";

    // Normalization
    std::cout << "\n--- Normalization ---\n";
    Quaternion qNorm = q.Normalized();
    std::cout << "q.Normalized(): " << qNorm << "\n";
    std::cout << "||Normalized()|| = " << qNorm.Norm() << " (should be 1)\n";

    // Unit checks
    std::cout << "\n--- Unit Checks ---\n";
    std::cout << "q.IsIdentity(): " << (q.IsIdentity() ? "yes" : "no") << "\n";
    std::cout << "qNorm.IsIdentity(): " << (qNorm.IsIdentity() ? "yes" : "no") << "\n";
    std::cout << "Identity().IsIdentity(): " << (Quaternion::Identity().IsIdentity() ? "yes" : "no") << "\n";
}

///////////////////////////////////////////////////////////////////////////////////////////
///                        3D Rotations with Quaternions                                ///
///////////////////////////////////////////////////////////////////////////////////////////

void Demo_Quaternion_Rotations()
{
    std::cout << "\n=== 3D Rotations ===\n\n";

    std::cout << "Quaternions excel at representing 3D rotations:\n";
    std::cout << "  - No gimbal lock (unlike Euler angles)\n";
    std::cout << "  - Smooth interpolation (SLERP)\n";
    std::cout << "  - Compact representation (4 numbers)\n\n";

    // 90° rotation around Z-axis
    std::cout << "--- Rotating Vectors ---\n";
    auto q = Quaternion::FromAxisAngle(Vec3Cart(0, 0, 1), Constants::PI / 2);
    std::cout << "Rotation: 90° around Z-axis\n";
    std::cout << "Quaternion: " << q << "\n\n";

    Vec3Cart vx(1, 0, 0);
    Vec3Cart vy(0, 1, 0);
    Vec3Cart vz(0, 0, 1);

    std::cout << "X-axis " << vx << " → " << q.Rotate(vx) << " (becomes Y-axis)\n";
    std::cout << "Y-axis " << vy << " → " << q.Rotate(vy) << " (becomes -X-axis)\n";
    std::cout << "Z-axis " << vz << " → " << q.Rotate(vz) << " (unchanged - rotation axis)\n";

    // Extract rotation parameters
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
}

///////////////////////////////////////////////////////////////////////////////////////////
///                        Combining Rotations                                          ///
///////////////////////////////////////////////////////////////////////////////////////////

void Demo_Quaternion_Combining()
{
    std::cout << "\n=== Combining Rotations ===\n\n";

    // Two sequential rotations
    auto rotX = Quaternion::FromAxisAngle(Vec3Cart(1, 0, 0), Constants::PI / 2);  // 90° around X
    auto rotZ = Quaternion::FromAxisAngle(Vec3Cart(0, 0, 1), Constants::PI / 2);  // 90° around Z

    std::cout << "rotX: 90° around X-axis\n";
    std::cout << "rotZ: 90° around Z-axis\n\n";

    // Combined rotation (note order matters!)
    // q1 * q2 means: apply q2 first, then q1
    auto combined1 = rotX * rotZ;  // Z first, then X
    auto combined2 = rotZ * rotX;  // X first, then Z

    std::cout << "--- Order Matters! ---\n";
    Vec3Cart test(1, 0, 0);
    std::cout << "Test vector: " << test << "\n\n";

    std::cout << "rotZ then rotX (rotX * rotZ):\n";
    std::cout << "  Result: " << combined1.Rotate(test) << "\n";

    std::cout << "\nrotX then rotZ (rotZ * rotX):\n";
    std::cout << "  Result: " << combined2.Rotate(test) << "\n";

    // Inverse rotation
    std::cout << "\n--- Inverse Rotation ---\n";
    Vec3Cart original(1, 2, 3);
    Vec3Cart rotated = rotX.Rotate(original);
    Vec3Cart back = rotX.Inverse().Rotate(rotated);

    std::cout << "Original:    " << original << "\n";
    std::cout << "After rotX:  " << rotated << "\n";
    std::cout << "After rotX⁻¹: " << back << " (back to original)\n";
}

///////////////////////////////////////////////////////////////////////////////////////////
///                        Interpolation (SLERP)                                        ///
///////////////////////////////////////////////////////////////////////////////////////////

void Demo_Quaternion_Interpolation()
{
    std::cout << "\n=== Interpolation (SLERP) ===\n\n";

    std::cout << "SLERP: Spherical Linear Interpolation\n";
    std::cout << "Smoothly interpolates between two orientations.\n\n";

    auto q_start = Quaternion::Identity();
    auto q_end = Quaternion::FromAxisAngle(Vec3Cart(0, 0, 1), Constants::PI);  // 180° around Z

    std::cout << "Start: " << q_start << " (no rotation)\n";
    std::cout << "End:   " << q_end << " (180° around Z)\n\n";

    std::cout << "Interpolated rotations applied to X-axis:\n";
    Vec3Cart xAxis(1, 0, 0);
    
    for (double t = 0.0; t <= 1.0; t += 0.25) {
        auto q_interp = Quaternion::Slerp(q_start, q_end, t);
        Vec3Cart result = q_interp.Rotate(xAxis);
        std::cout << "  t=" << t << ": " << result << "\n";
    }

    // Dot product for angle between quaternions
    std::cout << "\n--- Quaternion Dot Product ---\n";
    std::cout << "q_start.Dot(q_end) = " << q_start.Dot(q_end) << "\n";
    std::cout << "(cos of half the angle between orientations)\n";
}

///////////////////////////////////////////////////////////////////////////////////////////
///                        Practical Applications                                       ///
///////////////////////////////////////////////////////////////////////////////////////////

void Demo_Quaternion_Applications()
{
    std::cout << "\n=== Practical Applications ===\n\n";

    // Game dev: Camera orientation
    std::cout << "--- Game Dev: Camera Orientation ---\n";
    double yaw = 30.0 * Constants::PI / 180.0;    // Look left
    double pitch = -15.0 * Constants::PI / 180.0; // Look up
    
    auto camera_rot = Quaternion::FromEulerZYX(yaw, pitch, 0.0);
    Vec3Cart forward(0, 0, -1);  // Camera default forward
    Vec3Cart look_dir = camera_rot.Rotate(forward);
    
    std::cout << "Camera yaw=30°, pitch=-15° (looking left and up)\n";
    std::cout << "Look direction: " << look_dir << "\n";

    // Aerospace: Aircraft orientation
    std::cout << "\n--- Aerospace: Aircraft Orientation ---\n";
    double roll = 15.0 * Constants::PI / 180.0;
    double aircraft_pitch = 5.0 * Constants::PI / 180.0;
    double heading = 45.0 * Constants::PI / 180.0;
    
    auto aircraft_rot = Quaternion::FromEulerZYX(heading, aircraft_pitch, roll);
    Vec3Cart nose(1, 0, 0);  // Aircraft nose direction
    Vec3Cart nose_dir = aircraft_rot.Rotate(nose);
    
    std::cout << "Aircraft: roll=15°, pitch=5°, heading=45°\n";
    std::cout << "Nose direction: " << nose_dir << "\n";

    // Robotics: End effector orientation
    std::cout << "\n--- Robotics: Smooth Trajectory ---\n";
    auto start_orientation = Quaternion::Identity();
    auto target_orientation = Quaternion::FromAxisAngle(Vec3Cart(0, 1, 0), Constants::PI / 4);
    
    std::cout << "Interpolating robot arm from start to target (45° tilt):\n";
    for (int step = 0; step <= 4; step++) {
        double t = step / 4.0;
        auto current = Quaternion::Slerp(start_orientation, target_orientation, t);
        Vec3Cart tool_dir = current.Rotate(Vec3Cart(0, 0, 1));
        std::cout << "  Step " << step << " (t=" << t << "): tool points " << tool_dir << "\n";
    }
}

///////////////////////////////////////////////////////////////////////////////////////////
///                        Main Demo Entry Point                                        ///
///////////////////////////////////////////////////////////////////////////////////////////

void Demo_Quaternion()
{
    std::cout << "\n";
    std::cout << "***********************************************************************\n";
    std::cout << "****                    QUATERNIONS IN MML                         ****\n";
    std::cout << "****           3D Rotations Without Gimbal Lock                    ****\n";
    std::cout << "***********************************************************************\n";

    Demo_Quaternion_Basics();
    Demo_Quaternion_Components();
    Demo_Quaternion_Arithmetic();
    Demo_Quaternion_Operations();
    Demo_Quaternion_Rotations();
    Demo_Quaternion_Combining();
    Demo_Quaternion_Interpolation();
    Demo_Quaternion_Applications();
}
