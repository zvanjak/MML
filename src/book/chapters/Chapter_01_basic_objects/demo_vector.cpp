///////////////////////////////////////////////////////////////////////////////////////////
///  File:        demo_vector.cpp
///  Description: Comprehensive Vector and VectorN demonstrations for Chapter 01
///               Shows construction, operations, norms, and practical applications
///////////////////////////////////////////////////////////////////////////////////////////

#ifdef MML_USE_SINGLE_HEADER
#include "MML.h"
#else
#include "MMLBase.h"

#include "mml/base/Vector/Vector.h"
#include "mml/base/Vector/VectorN.h"
#include "mml/base/Vector/VectorTypes.h"
#include "mml/base/BaseUtils.h"
#endif

using namespace MML;

///////////////////////////////////////////////////////////////////////////////////////////
///                        Vector Construction Methods                                  ///
///////////////////////////////////////////////////////////////////////////////////////////

void Demo_Vector_Construction()
{
    std::cout << "\n=== Vector Construction ===\n\n";

    // Dynamic Vector<T> - size determined at runtime
    std::cout << "--- Dynamic Vectors (Vector<T>) ---\n";
    
    Vector<Real> v1(5);                           // 5 elements, uninitialized
    Vector<Real> v2(3, 3.14159);                  // 3 elements, all = π
    Vector<Real> v3({ 1.5, -2.1, 0.48 });         // Initializer list
    Vector<Real> v4(v3);                          // Copy constructor
    Vector<Real> v5 = v2;                         // Assignment initialization

    std::cout << "v2 (3 elements, all π):  " << v2 << "\n";
    std::cout << "v3 (initializer list):   " << v3 << "\n";

    // From STL container and C array
    std::vector<Real> std_vec{ -1.0, 5.0, -2.0 };
    double arr[4] = { 2.0, 4.0, 6.0, 8.0 };

    Vector<Real> v6(std_vec);                     // From std::vector
    Vector<Real> v7(4, arr);                      // From C array

    std::cout << "v6 (from std::vector):   " << v6 << "\n";
    std::cout << "v7 (from C array):       " << v7 << "\n";

    // Fixed-size VectorN<T, N> - compile-time size
    std::cout << "\n--- Fixed-Size Vectors (VectorN<T, N>) ---\n";
    
    VectorN<double, 3> vn1;                       // Default (zeros)
    VectorN<double, 3> vn2(1.5);                  // All elements = 1.5
    VectorN<double, 3> vn3{1.0, 2.0, 3.0};        // Initializer list
    Vec3 vn4 = vn3;                               // Using type alias

    std::cout << "vn1 (default):           " << vn1 << "\n";
    std::cout << "vn2 (all 1.5):           " << vn2 << "\n";
    std::cout << "vn3 (1, 2, 3):           " << vn3 << "\n";

    // Complex vectors
    std::cout << "\n--- Complex Vectors ---\n";
    
    Vector<Complex> vc1({ 1.0, 2.0, 3.0 });       // Real values → Complex
    VecC vc2({ Complex(1,1), Complex(-1,2), Complex(2,-0.5) });

    std::cout << "vc1 (from reals):        " << vc1 << "\n";
    std::cout << "vc2 (complex values):    " << vc2 << "\n";
}

///////////////////////////////////////////////////////////////////////////////////////////
///                        Vector Norms and Metrics                                     ///
///////////////////////////////////////////////////////////////////////////////////////////

void Demo_Vector_Norms()
{
    std::cout << "\n=== Vector Norms ===\n\n";

    VectorN<double, 3> v{3.0, -4.0, 0.0};
    std::cout << "Vector v = " << v << "\n\n";

    // L2 norm (Euclidean) - most common in physics
    double normL2 = v.NormL2();
    std::cout << "NormL2 (Euclidean):  " << normL2 << "  (√(9+16+0) = 5)\n";

    // L1 norm (Manhattan) - used in optimization, robust statistics
    double normL1 = v.NormL1();
    std::cout << "NormL1 (Manhattan):  " << normL1 << "  (|3|+|-4|+|0| = 7)\n";

    // L∞ norm (Maximum) - worst-case analysis
    double normLinf = v.NormLInf();
    std::cout << "NormLInf (Maximum):  " << normLinf << "  (max(|3|,|4|,|0|) = 4)\n";

    // Unit vectors (normalized)
    std::cout << "\n--- Normalization ---\n";
    VectorN<double, 3> unit = v.Normalized();
    std::cout << "v.Normalized() = " << unit << "\n";
    std::cout << "  Length after: " << unit.NormL2() << " (should be 1.0)\n";

    // Practical: Physics unit vectors
    std::cout << "\n--- Physics: Direction Vectors ---\n";
    VectorN<double, 3> force{10.0, 20.0, 15.0};
    VectorN<double, 3> direction = force.Normalized();
    double magnitude = force.NormL2();
    
    std::cout << "Force vector:     " << force << "\n";
    std::cout << "Magnitude:        " << magnitude << " N\n";
    std::cout << "Direction (unit): " << direction << "\n";
}

///////////////////////////////////////////////////////////////////////////////////////////
///                        Vector Products and Angles                                   ///
///////////////////////////////////////////////////////////////////////////////////////////

void Demo_Vector_Products()
{
    std::cout << "\n=== Vector Products ===\n\n";

    Vec3 a{1.0, 0.0, 0.0};  // X-axis
    Vec3 b{0.0, 1.0, 0.0};  // Y-axis
    Vec3 c{1.0, 1.0, 0.0};  // 45° in XY plane

    std::cout << "a = " << a << " (X-axis)\n";
    std::cout << "b = " << b << " (Y-axis)\n";
    std::cout << "c = " << c << " (45° in XY)\n\n";

    // Dot product (scalar product)
    std::cout << "--- Dot Product (Scalar Product) ---\n";
    Real dot_ab = Utils::ScalarProduct(a, b);
    Real dot_ac = Utils::ScalarProduct(a, c);
    std::cout << "a · b = " << dot_ab << " (perpendicular → 0)\n";
    std::cout << "a · c = " << dot_ac << " (1·1 + 0·1 + 0·0 = 1)\n";

    // Angle between vectors
    std::cout << "\n--- Angle Between Vectors ---\n";
    Real angle_ab = Utils::VectorsAngle(a, b);
    Real angle_ac = Utils::VectorsAngle(a, c);
    std::cout << "Angle(a, b) = " << angle_ab << " rad = " 
              << (angle_ab * 180.0 / Constants::PI) << "°\n";
    std::cout << "Angle(a, c) = " << angle_ac << " rad = " 
              << (angle_ac * 180.0 / Constants::PI) << "°\n";

    // Cross product (vector product) - 3D only
    std::cout << "\n--- Cross Product (Vector Product) ---\n";
    Vector3Cartesian v1{1.0, 0.0, 0.0};
    Vector3Cartesian v2{0.0, 1.0, 0.0};
    auto cross = VectorProduct(v1, v2);
    std::cout << "X × Y = " << cross << " (Z-axis, right-hand rule)\n";

    // Triple product - volume of parallelepiped
    std::cout << "\n--- Physics: Torque Calculation ---\n";
    Vector3Cartesian r{0.5, 0.0, 0.0};      // Position vector (0.5m along X)
    Vector3Cartesian F{0.0, 10.0, 0.0};    // Force (10N along Y)
    auto torque = VectorProduct(r, F);
    std::cout << "r = " << r << " m\n";
    std::cout << "F = " << F << " N\n";
    std::cout << "τ = r × F = " << torque << " N·m\n";
}

///////////////////////////////////////////////////////////////////////////////////////////
///                        Vector Arithmetic Operations                                 ///
///////////////////////////////////////////////////////////////////////////////////////////

void Demo_Vector_Arithmetic()
{
    std::cout << "\n=== Vector Arithmetic ===\n\n";

    VectorN<double, 3> a{1.0, 2.0, 3.0};
    VectorN<double, 3> b{4.0, 5.0, 6.0};

    std::cout << "a = " << a << "\n";
    std::cout << "b = " << b << "\n\n";

    // Basic operations
    std::cout << "--- Basic Operations ---\n";
    std::cout << "a + b = " << (a + b) << "\n";
    std::cout << "a - b = " << (a - b) << "\n";
    std::cout << "-a    = " << (-a) << "\n";

    // Scalar multiplication
    std::cout << "\n--- Scalar Operations ---\n";
    std::cout << "a * 2.0 = " << (a * 2.0) << "\n";
    std::cout << "3.0 * b = " << (3.0 * b) << "\n";
    std::cout << "a / 2.0 = " << (a / 2.0) << "\n";

    // Compound assignment
    std::cout << "\n--- Compound Assignment ---\n";
    VectorN<double, 3> c = a;
    c += b;
    std::cout << "c = a; c += b → " << c << "\n";
    c *= 0.5;
    std::cout << "c *= 0.5     → " << c << "\n";

    // Distance calculation
    std::cout << "\n--- Distance Between Points ---\n";
    VectorN<double, 3> p1{1.0, 2.0, 3.0};
    VectorN<double, 3> p2{4.0, 6.0, 3.0};
    double dist = (p2 - p1).NormL2();
    std::cout << "p1 = " << p1 << "\n";
    std::cout << "p2 = " << p2 << "\n";
    std::cout << "Distance = " << dist << "\n";
}

///////////////////////////////////////////////////////////////////////////////////////////
///                        Vector Comparison and Equality                               ///
///////////////////////////////////////////////////////////////////////////////////////////

void Demo_Vector_Comparison()
{
    std::cout << "\n=== Vector Comparison ===\n\n";

    VectorN<double, 3> a{1.0, 2.0, 3.0};
    VectorN<double, 3> b{1.0, 2.0, 3.0};
    VectorN<double, 3> c{1.0, 2.0, 3.0001};

    std::cout << "a = " << a << "\n";
    std::cout << "b = " << b << "\n";
    std::cout << "c = " << c << " (slightly different)\n\n";

    // Exact equality
    std::cout << "--- Exact Equality ---\n";
    std::cout << "a == b: " << (a == b ? "true" : "false") << "\n";
    std::cout << "a == c: " << (a == c ? "true" : "false") << "\n";

    // Approximate equality (essential for numerical computing!)
    std::cout << "\n--- Approximate Equality (Tolerances) ---\n";
    std::cout << "a.IsEqualTo(c, 0.001):  " << (a.IsEqualTo(c, 0.001) ? "true" : "false") << "\n";
    std::cout << "a.IsEqualTo(c, 0.0001): " << (a.IsEqualTo(c, 0.0001) ? "true" : "false") << "\n";

    // Null vector check
    std::cout << "\n--- Special Checks ---\n";
    VectorN<double, 3> zero{0.0, 0.0, 0.0};
    std::cout << "zero.isZero(): " << (zero.isZero() ? "true" : "false") << "\n";
    std::cout << "a.isZero():    " << (a.isZero() ? "true" : "false") << "\n";
}

///////////////////////////////////////////////////////////////////////////////////////////
///                        Complex Vector Operations                                    ///
///////////////////////////////////////////////////////////////////////////////////////////

void Demo_Complex_Vectors()
{
    std::cout << "\n=== Complex Vectors ===\n\n";

    // Quantum mechanics: State vectors
    std::cout << "--- Quantum Mechanics: Qubit State ---\n";
    Complex alpha(1.0/std::sqrt(2.0), 0.0);    // |0⟩ coefficient
    Complex beta(0.0, 1.0/std::sqrt(2.0));     // |1⟩ coefficient (with phase)
    
    Vector<Complex> psi({alpha, beta});
    std::cout << "|ψ⟩ = " << psi << "\n";

    // Normalization check (probability must sum to 1)
    Complex norm_sq = Utils::ScalarProduct(psi, psi);
    std::cout << "⟨ψ|ψ⟩ = " << norm_sq << " (should be 1.0)\n";

    // Signal processing: Phasor representation
    std::cout << "\n--- Electrical Engineering: Phasors ---\n";
    double V_mag = 120.0;       // 120V RMS
    double phase = 30.0 * Constants::PI / 180.0;  // 30° phase
    Complex V_phasor = std::polar(V_mag, phase);
    
    Vector<Complex> circuit({V_phasor, V_phasor * Complex(0.8, 0.6)});
    std::cout << "Voltage phasor: " << V_phasor << "\n";
    std::cout << "Circuit voltages: " << circuit << "\n";
}

///////////////////////////////////////////////////////////////////////////////////////////
///                        I/O and Formatting                                           ///
///////////////////////////////////////////////////////////////////////////////////////////

void Demo_Vector_IO()
{
    std::cout << "\n=== Vector I/O and Formatting ===\n\n";

    VectorN<double, 3> v{1.23456789, 2.718281828, 3.14159265};
    std::cout << "v = " << v << " (default format)\n";

    // Custom precision with to_string
    Vector<Real> v_dyn({1.23456789, 2.718281828, 3.14159265});
    std::cout << "v.to_string(10, 3): " << v_dyn.to_string(10, 3) << "\n";
    std::cout << "v.to_string(12, 6): " << v_dyn.to_string(12, 6) << "\n";

    // Print method with width and precision
    std::cout << "v.Print(cout, 10, 4): ";
    v_dyn.Print(std::cout, 10, 4);
    std::cout << "\n";
}

///////////////////////////////////////////////////////////////////////////////////////////
///                        Main Demo Entry Point                                        ///
///////////////////////////////////////////////////////////////////////////////////////////

void Demo_Vector()
{
    std::cout << "\n";
    std::cout << "***********************************************************************\n";
    std::cout << "****                      VECTORS IN MML                           ****\n";
    std::cout << "****          Dynamic (Vector<T>) and Fixed (VectorN<T,N>)         ****\n";
    std::cout << "***********************************************************************\n";

    Demo_Vector_Construction();
    Demo_Vector_Norms();
    Demo_Vector_Products();
    Demo_Vector_Arithmetic();
    Demo_Vector_Comparison();
    Demo_Complex_Vectors();
    Demo_Vector_IO();
}