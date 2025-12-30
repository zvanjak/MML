///////////////////////////////////////////////////////////////////////////////////////////
///                         MinimalMathLibrary (MML)                                  ///
///                                                                                   ///
///  File:        docs_demo_vectorn.cpp                                               ///
///  Description: Documentation examples for VectorN<T,N> class                       ///
///               These examples are referenced from docs/base/VectorN.md             ///
///                                                                                   ///
///  Usage:       Run MML_DocsDemo application to execute these examples              ///
///////////////////////////////////////////////////////////////////////////////////////////

#ifdef MML_USE_SINGLE_HEADER
#include "MML.h"
#else
#include "MMLBase.h"

#include "base/VectorN.h"
#include "base/BaseUtils.h"
#endif

#include <iostream>
#include <iomanip>

using namespace MML;

///////////////////////////////////////////////////////////////////////////////////////////
///                           Construction Examples                                     ///
///////////////////////////////////////////////////////////////////////////////////////////

void Demo_VectorN_Construction()
{
    std::cout << "\n=== VectorN Construction Examples ===\n\n";

    // Default construction (zero-initialized)
    VectorN<double, 3> v1;
    std::cout << "Default VectorN<double,3>(): " << v1 << std::endl;

    // Value construction (all same value)
    VectorN<double, 3> v2(5.0);
    std::cout << "VectorN<double,3>(5.0):      " << v2 << std::endl;

    // From initializer list
    VectorN<double, 3> v3{1.0, 2.0, 3.0};
    std::cout << "VectorN{1, 2, 3}:            " << v3 << std::endl;

    // From std::vector (takes first N elements)
    std::vector<double> data = {10, 20, 30, 40, 50};
    VectorN<double, 3> v4(data);
    std::cout << "From std::vector{10..50}:    " << v4 << " (first 3 elements)" << std::endl;

    // From C array
    double arr[] = {100.0, 200.0, 300.0};
    VectorN<double, 3> v5(arr);
    std::cout << "From C array:                " << v5 << std::endl;

    // Unit vectors
    std::cout << "\nUnit vectors:" << std::endl;
    auto e0 = VectorN<double, 3>::GetUnitVector(0);
    auto e1 = VectorN<double, 3>::GetUnitVector(1);
    auto e2 = VectorN<double, 3>::GetUnitVector(2);
    std::cout << "  e0 = GetUnitVector(0): " << e0 << std::endl;
    std::cout << "  e1 = GetUnitVector(1): " << e1 << std::endl;
    std::cout << "  e2 = GetUnitVector(2): " << e2 << std::endl;

    // Different sizes
    std::cout << "\nDifferent sizes:" << std::endl;
    VectorN<double, 2> v2d{1, 2};
    VectorN<double, 4> v4d{1, 2, 3, 4};
    VectorN<double, 5> v5d{1, 2, 3, 4, 5};
    std::cout << "  2D: " << v2d << std::endl;
    std::cout << "  4D: " << v4d << std::endl;
    std::cout << "  5D: " << v5d << std::endl;
}

///////////////////////////////////////////////////////////////////////////////////////////
///                        Access & Properties Examples                                 ///
///////////////////////////////////////////////////////////////////////////////////////////

void Demo_VectorN_Access()
{
    std::cout << "\n=== VectorN Access & Properties Examples ===\n\n";

    VectorN<double, 4> v{1.5, 2.5, 3.5, 4.5};
    std::cout << "Vector: " << v << std::endl;

    // Element access with operator[]
    std::cout << "\nElement access with operator[]:" << std::endl;
    std::cout << "  v[0] = " << v[0] << std::endl;
    std::cout << "  v[1] = " << v[1] << std::endl;
    std::cout << "  v[2] = " << v[2] << std::endl;
    std::cout << "  v[3] = " << v[3] << std::endl;

    // Modify element
    v[2] = 99.0;
    std::cout << "\nAfter v[2] = 99.0: " << v << std::endl;
    v[2] = 3.5;  // Restore

    // Safe access with at() - bounds checked
    std::cout << "\nSafe access with at():" << std::endl;
    std::cout << "  v.at(1) = " << v.at(1) << std::endl;
    
    try {
        double bad = v.at(10);  // Out of bounds!
        (void)bad;
    } catch (const VectorDimensionError& e) {
        std::cout << "  Caught exception for v.at(10): " << e.what() << std::endl;
    }

    // Size (compile-time constant)
    std::cout << "\nSize and state:" << std::endl;
    std::cout << "  v.size() = " << v.size() << " (always N=4)" << std::endl;

    // IsNullVec
    VectorN<double, 3> zero;
    VectorN<double, 3> nonzero{0, 0.001, 0};
    std::cout << "  VectorN<3>().IsNullVec() = " << (zero.IsNullVec() ? "true" : "false") << std::endl;
    std::cout << "  {0, 0.001, 0}.IsNullVec() = " << (nonzero.IsNullVec() ? "true" : "false") << std::endl;

    // Clear
    VectorN<double, 3> toClear{1, 2, 3};
    std::cout << "\nBefore clear(): " << toClear << std::endl;
    toClear.clear();
    std::cout << "After clear():  " << toClear << std::endl;
}

///////////////////////////////////////////////////////////////////////////////////////////
///                        Arithmetic Operations Examples                               ///
///////////////////////////////////////////////////////////////////////////////////////////

void Demo_VectorN_Arithmetic()
{
    std::cout << "\n=== VectorN Arithmetic Operations ===\n\n";

    VectorN<double, 3> a{1, 2, 3};
    VectorN<double, 3> b{4, 5, 6};
    
    std::cout << "a = " << a << std::endl;
    std::cout << "b = " << b << std::endl;

    // Unary minus
    auto negA = -a;
    std::cout << "\n-a = " << negA << std::endl;

    // Addition and subtraction
    auto sum = a + b;
    auto diff = a - b;
    std::cout << "\na + b = " << sum << std::endl;
    std::cout << "a - b = " << diff << std::endl;

    // Compound assignment
    VectorN<double, 3> c = a;
    c += b;
    std::cout << "\nAfter c = a; c += b: " << c << std::endl;
    c -= b;
    std::cout << "After c -= b:        " << c << std::endl;

    // Scalar multiplication and division
    auto scaled1 = a * 2.0;
    auto scaled2 = 3.0 * a;
    auto divided = a / 2.0;
    
    std::cout << "\na * 2.0 = " << scaled1 << std::endl;
    std::cout << "3.0 * a = " << scaled2 << std::endl;
    std::cout << "a / 2.0 = " << divided << std::endl;

    // Compound scalar operations
    VectorN<double, 3> d = a;
    d *= 10.0;
    std::cout << "\nAfter d = a; d *= 10: " << d << std::endl;
    d /= 5.0;
    std::cout << "After d /= 5:         " << d << std::endl;
}

///////////////////////////////////////////////////////////////////////////////////////////
///                        Vector Operations Examples                                   ///
///////////////////////////////////////////////////////////////////////////////////////////

void Demo_VectorN_Operations()
{
    std::cout << "\n=== VectorN Operations ===\n\n";

    // Norms
    VectorN<double, 2> v{3, 4};
    std::cout << "Vector for norms: " << v << std::endl;
    
    double normL2 = v.NormL2();      // sqrt(9 + 16) = 5
    double normL1 = v.NormL1();      // |3| + |4| = 7
    double normLinf = v.NormLInf();  // max(|3|, |4|) = 4
    
    std::cout << "  NormL2 (Euclidean): " << normL2 << " (sqrt(3² + 4²) = 5)" << std::endl;
    std::cout << "  NormL1 (Manhattan): " << normL1 << " (|3| + |4| = 7)" << std::endl;
    std::cout << "  NormLInf (Max):     " << normLinf << " (max(|3|, |4|) = 4)" << std::endl;

    // Normalization - VectorN HAS Normalized() method!
    std::cout << "\nNormalization (VectorN has Normalized method):" << std::endl;
    VectorN<double, 2> toNorm{3, 4};
    auto unit = toNorm.Normalized();
    
    std::cout << "  Original:   " << toNorm << " (length = " << toNorm.NormL2() << ")" << std::endl;
    std::cout << "  Normalized: " << unit << " (length = " << unit.NormL2() << ")" << std::endl;

    // Dot product (scalar product) using Utils
    std::cout << "\nDot product (using Utils::ScalarProduct):" << std::endl;
    Vec3 a{1, 2, 3};
    Vec3 b{4, 5, 6};
    
    Real dot = Utils::ScalarProduct(a, b);  // 1*4 + 2*5 + 3*6 = 32
    std::cout << "  a = " << a << std::endl;
    std::cout << "  b = " << b << std::endl;
    std::cout << "  Utils::ScalarProduct(a, b) = " << dot << " (1*4 + 2*5 + 3*6 = 32)" << std::endl;

    // Angle between vectors
    std::cout << "\nAngle between vectors:" << std::endl;
    Vec3 v1{1, 0, 0};
    Vec3 v2{1, 1, 0};
    
    Real angle = Utils::VectorsAngle(v1, v2);
    std::cout << "  v1 = " << v1 << std::endl;
    std::cout << "  v2 = " << v2 << std::endl;
    std::cout << "  Angle = " << angle << " rad = " << (angle * 180.0 / Constants::PI) << " degrees" << std::endl;

    // Distance calculation
    std::cout << "\nDistance calculation:" << std::endl;
    Vec3 p1{1, 2, 3};
    Vec3 p2{4, 6, 3};
    
    double dist = (p1 - p2).NormL2();
    std::cout << "  p1 = " << p1 << std::endl;
    std::cout << "  p2 = " << p2 << std::endl;
    std::cout << "  Distance = (p1 - p2).NormL2() = " << dist << std::endl;
}

///////////////////////////////////////////////////////////////////////////////////////////
///                        Comparison & Equality Examples                               ///
///////////////////////////////////////////////////////////////////////////////////////////

void Demo_VectorN_Equality()
{
    std::cout << "\n=== VectorN Comparison & Equality ===\n\n";

    VectorN<double, 3> a{1.0, 2.0, 3.0};
    VectorN<double, 3> b{1.0, 2.0, 3.0};
    VectorN<double, 3> c{1.0, 2.0, 3.0001};
    
    std::cout << "a = " << a << std::endl;
    std::cout << "b = " << b << std::endl;
    std::cout << "c = " << c << std::endl;

    // Exact equality
    std::cout << "\nExact equality (operator==):" << std::endl;
    std::cout << "  a == b: " << (a == b ? "true" : "false") << std::endl;
    std::cout << "  a == c: " << (a == c ? "true" : "false") << std::endl;
    std::cout << "  a != c: " << (a != c ? "true" : "false") << std::endl;

    // Approximate equality with tolerance
    std::cout << "\nApproximate equality:" << std::endl;
    
    // Static method
    std::cout << "  VectorN::AreEqual(a, c, 0.001):  " 
              << (VectorN<double,3>::AreEqual(a, c, 0.001) ? "true" : "false") << std::endl;
    std::cout << "  VectorN::AreEqual(a, c, 0.0001): " 
              << (VectorN<double,3>::AreEqual(a, c, 0.0001) ? "true" : "false") << std::endl;

    // Member method
    std::cout << "  a.IsEqualTo(c, 0.001):  " << (a.IsEqualTo(c, 0.001) ? "true" : "false") << std::endl;
    std::cout << "  a.IsEqualTo(c, 0.0001): " << (a.IsEqualTo(c, 0.0001) ? "true" : "false") << std::endl;

    // Very close values
    VectorN<double, 3> d{1.0 + 1e-10, 2.0 + 1e-10, 3.0 + 1e-10};
    std::cout << "\nVery close values:" << std::endl;
    std::cout << "  a.IsEqualTo(d, 1e-9): " << (a.IsEqualTo(d, 1e-9) ? "true" : "false") << std::endl;
}

///////////////////////////////////////////////////////////////////////////////////////////
///                           I/O and Printing Examples                                 ///
///////////////////////////////////////////////////////////////////////////////////////////

void Demo_VectorN_IO()
{
    std::cout << "\n=== VectorN I/O and Printing Examples ===\n\n";

    VectorN<double, 3> v{1.23456789, 2.718281828, 3.14159265};

    // Default stream output
    std::cout << "Default output (operator<<):" << std::endl;
    std::cout << "  " << v << std::endl;

    // Custom width and precision with Print()
    std::cout << "\nCustom formatting with Print(stream, width, precision):" << std::endl;
    std::cout << "  Print(cout, 12, 6): ";
    v.Print(std::cout, 12, 6);
    std::cout << std::endl;

    std::cout << "  Print(cout, 8, 2):  ";
    v.Print(std::cout, 8, 2);
    std::cout << std::endl;

    // Print with zero threshold
    VectorN<double, 3> small{1.5, 0.0001, 3.14};
    std::cout << "\nPrint with zero threshold:" << std::endl;
    std::cout << "  Original:                  " << small << std::endl;
    std::cout << "  Print(cout, 8, 3, 0.001):  ";
    small.Print(std::cout, 8, 3, 0.001);
    std::cout << " (values < 0.001 shown as 0)" << std::endl;

    // PrintLine with label
    std::cout << "\nPrintLine with label:" << std::endl;
    v.PrintLine(std::cout, "  My vector: ", 10, 4);

    // to_string
    std::cout << "\nString conversion with to_string():" << std::endl;
    std::string str1 = v.to_string(10, 3);
    std::string str2 = v.to_string(6, 1);
    std::cout << "  to_string(10, 3): " << str1 << std::endl;
    std::cout << "  to_string(6, 1):  " << str2 << std::endl;
}

///////////////////////////////////////////////////////////////////////////////////////////
///                           Type Aliases Examples                                     ///
///////////////////////////////////////////////////////////////////////////////////////////

void Demo_VectorN_TypeAliases()
{
    std::cout << "\n=== VectorN Type Aliases Examples ===\n\n";

    // Generic Real vectors (Vec2, Vec3, Vec4)
    std::cout << "Generic Real vectors (Vec2, Vec3, Vec4):" << std::endl;
    Vec2 v2{1.0, 2.0};
    Vec3 v3{1.0, 2.0, 3.0};
    Vec4 v4{1.0, 2.0, 3.0, 4.0};
    std::cout << "  Vec2: " << v2 << std::endl;
    std::cout << "  Vec3: " << v3 << std::endl;
    std::cout << "  Vec4: " << v4 << std::endl;

    // Double vectors
    std::cout << "\nDouble vectors (Vec2Dbl/Vec2D, Vec3Dbl/Vec3D, Vec4Dbl/Vec4D):" << std::endl;
    Vec3Dbl vd3{1.5, 2.5, 3.5};
    Vec3D   vd3_alt{1.5, 2.5, 3.5};  // Same as Vec3Dbl
    std::cout << "  Vec3Dbl: " << vd3 << std::endl;
    std::cout << "  Vec3D:   " << vd3_alt << " (same type)" << std::endl;

    // Float vectors
    std::cout << "\nFloat vectors (Vec2Flt/Vec2F, Vec3Flt/Vec3F, Vec4Flt/Vec4F):" << std::endl;
    Vec3Flt vf3{1.5f, 2.5f, 3.5f};
    Vec3F   vf3_alt{1.5f, 2.5f, 3.5f};
    std::cout << "  Vec3Flt: " << vf3 << std::endl;
    std::cout << "  Vec3F:   " << vf3_alt << " (same type)" << std::endl;

    // Complex vectors
    std::cout << "\nComplex vectors (Vec2Complex/Vec2C, Vec3Complex/Vec3C, Vec4Complex/Vec4C):" << std::endl;
    Vec3Complex vc3{Complex(1,1), Complex(2,-1), Complex(0,3)};
    Vec3C       vc3_alt{Complex(1,1), Complex(2,-1), Complex(0,3)};
    std::cout << "  Vec3Complex: " << vc3 << std::endl;
    std::cout << "  Vec3C:       " << vc3_alt << " (same type)" << std::endl;

    // Unit vector operations with type aliases
    std::cout << "\nUnit vector operations with type aliases:" << std::endl;
    Vec3 i = Vec3::GetUnitVector(0);  // [1, 0, 0]
    Vec3 j = Vec3::GetUnitVector(1);  // [0, 1, 0]
    Vec3 k = Vec3::GetUnitVector(2);  // [0, 0, 1]

    Vec3 combined = 3.0*i + 4.0*j + 0.0*k;  // [3, 4, 0]
    std::cout << "  3*i + 4*j + 0*k = " << combined << std::endl;
    std::cout << "  Length: " << combined.NormL2() << std::endl;  // 5.0
}

///////////////////////////////////////////////////////////////////////////////////////////
///                           Complex Vector Examples                                   ///
///////////////////////////////////////////////////////////////////////////////////////////

void Demo_VectorN_Complex()
{
    std::cout << "\n=== VectorN Complex Vector Examples ===\n\n";

    // Create complex vectors
    Vec3C v1{Complex(1,1), Complex(2,-1), Complex(0,3)};
    Vec3C v2{Complex(2,0), Complex(1,1), Complex(-1,1)};
    
    std::cout << "v1 = " << v1 << std::endl;
    std::cout << "v2 = " << v2 << std::endl;

    // Basic operations
    Vec3C sum = v1 + v2;
    Vec3C diff = v1 - v2;
    
    std::cout << "\nv1 + v2 = " << sum << std::endl;
    std::cout << "v1 - v2 = " << diff << std::endl;

    // Scalar multiplication (complex scalar)
    Vec3C scaled = v1 * Complex(0, 1);  // Multiply by i
    std::cout << "\nv1 * i = " << scaled << std::endl;

    // Complex dot product (conjugate) using Utils
    Complex dot = Utils::ScalarProduct(v1, v2);
    std::cout << "\nConjugate dot product:" << std::endl;
    std::cout << "  Utils::ScalarProduct(v1, v2) = " << dot << std::endl;
}

///////////////////////////////////////////////////////////////////////////////////////////
///                           Physics Example                                           ///
///////////////////////////////////////////////////////////////////////////////////////////

void Demo_VectorN_Physics()
{
    std::cout << "\n=== VectorN Physics Example ===\n\n";

    // Position and velocity simulation
    Vec3 position{0.0, 0.0, 0.0};
    Vec3 velocity{1.0, 0.5, -0.2};
    Real dt = 0.1;  // Time step

    std::cout << "Simple kinematics simulation:" << std::endl;
    std::cout << "  Initial position: " << position << std::endl;
    std::cout << "  Velocity:         " << velocity << std::endl;
    std::cout << "  Time step:        " << dt << std::endl;

    std::cout << "\nSimulation:" << std::endl;
    for (int step = 0; step <= 5; step++) {
        std::cout << "  t=" << std::fixed << std::setprecision(1) << step * dt 
                  << ": pos=" << position << std::endl;
        position += velocity * dt;
    }

    // 3D rotation concept (simplified)
    std::cout << "\nUnit vector basis:" << std::endl;
    Vec3 ex = Vec3::GetUnitVector(0);
    Vec3 ey = Vec3::GetUnitVector(1);
    Vec3 ez = Vec3::GetUnitVector(2);
    
    std::cout << "  ex (right):   " << ex << std::endl;
    std::cout << "  ey (up):      " << ey << std::endl;
    std::cout << "  ez (forward): " << ez << std::endl;
}

///////////////////////////////////////////////////////////////////////////////////////////
///                           Main Entry Point                                          ///
///////////////////////////////////////////////////////////////////////////////////////////

void Docs_Demo_VectorN()
{
    std::cout << "╔════════════════════════════════════════════════════════════════════╗" << std::endl;
    std::cout << "║              VectorN<T,N> Documentation Examples                   ║" << std::endl;
    std::cout << "║                  (from docs/base/VectorN.md)                       ║" << std::endl;
    std::cout << "╚════════════════════════════════════════════════════════════════════╝" << std::endl;

    Demo_VectorN_Construction();
    Demo_VectorN_Access();
    Demo_VectorN_Arithmetic();
    Demo_VectorN_Operations();
    Demo_VectorN_Equality();
    Demo_VectorN_IO();
    Demo_VectorN_TypeAliases();
    Demo_VectorN_Complex();
    Demo_VectorN_Physics();

    std::cout << "\n=== All VectorN Examples Complete ===\n" << std::endl;
}
