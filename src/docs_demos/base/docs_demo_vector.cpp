///////////////////////////////////////////////////////////////////////////////////////////
///                         MinimalMathLibrary (MML)                                  ///
///                                                                                   ///
///  File:        docs_demo_vector.cpp                                                ///
///  Description: Documentation examples for Vector<T> class                          ///
///               These examples are referenced from docs/base/Vector.md              ///
///                                                                                   ///
///  Usage:       Run MML_DocsDemo application to execute these examples              ///
///////////////////////////////////////////////////////////////////////////////////////////

#ifdef MML_USE_SINGLE_HEADER
#include "MML.h"
#else
#include "MMLBase.h"

#include "base/Vector.h"
#include "base/VectorN.h"
#include "base/VectorTypes.h"
#include "base/BaseUtils.h"
#endif

#include <iostream>
#include <iomanip>

using namespace MML;

///////////////////////////////////////////////////////////////////////////////////////////
///                           Construction Examples                                     ///
///////////////////////////////////////////////////////////////////////////////////////////

void Demo_Vector_Construction()
{
    std::cout << "\n=== Vector Construction Examples ===\n\n";

    // Empty vector
    Vector<double> empty;
    std::cout << "Empty vector size: " << empty.size() << std::endl;

    // Size-based construction
    Vector<double> zeros(5);            // Size 5, zero-initialized
    Vector<double> filled(5, 3.14);     // Size 5, all values = 3.14
    
    std::cout << "Zeros(5):  " << zeros << std::endl;
    std::cout << "Filled(5): " << filled << std::endl;

    // From C array
    double data[] = {1.0, 2.0, 3.0, 4.0, 5.0};
    Vector<double> fromArray(5, data);
    std::cout << "From array: " << fromArray << std::endl;

    // From std::vector
    std::vector<double> stdvec = {10.0, 20.0, 30.0};
    Vector<double> fromStdVec(stdvec);
    std::cout << "From std::vector: " << fromStdVec << std::endl;

    // From initializer list
    Vector<double> fromInit({1.1, 2.2, 3.3, 4.4});
    std::cout << "From initializer list: " << fromInit << std::endl;

    // Unit vectors
    Vector<double> e1 = Vector<double>::GetUnitVector(4, 0);  // [1, 0, 0, 0]
    Vector<double> e2 = Vector<double>::GetUnitVector(4, 1);  // [0, 1, 0, 0]
    Vector<double> e3 = Vector<double>::GetUnitVector(4, 2);  // [0, 0, 1, 0]
    Vector<double> e4 = Vector<double>::GetUnitVector(4, 3);  // [0, 0, 0, 1]
    
    std::cout << "e1 = " << e1 << std::endl;
    std::cout << "e2 = " << e2 << std::endl;
    std::cout << "e3 = " << e3 << std::endl;
    std::cout << "e4 = " << e4 << std::endl;

    // Copy and move semantics
    Vector<double> original({1, 2, 3, 4, 5});
    Vector<double> copied(original);                // Copy constructor
    Vector<double> moved = std::move(copied);       // Move constructor
    
    std::cout << "Original: " << original << std::endl;
    std::cout << "Moved:    " << moved << std::endl;
    std::cout << "Copied size after move: " << copied.size() << " (moved-from state)" << std::endl;
}

///////////////////////////////////////////////////////////////////////////////////////////
///                        Access & Iteration Examples                                  ///
///////////////////////////////////////////////////////////////////////////////////////////

void Demo_Vector_Access()
{
    std::cout << "\n=== Vector Access & Iteration Examples ===\n\n";

    Vector<double> v({1.5, 2.5, 3.5, 4.5, 5.5});
    std::cout << "Vector: " << v << std::endl;

    // Element access - operator[]
    std::cout << "\nElement access with operator[]:" << std::endl;
    std::cout << "  v[0] = " << v[0] << std::endl;
    std::cout << "  v[2] = " << v[2] << std::endl;
    std::cout << "  v[4] = " << v[4] << std::endl;

    // Element access - at() with bounds checking
    std::cout << "\nElement access with at() (bounds-checked):" << std::endl;
    std::cout << "  v.at(1) = " << v.at(1) << std::endl;
    
    // Demonstrate exception handling
    try {
        double bad = v.at(10);  // Out of bounds!
        (void)bad;  // Suppress unused warning
    } catch (const VectorDimensionError& e) {
        std::cout << "  Caught exception for v.at(10): " << e.what() << std::endl;
    }

    // Front and back
    std::cout << "\nFront and back:" << std::endl;
    std::cout << "  v.front() = " << v.front() << std::endl;
    std::cout << "  v.back()  = " << v.back() << std::endl;

    // Size and state
    std::cout << "\nSize and state:" << std::endl;
    std::cout << "  v.size()    = " << v.size() << std::endl;
    std::cout << "  v.isEmpty() = " << (v.isEmpty() ? "true" : "false") << std::endl;

    // Range-based iteration
    std::cout << "\nRange-based iteration:" << std::endl;
    std::cout << "  Elements: ";
    for (const double& x : v) {
        std::cout << x << " ";
    }
    std::cout << std::endl;

    // Modifying through iteration
    Vector<double> toModify({1, 2, 3, 4, 5});
    std::cout << "\nBefore modification: " << toModify << std::endl;
    for (double& x : toModify) {
        x *= 2;  // Double each element
    }
    std::cout << "After x *= 2:        " << toModify << std::endl;

    // Iterator-based access
    std::cout << "\nIterator-based access:" << std::endl;
    std::cout << "  Sum via iterators: ";
    double sum = 0.0;
    for (auto it = v.begin(); it != v.end(); ++it) {
        sum += *it;
    }
    std::cout << sum << std::endl;
}

///////////////////////////////////////////////////////////////////////////////////////////
///                        Arithmetic Operations Examples                               ///
///////////////////////////////////////////////////////////////////////////////////////////

void Demo_Vector_Arithmetic()
{
    std::cout << "\n=== Vector Arithmetic Operations ===\n\n";

    Vector<double> a({1, 2, 3, 4, 5});
    Vector<double> b({5, 4, 3, 2, 1});
    
    std::cout << "a = " << a << std::endl;
    std::cout << "b = " << b << std::endl;

    // Unary minus
    Vector<double> negA = -a;
    std::cout << "\n-a = " << negA << std::endl;

    // Addition and subtraction
    Vector<double> sum = a + b;
    Vector<double> diff = a - b;
    std::cout << "\na + b = " << sum << std::endl;
    std::cout << "a - b = " << diff << std::endl;

    // Compound assignment
    Vector<double> c = a;
    c += b;
    std::cout << "\nAfter c = a; c += b: " << c << std::endl;
    c -= b;
    std::cout << "After c -= b:        " << c << std::endl;

    // Scalar multiplication and division
    Vector<double> scaled1 = a * 2.0;
    Vector<double> scaled2 = 3.0 * a;
    Vector<double> divided = a / 2.0;
    
    std::cout << "\na * 2.0 = " << scaled1 << std::endl;
    std::cout << "3.0 * a = " << scaled2 << std::endl;
    std::cout << "a / 2.0 = " << divided << std::endl;

    // Compound scalar operations
    Vector<double> d = a;
    d *= 10.0;
    std::cout << "\nAfter d = a; d *= 10: " << d << std::endl;
    d /= 5.0;
    std::cout << "After d /= 5:         " << d << std::endl;
}

///////////////////////////////////////////////////////////////////////////////////////////
///                        Vector Operations Examples                                   ///
///////////////////////////////////////////////////////////////////////////////////////////

void Demo_Vector_Operations()
{
    std::cout << "\n=== Vector Operations ===\n\n";

    // Dot product (scalar product)
    Vector<double> a({1, 2, 3});
    Vector<double> b({4, 5, 6});
    
    std::cout << "a = " << a << std::endl;
    std::cout << "b = " << b << std::endl;
    
    // Use free function Utils::ScalarProduct
    double dot = Utils::ScalarProduct(a, b);  // 1*4 + 2*5 + 3*6 = 32
    std::cout << "\nDot product: Utils::ScalarProduct(a, b) = " << dot << std::endl;
    std::cout << "  (computed as 1*4 + 2*5 + 3*6 = 4 + 10 + 18 = 32)" << std::endl;

    // Norms
    Vector<double> v({3, 4});
    std::cout << "\nVector for norms: " << v << std::endl;
    
    double normL2 = v.NormL2();      // sqrt(9 + 16) = 5
    double normL1 = v.NormL1();      // |3| + |4| = 7
    double normLinf = v.NormLInf();  // max(|3|, |4|) = 4
    
    std::cout << "  NormL2 (Euclidean):   " << normL2 << " (sqrt(3² + 4²) = 5)" << std::endl;
    std::cout << "  NormL1 (Manhattan):   " << normL1 << " (|3| + |4| = 7)" << std::endl;
    std::cout << "  NormLInf (Max abs):   " << normLinf << " (max(|3|, |4|) = 4)" << std::endl;

    // Manual normalization (Vector<T> has no Normalize method)
    std::cout << "\nManual normalization:" << std::endl;
    Vector<double> toNorm({3, 4});
    Real len = toNorm.NormL2();
    Vector<double> unit = toNorm / len;
    
    std::cout << "  Original:   " << toNorm << " (length = " << len << ")" << std::endl;
    std::cout << "  Normalized: " << unit << " (length = " << unit.NormL2() << ")" << std::endl;

    // Distance calculation (manual)
    std::cout << "\nDistance calculation:" << std::endl;
    Vector<double> p1({1, 2});
    Vector<double> p2({4, 6});
    
    double dist = (p1 - p2).NormL2();  // sqrt((4-1)² + (6-2)²) = 5
    std::cout << "  p1 = " << p1 << std::endl;
    std::cout << "  p2 = " << p2 << std::endl;
    std::cout << "  Distance = (p1 - p2).NormL2() = " << dist << std::endl;

    // Cross product - note: requires Vector3Cartesian, not Vector<T>
    std::cout << "\nCross product (using Vector3Cartesian):" << std::endl;
    Vector3Cartesian i(1, 0, 0);
    Vector3Cartesian j(0, 1, 0);
    Vector3Cartesian k = VectorProduct(i, j);
    
    std::cout << "  i = " << i << std::endl;
    std::cout << "  j = " << j << std::endl;
    std::cout << "  i × j = VectorProduct(i, j) = " << k << std::endl;

    // Angle between vectors
    std::cout << "\nAngle between vectors:" << std::endl;
    Vector<double> v1({1, 0, 0});
    Vector<double> v2({1, 1, 0});
    
    Real angle = Utils::VectorsAngle(v1, v2);
    std::cout << "  v1 = " << v1 << std::endl;
    std::cout << "  v2 = " << v2 << std::endl;
    std::cout << "  Angle = " << angle << " rad = " << (angle * 180.0 / Constants::PI) << " degrees" << std::endl;
}

///////////////////////////////////////////////////////////////////////////////////////////
///                           Modification Examples                                     ///
///////////////////////////////////////////////////////////////////////////////////////////

void Demo_Vector_Modification()
{
    std::cout << "\n=== Vector Modification Examples ===\n\n";

    // Dynamic building with push_back
    std::cout << "Building vector dynamically:" << std::endl;
    Vector<double> dynamic;
    std::cout << "  Initial size: " << dynamic.size() << std::endl;
    
    for (int i = 1; i <= 5; i++) {
        dynamic.push_back(i * 1.5);
        std::cout << "  After push_back(" << i * 1.5 << "): " << dynamic << std::endl;
    }

    // Insert
    std::cout << "\nInsert operations:" << std::endl;
    Vector<double> v({10, 20, 30, 40, 50});
    std::cout << "  Original: " << v << std::endl;
    
    v.insert(2, 25);  // Insert 25 at position 2
    std::cout << "  After insert(2, 25): " << v << std::endl;
    
    v.insert(0, 5);   // Insert at beginning
    std::cout << "  After insert(0, 5):  " << v << std::endl;

    // Erase
    std::cout << "\nErase operations:" << std::endl;
    Vector<double> w({1, 2, 3, 4, 5, 6, 7, 8, 9, 10});
    std::cout << "  Original: " << w << std::endl;
    
    w.erase(5);  // Remove element at position 5
    std::cout << "  After erase(5):    " << w << std::endl;
    
    w.erase(0, 3);  // Remove range [0, 3)
    std::cout << "  After erase(0, 3): " << w << std::endl;

    // Resize
    std::cout << "\nResize operations:" << std::endl;
    Vector<double> r({1, 2, 3, 4, 5, 6, 7, 8, 9, 10});
    std::cout << "  Original:                " << r << std::endl;
    
    r.Resize(5, true);  // Shrink, preserve elements
    std::cout << "  After Resize(5, true):   " << r << std::endl;
    
    r.Resize(8, true);  // Grow, preserve elements (new ones are zero)
    std::cout << "  After Resize(8, true):   " << r << std::endl;
    
    r.Resize(3);  // Shrink without preserve flag
    std::cout << "  After Resize(3):         " << r << std::endl;

    // Clear
    std::cout << "\nClear operation:" << std::endl;
    Vector<double> toClear({1, 2, 3, 4, 5});
    std::cout << "  Before Clear(): " << toClear << " (size=" << toClear.size() << ")" << std::endl;
    toClear.Clear();
    std::cout << "  After Clear():  size=" << toClear.size() << ", isEmpty=" << (toClear.isEmpty() ? "true" : "false") << std::endl;
}

///////////////////////////////////////////////////////////////////////////////////////////
///                           I/O and Printing Examples                                 ///
///////////////////////////////////////////////////////////////////////////////////////////

void Demo_Vector_IO()
{
    std::cout << "\n=== Vector I/O and Printing Examples ===\n\n";

    Vector<double> v({1.23456789, 2.718281828, 3.14159265, 42.0, 0.001});

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

    // PrintLine with label
    std::cout << "\nPrintLine with label:" << std::endl;
    v.PrintLine(std::cout, "  My vector: ", 10, 4);

    // to_string
    std::cout << "\nString conversion with to_string():" << std::endl;
    std::string str1 = v.to_string(10, 3);
    std::string str2 = v.to_string(6, 1);
    std::cout << "  to_string(10, 3): " << str1 << std::endl;
    std::cout << "  to_string(6, 1):  " << str2 << std::endl;

    // PrintCol for column output
    std::cout << "\nColumn output with PrintCol():" << std::endl;
    Vector<double> small({1.5, 2.5, 3.5});
    small.PrintCol(std::cout, 8, 3);
}

///////////////////////////////////////////////////////////////////////////////////////////
///                           Complex Vector Examples                                   ///
///////////////////////////////////////////////////////////////////////////////////////////

void Demo_Vector_Complex()
{
    std::cout << "\n=== Complex Vector Examples ===\n\n";

    // Create complex vectors
    VectorComplex v1({Complex(1,1), Complex(2,-1), Complex(0,3)});
    VectorComplex v2({Complex(2,0), Complex(1,1), Complex(-1,1)});
    
    std::cout << "v1 = " << v1 << std::endl;
    std::cout << "v2 = " << v2 << std::endl;

    // Basic operations
    VectorComplex sum = v1 + v2;
    VectorComplex diff = v1 - v2;
    
    std::cout << "\nv1 + v2 = " << sum << std::endl;
    std::cout << "v1 - v2 = " << diff << std::endl;

    // Scalar multiplication (complex scalar)
    VectorComplex scaled = v1 * Complex(0, 1);  // Multiply by i
    std::cout << "\nv1 * i = " << scaled << std::endl;

    // Complex dot product (conjugate)
    Complex dot = Utils::ScalarProduct(v1, v2);
    std::cout << "\nConjugate dot product:" << std::endl;
    std::cout << "  Utils::ScalarProduct(v1, v2) = " << dot << std::endl;

    // Note: NormL2() on VectorComplex requires proper complex norm implementation
    // The complex L2 norm should use |z|² = z * conj(z), not z * z
    // For now, compute manually:
    Real norm = 0.0;
    for (const auto& z : v1) {
        norm += std::norm(z);  // |z|² = real² + imag²
    }
    norm = std::sqrt(norm);
    std::cout << "\nManual NormL2(v1) = " << norm << std::endl;
}

///////////////////////////////////////////////////////////////////////////////////////////
///                           Equality Testing Examples                                 ///
///////////////////////////////////////////////////////////////////////////////////////////

void Demo_Vector_Equality()
{
    std::cout << "\n=== Vector Equality Testing ===\n\n";

    Vector<double> a({1.0, 2.0, 3.0});
    Vector<double> b({1.0, 2.0, 3.0});
    Vector<double> c({1.0, 2.0, 3.001});
    
    std::cout << "a = " << a << std::endl;
    std::cout << "b = " << b << std::endl;
    std::cout << "c = " << c << std::endl;

    // Exact equality
    std::cout << "\nExact equality (operator==):" << std::endl;
    std::cout << "  a == b: " << (a == b ? "true" : "false") << std::endl;
    std::cout << "  a == c: " << (a == c ? "true" : "false") << std::endl;
    std::cout << "  a != c: " << (a != c ? "true" : "false") << std::endl;

    // Approximate equality with tolerance
    std::cout << "\nApproximate equality (IsEqualTo with tolerance):" << std::endl;
    std::cout << "  a.IsEqualTo(c, 0.01):   " << (a.IsEqualTo(c, 0.01) ? "true" : "false") << std::endl;
    std::cout << "  a.IsEqualTo(c, 0.001):  " << (a.IsEqualTo(c, 0.001) ? "true" : "false") << std::endl;
    std::cout << "  a.IsEqualTo(c, 0.0001): " << (a.IsEqualTo(c, 0.0001) ? "true" : "false") << std::endl;

    // Null vector check
    std::cout << "\nNull vector check:" << std::endl;
    Vector<double> zero(5);  // All zeros
    Vector<double> nonzero({0, 0, 0.0001, 0, 0});
    std::cout << "  Vector(5).IsNullVec():     " << (zero.IsNullVec() ? "true" : "false") << std::endl;
    std::cout << "  {0,0,0.0001,0,0}.IsNullVec(): " << (nonzero.IsNullVec() ? "true" : "false") << std::endl;
}

///////////////////////////////////////////////////////////////////////////////////////////
///                           Main Entry Point                                          ///
///////////////////////////////////////////////////////////////////////////////////////////

void Docs_Demo_Vector()
{
    std::cout << "╔════════════════════════════════════════════════════════════════════╗" << std::endl;
    std::cout << "║               Vector<T> Documentation Examples                     ║" << std::endl;
    std::cout << "║                  (from docs/base/Vector.md)                        ║" << std::endl;
    std::cout << "╚════════════════════════════════════════════════════════════════════╝" << std::endl;

    Demo_Vector_Construction();
    Demo_Vector_Access();
    Demo_Vector_Arithmetic();
    Demo_Vector_Operations();
    Demo_Vector_Modification();
    Demo_Vector_IO();
    Demo_Vector_Complex();
    Demo_Vector_Equality();

    std::cout << "\n=== All Vector Examples Complete ===\n" << std::endl;
}