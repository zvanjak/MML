///////////////////////////////////////////////////////////////////////////////////////////
// MML Documentation Demo: BaseUtils
// Demonstrates utility functions from BaseUtils.h
///////////////////////////////////////////////////////////////////////////////////////////

#ifdef MML_USE_SINGLE_HEADER
#include "MML.h"
#else
#include "MMLBase.h"

#include "base/Vector.h"
#include "base/VectorN.h"
#include "base/Matrix.h"
#include "base/MatrixNM.h"
#include "base/BaseUtils.h"
#endif

#include <iostream>
#include <iomanip>

using namespace MML;

void Demo_AngleUtilities()
{
    std::cout << "=== ANGLE UTILITIES ===" << std::endl;
    
    // Degree-Radian conversion
    Real rad = Utils::DegToRad(180.0);
    std::cout << "DegToRad(180.0) = " << rad << " (PI = " << Constants::PI << ")" << std::endl;
    
    Real deg = Utils::RadToDeg(Constants::PI);
    std::cout << "RadToDeg(PI) = " << deg << std::endl;
    
    // Explicit angle format (degrees, minutes, seconds)
    Real angleDeg = 45.75;  // 45.75° = 45° 45' 0"
    Real d, m, s;
    Utils::AngleDegToExplicit(angleDeg, d, m, s);
    std::cout << "\nAngleDegToExplicit(" << angleDeg << "°):" << std::endl;
    std::cout << "  " << d << "° " << m << "' " << s << "\"" << std::endl;
    
    // Convert back from deg-min-sec
    Real angle = Utils::ExplicitToAngleDeg(45, 45, 0);
    std::cout << "ExplicitToAngleDeg(45, 45, 0) = " << angle << "°" << std::endl;
    
    // Angle normalization
    Real bigAngle = 5.0 * Constants::PI;  // 900°
    Real normalized = Utils::AngleTo2PiRange(bigAngle);
    std::cout << "\nAngleTo2PiRange(5*PI) = " << normalized << " rad" << std::endl;
    
    Real negAngle = -0.5 * Constants::PI;
    Real normPi = Utils::AngleToPiPiRange(negAngle);
    std::cout << "AngleToPiPiRange(-PI/2) = " << normPi << " rad" << std::endl;
}

void Demo_ComplexUtilities()
{
    std::cout << "\n=== COMPLEX NUMBER UTILITIES ===" << std::endl;
    
    Complex a(1.0, 2.0);
    Complex b(1.0000001, 2.0000001);
    
    // Check equality with default tolerance
    bool equal1 = Utils::AreEqual(a, b);
    std::cout << "AreEqual(" << a << ", " << b << ") = " << (equal1 ? "true" : "false") << std::endl;
    
    // Custom tolerance
    bool equal2 = Utils::AreEqual(a, b, 1e-10);
    std::cout << "AreEqual with eps=1e-10: " << (equal2 ? "true" : "false") << std::endl;
    
    // Absolute comparison
    bool equalAbs = Utils::AreEqualAbs(a, b);
    std::cout << "AreEqualAbs: " << (equalAbs ? "true" : "false") << std::endl;
    
    // Vector of complex numbers
    Vector<Complex> v1(std::vector<Complex>{Complex(1,0), Complex(0,1)});
    Vector<Complex> v2(std::vector<Complex>{Complex(1,0), Complex(0,1.000001)});
    bool vecEqual = Utils::AreEqual(v1, v2);
    std::cout << "\nComplex vector equality: " << (vecEqual ? "true" : "false") << std::endl;
}

void Demo_VectorUtilities()
{
    std::cout << "\n=== VECTOR UTILITIES ===" << std::endl;
    
    // Vector projections
    Vector<Real> v{3, 4, 0};
    Vector<Real> onto{1, 0, 0};  // x-axis
    
    Vector<Real> parallel = Utils::VectorProjectionParallelTo(v, onto);
    Vector<Real> perp = Utils::VectorProjectionPerpendicularTo(v, onto);
    
    std::cout << "Vector v = " << v << std::endl;
    std::cout << "Projection parallel to x-axis: " << parallel << std::endl;
    std::cout << "Projection perpendicular to x-axis: " << perp << std::endl;
    
    // Scalar product
    Vector<Real> a{1, 2, 3};
    Vector<Real> b{4, 5, 6};
    Real dot = Utils::ScalarProduct(a, b);
    std::cout << "\nScalarProduct(" << a << ", " << b << ") = " << dot << std::endl;
    
    // Angle between vectors
    Real angle = Utils::VectorsAngle(a, b);
    std::cout << "Angle between vectors: " << Utils::RadToDeg(angle) << "°" << std::endl;
    
    // Complex vectors
    Vector<Complex> ca(std::vector<Complex>{{1,0}, {0,1}});
    Vector<Complex> cb(std::vector<Complex>{{1,0}, {0,-1}});
    Complex complexDot = Utils::ScalarProduct(ca, cb);
    std::cout << "\nComplex scalar product: " << complexDot << std::endl;
    
    // Outer product (tensor product)
    Vector<Real> u{1, 2};
    Vector<Real> w{3, 4};
    Matrix<Real> outer = Utils::OuterProduct(u, w);
    std::cout << "\nOuter product of " << u << " and " << w << ":" << std::endl;
    std::cout << outer << std::endl;
    
    // Mixed Complex-Real operations
    Vector<Complex> vc(std::vector<Complex>{{1,2}, {3,4}});
    Vector<Real> vr{1, 1};
    Vector<Complex> sum = Utils::AddVec(vc, vr);
    std::cout << "AddVec(complex, real): " << sum << std::endl;
}

void Demo_MatrixConstruction()
{
    std::cout << "\n=== MATRIX CONSTRUCTION FROM VECTORS ===" << std::endl;
    
    Vector<Real> v{1, 2, 3};
    
    // Row matrix
    Matrix<Real> row = Utils::RowMatrixFromVector(v);
    std::cout << "RowMatrixFromVector([1,2,3]):" << std::endl;
    std::cout << row << std::endl;
    
    // Column matrix
    Matrix<Real> col = Utils::ColumnMatrixFromVector(v);
    std::cout << "ColumnMatrixFromVector([1,2,3]):" << std::endl;
    std::cout << col << std::endl;
    
    // Diagonal matrix
    Matrix<Real> diag = Utils::DiagonalMatrixFromVector(v);
    std::cout << "DiagonalMatrixFromVector([1,2,3]):" << std::endl;
    std::cout << diag << std::endl;
    
    // Matrix from vectors in rows/columns
    Vector<Real> row1{1, 2};
    Vector<Real> row2{3, 4};
    Vector<Real> row3{5, 6};
    std::vector<Vector<Real>> rows;
    rows.push_back(row1);
    rows.push_back(row2);
    rows.push_back(row3);
    Matrix<Real> fromRows = Utils::MatrixFromVectorsInRows(rows);
    std::cout << "MatrixFromVectorsInRows:" << std::endl;
    std::cout << fromRows << std::endl;
}

void Demo_MatrixOperations()
{
    std::cout << "\n=== MATRIX OPERATIONS ===" << std::endl;
    
    // Matrix decomposition
    Matrix<Real> A(3, 3);
    A[0][0] = 1; A[0][1] = 2; A[0][2] = 3;
    A[1][0] = 4; A[1][1] = 5; A[1][2] = 6;
    A[2][0] = 7; A[2][1] = 8; A[2][2] = 9;
    
    Matrix<Real> symmetric, antisymmetric;
    Utils::MatrixDecomposeToSymAntisym(A, symmetric, antisymmetric);
    
    std::cout << "Original matrix A:" << std::endl;
    std::cout << A << std::endl;
    std::cout << "Symmetric part:" << std::endl;
    std::cout << symmetric << std::endl;
    std::cout << "Antisymmetric part:" << std::endl;
    std::cout << antisymmetric << std::endl;
    
    // Commutators
    Matrix<Real> M1(2, 2), M2(2, 2);
    M1[0][0] = 1; M1[0][1] = 2; M1[1][0] = 3; M1[1][1] = 4;
    M2[0][0] = 5; M2[0][1] = 6; M2[1][0] = 7; M2[1][1] = 8;
    
    Matrix<Real> comm = Utils::Commutator(M1, M2);
    std::cout << "Commutator [A, B] = AB - BA:" << std::endl;
    std::cout << comm << std::endl;
    
    Matrix<Real> antiComm = Utils::AntiCommutator(M1, M2);
    std::cout << "Anti-commutator {A, B} = AB + BA:" << std::endl;
    std::cout << antiComm << std::endl;
}

void Demo_MatrixFunctions()
{
    std::cout << "\n=== MATRIX FUNCTIONS ===" << std::endl;
    
    // Matrix exponential
    Matrix<Real> omega(2, 2);
    omega[0][0] = 0; omega[0][1] = -1;
    omega[1][0] = 1; omega[1][1] = 0;
    
    std::cout << "Skew-symmetric matrix (infinitesimal rotation):" << std::endl;
    std::cout << omega << std::endl;
    
    Matrix<Real> expA = Utils::Exp(omega * (Constants::PI / 4), 20);
    std::cout << "exp(ω * π/4) (45° rotation matrix):" << std::endl;
    std::cout << expA << std::endl;
    
    // Verify it's a rotation
    std::cout << "Is result orthogonal? " << (Utils::IsOrthogonal(expA) ? "Yes" : "No") << std::endl;
}

void Demo_MatrixProperties()
{
    std::cout << "\n=== MATRIX PROPERTIES ===" << std::endl;
    
    // Orthogonal matrix check
    Matrix<Real> R(2, 2);
    Real theta = Constants::PI / 6;  // 30 degrees
    R[0][0] = std::cos(theta); R[0][1] = -std::sin(theta);
    R[1][0] = std::sin(theta); R[1][1] = std::cos(theta);
    
    std::cout << "Rotation matrix R (30°):" << std::endl;
    std::cout << R << std::endl;
    std::cout << "IsOrthogonal(R) = " << (Utils::IsOrthogonal(R) ? "true" : "false") << std::endl;
    
    // Hermitian matrix check
    Matrix<Complex> H(2, 2);
    H[0][0] = Complex(1, 0);
    H[0][1] = Complex(0, 1);
    H[1][0] = Complex(0, -1);  // Conjugate of H(0,1)
    H[1][1] = Complex(2, 0);
    
    std::cout << "\nHermitian matrix H:" << std::endl;
    std::cout << H << std::endl;
    std::cout << "IsHermitian(H) = " << (Utils::IsHermitian(H) ? "true" : "false") << std::endl;
    
    // Complex matrix operations
    Matrix<Real> realPart = Utils::GetRealPart(H);
    Matrix<Real> imagPart = Utils::GetImagPart(H);
    std::cout << "\nReal part:" << std::endl << realPart << std::endl;
    std::cout << "Imaginary part:" << std::endl << imagPart << std::endl;
    
    Matrix<Complex> Hdagger = Utils::GetConjugateTranspose(H);
    std::cout << "Conjugate transpose H†:" << std::endl << Hdagger << std::endl;
}

void Demo_PhysicsConstants()
{
    std::cout << "\n=== PHYSICS CONSTANTS ===" << std::endl;
    
    // Pauli matrices
    std::cout << "Pauli matrices σ_x, σ_y, σ_z:" << std::endl;
    std::cout << "σ_x:" << std::endl << Utils::Pauli[0] << std::endl;
    std::cout << "σ_y:" << std::endl << Utils::Pauli[1] << std::endl;
    std::cout << "σ_z:" << std::endl << Utils::Pauli[2] << std::endl;
    
    // Dirac gamma matrices (just show first one as example)
    std::cout << "Dirac γ⁰ matrix:" << std::endl;
    std::cout << Utils::DiracGamma[0] << std::endl;
    
    std::cout << "Dirac γ⁵ matrix:" << std::endl;
    std::cout << Utils::DiracGamma5 << std::endl;
}

void Demo_LeviCivita()
{
    std::cout << "\n=== LEVI-CIVITA SYMBOL ===" << std::endl;
    
    // 3D Levi-Civita
    std::cout << "3D Levi-Civita ε_{ijk}:" << std::endl;
    std::cout << "ε_{123} = " << Utils::LeviCivita(1, 2, 3) << " (even permutation)" << std::endl;
    std::cout << "ε_{231} = " << Utils::LeviCivita(2, 3, 1) << " (even permutation)" << std::endl;
    std::cout << "ε_{321} = " << Utils::LeviCivita(3, 2, 1) << " (odd permutation)" << std::endl;
    std::cout << "ε_{112} = " << Utils::LeviCivita(1, 1, 2) << " (repeated index)" << std::endl;
    
    // 4D Levi-Civita
    std::cout << "\n4D Levi-Civita ε_{ijkl}:" << std::endl;
    std::cout << "ε_{1234} = " << Utils::LeviCivita(1, 2, 3, 4) << " (identity)" << std::endl;
    std::cout << "ε_{2134} = " << Utils::LeviCivita(2, 1, 3, 4) << " (one swap)" << std::endl;
    std::cout << "ε_{4321} = " << Utils::LeviCivita(4, 3, 2, 1) << " (reverse)" << std::endl;
}

void Docs_Demo_BaseUtils()
{
    std::cout << "\n" << std::string(70, '=') << std::endl;
    std::cout << "             BASEUTILS - UTILITY FUNCTIONS DEMO" << std::endl;
    std::cout << std::string(70, '=') << std::endl;
    
    Demo_AngleUtilities();
    Demo_ComplexUtilities();
    Demo_VectorUtilities();
    Demo_MatrixConstruction();
    Demo_MatrixOperations();
    Demo_MatrixFunctions();
    Demo_MatrixProperties();
    Demo_PhysicsConstants();
    Demo_LeviCivita();
    
    std::cout << "\n" << std::string(70, '=') << std::endl;
    std::cout << "           END OF BASEUTILS DEMO" << std::endl;
    std::cout << std::string(70, '=') << std::endl;
}
