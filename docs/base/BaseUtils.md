# BaseUtils - Utility Functions

**File**: `mml/base/BaseUtils.h`

Collection of utility functions for angles, complex numbers, vectors, matrices, and physics constants.

## Table of Contents
- [Levi-Civita Symbol](#levi-civita-symbol)
- [Angle Utilities](#angle-utilities)
- [Complex Number Utilities](#complex-number-utilities)
- [Vector Utilities](#vector-utilities)
- [Matrix Utilities](#matrix-utilities)
- [Matrix Properties](#matrix-properties)
- [Physics Constants](#physics-constants)

---

## Levi-Civita Symbol

### 3D Levi-Civita
```cpp
using namespace MML::Utils;

// ε_{ijk} for i,j,k ∈ {1,2,3}
int eps123 = LeviCivita(1, 2, 3);  // +1 (even permutation)
int eps231 = LeviCivita(2, 3, 1);  // +1 (even permutation)
int eps321 = LeviCivita(3, 2, 1);  // -1 (odd permutation)
int eps112 = LeviCivita(1, 1, 2);  //  0 (repeated index)
```

### 4D Levi-Civita
```cpp
// ε_{ijkl} for i,j,k,l ∈ {1,2,3,4}
int eps1234 = LeviCivita(1, 2, 3, 4);  // +1 (identity permutation)
int eps2134 = LeviCivita(2, 1, 3, 4);  // -1 (one transposition)
int eps4321 = LeviCivita(4, 3, 2, 1);  // +1 (even number of swaps)
```

---

## Angle Utilities

### Degree-Radian Conversion
```cpp
using namespace MML::Utils;

// Convert degrees to radians
Real rad = DegToRad(180.0);  // π

// Convert radians to degrees
Real deg = RadToDeg(Constants::PI);  // 180.0
```

### Explicit Angle Format (Degrees, Minutes, Seconds)
```cpp
// Convert angle to deg-min-sec
Real angleDeg = 45.75;  // 45.75° = 45° 45' 0"
Real d, m, s;
AngleDegToExplicit(angleDeg, d, m, s);
// d = 45, m = 45, s = 0

// Convert from deg-min-sec
Real angle = ExplicitToAngleDeg(45, 45, 0);  // 45.75°

// Radian versions
Real angleRad = Constants::PI / 4;
AngleRadToExplicit(angleRad, d, m, s);
Real reconstructed = ExplicitToAngleRad(d, m, s);
```

### Angle Normalization
```cpp
// Normalize angle to [0, 2π) range
Real bigAngle = 5 * Constants::PI;  // 900°
Real normalized = AngleTo2PiRange(bigAngle);  // ~π

// Normalize angle to [-π, π) range
Real negAngle = 3 * Constants::PI;  // 540°
Real symNorm = AngleToPiPiRange(negAngle);  // ~-π
```

---

## Complex Number Utilities

### Equality with Tolerance
```cpp
using namespace MML::Utils;

Complex a(1.0, 2.0);
Complex b(1.0000001, 2.0000001);

// Check equality with default tolerance
bool equal = AreEqual(a, b);  // true (within default eps)

// Custom tolerance
bool equal2 = AreEqual(a, b, 1e-10);  // false

// Absolute value comparison
bool equalAbs = AreEqualAbs(a, b);  // uses |a - b| comparison

// Vector of complex numbers
Vector<Complex> v1 = {Complex(1,0), Complex(0,1)};
Vector<Complex> v2 = {Complex(1,0), Complex(0,1.000001)};
bool vec_equal = AreEqual(v1, v2);
bool vec_equal_abs = AreEqualAbs(v1, v2);  // absolute comparison
```

---

## Vector Utilities

### Vector Projections
```cpp
using namespace MML::Utils;

Vector<Real> v{3, 4, 0};
Vector<Real> onto{1, 0, 0};  // x-axis

// Project v parallel to 'onto'
Vector<Real> parallel = VectorProjectionParallelTo(v, onto);
// Result: {3, 0, 0}

// Project v perpendicular to 'onto'
Vector<Real> perp = VectorProjectionPerpendicularTo(v, onto);
// Result: {0, 4, 0}

// Verify: v = parallel + perp
```

### Scalar and Outer Products
```cpp
// Scalar product (dot product)
Vector<Real> a{1, 2, 3};
Vector<Real> b{4, 5, 6};
Real dot = ScalarProduct(a, b);  // 1*4 + 2*5 + 3*6 = 32

// Angle between vectors
Real angle = VectorsAngle(a, b);  // radians

// Complex vectors (uses conjugate in dot product)
Vector<Complex> ca{{1,0}, {0,1}};
Vector<Complex> cb{{1,0}, {0,-1}};
Complex complex_dot = ScalarProduct(ca, cb);  // ⟨ca|cb⟩

// Outer product (tensor product)
Vector<Real> u{1, 2};
Vector<Real> v{3, 4};
Matrix<Real> outer = OuterProduct(u, v);
// Result: [[3, 4],
//          [6, 8]]
```

### Mixed Complex-Real Vector Operations
```cpp
Vector<Complex> vc{{1,2}, {3,4}};
Vector<Real> vr{1, 1};

// Add complex and real vectors
Vector<Complex> sum1 = AddVec(vc, vr);
Vector<Complex> sum2 = AddVec(vr, vc);

// Subtract
Vector<Complex> diff1 = SubVec(vc, vr);
Vector<Complex> diff2 = SubVec(vr, vc);

// Scalar multiplication
Vector<Complex> scaled = MulVec(Complex(2, 1), vr);  // Complex * Real vector
```

---

## Matrix Utilities

### Matrix Construction from Vectors
```cpp
using namespace MML::Utils;

Vector<Real> v{1, 2, 3};

// Row matrix [1 2 3]
Matrix<Real> row = RowMatrixFromVector(v);
// 1x3 matrix

// Column matrix [1]
//               [2]
//               [3]
Matrix<Real> col = ColumnMatrixFromVector(v);
// 3x1 matrix

// Diagonal matrix [1 0 0]
//                 [0 2 0]
//                 [0 0 3]
Matrix<Real> diag = DiagonalMatrixFromVector(v);
// 3x3 matrix

// Matrix from list of vectors
std::vector<Vector<Real>> rows = {{1, 2}, {3, 4}};
Matrix<Real> fromRows = MatrixFromVectorsInRows(rows);
// [[1 2]
//  [3 4]]

Matrix<Real> fromCols = MatrixFromVectorsInColumns(rows);
// [[1 3]
//  [2 4]]
```

### Matrix Decomposition (Symmetric/Antisymmetric)
```cpp
Matrix<Real> A(3, 3);
// ... fill A ...

Matrix<Real> symmetric, antisymmetric;
MatrixDecomposeToSymAntisym(A, symmetric, antisymmetric);

// A = symmetric + antisymmetric
// symmetric^T = symmetric
// antisymmetric^T = -antisymmetric
```

### Commutators
```cpp
Matrix<Real> A(2, 2), B(2, 2);
// ... initialize ...

// Commutator: [A, B] = AB - BA
Matrix<Real> comm = Commutator(A, B);

// Anti-commutator: {A, B} = AB + BA
Matrix<Real> anti_comm = AntiCommutator(A, B);
```

### Matrix Functions
```cpp
// Matrix exponential: e^A = I + A + A²/2! + A³/3! + ...
Matrix<Real> A(2, 2);
A(0,0) = 0; A(0,1) = -1;
A(1,0) = 1; A(1,1) = 0;

// Compute exp(A) with 10 terms
Matrix<Real> expA = Exp(A, 10);

// For rotation matrices, exp of skew-symmetric gives rotation

// Matrix Sin and Cos (via Taylor series)
Matrix<Real> sinA = Sin(A, 10);  // sin(A)
Matrix<Real> cosA = Cos(A, 10);  // cos(A)
```

---

## Matrix Properties

### Orthogonal Matrices
```cpp
using namespace MML::Utils;

Matrix<Real> Q(3, 3);
// ... fill with rotation matrix ...

// Check if Q is orthogonal (Q^T * Q = I)
bool is_orthogonal = IsOrthogonal(Q);

// Custom tolerance
bool is_orthogonal2 = IsOrthogonal(Q, 1e-10);
```

### Unitary Matrices
```cpp
Matrix<Complex> U(2, 2);
// ... fill with unitary matrix ...

// Check if U is unitary (U† * U = I)
bool is_unitary = IsUnitary(U);

// Custom tolerance
bool is_unitary2 = IsUnitary(U, 1e-10);
```

### Hermitian Matrices
```cpp
Matrix<Complex> H(2, 2);
H(0,0) = Complex(1, 0);
H(0,1) = Complex(0, 1);
H(1,0) = Complex(0, -1);  // Conjugate of H(0,1)
H(1,1) = Complex(2, 0);

// Check if H = H† (conjugate transpose)
bool is_hermitian = IsHermitian(H);
```

### Complex Matrix Operations
```cpp
Matrix<Complex> C(2, 2);
C(0,0) = Complex(1, 2);
C(0,1) = Complex(3, 4);
// ...

// Extract real part
Matrix<Real> real_part = GetRealPart(C);

// Extract imaginary part
Matrix<Real> imag_part = GetImagPart(C);

// Conjugate transpose (Hermitian adjoint)
Matrix<Complex> C_dagger = GetConjugateTranspose(C);

// Convert real matrix to complex
Matrix<Real> R(2, 2);
Matrix<Complex> Rc = CmplxMatFromRealMat(R);

// Check if complex matrix is purely real
bool isReal = IsComplexMatReal(C);
```

### Mixed Complex-Real Matrix Operations
```cpp
Matrix<Complex> Mc(2, 2);
Matrix<Real> Mr(2, 2);

// Addition
Matrix<Complex> sum1 = AddMat(Mc, Mr);
Matrix<Complex> sum2 = AddMat(Mr, Mc);

// Subtraction
Matrix<Complex> diff1 = SubMat(Mc, Mr);
Matrix<Complex> diff2 = SubMat(Mr, Mc);

// Multiplication
Matrix<Complex> prod = MulMat(Mc, Mr);

// Scalar multiplication
Matrix<Complex> scaled = MulMat(Complex(2, 1), Mr);

// Matrix-vector multiplication
Vector<Complex> vc(2);
Vector<Real> vr(2);
Vector<Complex> result1 = MulMatVec(Mr, vc);  // Real matrix * Complex vector
Vector<Complex> result2 = MulMatVec(Mc, vr);  // Complex matrix * Real vector
Vector<Complex> result3 = MulVecMat(vc, Mr);  // Complex vector * Real matrix
```

---

## Physics Constants

### Pauli Matrices
```cpp
using namespace MML::Utils;

// σ_x, σ_y, σ_z
const MatrixNM<Complex, 2, 2>& sigma_x = Pauli[0];
const MatrixNM<Complex, 2, 2>& sigma_y = Pauli[1];
const MatrixNM<Complex, 2, 2>& sigma_z = Pauli[2];

// σ_x = [0 1]
//       [1 0]
// σ_y = [0 -i]
//       [i  0]
// σ_z = [1  0]
//       [0 -1]
```

### Dirac Gamma Matrices
```cpp
const MatrixNM<Complex, 4, 4>& gamma_0 = DiracGamma[0];
const MatrixNM<Complex, 4, 4>& gamma_1 = DiracGamma[1];
const MatrixNM<Complex, 4, 4>& gamma_2 = DiracGamma[2];
const MatrixNM<Complex, 4, 4>& gamma_3 = DiracGamma[3];

// γ⁵ = iγ⁰γ¹γ²γ³
const MatrixNM<Complex, 4, 4>& gamma_5 = DiracGamma5;

// Standard representation (Dirac basis)
```

---

## Examples

### Example 1: Vector Decomposition
```cpp
using namespace MML::Utils;

// Decompose velocity into radial and tangential components
Vector<Real> velocity{3, 4, 0};
Vector<Real> radial_direction{1, 0, 0};

auto radial_component = VectorProjectionParallelTo(velocity, radial_direction);
auto tangential_component = VectorProjectionPerpendicularTo(velocity, radial_direction);

std::cout << "Radial: " << radial_component << std::endl;
std::cout << "Tangential: " << tangential_component << std::endl;
```

### Example 2: Matrix Exponential for Rotation
```cpp
using namespace MML::Utils;

// Rotation matrix from exponential of skew-symmetric matrix
Matrix<Real> omega(3, 3);
omega(0,0) = 0;  omega(0,1) = 0;  omega(0,2) = 1;   // ω_z = 1
omega(1,0) = 0;  omega(1,1) = 0;  omega(1,2) = 0;
omega(2,0) = -1; omega(2,1) = 0;  omega(2,2) = 0;

// R = exp(ω * θ) rotates by angle θ around z-axis
Matrix<Real> R = Exp(omega * Constants::PI/4, 20);  // 45° rotation
```

### Example 3: Quantum Mechanics - Pauli Spin
```cpp
using namespace MML::Utils;

// Spin-1/2 state |+⟩ along z-axis
Vector<Complex> spin_up{{1,0}, {0,0}};

// Measure spin along x-axis
const auto& sigma_x = Pauli[0];
Matrix<Complex> sigma_x_mat = // convert to dynamic matrix if needed

// Expectation value ⟨σ_x⟩ = ⟨ψ|σ_x|ψ⟩
```

---

## See Also
- [Vector.md](Vector.md), [VectorN.md](VectorN.md) - Vector types
- [Matrix.md](Matrix.md), [MatrixNM.md](MatrixNM.md) - Matrix types
- [Linear_equations_solvers.md](../core/Linear_equations_solvers.md) - Matrix decompositions
- [Coordinate_transformations.md](../core/Coordinate_transformations.md) - Using rotation matrices

---

## Runnable Examples

| Example | Source File | Description |
|---------|------------|-------------|
| BaseUtils Demo | [docs_demo_base_utils.cpp](../../src/docs_demos/docs_demo_base_utils.cpp) | Comprehensive demo of all utility functions |

**Build and Run:**
```bash
cmake --build build --target MML_DocsApp
./build/src/docs_demos/Release/MML_DocsApp
```





