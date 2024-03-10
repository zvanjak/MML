# BaseUtils

## Angle calculations

~~~C++
static Real DegToRad(Real angleDeg) { return angleDeg * Constants::PI / 180.0; }
static Real RadToDeg(Real angleRad) { return angleRad * 180.0 / Constants::PI; }

static void AngleDegToExplicit(Real angle, Real& deg, Real& min, Real& sec);
static void AngleRadToExplicit(Real angleRad, Real& deg, Real& min, Real& sec) { AngleDegToExplicit(angleRad * 180.0 / Constants::PI, deg, min, sec); }
static Real ExplicitToAngleDeg(Real deg, Real min, Real sec) { return deg + min / 60.0 + sec / 3600.0; }
static Real ExplicitToAngleRad(Real deg, Real min, Real sec) { return ExplicitToAngleDeg(deg, min, sec) * Constants::PI / 180.0; }

~~~

## Complex helpers     

~~~c++
static bool AreEqual(const Complex &a, const Complex &b, double eps=Defaults::ComplexEqualityPrecision)
static bool AreEqual(const Vector<Complex> &a, const Vector<Complex> &b, double eps=Defaults::ComplexEqualityPrecision)
~~~

## Vector helpers

~~~c++
static Vector<Real> VectorProjectionParallelTo(const Vector<Real> &orig, const Vector<Real> &b);
static Vector<Real> VectorProjectionPerpendicularTo(const Vector<Real> &orig, const Vector<Real> &b);

static Real ScalarProduct(const Vector<Real>& a, const Vector<Real>& b);
static Complex ScalarProduct(const Vector<Complex>& a, const Vector<Complex>& b);
		
template<class Type> static Matrix<Type> OuterProduct(const Vector<Type> &a, const Vector<Type> &b);
~~~

## Vector\<Complex> - Vector\<Real> operations

~~~c++
static Vector<Complex> AddVec(const Vector<Complex> &a, const Vector<Real> &b);
static Vector<Complex> AddVec(const Vector<Real> &a, const Vector<Complex> &b);
static Vector<Complex> SubVec(const Vector<Complex> &a, const Vector<Real> &b);
static Vector<Complex> SubVec(const Vector<Real> &a, const Vector<Complex> &b);
~~~

## Matrix helpers

~~~c++
const static inline MatrixNM<Complex, 2, 2> Pauli[];
const static inline MatrixNM<Complex, 4, 4> DiracGamma[];

/////////////////////             Creating Matrix from Vector             ///////////////////
template<class Type> static Matrix<Type> RowMatrixFromVector(const Vector<Type> &b);
template<class Type> static Matrix<Type> ColumnMatrixFromVector(const Vector<Type> &b);
template<class Type> static Matrix<Type> DiagonalMatrixFromVector(const Vector<Type> &b);

template<class Type> static Matrix<Type> Commutator(const Matrix<Type> &a, const Matrix<Type> &b);
template<class Type> static Matrix<Type> AntiCommutator(const Matrix<Type> &a, const Matrix<Type> &b);

template<class Type> static void MatrixDecompose(const Matrix<Type> &orig, Matrix<Type> &outSym, Matrix<Type> &outAntiSym);
~~~

## Matrix functions

~~~c++
template<class Type> static Matrix<Type> Exp(const Matrix<Type> &a, int n = 10);
~~~

## Real matrix helpers

~~~c++
static bool IsOrthogonal(const Matrix<Real> &mat, double eps = Defaults::IsMatrixOrthogonalPrecision);
~~~

## Complex matrix helpers

~~~c++
static Matrix<Real> GetRealPart(const Matrix<Complex>& a);
static Matrix<Real> GetImagPart(const Matrix<Complex>& a);

static Matrix<Complex> GetConjugateTranspose(const Matrix<Complex> &mat)
static Matrix<Complex> CmplxMatFromRealMat(const Matrix<Real> &mat)
static bool IsComplexMatReal(const Matrix<Complex> &mat)
static bool IsHermitian(const Matrix<Complex> &mat)
static bool IsUnitary(const Matrix<Complex> &mat)
~~~

## Enabling Matrix<Complex> - Matrix<Real>  operations

~~~c++
static Matrix<Complex> AddMat(const Matrix<Complex>& a, const Matrix<Real>& b);
static Matrix<Complex> AddMat(const Matrix<Real>& a, const Matrix<Complex>& b);

static Matrix<Complex> SubMat(const Matrix<Complex>& a, const Matrix<Real>& b);
static Matrix<Complex> SubMat(const Matrix<Real>& a, const Matrix<Complex>& b);

static Matrix<Complex> MulMat(const Complex& a, const Matrix<Real>& b);
static Matrix<Complex> MulMat(const Matrix<Complex>& a, const Matrix<Real>& b);
static Matrix<Complex> MulMat(const Matrix<Real>& a, const Matrix<Complex>& b);

static Vector<Complex> MulMatVec(const Matrix<Real>& a, const Vector<Complex>& b);
static Vector<Complex> MulMatVec(const Matrix<Complex>& a, const Vector<Real>& b);

static Vector<Complex> MulVecMat(const Vector<Complex>& a, const Matrix<Real>& b);
static Vector<Complex> MulVecMat(const Vector<Real>& a, const Matrix<Complex>& b);
~~~





