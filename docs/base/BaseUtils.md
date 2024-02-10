# BaseUtils

## Levi-Civitta

## Complex helpers     

~~~c++
static bool AreEqual(const Complex &a, const Complex &b, double eps=Defaults::ComplexEqualityPrecision)
static bool AreEqual(const Vector<Complex> &a, const Vector<Complex> &b, double eps=Defaults::ComplexEqualityPrecision)
~~~

## Vector helpers

~~~c++
static Vector<Real> VectorProjectionParallelTo(const Vector<Real> &orig, const Vector<Real> &b);
static Vector<Real> VectorProjectionPerpendicularTo(const Vector<Real> &orig, const Vector<Real> &b);
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
template<class Type> static Matrix<Type> RowMatrixFromVector(const Vector<Type> &b)
template<class Type> static Matrix<Type> ColumnMatrixFromVector(const Vector<Type> &b)
template<class Type> static Matrix<Type> DiagonalMatrixFromVector(const Vector<Type> &b)

template<class Type> static Matrix<Type> Commutator(const Matrix<Type> &a, const Matrix<Type> &b)
template<class Type> static Matrix<Type> AntiCommutator(const Matrix<Type> &a, const Matrix<Type> &b)

template<class Type> static void MatrixDecompose(const Matrix<Type> &orig, Matrix<Type> &outSym, Matrix<Type> &outAntiSym)
~~~

## Matrix functions

~~~c++
template<class Type> static Matrix<Type> Exp(const Matrix<Type> &a, int n = 10)
~~~

## Real matrix helpers

~~~c++
static bool IsOrthogonal(const Matrix<Real> &mat, double eps = Defaults::IsMatrixOrthogonalPrecision)
~~~

## Complex matrix helpers

~~~c++
static Matrix<Complex> GetConjugateTranspose(const Matrix<Complex> &mat)
static Matrix<Complex> CmplxMatFromRealMat(const Matrix<Real> &mat)
static bool IsComplexMatReal(const Matrix<Complex> &mat)
static bool IsHermitian(const Matrix<Complex> &mat)
static bool IsUnitary(const Matrix<Complex> &mat)
~~~

## enabling Matrix<Complex> - Matrix<Real>  operations




