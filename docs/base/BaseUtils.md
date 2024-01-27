# BaseUtils

## Levi-Civitta

## Complex helpers              
static bool AreEqual(const Complex &a, const Complex &b, double eps=Defaults::ComplexEqualityPrecision)
static bool AreEqual(const Vector<Complex> &a, const Vector<Complex> &b, double eps=Defaults::ComplexEqualityPrecision)

## Vector helpers
static Vector<Real> VectorProjectionParallelTo(const Vector<Real> &orig, const Vector<Real> &b)
static Vector<Real> VectorProjectionPerpendicularTo(const Vector<Real> &orig, const Vector<Real> &b)
template<class _Type>
static Matrix<_Type> OuterProduct(const Vector<_Type> &a, const Vector<_Type> &b)

## Vector\<Complex> - Vector<Real> operations
static Vector<Complex> AddVec(const Vector<Complex> &a, const Vector<Real> &b)
static Vector<Complex> AddVec(const Vector<Real> &a, const Vector<Complex> &b)
static Vector<Complex> SubVec(const Vector<Complex> &a, const Vector<Real> &b)
static Vector<Complex> SubVec(const Vector<Real> &a, const Vector<Complex> &b)


namespace MatrixUtils
const static inline MatrixNM<Complex, 2, 2> Pauli[]

const static inline MatrixNM<Complex, 4, 4> DiracGamma[]

/////////////////////             Creating Matrix from Vector             ///////////////////
template<class _Type>
static Matrix<_Type> RowMatrixFromVector(const Vector<_Type> &b)
template<class _Type>
static Matrix<_Type> ColumnMatrixFromVector(const Vector<_Type> &b)
template<class _Type>
static Matrix<_Type> DiagonalMatrixFromVector(const Vector<_Type> &b)

///////////////////////             Matrix helpers               ////////////////////
template<class _Type>
static Matrix<_Type> Commutator(const Matrix<_Type> &a, const Matrix<_Type> &b)
template<class _Type>
static Matrix<_Type> AntiCommutator(const Matrix<_Type> &a, const Matrix<_Type> &b)
template<class _Type>
static void MatrixDecompose(const Matrix<_Type> &orig, Matrix<_Type> &outSym, Matrix<_Type> &outAntiSym)

///////////////////////             Matrix functions               ////////////////////
template<class _Type>
static Matrix<_Type> Exp(const Matrix<_Type> &a, int n = 10)

///////////////////////             Real matrix helpers               ////////////////////
static bool IsOrthogonal(const Matrix<Real> &mat, double eps = Defaults::IsMatrixOrthogonalPrecision)

//////////////////////             Complex matrix helpers               ///////////////////
static Matrix<Complex> GetConjugateTranspose(const Matrix<Complex> &mat)
static Matrix<Complex> CmplxMatFromRealMat(const Matrix<Real> &mat)
static bool IsComplexMatReal(const Matrix<Complex> &mat)
static bool IsHermitian(const Matrix<Complex> &mat)
static bool IsUnitary(const Matrix<Complex> &mat)

//////////////////////         enabling Matrix<Complex> - Matrix<Real>  operations         ///////////////////////



