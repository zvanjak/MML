# MatrixNM class

Class representing matrix of statically define type and size (as template parameters).

## Class definition
~~~ c++
template <class Type, int N, int M>
class MatrixNM
{
public:
	Type _vals[N][M] = { {0} };

public:
	//////////////////////////             Constructors           /////////////////////////
	MatrixNM() {}
	MatrixNM(std::initializer_list<Type> values)
	MatrixNM(const MatrixNM& m)
	MatrixNM(const Type& m)        // initialize as diagonal matrix

	typedef Type value_type;      // make T available externally

	////////////////////////            Standard stuff             ////////////////////////
	int RowNum() const { return N; }
	int ColNum() const { return M; }

	static MatrixNM GetUnitMatrix()
	void   MakeUnitMatrix(void)

	MatrixNM GetLower(bool includeDiagonal = true) const
	MatrixNM GetUpper(bool includeDiagonal = true) const

	/////////////////////          Vector-Matrix conversion           /////////////////////
	static MatrixNM<Type, 1, M>  RowMatrixFromVector(const VectorN<Type, M>& b)
	static MatrixNM<Type, N, 1>  ColumnMatrixFromVector(const VectorN<Type, N>& b)
	static MatrixNM<Type, N, N> DiagonalMatrixFromVector(const VectorN<Type, N>& b)
	static VectorN<Type, M> VectorFromRow(const MatrixNM<Type, N, M>& a, int rowInd)
	static VectorN<Type, N> VectorFromColumn(const MatrixNM<Type, N, M>& a, int colInd)
	static VectorN<Type, N> VectorFromDiagonal(const MatrixNM<Type, N, N>& a)

	/////////////////////            Assignment operators             ////////////////////
	MatrixNM& operator=(const MatrixNM& m)
	MatrixNM& operator=(const Type& m)

	////////////////////            Access operators             ///////////////////////
	Type* operator[](int i) { return _vals[i]; }
	const Type* operator[](const int i) const { return _vals[i]; }

	Type  operator()(int i, int j) const { return _vals[i][j]; }
	Type& operator()(int i, int j) { return _vals[i][j]; }

	// version with checking bounds
	Type  ElemAt(int i, int j) const
	Type& ElemAt(int i, int j)

	////////////////////            Arithmetic operators             ////////////////////
	MatrixNM operator-()            // unary minus
	MatrixNM operator+(const MatrixNM& b) const
	MatrixNM operator-(const MatrixNM& b) const
	template<int K>
	MatrixNM<Type, N, K>  operator*(const MatrixNM<Type, M, K>& b) const

	friend MatrixNM operator*(const MatrixNM& a, Type b)
	friend MatrixNM operator/(const MatrixNM& a, Type b)
	friend MatrixNM operator*(Type a, const MatrixNM& b)

	friend VectorN<Type, N> operator*(const MatrixNM<Type, N, M>& a, const VectorN<Type, M>& b)
	friend VectorN<Type, M> operator*(const VectorN<Type, N>& a, const MatrixNM<Type, N, M>& b)

	///////////////////////            Equality operations             //////////////////////
	bool operator==(const MatrixNM& b) const
	bool operator!=(const MatrixNM& b) const

	bool IsEqual(const MatrixNM& b, Type eps = Defaults::MatrixEqualityPrecision) const
	bool AreEqual(const MatrixNM& a, const MatrixNM& b, Type eps = Defaults::MatrixEqualityPrecision) const

	///////////////////            Trace, Inverse & Transpose             ///////////////////
	Type   Trace() const

	void Invert()
	MatrixNM GetInverse() const

	void Transpose()
	MatrixNM<Type, M, N> GetTranspose() const

	///////////////////////////               I/O                 ///////////////////////////
	std::string to_string(int width, int precision) const
	void Print(std::ostream& stream, int width, int precision) const
	friend std::ostream& operator<<(std::ostream& stream, const MatrixNM& a)
};
~~~

## MatrixNM initialization
~~~ c++
MatrixNM<Real, 2, 2> a;
MatrixNM<Real, 2, 2> b({ 1.0, 0.0, 0.0, 1.0 });
MatrixNM<Real, 2, 2> c(b);
MatrixNM<Real, 2, 2> d = c;
auto e = MatrixNM<Real, 3, 3>::GetUnitMatrix();

std::cout << "a = " << a << std::endl;
std::cout << "b = " << b << std::endl;
std::cout << "c = " << c << std::endl;
std::cout << "d = " << d << std::endl;
std::cout << "e = " << e << std::endl;

	/* OUTPUT
a = Rows: 2  Cols: 2
[          0,          0,  ]
[          0,          0,  ]

b = Rows: 2  Cols: 2
[          1,          0,  ]
[          0,          1,  ]

c = Rows: 2  Cols: 2
[          1,          0,  ]
[          0,          1,  ]

d = Rows: 2  Cols: 2
[          1,          0,  ]
[          0,          1,  ]

e = Rows: 3  Cols: 3
[          1,          0,          0,  ]
[          0,          1,          0,  ]
[          0,          0,          1,  ]
	*/
~~~

## MatrixNM - VectorN init operations
~~~ c++

VectorN<Real, 3> a({ 1.0, 1.0, 1.0 });
MatrixNM<Real, 1, 3> matA = MatrixNM<Real, 1, 3>::RowMatrixFromVector(a);
auto matAauto = MatrixNM<Real, 1, 3>::RowMatrixFromVector(a);
MatrixNM<Real, 3, 1> matB = MatrixNM<Real, 3, 1>::ColumnMatrixFromVector(a);
auto matBauto = MatrixNM<Real, 3, 1>::ColumnMatrixFromVector(a);

std::cout << "Vector a = " << a << std::endl;
std::cout << "Matrix matA = Matrix::RowMatrixFromVector(a);\nmatA = " << matA << std::endl;
std::cout << "Matrix matB = Matrix::ColMatrixFromVector(a);\nmatB = " << matB << std::endl;

MatrixNM<Real, 2, 2> m1({ 1.0, -1.0, 1.5, 3.0 });
VectorN<Real, 2> vecRow = MatrixNM<Real, 2, 2>::VectorFromRow(m1, 0);
VectorN<Real, 2> vecCol = MatrixNM<Real, 2, 2>::VectorFromColumn(m1, 0);
VectorN<Real, 2> vecDiag = MatrixNM<Real, 2, 2>::VectorFromDiagonal(m1);

std::cout << "Matrix m1 = " << m1 << std::endl;
std::cout << "Vector vecRow = Matrix::VectorFromRow(a,0)     = " << vecRow << std::endl;
std::cout << "Vector vecCol = Matrix::VectorFromColumn(a, 0) = " << vecCol << std::endl;
std::cout << "Vector vecCol = Matrix::VectorFromDiagonal(a)  = " << vecDiag << std::endl;

/* OUTPUT
Vector a = [              1,               1,               1]
Matrix matA = Matrix::RowMatrixFromVector(a);
matA = Rows: 1  Cols: 3
[          1,          1,          1,  ]

Matrix matB = Matrix::ColMatrixFromVector(a);
matB = Rows: 3  Cols: 1
[          1,  ]
[          0,  ]
[          0,  ]

Matrix m1 = Rows: 2  Cols: 2
[          1,         -1,  ]
[        1.5,          3,  ]

Vector vecRow = Matrix::VectorFromRow(a,0)     = [              1,              -1]
Vector vecCol = Matrix::VectorFromColumn(a, 0) = [              1,             1.5]
Vector vecCol = Matrix::VectorFromDiagonal(a)  = [              1,               3]
*/
~~~

## MatrixNM Basic operations
~~~ c++
MatrixNM<Real, 2, 2> m1({ 1.0, -1.0, 1.5, 3.0 }), m2;
m2.MakeUnitMatrix();

std::cout << "m1 = " << m1 << std::endl;
std::cout << "m2 = " << m2 << std::endl;

std::cout << "m1 + m2 = " << m1 + m2 << std::endl;
std::cout << "m1 - m2 = " << m1 - m2 << std::endl;
std::cout << "m1 * m2 = " << m1 * m2 << std::endl;
std::cout << "2.0 * m1 = " << 2.0 * m1 << std::endl;
std::cout << "m1 * 2.0 = " << m1 * 2.0 << std::endl;
std::cout << "m1 / 2.0 = " << m1 / 2.0 << std::endl;

/* OUTPUT
m1 = Rows: 2  Cols: 2
[          1,         -1,  ]
[        1.5,          3,  ]

m2 = Rows: 2  Cols: 2
[          1,          0,  ]
[          0,          1,  ]

m1 + m2 = Rows: 2  Cols: 2
[          2,         -1,  ]
[        1.5,          4,  ]

m1 - m2 = Rows: 2  Cols: 2
[          0,         -1,  ]
[        1.5,          2,  ]

m1 * m2 = Rows: 2  Cols: 2
[          1,         -1,  ]
[        1.5,          3,  ]

2.0 * m1 = Rows: 2  Cols: 2
[          2,         -2,  ]
[          3,          6,  ]

m1 * 2.0 = Rows: 2  Cols: 2
[          2,         -2,  ]
[          3,          6,  ]

m1 / 2.0 = Rows: 2  Cols: 2
[        0.5,       -0.5,  ]
[       0.75,        1.5,  ]
*/
~~~

## MatrixNM-VectorN multiplication
~~~ c++
VectorN<Real, 2> v1({ 1.0, 2.0 });
MatrixNM<Real, 2, 2> m1({ 1.0, -1.0, 1.5, 3.0 });

std::cout << "v1 = " << v1 << std::endl;
std::cout << "m1 = " << m1 << std::endl;

std::cout << "v1 * m1 = " << v1 * m1 << std::endl;
std::cout << "m1 * v1 = " << m1 * v1 << std::endl;

/* OUTPUT
v1 = [              1,               2]
m1 = Rows: 2  Cols: 2
[          1,         -1,  ]
[        1.5,          3,  ]

v1 * m1 = [              4,               5]
m1 * v1 = [             -1,             7.5]
*/
~~~

## MatrixNM-MatrixNM multiplication
~~~ c++
MatrixNM<Real, 1, 3> m3{ 1.0, -2.0, 3.0 };
MatrixNM<Real, 3, 4> m4{ 1.0, 0.0, 0.0, 0.0,
												0.0, 1.0, 0.0, 0.0,
												0.0, 0.0, 1.0, 1.0 };
MatrixNM<Real, 1, 4> m5 = m3 * m4;

std::cout << "m3 = " << m3 << std::endl;
std::cout << "m4 = " << m4 << std::endl;
std::cout << "m3 * m4 = " << m5 << std::endl;

/* OUTPUT
m3 = Rows: 1  Cols: 3
[          1,         -2,          3,  ]

m4 = Rows: 3  Cols: 4
[          1,          0,          0,          0,  ]
[          0,          1,          0,          0,  ]
[          0,          0,          1,          1,  ]

m3 * m4 = Rows: 1  Cols: 4
[          1,         -2,          3,          3,  ]
*/
~~~

## MatrixNM Invert
~~~ c++
MatrixNM<Real, 2, 2> m1({ 1.0, -1.0, 1.5, 3.0 });

std::cout << "m1 = " << m1 << std::endl;
auto m2 = m1.GetInverse();

std::cout << "m2 (inv) = " << m2 << std::endl;

auto munit = m1 * m2;
std::cout << "m1 * m2 = " << munit << std::endl;

/* OUTPUT
m1 = Rows: 2  Cols: 2
[          1,         -1,  ]
[        1.5,          3,  ]

m2 (inv) = Rows: 2  Cols: 2
[      0.667,      0.222,  ]
[     -0.333,      0.222,  ]

m1 * m2 = Rows: 2  Cols: 2
[          1,          0,  ]
[          0,          1,  ]
*/
~~~

## MatrixNM Transpose
~~~ c++
MatrixNM<Real, 2, 2> m1({ 1.0, -1.0, 1.5, 3.0 });
std::cout << "m1 = " << m1 << std::endl;

auto m2 = m1.GetTranspose();
std::cout << "m2 (transp) = " << m2 << std::endl;

m1.Transpose();
std::cout << "Transposing m1 (in place) = " << m1 << std::endl;

/* OUTPUT
m1 = Rows: 2  Cols: 2
[          1,         -1,  ]
[        1.5,          3,  ]

m2 (transp) = Rows: 2  Cols: 2
[          1,        1.5,  ]
[         -1,          3,  ]

Transposing m1 (in place) = Rows: 2  Cols: 2
[          1,        1.5,  ]
[         -1,          3,  ]
*/
~~~

## Defined typedefs for easier use

MatrixNM class
~~~ c++
    typedef MatrixNM<float, 2, 2> Matrix22Flt;
    typedef MatrixNM<float, 3, 3> Matrix33Flt;
    typedef MatrixNM<float, 4, 4> Matrix44Flt;

    typedef MatrixNM<double, 2, 2> Matrix22Dbl;
    typedef MatrixNM<double, 3, 3> Matrix33Dbl;
    typedef MatrixNM<double, 4, 4> Matrix44Dbl;

    typedef MatrixNM<Complex, 2, 2> Matrix22Complex;
    typedef MatrixNM<Complex, 3, 3> Matrix33Complex;
    typedef MatrixNM<Complex, 4, 4> Matrix44Complex;

    typedef MatrixNM<float, 2, 2> Mat22F;
    typedef MatrixNM<float, 3, 3> Mat33F;
    typedef MatrixNM<float, 4, 4> Mat44F;

    typedef MatrixNM<double, 2, 2> Mat22D;
    typedef MatrixNM<double, 3, 3> Mat33D;
    typedef MatrixNM<double, 4, 4> Mat44D;

    typedef MatrixNM<Complex, 2, 2> Mat22C;
    typedef MatrixNM<Complex, 3, 3> Mat33C;
    typedef MatrixNM<Complex, 4, 4> Mat44C;
~~~