#ifdef MML_USE_SINGLE_HEADER
#include "MML.h"
#else
#include "MMLBase.h"

#include "mml/base/BaseUtils.h"
#include "mml/base/Matrix.h"
#include "mml/base/MatrixNM.h"
#include "mml/base/Polynom.h"

#include "mml/core/MatrixUtils.h"
#endif

using namespace MML;

//void FaddeevAlg(const Matrix<Real>& A, RealPolynom& outCharPoly, Real& outDet, Matrix<Real>& outInv)
//{
//	int N = A.RowNum();
//
//	Matrix<Real> nullMat(N, N);		// null matrix
//	auto identMat = Matrix<Real>::GetUnitMatrix(N);
//
//	Vector<Matrix<Real>> M(N + 1, nullMat);
//
//	outCharPoly.SetDegree(N);
//
//	M[0] = nullMat;
//	outCharPoly[N] = 1.0;
//	for (int k = 1; k <= A.RowNum(); k++)
//	{
//		M[k] = A * M[k - 1] + outCharPoly[N - k + 1] * identMat;
//		outCharPoly[N - k] = -1.0 / k * (A * M[k]).Trace();
//	}
//
//	outDet = std::pow(-1.0, N) * outCharPoly[0];
//	outInv = M[N] / outDet; // inverse of A
//}

// initialization

// operations

// matrix - vector

// matrix - matrix

// base utils

void Matrix_init()
{
	Matrix<Real> a, b(2, 2), c(2, 2, { 1.0, -2.0, 3.75, -1.0 });     
	Matrix<Real> d(c), e = c;
	
	// different types
	MatrixInt      mat_int(2, 2, { 1, 2, 3, 4 });
	MatrixComplex  mat_cmplx(2, 2, { Complex(1,1), Complex(-1,2),
																	 Complex(2,-0.5), Complex(1,1) });

	// initializing matrix with submatrix from existing matrix
	Matrix<Real> f(5, 5, { 1,  2,  3,  4,  5,
												 6,  7,  8,  9,  10,
												 11, 12, 13, 14, 15,
												 16, 17, 18, 19, 20,
												 21, 22, 23, 24, 25 });
	Matrix<Real> g(f, 3, 1, 2, 2);	

	// or using GetSubmatrix function()
	g = f.GetSubmatrix(3, 1, 2, 2);
	std::cout << "Matrix g: " << g << std::endl;

	Matrix<Real>    unit_r = Matrix<Real>::GetUnitMatrix(3);
	Matrix<Complex> unit_c = Matrix<Complex>::GetUnitMatrix(2);
	
	Vector<Real> vec_a({ 1, 2, 3, 4, 5 });
	Matrix<Real> diag = Utils::DiagonalMatrixFromVector(vec_a);
	Matrix<Real> matA = Utils::RowMatrixFromVector<Real>(vec_a);
	Matrix<Real> matB = Utils::ColumnMatrixFromVector<Real>(vec_a);
}

void Matrix_operations()
{
	Matrix<Real> m1(3, 2, { 1.0, -1.0, 
													2.0, 0.5,
													1.5, 3.0 });
	Matrix<Real> m2(2, 4, {-2.5, 3.0, 1.0, -1.0, 
											0.5,-2.0, 5.0, -3 });
	
	Vector<Real> v1({ 1.0, 2.0 }), v2({ 5.0, -1.0, 2.0, 3.0 });

	// scalar results
	Real n1 = m1.NormL1();		// L1 norm of matrix m1
	// Real tr = m1.Trace();			// trace only works on square matrices

	// vector results
	Vector<Real>  v = m1 * v1;			// matrix - vector multiplication (3x2) * (2x1)	// matrix results
	Matrix<Real> m3 = m1 * m2;			// matrix - matrix multiplication
	Matrix<Real> m4 = 2.5 * m1 + m2.GetSubmatrix(0, 0, 2, 3).GetTranspose() / 1.5; 
}

void FaddeevAlgExample()
{
	//Matrix<Real> A(4, 4, { 3, 2, 5, 5, 
	//											 1, 4, 4, 3,
	//											 3, 1, 9, 2,
	//											 4, 6, 2, 7});

	Matrix<Real> A(3, 3, { 3, 1, 5, 
												 3, 3, 1, 
												 4, 6, 4 });

	int N = A.RowNum();
	Matrix<Real> M0(N, N);		// null matrix
	Matrix<Real> I = Matrix<Real>::GetUnitMatrix(N);

	Vector<Matrix<Real>> M(N + 1, M0);
	PolynomRealFunc c(N);

	M[0] = M0;
	c[N] = 1.0;
	std::cout << "M[0] = " << M[0] << std::endl;
	std::cout << "c[" << N << "] = " << c[N] << std::endl;

	for (int k = 1; k <= A.RowNum(); k++)
	{
		M[k] = A * M[k - 1] + c[N - k + 1] * I;
		c[N - k] = -1.0 / k  * (A * M[k]).Trace();

		std::cout << "M[" << k << "] = " << M[k] << std::endl;
		std::cout << "c[" << N - k << "] = " << c[N - k] << std::endl;
	}
	Real det = std::pow(-1.0, N) * c[0];
	Matrix<Real> Ainv = M[N] / det; // inverse of A

	std::cout << "Characteristic polynomial: " << c << std::endl;	
	std::cout << "Determinant of A: " << det << std::endl;
	std::cout << "Inverse of A: " << Ainv << std::endl;

	std::cout << "A * Ainv = " << A * Ainv << std::endl;

	// call FaddeevAlg function to compute the inverse of A
	Real detValue;
	PolynomRealFunc charPoly;
	Matrix<Real> Ainv2;

	Utils::FaddeevAlg(A, charPoly, detValue, Ainv2);
	
	std::cout << "FaddeevAlg:" << std::endl;
	std::cout << "Characteristic polynomial: " << charPoly << std::endl;
	std::cout << "Determinant of A: " << detValue << std::endl;
	std::cout << "Inverse of A: " << Ainv2 << std::endl;
	std::cout << "A * Ainv2 = " << A * Ainv2 << std::endl;
}

void Demo_Matrix()
{

	std::cout << "***********************************************************************" << std::endl;
	std::cout << "****                         EXAMPLE 1 - Matrix                    ****" << std::endl;
	std::cout << "***********************************************************************" << std::endl;

	Matrix_init();
	Matrix_operations();
	FaddeevAlgExample();

	// working with views
	MML::Matrix<double> mat(5, 5, 1.0);
	// auto sub = mat.block(1, 1, 3, 3); // 3x3 view starting at (1,1)
	// sub(0, 0) = 42; // modifies mat(1,1)

	// auto row2 = mat.row?view(2);
	// row2(0, 3) = 99; // modifies mat(2,3)

	// auto col1 = mat.col_view(1);
	// col1(2, 0) = 77; // modifies mat(2,1)
	
	std::cout << "\nMatrix after view modifications:\n" << mat << std::endl;
}