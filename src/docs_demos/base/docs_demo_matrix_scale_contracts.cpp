#ifdef MML_USE_SINGLE_HEADER
#include <MML.h>
#else
#include <mml/base/Matrix/Matrix.h>
#include <mml/base/Matrix/MatrixSym.h>
#include <mml/base/Matrix/MatrixTriDiag.h>
#include <mml/base/Matrix/MatrixBandDiag.h>
#include <mml/base/SparseMatrix/SparseMatrixCSR.h>
#endif

#include <iostream>
#include <limits>

void Docs_Demo_Matrix_ScaleContracts()
{
	using namespace MML;
	const Real nan = std::numeric_limits<Real>::quiet_NaN();
	const TridiagonalMatrix<Real> zero(1, {0}), invalid(1, {nan}), one(1, {1});
	const bool finiteOnly = invalid.IsEqualTo(zero);
	const bool invalidTolerance = zero.IsEqualTo(zero, nan);
	const bool inclusiveBoundary = zero.IsEqualTo(one, 1);
	if (finiteOnly || invalidTolerance || !inclusiveBoundary ||
		zero.IsEqualTo(zero, std::numeric_limits<Real>::infinity()))
		throw MatrixNumericalError("Tridiagonal equality contract demo failed");

	const Matrix<Real> nonzeroDiagonal(1, 1, {2});
	const Matrix<Complex> imaginaryDiagonal(1, 1, {Complex(0, 2)});
	if (nonzeroDiagonal.IsAntiSymmetric() || nonzeroDiagonal.IsSkewHermitian() ||
		!imaginaryDiagonal.IsSkewHermitian() || imaginaryDiagonal.IsAntiSymmetric())
		throw MatrixNumericalError("Matrix diagonal predicate demo failed");

	const Real scale = std::numeric_limits<Real>::min() * 64;
	const TridiagonalMatrix<Real> system(3,
		{4 * scale, -scale, -scale, 4 * scale, -scale, -scale, 4 * scale});
	const Vector<Real> rhs{2 * scale, 4 * scale, 10 * scale};
	const auto solution = system.Solve(rhs);
	for (int i = 0; i < 3; ++i)
		if (!IsNearlyEqual(solution[i], Real(i + 1), Precision::EqualityTolerance))
			throw MatrixNumericalError("Scaled tridiagonal solve demo failed");

	const Real large = std::numeric_limits<Real>::max() / 64;
	const MatrixSym<Real> symmetric(2, {3 * large, 4 * large, 0});
	Matrix<Real> bandData(2, 3);
	bandData(0, 1) = 3 * large;
	bandData(0, 2) = 4 * large;
	const BandDiagonalMatrix band(2, 1, 1, bandData);
	const SparseMatrix::SparseMatrixCSR<Real> sparse(2, 2,
		{3 * large, 4 * large}, {0, 1}, {0, 2, 2});
	if (!std::isfinite(symmetric.NormFrobenius()) ||
		!IsNearlyEqual(band.NormFrobenius() / large, Real{5}, Precision::EqualityTolerance) ||
		!IsNearlyEqual(sparse.NormFrobenius() / large, Real{5}, Precision::EqualityTolerance))
		throw MatrixNumericalError("Scale-safe Frobenius norm demo failed");
	std::cout << "\nMatrix scale contracts: finite predicates, scaled Thomas solve, "
		"and structured/sparse Frobenius norms verified.\n";
}
