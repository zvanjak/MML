#if !defined MML_MATRIX_INVERSION_H
#define MML_MATRIX_INVERSION_H

#include <mml/core/LinAlgEqSolvers/LinAlgDirect.h>
#include <mml/base/Matrix/MatrixNM.h>
#include <mml/base/Matrix/MatrixSym.h>
#include <mml/base/Matrix/MatrixTriDiag.h>

namespace MML::MatrixAlg {
	template<MMLScalar Scalar>
	Matrix<Scalar> Inverse(const Matrix<Scalar>& matrix) {
		if (matrix.rows() == 0 || matrix.rows() != matrix.cols())
			throw MatrixDimensionError("Inverse requires a nonempty square matrix", matrix.rows(), matrix.cols(), -1, -1);
		LUSolver<Scalar> solver(matrix);
		Matrix<Scalar> result;
		solver.inverse(result);
		return result;
	}

	template<MMLScalar Scalar, int Rows, int Cols>
	MatrixNM<Scalar, Rows, Cols> Inverse(const MatrixNM<Scalar, Rows, Cols>& matrix) {
		if constexpr (Rows != Cols)
			throw MatrixDimensionError("Inverse requires a square matrix", Rows, Cols, -1, -1);
		Matrix<Scalar> dynamic(Rows, Cols);
		for (int i = 0; i < Rows; ++i)
			for (int j = 0; j < Cols; ++j)
				dynamic(i, j) = matrix(i, j);
		const auto inverse = Inverse(dynamic);
		MatrixNM<Scalar, Rows, Cols> result;
		for (int i = 0; i < Rows; ++i)
			for (int j = 0; j < Cols; ++j)
				result(i, j) = inverse(i, j);
		return result;
	}

	template<MMLScalar Scalar>
	Matrix<Scalar> Inverse(const MatrixSym<Scalar>& matrix) {
		return Inverse(matrix.ToFullMatrix());
	}

	template<MMLScalar Scalar>
	Matrix<Scalar> Inverse(const TridiagonalMatrix<Scalar>& matrix) {
		return Inverse(matrix.ToFullMatrix());
	}

	template<MMLScalar Scalar>
	void Invert(Matrix<Scalar>& matrix) {
		matrix = Inverse(matrix);
	}

	template<MMLScalar Scalar, int Rows, int Cols>
	void Invert(MatrixNM<Scalar, Rows, Cols>& matrix) {
		matrix = Inverse(matrix);
	}
}

#endif // MML_MATRIX_INVERSION_H
