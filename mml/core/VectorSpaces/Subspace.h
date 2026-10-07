///////////////////////////////////////////////////////////////////////////////////////////
///                         MinimalMathLibrary (MML)                                  ///
///                                                                                   ///
///  File:        VectorSpaces/Subspace.h                                             ///
///  Description: Dynamic finite-dimensional subspace wrapper                         ///
///                                                                                   ///
///  Copyright:   (c) 2024-2026 Zvonimir Vanjak                                       ///
///  License:     MIT License (see LICENSE.md)                                         ///
///                                                                                   ///
///////////////////////////////////////////////////////////////////////////////////////////
#if !defined MML_VECTOR_SPACES_SUBSPACE_H
#define MML_VECTOR_SPACES_SUBSPACE_H

#include <mml/MMLBase.h>

#include <mml/base/Matrix/Matrix.h>
#include <mml/base/Vector/Vector.h>
#include <mml/core/LinAlgEqSolvers/LinAlgSVD.h>

namespace MML::VectorSpaces
{
	template<class Scalar>
	class Subspace
	{
		int _ambientDimension;
		Matrix<Scalar> _orthonormalBasis;

		static Scalar Dot(const Vector<Scalar>& a, const Vector<Scalar>& b)
		{
			if (a.size() != b.size())
				throw VectorDimensionError("Subspace::Dot - vectors must have same size", a.size(), b.size());

			Scalar result{};
			for (int i = 0; i < a.size(); i++)
				result += a[i] * b[i];
			return result;
		}

		static Vector<Scalar> ZeroVector(int size)
		{
			return Vector<Scalar>(size, Scalar{});
		}

		static Matrix<Scalar> ConcatenateColumns(const Matrix<Scalar>& left, const Matrix<Scalar>& right, int ambientDimension)
		{
			int leftCols = left.rows() == 0 ? 0 : left.cols();
			int rightCols = right.rows() == 0 ? 0 : right.cols();
			Matrix<Scalar> result(ambientDimension, leftCols + rightCols);

			for (int j = 0; j < leftCols; j++)
				for (int i = 0; i < ambientDimension; i++)
					result(i, j) = left(i, j);

			for (int j = 0; j < rightCols; j++)
				for (int i = 0; i < ambientDimension; i++)
					result(i, leftCols + j) = right(i, j);

			return result;
		}

		static Matrix<Scalar> NullSpaceBasis(const Matrix<Scalar>& matrix, Real tolerance) requires MMLReal<Scalar>
		{
			SVDecompositionSolver<Scalar> decomposition(matrix);
			return decomposition.Nullspace(tolerance);
		}

		static Matrix<Scalar> ColumnSpaceBasis(const Matrix<Scalar>& matrix, Real tolerance) requires MMLReal<Scalar>
		{
			SVDecompositionSolver<Scalar> decomposition(matrix);
			return decomposition.Range(tolerance);
		}

	public:
		Subspace() : _ambientDimension(0), _orthonormalBasis() { }

		Subspace(int ambientDimension, const Matrix<Scalar>& orthonormalBasisColumns)
			: _ambientDimension(ambientDimension), _orthonormalBasis(orthonormalBasisColumns)
		{
			if (ambientDimension < 0)
				throw ArgumentError("Subspace ambient dimension must be non-negative");

			if (_orthonormalBasis.rows() != 0 && _orthonormalBasis.rows() != _ambientDimension)
				throw MatrixDimensionError("Subspace basis row count must match ambient dimension",
					_orthonormalBasis.rows(), _orthonormalBasis.cols(), _ambientDimension, -1);
		}

		int dimension() const noexcept { return _orthonormalBasis.rows() == 0 ? 0 : _orthonormalBasis.cols(); }
		int ambientDimension() const noexcept { return _ambientDimension; }

		const Matrix<Scalar>& basis() const noexcept { return _orthonormalBasis; }
		const Matrix<Scalar>& orthonormalBasis() const noexcept { return _orthonormalBasis; }

		bool contains(const Vector<Scalar>& vector, Real tolerance = Defaults::VectorIsEqualTolerance) const
		{
			return residual(vector).NormL2() <= tolerance;
		}

		Vector<Scalar> project(const Vector<Scalar>& vector) const
		{
			if (vector.size() != _ambientDimension)
				throw VectorDimensionError("Subspace::project - vector size must match ambient dimension",
					_ambientDimension, vector.size());

			Vector<Scalar> result = ZeroVector(_ambientDimension);
			for (int j = 0; j < dimension(); j++) {
				Vector<Scalar> basisVector = _orthonormalBasis.VectorFromColumn(j);
				Scalar coefficient = Dot(basisVector, vector);
				for (int i = 0; i < _ambientDimension; i++)
					result[i] += coefficient * basisVector[i];
			}
			return result;
		}

		Vector<Scalar> residual(const Vector<Scalar>& vector) const
		{
			return vector - project(vector);
		}

		Real distance(const Vector<Scalar>& vector) const
		{
			return residual(vector).NormL2();
		}

		Subspace orthogonalComplement(Real tolerance = -1.0) const requires MMLReal<Scalar>
		{
			if (dimension() == 0)
				return Subspace(_ambientDimension, Matrix<Scalar>::Identity(_ambientDimension));

			Matrix<Scalar> complement = NullSpaceBasis(_orthonormalBasis.transpose(), tolerance);
			return Subspace(_ambientDimension, complement);
		}

		Subspace sum(const Subspace& other, Real tolerance = 1e-12) const requires MMLReal<Scalar>
		{
			if (_ambientDimension != other._ambientDimension)
				throw MatrixDimensionError("Subspace::sum - ambient dimensions must match",
					_ambientDimension, other._ambientDimension, -1, -1);

			Matrix<Scalar> combined = ConcatenateColumns(_orthonormalBasis, other._orthonormalBasis, _ambientDimension);
			Matrix<Scalar> orthonormal = ColumnSpaceBasis(combined, tolerance);
			return Subspace(_ambientDimension, orthonormal);
		}

		static Subspace NullSpace(const Matrix<Scalar>& matrix, Real tolerance = -1.0) requires MMLReal<Scalar>
		{
			return Subspace(matrix.cols(), NullSpaceBasis(matrix, tolerance));
		}

		static Subspace ColumnSpace(const Matrix<Scalar>& matrix, Real tolerance = -1.0) requires MMLReal<Scalar>
		{
			return Subspace(matrix.rows(), ColumnSpaceBasis(matrix, tolerance));
		}

		static Subspace RowSpace(const Matrix<Scalar>& matrix, Real tolerance = -1.0) requires MMLReal<Scalar>
		{
			return Subspace(matrix.cols(), ColumnSpaceBasis(matrix.transpose(), tolerance));
		}

		static Subspace LeftNullSpace(const Matrix<Scalar>& matrix, Real tolerance = -1.0) requires MMLReal<Scalar>
		{
			return Subspace(matrix.rows(), NullSpaceBasis(matrix.transpose(), tolerance));
		}
	};
}

#endif // MML_VECTOR_SPACES_SUBSPACE_H