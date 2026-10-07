///////////////////////////////////////////////////////////////////////////////////////////
///                         MinimalMathLibrary (MML)                                  ///
///                                                                                   ///
///  File:        VectorSpaces/InnerProductSpace.h                                    ///
///  Description: Fixed-size inner-product spaces and orthogonal operations           ///
///                                                                                   ///
///  Copyright:   (c) 2024-2026 Zvonimir Vanjak                                       ///
///  License:     MIT License (see LICENSE.md)                                         ///
///                                                                                   ///
///////////////////////////////////////////////////////////////////////////////////////////
#if !defined MML_VECTOR_SPACES_INNER_PRODUCT_SPACE_H
#define MML_VECTOR_SPACES_INNER_PRODUCT_SPACE_H

#include <mml/MMLBase.h>

#include <mml/base/DifferentialGeometry/Metric.h>
#include <mml/base/Matrix/MatrixNM.h>
#include <mml/base/Vector/Vector.h>
#include <mml/base/Vector/VectorN.h>
#include <mml/core/LinAlgEqSolvers/LinAlgDirect.h>
#include <mml/core/VectorSpaces/Basis.h>
#include <mml/core/VectorSpaces/LinearMap.h>
#include <mml/core/VectorSpaces/Subspace.h>

#include <type_traits>

namespace MML::VectorSpaces
{
	template<class Scalar, int N>
		requires MMLArithmetic<Scalar>
	class InnerProductSpace
	{
		MatrixNM<Scalar, N, N> _gram;
		MatrixNM<Scalar, N, N> _inverseGram;

		void validate(Real tolerance)
		{
			if (!isSymmetric(tolerance))
				throw MatrixNumericalError("InnerProductSpace Gram matrix must be symmetric");
			if (!isPositiveDefinite(tolerance))
				throw MatrixNumericalError("InnerProductSpace Gram matrix must be positive definite; indefinite metrics are represented by Metric, not InnerProductSpace");
			_inverseGram = _gram.GetInverse();
		}

		static Vector<Scalar> ZeroVector(int size)
		{
			return Vector<Scalar>(size, Scalar{});
		}

		Matrix<Scalar> dynamicGram() const
		{
			Matrix<Scalar> result(N, N);
			for (int i = 0; i < N; i++)
				for (int j = 0; j < N; j++)
					result(i, j) = _gram(i, j);
			return result;
		}

		void ensureAmbientDimension(int ambientDimension) const
		{
			if (ambientDimension != N)
				throw MatrixDimensionError("InnerProductSpace ambient dimension mismatch", N, ambientDimension, -1, -1);
		}

	public:
		using scalar_type = Scalar;
		static constexpr int Dimension = N;

		InnerProductSpace()
			: _gram(MatrixNM<Scalar, N, N>::Identity()), _inverseGram(MatrixNM<Scalar, N, N>::Identity())
		{
		}

		explicit InnerProductSpace(const MatrixNM<Scalar, N, N>& gram, Real tolerance = Defaults::MatrixIsEqualTolerance)
			: _gram(gram), _inverseGram()
		{
			validate(tolerance);
		}

		template<class FrameTag>
		explicit InnerProductSpace(const Metric<N, FrameTag>& metric, Real tolerance = Defaults::MatrixIsEqualTolerance)
			requires std::same_as<Scalar, Real>
			: InnerProductSpace(metric.covariantComponents(), tolerance)
		{
		}

		static InnerProductSpace Euclidean()
		{
			return InnerProductSpace();
		}

		const MatrixNM<Scalar, N, N>& gramMatrix() const noexcept { return _gram; }
		const MatrixNM<Scalar, N, N>& inverseGramMatrix() const noexcept { return _inverseGram; }

		bool isSymmetric(Real tolerance = Defaults::MatrixIsEqualTolerance) const
		{
			for (int i = 0; i < N; i++)
				for (int j = i + 1; j < N; j++)
					if (std::abs(_gram(i, j) - _gram(j, i)) > Scalar(tolerance))
						return false;
			return true;
		}

		bool isPositiveDefinite(Real tolerance = Defaults::MatrixIsEqualTolerance) const
		{
			long double lower[N][N] = {};
			for (int i = 0; i < N; i++) {
				for (int j = 0; j <= i; j++) {
					long double sum = static_cast<long double>(_gram(i, j));
					for (int k = 0; k < j; k++)
						sum -= lower[i][k] * lower[j][k];

					if (i == j) {
						if (sum <= static_cast<long double>(tolerance))
							return false;
						lower[i][j] = std::sqrt(sum);
					}
					else {
						lower[i][j] = sum / lower[j][j];
					}
				}
			}
			return true;
		}

		Scalar inner(const VectorN<Scalar, N>& a, const VectorN<Scalar, N>& b) const
		{
			VectorN<Scalar, N> gramB = _gram * b;
			Scalar result{};
			for (int i = 0; i < N; i++)
				result += a[i] * gramB[i];
			return result;
		}

		Real norm(const VectorN<Scalar, N>& vector) const
		{
			Scalar squaredNorm = inner(vector, vector);
			if (squaredNorm < Scalar{})
				throw MatrixNumericalError("InnerProductSpace::norm - negative squared norm; Gram matrix is not positive definite within tolerance");
			return std::sqrt(static_cast<Real>(squaredNorm));
		}

		Real distance(const VectorN<Scalar, N>& a, const VectorN<Scalar, N>& b) const
		{
			return norm(a - b);
		}

		MatrixNM<Scalar, N, N> gramMatrixInBasis(const Basis<Scalar, N>& basis) const
		{
			return basis.matrix().transpose() * _gram * basis.matrix();
		}

		Basis<Scalar, N> orthonormalize(const Basis<Scalar, N>& basis, Real tolerance = Defaults::MatrixIsEqualTolerance) const
		{
			MatrixNM<Scalar, N, N> orthonormalColumns;

			for (int col = 0; col < N; col++) {
				VectorN<Scalar, N> vector;
				for (int row = 0; row < N; row++)
					vector[row] = basis.matrix()(row, col);

				for (int prev = 0; prev < col; prev++) {
					VectorN<Scalar, N> q;
					for (int row = 0; row < N; row++)
						q[row] = orthonormalColumns(row, prev);
					Scalar coefficient = inner(q, vector);
					vector = vector - coefficient * q;
				}

				Real length = norm(vector);
				if (length <= tolerance)
					throw SingularMatrixError("InnerProductSpace::orthonormalize - basis vectors are linearly dependent");

				for (int row = 0; row < N; row++)
					orthonormalColumns(row, col) = vector[row] / Scalar(length);
			}

			return Basis<Scalar, N>(orthonormalColumns, tolerance);
		}

		Vector<Scalar> orthogonalProjection(const Subspace<Scalar>& subspace, const Vector<Scalar>& vector) const
		{
			ensureAmbientDimension(subspace.ambientDimension());
			if (vector.size() != N)
				throw VectorDimensionError("InnerProductSpace::orthogonalProjection - vector size must match space dimension", N, vector.size());
			if (subspace.dimension() == 0)
				return ZeroVector(N);

			const Matrix<Scalar>& basis = subspace.orthonormalBasis();
			Matrix<Scalar> gram = dynamicGram();
			Matrix<Scalar> normal(subspace.dimension(), subspace.dimension());
			Vector<Scalar> rhs(subspace.dimension());

			for (int i = 0; i < subspace.dimension(); i++) {
				for (int row = 0; row < N; row++)
					for (int col = 0; col < N; col++)
						rhs[i] += basis(row, i) * gram(row, col) * vector[col];

				for (int j = 0; j < subspace.dimension(); j++)
					for (int row = 0; row < N; row++)
						for (int col = 0; col < N; col++)
							normal(i, j) += basis(row, i) * gram(row, col) * basis(col, j);
			}

			LUSolver<Scalar> solver(normal);
			Vector<Scalar> coefficients = solver.Solve(rhs);
			Vector<Scalar> projection = ZeroVector(N);
			for (int j = 0; j < subspace.dimension(); j++)
				for (int i = 0; i < N; i++)
					projection[i] += coefficients[j] * basis(i, j);
			return projection;
		}

		Subspace<Scalar> orthogonalComplement(const Subspace<Scalar>& subspace, Real tolerance = -1.0) const
		{
			ensureAmbientDimension(subspace.ambientDimension());
			if (subspace.dimension() == 0)
				return Subspace<Scalar>(N, Matrix<Scalar>::Identity(N));

			const Matrix<Scalar>& basis = subspace.orthonormalBasis();
			Matrix<Scalar> gram = dynamicGram();
			Matrix<Scalar> constraints(subspace.dimension(), N);
			for (int i = 0; i < subspace.dimension(); i++)
				for (int j = 0; j < N; j++)
					for (int k = 0; k < N; k++)
						constraints(i, j) += basis(k, i) * gram(k, j);

			return Subspace<Scalar>::NullSpace(constraints, tolerance);
		}

		template<int CodomainN>
		LinearMap<Scalar, CodomainN, N> adjoint(const LinearMap<Scalar, N, CodomainN>& map,
			const InnerProductSpace<Scalar, CodomainN>& codomainSpace) const
		{
			return LinearMap<Scalar, CodomainN, N>(_inverseGram * map.matrix().transpose() * codomainSpace.gramMatrix());
		}
	};

	template<class Scalar, int DomainN, int CodomainN>
	LinearMap<Scalar, CodomainN, DomainN> Adjoint(const LinearMap<Scalar, DomainN, CodomainN>& map,
		const InnerProductSpace<Scalar, DomainN>& domainSpace,
		const InnerProductSpace<Scalar, CodomainN>& codomainSpace)
	{
		return domainSpace.adjoint(map, codomainSpace);
	}

	template<class Scalar, int DomainN, int CodomainN>
	LinearMap<Scalar, CodomainN, DomainN> adjoint(const LinearMap<Scalar, DomainN, CodomainN>& map,
		const InnerProductSpace<Scalar, DomainN>& domainSpace,
		const InnerProductSpace<Scalar, CodomainN>& codomainSpace)
	{
		return Adjoint(map, domainSpace, codomainSpace);
	}
}

#endif // MML_VECTOR_SPACES_INNER_PRODUCT_SPACE_H