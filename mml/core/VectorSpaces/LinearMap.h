///////////////////////////////////////////////////////////////////////////////////////////
///                         MinimalMathLibrary (MML)                                  ///
///                                                                                   ///
///  File:        VectorSpaces/LinearMap.h                                            ///
///  Description: Fixed-size linear maps with explicit domain and codomain            ///
///                                                                                   ///
///  Copyright:   (c) 2024-2026 Zvonimir Vanjak                                       ///
///  License:     MIT License (see LICENSE.md)                                         ///
///                                                                                   ///
///////////////////////////////////////////////////////////////////////////////////////////
#if !defined MML_VECTOR_SPACES_LINEAR_MAP_H
#define MML_VECTOR_SPACES_LINEAR_MAP_H

#include <mml/MMLBase.h>

#include <mml/base/Matrix/MatrixNM.h>
#include <mml/base/Vector/VectorN.h>
#include <mml/core/LinAlgEqSolvers/LinAlgSVD.h>
#include <mml/core/VectorSpaces/Basis.h>
#include <mml/core/VectorSpaces/Subspace.h>

#include <type_traits>

namespace MML::VectorSpaces
{
	namespace Detail
	{
		template<class Scalar, int Rows, int Cols>
		Matrix<Scalar> ToDynamicMatrix(const MatrixNM<Scalar, Rows, Cols>& matrix)
		{
			Matrix<Scalar> result(Rows, Cols);
			for (int i = 0; i < Rows; i++)
				for (int j = 0; j < Cols; j++)
					result(i, j) = matrix(i, j);
			return result;
		}
	}

	template<class Scalar, int DomainN, int CodomainN>
	class LinearMap
	{
		MatrixNM<Scalar, CodomainN, DomainN> _matrix;

	public:
		using scalar_type = Scalar;
		static constexpr int DomainDimension = DomainN;
		static constexpr int CodomainDimension = CodomainN;

		LinearMap() : _matrix() { }

		explicit LinearMap(const MatrixNM<Scalar, CodomainN, DomainN>& matrix)
			: _matrix(matrix)
		{
		}

		static LinearMap Zero()
		{
			return LinearMap();
		}

		static LinearMap Identity() requires (DomainN == CodomainN)
		{
			return LinearMap(MatrixNM<Scalar, CodomainN, DomainN>::Identity());
		}

		const MatrixNM<Scalar, CodomainN, DomainN>& matrix() const noexcept { return _matrix; }

		VectorN<Scalar, CodomainN> apply(const VectorN<Scalar, DomainN>& vector) const
		{
			return _matrix * vector;
		}

		VectorN<Scalar, CodomainN> operator()(const VectorN<Scalar, DomainN>& vector) const
		{
			return apply(vector);
		}

		int rank(Real tolerance = -1.0) const requires MMLReal<Scalar>
		{
			SVDecompositionSolver<Scalar> decomposition(Detail::ToDynamicMatrix(_matrix));
			return decomposition.Rank(tolerance);
		}

		int nullity(Real tolerance = -1.0) const requires MMLReal<Scalar>
		{
			return DomainN - rank(tolerance);
		}

		Subspace<Scalar> kernel(Real tolerance = -1.0) const requires MMLReal<Scalar>
		{
			return Subspace<Scalar>::NullSpace(Detail::ToDynamicMatrix(_matrix), tolerance);
		}

		Subspace<Scalar> image(Real tolerance = -1.0) const requires MMLReal<Scalar>
		{
			return Subspace<Scalar>::ColumnSpace(Detail::ToDynamicMatrix(_matrix), tolerance);
		}

		LinearMap inverse() const requires (DomainN == CodomainN)
		{
			return LinearMap(_matrix.GetInverse());
		}

		MatrixNM<Scalar, CodomainN, DomainN> matrixInBases(const Basis<Scalar, DomainN>& domainBasis,
			const Basis<Scalar, CodomainN>& codomainBasis) const
		{
			return codomainBasis.inverseMatrix() * _matrix * domainBasis.matrix();
		}
	};

	template<class Scalar, int DomainN, int MiddleN, int CodomainN>
	LinearMap<Scalar, DomainN, CodomainN> Compose(const LinearMap<Scalar, MiddleN, CodomainN>& g,
		const LinearMap<Scalar, DomainN, MiddleN>& f)
	{
		return LinearMap<Scalar, DomainN, CodomainN>(g.matrix() * f.matrix());
	}

	template<class Scalar, int DomainN, int MiddleN, int CodomainN>
	LinearMap<Scalar, DomainN, CodomainN> compose(const LinearMap<Scalar, MiddleN, CodomainN>& g,
		const LinearMap<Scalar, DomainN, MiddleN>& f)
	{
		return Compose(g, f);
	}
}

#endif // MML_VECTOR_SPACES_LINEAR_MAP_H