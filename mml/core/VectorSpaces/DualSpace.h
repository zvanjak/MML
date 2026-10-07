///////////////////////////////////////////////////////////////////////////////////////////
///                         MinimalMathLibrary (MML)                                  ///
///                                                                                   ///
///  File:        VectorSpaces/DualSpace.h                                            ///
///  Description: Dual spaces, vectors-in-space, and covector bridge APIs             ///
///                                                                                   ///
///  Copyright:   (c) 2024-2026 Zvonimir Vanjak                                       ///
///  License:     MIT License (see LICENSE.md)                                         ///
///                                                                                   ///
///////////////////////////////////////////////////////////////////////////////////////////
#if !defined MML_VECTOR_SPACES_DUAL_SPACE_H
#define MML_VECTOR_SPACES_DUAL_SPACE_H

#include <mml/MMLBase.h>

#include <mml/base/Tensor/Rank1Tensor.h>
#include <mml/core/VectorSpaces/Basis.h>
#include <mml/core/VectorSpaces/LinearMap.h>
#include <mml/core/VectorSpaces/Subspace.h>

namespace MML::VectorSpaces
{
	template<class Scalar, int N>
	class VectorInSpace
	{
		VectorN<Scalar, N> _coordinates;
		Basis<Scalar, N> _basis;

	public:
		VectorInSpace(const VectorN<Scalar, N>& coordinates, const Basis<Scalar, N>& basis = Basis<Scalar, N>::Standard())
			: _coordinates(coordinates), _basis(basis)
		{
		}

		const VectorN<Scalar, N>& coordinates() const noexcept { return _coordinates; }
		const Basis<Scalar, N>& basis() const noexcept { return _basis; }

		VectorN<Scalar, N> coordinatesInStandard() const
		{
			return _basis.coordinatesInStandard(_coordinates);
		}

		VectorN<Scalar, N> coordinatesIn(const Basis<Scalar, N>& targetBasis) const
		{
			return ChangeBasis(_coordinates, _basis, targetBasis);
		}
	};

	template<class Scalar, int N>
	class CovectorInSpace
	{
		VectorN<Scalar, N> _coordinates;
		Basis<Scalar, N> _basis;

	public:
		CovectorInSpace(const VectorN<Scalar, N>& coordinates, const Basis<Scalar, N>& basis = Basis<Scalar, N>::Standard())
			: _coordinates(coordinates), _basis(basis)
		{
		}

		const VectorN<Scalar, N>& coordinates() const noexcept { return _coordinates; }
		const Basis<Scalar, N>& basis() const noexcept { return _basis; }

		Scalar operator()(const VectorInSpace<Scalar, N>& vector) const
		{
			VectorN<Scalar, N> vectorCoords = vector.coordinatesIn(_basis);
			Scalar result{};
			for (int i = 0; i < N; i++)
				result += _coordinates[i] * vectorCoords[i];
			return result;
		}

		VectorN<Scalar, N> coordinatesInStandardDual() const
		{
			return _coordinates * _basis.inverseMatrix();
		}

		template<class FrameTag>
		Covector<N, FrameTag> toCovector() const
		{
			return Covector<N, FrameTag>(coordinatesInStandardDual());
		}
	};

	template<class Scalar, int N>
	class DualSpace
	{
		Basis<Scalar, N> _basis;

	public:
		explicit DualSpace(const Basis<Scalar, N>& basis = Basis<Scalar, N>::Standard())
			: _basis(basis)
		{
		}

		const Basis<Scalar, N>& primalBasis() const noexcept { return _basis; }

		CovectorInSpace<Scalar, N> covector(const VectorN<Scalar, N>& coordinates) const
		{
			return CovectorInSpace<Scalar, N>(coordinates, _basis);
		}

		CovectorInSpace<Scalar, N> dualBasisCovector(int index) const
		{
			if (index < 0 || index >= N)
				throw VectorDimensionError("DualSpace::dualBasisCovector - index out of bounds", N, index);

			VectorN<Scalar, N> coords;
			coords[index] = Scalar{ 1 };
			return covector(coords);
		}

		template<class FrameTag>
		CovectorInSpace<Scalar, N> fromCovector(const Covector<N, FrameTag>& covector) const
		{
			VectorN<Scalar, N> coords = covector.components() * _basis.matrix();
			return CovectorInSpace<Scalar, N>(coords, _basis);
		}

		template<int DomainN>
		CovectorInSpace<Scalar, DomainN> pullback(const LinearMap<Scalar, DomainN, N>& map,
			const CovectorInSpace<Scalar, N>& covector,
			const Basis<Scalar, DomainN>& domainBasis = Basis<Scalar, DomainN>::Standard()) const
		{
			MatrixNM<Scalar, N, DomainN> representation = map.matrixInBases(domainBasis, covector.basis());
			VectorN<Scalar, DomainN> pulledCoords = covector.coordinates() * representation;
			return CovectorInSpace<Scalar, DomainN>(pulledCoords, domainBasis);
		}

		Subspace<Scalar> annihilator(const Subspace<Scalar>& subspace, Real tolerance = -1.0) const
		{
			if (subspace.ambientDimension() != N)
				throw MatrixDimensionError("DualSpace::annihilator - subspace ambient dimension mismatch",
					subspace.ambientDimension(), N, -1, -1);

			Matrix<Scalar> constraints = subspace.orthonormalBasis().transpose();
			return Subspace<Scalar>::NullSpace(constraints, tolerance);
		}
	};

	template<class Scalar, int N>
	DualSpace<Scalar, N> dual(const Basis<Scalar, N>& basis = Basis<Scalar, N>::Standard())
	{
		return DualSpace<Scalar, N>(basis);
	}
}

#endif // MML_VECTOR_SPACES_DUAL_SPACE_H