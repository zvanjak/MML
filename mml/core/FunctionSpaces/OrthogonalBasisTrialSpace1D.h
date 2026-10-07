///////////////////////////////////////////////////////////////////////////////////////////
///                         MinimalMathLibrary (MML)                                  ///
///                                                                                   ///
///  File:        FunctionSpaces/OrthogonalBasisTrialSpace1D.h                        ///
///  Description: Trial-space wrapper for existing orthogonal bases                   ///
///                                                                                   ///
///  Copyright:   (c) 2024-2026 Zvonimir Vanjak                                       ///
///  License:     MIT License (see LICENSE.md)                                         ///
///                                                                                   ///
///////////////////////////////////////////////////////////////////////////////////////////
#if !defined MML_ORTHOGONAL_BASIS_TRIAL_SPACE_1D_H
#define MML_ORTHOGONAL_BASIS_TRIAL_SPACE_1D_H

#include <mml/core/FunctionSpaces/TrialSpace1D.h>
#include <mml/core/OrthogonalBasis.h>

#include <cmath>
#include <string>

namespace MML::FunctionSpaces
{
	class OrthogonalBasisFunctionSpace1D : public FunctionSpace1D
	{
		const OrthogonalBasis& _basis;

	public:
		explicit OrthogonalBasisFunctionSpace1D(const OrthogonalBasis& basis)
			: _basis(basis)
		{
			ValidateFiniteDomain(_basis.DomainMin(), _basis.DomainMax(), "OrthogonalBasisFunctionSpace1D");
		}

		const OrthogonalBasis& basis() const noexcept { return _basis; }

		Real domainMin() const noexcept override { return _basis.DomainMin(); }
		Real domainMax() const noexcept override { return _basis.DomainMax(); }
		Real weight(Real x) const override { return _basis.WeightFunction(x); }
	};

	class OrthogonalBasisTrialSpace1D : public TrialSpace1D
	{
		const OrthogonalBasis& _basis;
		OrthogonalBasisFunctionSpace1D _space;
		int _dimension;

	public:
		OrthogonalBasisTrialSpace1D(const OrthogonalBasis& basis, int dimension)
			: _basis(basis), _space(basis), _dimension(dimension)
		{
			if (_dimension <= 0)
				throw FunctionSpaceInputError("OrthogonalBasisTrialSpace1D: dimension must be positive");
		}

		const OrthogonalBasis& basis() const noexcept { return _basis; }
		const FunctionSpace1D& functionSpace() const noexcept override { return _space; }
		int dimension() const noexcept override { return _dimension; }

		Real basisValue(int basisIndex, Real x) const override
		{
			validateBasisIndex(basisIndex, "OrthogonalBasisTrialSpace1D::basisValue");
			return _basis.Evaluate(basisIndex, x);
		}
	};
}

#endif // MML_ORTHOGONAL_BASIS_TRIAL_SPACE_1D_H
