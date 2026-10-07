///////////////////////////////////////////////////////////////////////////////////////////
///                         MinimalMathLibrary (MML)                                  ///
///                                                                                   ///
///  File:        FunctionSpaces/FunctionExpansion1D.h                                ///
///  Description: One-dimensional finite function expansions                          ///
///                                                                                   ///
///  Copyright:   (c) 2024-2026 Zvonimir Vanjak                                       ///
///  License:     MIT License (see LICENSE.md)                                         ///
///                                                                                   ///
///////////////////////////////////////////////////////////////////////////////////////////
#if !defined MML_FUNCTION_EXPANSION_1D_H
#define MML_FUNCTION_EXPANSION_1D_H

#include <mml/core/FunctionSpaces/TrialSpace1D.h>

#include <mml/base/Vector/Vector.h>

namespace MML::FunctionSpaces
{
	class FunctionExpansion1D
	{
		const TrialSpace1D* _space;
		Vector<Real> _coefficients;

		void validate() const
		{
			if (_space == nullptr)
				throw FunctionSpaceInputError("FunctionExpansion1D: trial space cannot be null");
			if (_coefficients.size() != _space->dimension())
				throw FunctionSpaceInputError("FunctionExpansion1D: coefficient count must match trial-space dimension");
		}

	public:
		FunctionExpansion1D(const TrialSpace1D& space, const Vector<Real>& coefficients)
			: _space(&space), _coefficients(coefficients)
		{
			validate();
		}

		const TrialSpace1D& trialSpace() const noexcept { return *_space; }
		const FunctionSpace1D& functionSpace() const noexcept { return _space->functionSpace(); }
		const Vector<Real>& coefficients() const noexcept { return _coefficients; }
		int dimension() const noexcept { return _coefficients.size(); }

		Real evaluate(Real x) const
		{
			Real value = REAL(0.0);
			for (int index = 0; index < _coefficients.size(); ++index)
				value += _coefficients[index] * _space->basisValue(index, x);
			return value;
		}

		Real derivative(Real x, int order = 1) const
		{
			if (order < 0)
				throw FunctionSpaceInputError("FunctionExpansion1D::derivative: derivative order must be non-negative");
			if (order == 0)
				return evaluate(x);
			if (!_space->hasBasisDerivative(order))
				throw FunctionSpaceError("FunctionExpansion1D::derivative: trial space does not provide requested derivative order");

			Real value = REAL(0.0);
			for (int index = 0; index < _coefficients.size(); ++index)
				value += _coefficients[index] * _space->basisDerivative(index, order, x);
			return value;
		}
	};
}

#endif // MML_FUNCTION_EXPANSION_1D_H
