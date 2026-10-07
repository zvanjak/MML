///////////////////////////////////////////////////////////////////////////////////////////
///                         MinimalMathLibrary (MML)                                  ///
///                                                                                   ///
///  File:        FunctionSpaces/LinearDifferentialOperator1D.h                       ///
///  Description: Scalar linear differential operators in one dimension               ///
///                                                                                   ///
///  Copyright:   (c) 2024-2026 Zvonimir Vanjak                                       ///
///  License:     MIT License (see LICENSE.md)                                         ///
///                                                                                   ///
///////////////////////////////////////////////////////////////////////////////////////////
#if !defined MML_LINEAR_DIFFERENTIAL_OPERATOR_1D_H
#define MML_LINEAR_DIFFERENTIAL_OPERATOR_1D_H

#include <mml/core/FunctionSpaces/TrialSpace1D.h>

#include <functional>
#include <string>
#include <utility>
#include <vector>

namespace MML::FunctionSpaces
{
	class LinearDifferentialOperator1D
	{
		std::vector<std::function<Real(Real)>> _coefficients;

		void validateOrder(int order, const char* context) const
		{
			if (order < 0 || order >= static_cast<int>(_coefficients.size()))
				throw FunctionSpaceInputError(std::string(context) + ": derivative order out of range");
		}

	public:
		explicit LinearDifferentialOperator1D(std::vector<std::function<Real(Real)>> coefficients)
			: _coefficients(std::move(coefficients))
		{
			if (_coefficients.empty())
				throw FunctionSpaceInputError("LinearDifferentialOperator1D: at least one coefficient is required");
			for (const auto& coefficient : _coefficients)
				if (!coefficient)
					throw FunctionSpaceInputError("LinearDifferentialOperator1D: coefficient functions must be callable");
		}

		static LinearDifferentialOperator1D Multiplication(std::function<Real(Real)> a0)
		{
			return LinearDifferentialOperator1D({ std::move(a0) });
		}

		static LinearDifferentialOperator1D FirstOrder(std::function<Real(Real)> a1, std::function<Real(Real)> a0)
		{
			return LinearDifferentialOperator1D({ std::move(a0), std::move(a1) });
		}

		static LinearDifferentialOperator1D SecondOrder(std::function<Real(Real)> a2, std::function<Real(Real)> a1, std::function<Real(Real)> a0)
		{
			return LinearDifferentialOperator1D({ std::move(a0), std::move(a1), std::move(a2) });
		}

		int order() const noexcept { return static_cast<int>(_coefficients.size()) - 1; }

		Real coefficient(int derivativeOrder, Real x) const
		{
			validateOrder(derivativeOrder, "LinearDifferentialOperator1D::coefficient");
			return _coefficients[derivativeOrder](x);
		}

		Real applyToBasis(const TrialSpace1D& trialSpace, int basisIndex, Real x) const
		{
			Real value = REAL(0.0);
			for (int derivativeOrder = 0; derivativeOrder <= order(); ++derivativeOrder) {
				if (!trialSpace.hasBasisDerivative(derivativeOrder))
					throw FunctionSpaceError("LinearDifferentialOperator1D::applyToBasis: trial space does not provide required derivative order");
				value += coefficient(derivativeOrder, x) * trialSpace.basisDerivative(basisIndex, derivativeOrder, x);
			}
			return value;
		}
	};
}

#endif // MML_LINEAR_DIFFERENTIAL_OPERATOR_1D_H
