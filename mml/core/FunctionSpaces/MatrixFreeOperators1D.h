///////////////////////////////////////////////////////////////////////////////////////////
///                         MinimalMathLibrary (MML)                                  ///
///                                                                                   ///
///  File:        FunctionSpaces/MatrixFreeOperators1D.h                              ///
///  Description: Matrix-free finite-difference operators in one dimension            ///
///                                                                                   ///
///  Copyright:   (c) 2024-2026 Zvonimir Vanjak                                       ///
///  License:     MIT License (see LICENSE.md)                                         ///
///                                                                                   ///
///////////////////////////////////////////////////////////////////////////////////////////
#if !defined MML_FUNCTION_SPACES_MATRIX_FREE_OPERATORS_1D_H
#define MML_FUNCTION_SPACES_MATRIX_FREE_OPERATORS_1D_H

#include <mml/core/FunctionSpaces/LinearOperator.h>

namespace MML::FunctionSpaces
{
	class DirichletSecondDerivativeOperator1D : public LinearOperator
	{
		int _interiorPoints;
		Real _h;
		Real _leftBoundaryValue;
		Real _rightBoundaryValue;

	public:
		using LinearOperator::apply;

		DirichletSecondDerivativeOperator1D(int interiorPoints, Real h, Real leftBoundaryValue = REAL(0.0), Real rightBoundaryValue = REAL(0.0))
			: _interiorPoints(interiorPoints), _h(h), _leftBoundaryValue(leftBoundaryValue), _rightBoundaryValue(rightBoundaryValue)
		{
			if (_interiorPoints <= 0)
				throw FunctionSpaceInputError("DirichletSecondDerivativeOperator1D: interior point count must be positive");
			if (!std::isfinite(_h) || _h <= REAL(0.0))
				throw FunctionSpaceInputError("DirichletSecondDerivativeOperator1D: grid spacing must be positive and finite");
		}

		int rows() const noexcept override { return _interiorPoints; }
		int cols() const noexcept override { return _interiorPoints; }

		Real gridSpacing() const noexcept { return _h; }
		Real leftBoundaryValue() const noexcept { return _leftBoundaryValue; }
		Real rightBoundaryValue() const noexcept { return _rightBoundaryValue; }

		void apply(const Vector<Real>& input, Vector<Real>& output) const override
		{
			validateInputOutput(input, output, "DirichletSecondDerivativeOperator1D::apply");
			Real h2inv = REAL(1.0) / (_h * _h);
			for (int index = 0; index < _interiorPoints; ++index) {
				Real left = (index == 0) ? _leftBoundaryValue : input[index - 1];
				Real right = (index == _interiorPoints - 1) ? _rightBoundaryValue : input[index + 1];
				output[index] = (left - REAL(2.0) * input[index] + right) * h2inv;
			}
		}
	};
}

#endif // MML_FUNCTION_SPACES_MATRIX_FREE_OPERATORS_1D_H
