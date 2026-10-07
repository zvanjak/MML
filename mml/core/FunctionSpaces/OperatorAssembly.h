///////////////////////////////////////////////////////////////////////////////////////////
///                         MinimalMathLibrary (MML)                                  ///
///                                                                                   ///
///  File:        FunctionSpaces/OperatorAssembly.h                                   ///
///  Description: Dense operator assembly helpers for finite function spaces          ///
///                                                                                   ///
///  Copyright:   (c) 2024-2026 Zvonimir Vanjak                                       ///
///  License:     MIT License (see LICENSE.md)                                         ///
///                                                                                   ///
///////////////////////////////////////////////////////////////////////////////////////////
#if !defined MML_FUNCTION_SPACES_OPERATOR_ASSEMBLY_H
#define MML_FUNCTION_SPACES_OPERATOR_ASSEMBLY_H

#include <mml/core/FunctionSpaces/BoundaryCondition1D.h>
#include <mml/core/FunctionSpaces/ChebyshevCollocationSpace1D.h>
#include <mml/core/FunctionSpaces/LinearDifferentialOperator1D.h>

#include <mml/base/Matrix/Matrix.h>
#include <mml/base/Vector/Vector.h>
#include <mml/core/Integration/Integration1D.h>

namespace MML::FunctionSpaces
{
	inline Matrix<Real> BuildCollocationMatrix(const LinearDifferentialOperator1D& differentialOperator,
		const TrialSpace1D& trialSpace)
	{
		if (!trialSpace.hasNodes())
			throw FunctionSpaceInputError("BuildCollocationMatrix: trial space does not provide collocation nodes");
		if (trialSpace.nodeCount() != trialSpace.dimension())
			throw FunctionSpaceInputError("BuildCollocationMatrix: node count must match trial-space dimension");

		if (const auto* chebyshevSpace = dynamic_cast<const ChebyshevCollocationSpace1D*>(&trialSpace)) {
			if (differentialOperator.order() > 2)
				throw FunctionSpaceError("BuildCollocationMatrix: Chebyshev collocation currently supports operators up to order 2");

			Matrix<Real> firstDerivative = chebyshevSpace->firstDerivativeMatrix();
			Matrix<Real> secondDerivative = chebyshevSpace->secondDerivativeMatrix();
			Matrix<Real> basisValues(trialSpace.dimension(), trialSpace.dimension());
			Matrix<Real> firstBasisDerivatives(trialSpace.dimension(), trialSpace.dimension());
			Matrix<Real> secondBasisDerivatives(trialSpace.dimension(), trialSpace.dimension());

			for (int row = 0; row < trialSpace.dimension(); ++row)
				for (int column = 0; column < trialSpace.dimension(); ++column)
					basisValues(row, column) = trialSpace.basisValue(column, trialSpace.node(row));

			for (int row = 0; row < trialSpace.dimension(); ++row) {
				for (int column = 0; column < trialSpace.dimension(); ++column) {
					Real firstValue = REAL(0.0);
					Real secondValue = REAL(0.0);
					for (int inner = 0; inner < trialSpace.dimension(); ++inner) {
						firstValue += firstDerivative(row, inner) * basisValues(inner, column);
						secondValue += secondDerivative(row, inner) * basisValues(inner, column);
					}
					firstBasisDerivatives(row, column) = firstValue;
					secondBasisDerivatives(row, column) = secondValue;
				}
			}

			Matrix<Real> result(trialSpace.dimension(), trialSpace.dimension());
			for (int row = 0; row < trialSpace.dimension(); ++row) {
				Real x = trialSpace.node(row);
				for (int column = 0; column < trialSpace.dimension(); ++column) {
					Real value = differentialOperator.coefficient(0, x) * basisValues(row, column);
					if (differentialOperator.order() >= 1)
						value += differentialOperator.coefficient(1, x) * firstBasisDerivatives(row, column);
					if (differentialOperator.order() >= 2)
						value += differentialOperator.coefficient(2, x) * secondBasisDerivatives(row, column);
					result(row, column) = value;
				}
			}
			return result;
		}

		Matrix<Real> result(trialSpace.dimension(), trialSpace.dimension());
		for (int row = 0; row < trialSpace.dimension(); ++row) {
			Real x = trialSpace.node(row);
			for (int column = 0; column < trialSpace.dimension(); ++column)
				result(row, column) = differentialOperator.applyToBasis(trialSpace, column, x);
		}
		return result;
	}

	inline void ApplyBoundaryRow(Matrix<Real>& matrix,
		Vector<Real>& rhs,
		int rowIndex,
		const TrialSpace1D& trialSpace,
		const BoundaryCondition1D& condition)
	{
		if (matrix.rows() != trialSpace.dimension() || matrix.cols() != trialSpace.dimension())
			throw MatrixDimensionError("ApplyBoundaryRow: matrix dimensions must match trial-space dimension", matrix.rows(), matrix.cols(), trialSpace.dimension(), trialSpace.dimension());
		if (rhs.size() != trialSpace.dimension())
			throw VectorDimensionError("ApplyBoundaryRow: RHS dimension must match trial-space dimension", rhs.size(), trialSpace.dimension());
		if (rowIndex < 0 || rowIndex >= trialSpace.dimension())
			throw FunctionSpaceInputError("ApplyBoundaryRow: row index out of range");

		Vector<Real> boundaryRow = condition.row(trialSpace);
		for (int column = 0; column < trialSpace.dimension(); ++column)
			matrix(rowIndex, column) = boundaryRow[column];
		rhs[rowIndex] = condition.value;
	}

	inline Matrix<Real> BuildGalerkinStiffnessMatrix(const TrialSpace1D& trialSpace,
		const IRealFunction& p,
		const IRealFunction& q,
		Real eps = Precision::DefaultToleranceStrict)
	{
		if (!trialSpace.hasBasisDerivative(1))
			throw FunctionSpaceError("BuildGalerkinStiffnessMatrix: trial space must provide first derivatives");

		class Integrand : public IRealFunction
		{
			const TrialSpace1D& _trialSpace;
			const IRealFunction& _p;
			const IRealFunction& _q;
			int _testIndex;
			int _trialIndex;

		public:
			Integrand(const TrialSpace1D& trialSpace, const IRealFunction& p, const IRealFunction& q, int testIndex, int trialIndex)
				: _trialSpace(trialSpace), _p(p), _q(q), _testIndex(testIndex), _trialIndex(trialIndex) { }

			Real operator()(Real x) const override
			{
				Real testValue = _trialSpace.basisValue(_testIndex, x);
				Real trialValue = _trialSpace.basisValue(_trialIndex, x);
				Real testDerivative = _trialSpace.basisDerivative(_testIndex, 1, x);
				Real trialDerivative = _trialSpace.basisDerivative(_trialIndex, 1, x);
				return _p(x) * testDerivative * trialDerivative + _q(x) * testValue * trialValue;
			}
		};

		Matrix<Real> result(trialSpace.dimension(), trialSpace.dimension());
		for (int row = 0; row < trialSpace.dimension(); ++row) {
			for (int column = 0; column < trialSpace.dimension(); ++column) {
				Integrand integrand(trialSpace, p, q, row, column);
				result(row, column) = IntegrateTrap(integrand,
					trialSpace.functionSpace().domainMin(),
					trialSpace.functionSpace().domainMax(),
					eps).value;
			}
		}
		return result;
	}
}

#endif // MML_FUNCTION_SPACES_OPERATOR_ASSEMBLY_H
