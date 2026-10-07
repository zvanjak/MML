///////////////////////////////////////////////////////////////////////////////////////////
///                         MinimalMathLibrary (MML)                                  ///
///                                                                                   ///
///  File:        FunctionSpaces/Interpolation.h                                      ///
///  Description: Interpolation maps into finite function-space coordinates           ///
///                                                                                   ///
///  Copyright:   (c) 2024-2026 Zvonimir Vanjak                                       ///
///  License:     MIT License (see LICENSE.md)                                         ///
///                                                                                   ///
///////////////////////////////////////////////////////////////////////////////////////////
#if !defined MML_FUNCTION_SPACES_INTERPOLATION_H
#define MML_FUNCTION_SPACES_INTERPOLATION_H

#include <mml/core/FunctionSpaces/FunctionExpansion1D.h>
#include <mml/core/FunctionSpaces/FunctionSpaceResult.h>

#include <mml/base/Matrix/Matrix.h>
#include <mml/core/LinAlgEqSolvers/LinAlgDirect.h>

namespace MML::FunctionSpaces
{
	inline FunctionExpansion1D Interpolate(const IRealFunction& function, const TrialSpace1D& trialSpace)
	{
		if (!trialSpace.hasNodes())
			throw FunctionSpaceInputError("Interpolate: trial space does not provide interpolation nodes");
		if (trialSpace.nodeCount() != trialSpace.dimension())
			throw FunctionSpaceInputError("Interpolate: node count must match trial-space dimension");

		int dimension = trialSpace.dimension();
		Matrix<Real> basisValues(dimension, dimension);
		Vector<Real> sampledValues(dimension);

		for (int row = 0; row < dimension; ++row) {
			Real x = trialSpace.node(row);
			sampledValues[row] = function(x);
			for (int column = 0; column < dimension; ++column)
				basisValues(row, column) = trialSpace.basisValue(column, x);
		}

		LUSolver<Real> solver(basisValues);
		Vector<Real> coefficients = solver.Solve(sampledValues);
		return FunctionExpansion1D(trialSpace, coefficients);
	}

	inline FunctionSpaceOperationResult1D TryInterpolate(const IRealFunction& function, const TrialSpace1D& trialSpace)
	{
		try {
			return FunctionSpaceOperationResult1D::Success(Interpolate(function, trialSpace), "Interpolate");
		}
		catch (const FunctionSpaceInputError& e) {
			return FunctionSpaceOperationResult1D::Failure(AlgorithmStatus::InvalidInput, "Interpolate", e.what());
		}
		catch (const SingularMatrixError& e) {
			return FunctionSpaceOperationResult1D::Failure(AlgorithmStatus::SingularMatrix, "Interpolate", e.what());
		}
		catch (const std::exception& e) {
			return FunctionSpaceOperationResult1D::Failure(AlgorithmStatus::AlgorithmSpecificFailure, "Interpolate", e.what());
		}
	}
}

#endif // MML_FUNCTION_SPACES_INTERPOLATION_H
