///////////////////////////////////////////////////////////////////////////////////////////
///                         MinimalMathLibrary (MML)                                  ///
///                                                                                   ///
///  File:        FunctionSpaces/Projection.h                                         ///
///  Description: Projection maps into finite function-space coordinates              ///
///                                                                                   ///
///  Copyright:   (c) 2024-2026 Zvonimir Vanjak                                       ///
///  License:     MIT License (see LICENSE.md)                                         ///
///                                                                                   ///
///////////////////////////////////////////////////////////////////////////////////////////
#if !defined MML_FUNCTION_SPACES_PROJECTION_H
#define MML_FUNCTION_SPACES_PROJECTION_H

#include <mml/core/FunctionSpaces/FunctionExpansion1D.h>
#include <mml/core/FunctionSpaces/FunctionSpaceResult.h>
#include <mml/core/FunctionSpaces/OrthogonalBasisTrialSpace1D.h>

namespace MML::FunctionSpaces
{
	inline FunctionExpansion1D ProjectL2(const IRealFunction& function,
		const OrthogonalBasisTrialSpace1D& trialSpace,
		Real eps = Precision::DefaultToleranceStrict)
	{
		Vector<Real> coefficients(trialSpace.dimension());
		const OrthogonalBasis& basis = trialSpace.basis();
		for (int index = 0; index < trialSpace.dimension(); ++index)
			coefficients[index] = basis.ComputeCoefficient(function, index, eps);
		return FunctionExpansion1D(trialSpace, coefficients);
	}

	inline FunctionSpaceOperationResult1D TryProjectL2(const IRealFunction& function,
		const OrthogonalBasisTrialSpace1D& trialSpace,
		Real eps = Precision::DefaultToleranceStrict)
	{
		try {
			return FunctionSpaceOperationResult1D::Success(ProjectL2(function, trialSpace, eps), "ProjectL2");
		}
		catch (const FunctionSpaceInputError& e) {
			return FunctionSpaceOperationResult1D::Failure(AlgorithmStatus::InvalidInput, "ProjectL2", e.what());
		}
		catch (const std::exception& e) {
			return FunctionSpaceOperationResult1D::Failure(AlgorithmStatus::AlgorithmSpecificFailure, "ProjectL2", e.what());
		}
	}
}

#endif // MML_FUNCTION_SPACES_PROJECTION_H
