///////////////////////////////////////////////////////////////////////////////////////////
///                         MinimalMathLibrary (MML)                                  ///
///                                                                                   ///
///  File:        FunctionSpaces/FunctionSpaceResult.h                                ///
///  Description: Diagnostic result types for function-space operations               ///
///                                                                                   ///
///  Copyright:   (c) 2024-2026 Zvonimir Vanjak                                       ///
///  License:     MIT License (see LICENSE.md)                                         ///
///                                                                                   ///
///////////////////////////////////////////////////////////////////////////////////////////
#if !defined MML_FUNCTION_SPACE_RESULT_H
#define MML_FUNCTION_SPACE_RESULT_H

#include <mml/core/FunctionSpaces/FunctionExpansion1D.h>
#include <mml/base/AlgorithmTypes.h>

#include <optional>
#include <string>
#include <utility>

namespace MML::FunctionSpaces
{
	struct FunctionSpaceOperationResult1D
	{
		std::optional<FunctionExpansion1D> expansion;
		AlgorithmStatus status = AlgorithmStatus::Success;
		std::string algorithm_name;
		std::string error_message;

		bool success() const noexcept { return status == AlgorithmStatus::Success && expansion.has_value(); }

		static FunctionSpaceOperationResult1D Success(FunctionExpansion1D expansion, const std::string& algorithmName)
		{
			FunctionSpaceOperationResult1D result;
			result.expansion = std::move(expansion);
			result.status = AlgorithmStatus::Success;
			result.algorithm_name = algorithmName;
			return result;
		}

		static FunctionSpaceOperationResult1D Failure(AlgorithmStatus status, const std::string& algorithmName, const std::string& message)
		{
			FunctionSpaceOperationResult1D result;
			result.status = status;
			result.algorithm_name = algorithmName;
			result.error_message = message;
			return result;
		}
	};

	struct DenseBVPSolveResult1D
	{
		std::optional<FunctionExpansion1D> solution;
		Real residual_norm = REAL(0.0);
		AlgorithmStatus status = AlgorithmStatus::Success;
		std::string algorithm_name = "DenseBVPSolve1D";
		std::string error_message;
		int dimension = 0;
		int boundary_conditions_applied = 0;

		bool success() const noexcept { return status == AlgorithmStatus::Success && solution.has_value(); }

		static DenseBVPSolveResult1D Success(FunctionExpansion1D solution, Real residualNorm, int boundaryConditionCount)
		{
			DenseBVPSolveResult1D result;
			result.dimension = solution.dimension();
			result.solution = std::move(solution);
			result.residual_norm = residualNorm;
			result.status = AlgorithmStatus::Success;
			result.boundary_conditions_applied = boundaryConditionCount;
			return result;
		}

		static DenseBVPSolveResult1D Failure(AlgorithmStatus status, const std::string& message)
		{
			DenseBVPSolveResult1D result;
			result.status = status;
			result.error_message = message;
			return result;
		}
	};
}

#endif // MML_FUNCTION_SPACE_RESULT_H
