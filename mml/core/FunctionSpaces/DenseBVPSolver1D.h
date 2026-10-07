///////////////////////////////////////////////////////////////////////////////////////////
///                         MinimalMathLibrary (MML)                                  ///
///                                                                                   ///
///  File:        FunctionSpaces/DenseBVPSolver1D.h                                   ///
///  Description: Dense collocation BVP solver helpers in one dimension               ///
///                                                                                   ///
///  Copyright:   (c) 2024-2026 Zvonimir Vanjak                                       ///
///  License:     MIT License (see LICENSE.md)                                         ///
///                                                                                   ///
///////////////////////////////////////////////////////////////////////////////////////////
#if !defined MML_DENSE_BVP_SOLVER_1D_H
#define MML_DENSE_BVP_SOLVER_1D_H

#include <mml/core/FunctionSpaces/OperatorAssembly.h>
#include <mml/core/FunctionSpaces/FunctionSpaceResult.h>

#include <mml/core/LinAlgEqSolvers/LinAlgDirect.h>

namespace MML::FunctionSpaces
{
	inline DenseBVPSolveResult1D SolveDenseCollocationBVP(const LinearDifferentialOperator1D& differentialOperator,
		const IRealFunction& rhsFunction,
		const TrialSpace1D& trialSpace,
		const BoundaryConditions1D& boundaryConditions)
	{
		try {
			Matrix<Real> systemMatrix = BuildCollocationMatrix(differentialOperator, trialSpace);
			Vector<Real> rhs(trialSpace.dimension());
			for (int index = 0; index < trialSpace.dimension(); ++index)
				rhs[index] = rhsFunction(trialSpace.node(index));

			boundaryConditions.validate(trialSpace.functionSpace());
			for (int index = 0; index < boundaryConditions.size(); ++index)
				ApplyBoundaryRow(systemMatrix, rhs, index, trialSpace, boundaryConditions[index]);

			LUSolver<Real> solver(systemMatrix);
			Vector<Real> coefficients = solver.Solve(rhs);
			FunctionExpansion1D solution(trialSpace, coefficients);
			return DenseBVPSolveResult1D::Success(solution, REAL(0.0), boundaryConditions.size());
		}
		catch (const SingularMatrixError& e) {
			return DenseBVPSolveResult1D::Failure(AlgorithmStatus::SingularMatrix, e.what());
		}
		catch (const FunctionSpaceInputError& e) {
			return DenseBVPSolveResult1D::Failure(AlgorithmStatus::InvalidInput, e.what());
		}
		catch (const std::exception& e) {
			return DenseBVPSolveResult1D::Failure(AlgorithmStatus::AlgorithmSpecificFailure, e.what());
		}
	}
}

#endif // MML_DENSE_BVP_SOLVER_1D_H
