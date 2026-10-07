///////////////////////////////////////////////////////////////////////////////////////////
///                         MinimalMathLibrary (MML)                                  ///
///                                                                                   ///
///  File:        FixedPointAnalysis.h                                               ///
///  Description: Fixed-point finding and stability classification                   ///
///                                                                                   ///
///  Copyright:   (c) 2024-2026 Zvonimir Vanjak                                       ///
///  License:     MIT License (see LICENSE.md)                                         ///
///                                                                                   ///
///////////////////////////////////////////////////////////////////////////////////////////
#if !defined MML_FIXED_POINT_ANALYSIS_H
#define MML_FIXED_POINT_ANALYSIS_H

#include <mml/systems/DynamicalSystem/DynamicalAnalysisCommon.h>
#include <mml/systems/DynamicalSystem/DynamicalSystemTypes.h>
#include <mml/algorithms/EigenSystemSolvers.h>
#include <mml/core/LinAlgEqSolvers.h>

#include <algorithm>
#include <cmath>
#include <complex>
#include <string>
#include <vector>

namespace MML::Systems
{
	class FixedPointFinder
	{
	public:
		static FixedPoint<Real> Find(IDynamicalSystem& system, const Vector<Real>& initialGuess,
			Real tolerance = Precision::DefaultToleranceStrict, int maxIterations = 50)
		{
			const int dimension = system.getDim();
			FixedPoint<Real> result;
			result.location = initialGuess;
			result.jacobian.Resize(dimension, dimension);

			Vector<Real> state = initialGuess;
			Vector<Real> function(dimension), step(dimension);
			Matrix<Real> jacobian(dimension, dimension);

			for (int iteration = 0; iteration < maxIterations; ++iteration) {
				system.derivs(REAL(0.0), state, function);
				const Real residual = function.NormL2();
				if (residual < tolerance) {
					result.location = state;
					result.convergenceResidual = residual;
					result.iterations = iteration;
					system.jacobian(REAL(0.0), state, result.jacobian);
					ClassifyFixedPoint(result);
					return result;
				}

				system.jacobian(REAL(0.0), state, jacobian);
				try {
					LUSolver<Real> solver(jacobian);
					step = solver.Solve(function * REAL(-1.0));
				}
				catch (...) {
					result.type = FixedPointType::Unknown;
					result.convergenceResidual = residual;
					result.iterations = iteration;
					return result;
				}
				state = state + step;
			}

			system.derivs(REAL(0.0), state, function);
			result.location = state;
			result.convergenceResidual = function.NormL2();
			result.iterations = maxIterations;
			result.type = FixedPointType::Unknown;
			return result;
		}

		static std::vector<FixedPoint<Real>> FindMultiple(IDynamicalSystem& system,
			const std::vector<Vector<Real>>& initialGuesses,
			Real tolerance = Precision::FixedPointTolerance,
			Real uniqueTolerance = Precision::FixedPointUniquenessTolerance)
		{
			std::vector<FixedPoint<Real>> results;
			for (const auto& guess : initialGuesses) {
				auto fixedPoint = Find(system, guess, tolerance);
				if (fixedPoint.convergenceResidual >= tolerance)
					continue;

				bool isNew = true;
				for (const auto& existing : results) {
					if ((fixedPoint.location - existing.location).NormL2() < uniqueTolerance) {
						isNew = false;
						break;
					}
				}
				if (isNew)
					results.push_back(fixedPoint);
			}
			return results;
		}

		static void ClassifyFixedPoint(FixedPoint<Real>& fixedPoint)
		{
			const int dimension = fixedPoint.jacobian.rows();
			const auto eigenResult = EigenSolver::Solve(fixedPoint.jacobian);
			fixedPoint.eigenvalues.clear();
			for (const auto& eigenvalue : eigenResult.eigenvalues)
				fixedPoint.eigenvalues.emplace_back(eigenvalue.real, eigenvalue.imag);

			int positiveReal = 0;
			int negativeReal = 0;
			bool hasComplex = false;
			Real maxRealPart = -REAL(1e30);
			Real minRealPart = REAL(1e30);
			for (const auto& eigenvalue : fixedPoint.eigenvalues) {
				const Real realPart = eigenvalue.real();
				const Real imaginaryPart = eigenvalue.imag();
				maxRealPart = std::max(maxRealPart, realPart);
				minRealPart = std::min(minRealPart, realPart);
				if (std::abs(imaginaryPart) > Precision::DefaultToleranceStrict)
					hasComplex = true;
				if (realPart > Precision::DefaultToleranceStrict)
					++positiveReal;
				else if (realPart < -Precision::DefaultToleranceStrict)
					++negativeReal;
			}

			fixedPoint.isStable = maxRealPart < -Precision::DefaultToleranceStrict;
			if (hasComplex) {
				if (std::abs(maxRealPart) < Precision::DefaultToleranceStrict
					&& std::abs(minRealPart) < Precision::DefaultToleranceStrict)
					fixedPoint.type = FixedPointType::Center;
				else if (positiveReal > 0 && negativeReal > 0)
					fixedPoint.type = FixedPointType::SaddleFocus;
				else if (maxRealPart < -Precision::DefaultToleranceStrict)
					fixedPoint.type = FixedPointType::StableFocus;
				else
					fixedPoint.type = FixedPointType::UnstableFocus;
			}
			else if (positiveReal > 0 && negativeReal > 0)
				fixedPoint.type = FixedPointType::Saddle;
			else if (positiveReal == 0 && negativeReal == dimension)
				fixedPoint.type = FixedPointType::StableNode;
			else if (positiveReal == dimension && negativeReal == 0)
				fixedPoint.type = FixedPointType::UnstableNode;
			else
				fixedPoint.type = FixedPointType::Unknown;
		}
	};

	template<typename Type = Real>
	struct FindFixedPointResult : public EvaluationResultBase
	{
		FixedPoint<Type> fixed_point;
	};

	inline FindFixedPointResult<Real> FindFixedPointDetailed(IDynamicalSystem& system,
		const Vector<Real>& initialGuess, Real tolerance = Precision::DefaultToleranceStrict,
		int maxIterations = 50, const DynSysConfig& config = {})
	{
		return DynSysDetail::ExecuteDynSysDetailed<FindFixedPointResult<Real>>(
			"FixedPointFinder", config, [&](FindFixedPointResult<Real>& result) {
				result.fixed_point = FixedPointFinder::Find(system, initialGuess, tolerance, maxIterations);
				result.function_evaluations = result.fixed_point.iterations;
				if (result.fixed_point.convergenceResidual >= tolerance) {
					result.status = AlgorithmStatus::MaxIterationsExceeded;
					result.error_message = "Newton iteration did not converge (residual="
						+ std::to_string(result.fixed_point.convergenceResidual) + ")";
				}
			});
	}
}

#endif // MML_FIXED_POINT_ANALYSIS_H