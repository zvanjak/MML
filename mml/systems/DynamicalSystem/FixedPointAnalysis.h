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
#include <limits>
#include <string>
#include <vector>

namespace MML::Systems
{
	namespace FixedPointAnalysisDetail
	{
		template<class Values>
		bool IsFinite(const Values& values) {
			return std::all_of(values.begin(), values.end(), [](Real value) { return std::isfinite(value); });
		}

		inline void Fail(FixedPoint<Real>& result, AlgorithmStatus status, const std::string& message) {
			result.status = status;
			result.error_message = message;
			result.type = FixedPointType::Unknown;
			result.isStable = false;
			result.eigenvalues.clear();
		}
	}

	template<typename Type> struct FindFixedPointResult;
	FindFixedPointResult<Real> FindFixedPointDetailed(IDynamicalSystem&, const Vector<Real>&,
		Real, int, const DynSysConfig&);

	class FixedPointFinder
	{
		friend FindFixedPointResult<Real> FindFixedPointDetailed(IDynamicalSystem&, const Vector<Real>&,
			Real, int, const DynSysConfig&);
	public:
		static FixedPoint<Real> Find(IDynamicalSystem& system, const Vector<Real>& initialGuess,
			Real tolerance = Precision::DefaultToleranceStrict, int maxIterations = Precision::FixedPointMaxIterations)
		{
			int evaluations = 0;
			return FindWithWork(system, initialGuess, tolerance, maxIterations, evaluations);
		}

	private:
		static FixedPoint<Real> FindWithWork(IDynamicalSystem& system, const Vector<Real>& initialGuess,
			Real tolerance, int maxIterations, int& evaluations)
		{
			const int dimension = system.getDim();
			FixedPoint<Real> result;
			result.location = initialGuess;
			if (!system.isAutonomous() || dimension <= 0 || initialGuess.size() != dimension
				|| !FixedPointAnalysisDetail::IsFinite(initialGuess)
				|| !std::isfinite(tolerance) || tolerance <= Real(0) || maxIterations < 0) {
				FixedPointAnalysisDetail::Fail(result, AlgorithmStatus::InvalidInput,
					"FixedPointFinder requires an autonomous positive-dimensional system, a matching finite guess, "
					"finite positive tolerance and nonnegative iteration budget");
				return result;
			}
			result.jacobian = Matrix<Real>(dimension, dimension, NoEstimate);
			Vector<Real> function(dimension), step(dimension);
			for (;;) {
				std::fill(function.begin(), function.end(), NoEstimate);
				if (result.function_evaluations == std::numeric_limits<int>::max()) {
					FixedPointAnalysisDetail::Fail(result, AlgorithmStatus::ResourceLimitExceeded, "Fixed-point evaluation counter limit reached");
					return result;
				}
				++result.function_evaluations;
				++evaluations;
				system.derivs(Real(0), result.location, function);
				if (function.size() != dimension) {
					FixedPointAnalysisDetail::Fail(result, AlgorithmStatus::InvalidInput, "Fixed-point derivative output has the wrong dimension");
					return result;
				}
				if (!FixedPointAnalysisDetail::IsFinite(function)) {
					FixedPointAnalysisDetail::Fail(result, AlgorithmStatus::NumericalInstability, "Nonfinite fixed-point derivative output");
					return result;
				}
				result.convergenceResidual = function.NormL2();
				if (!std::isfinite(result.convergenceResidual)) {
					FixedPointAnalysisDetail::Fail(result, AlgorithmStatus::NumericalInstability, "Nonfinite fixed-point residual norm");
					return result;
				}
				const bool converged = result.convergenceResidual < tolerance;
				if (!converged && result.iterations == maxIterations) {
					FixedPointAnalysisDetail::Fail(result, AlgorithmStatus::MaxIterationsExceeded,
						"Newton iteration exhausted its update budget (residual=" + std::to_string(result.convergenceResidual) + ")");
					return result;
				}
				system.jacobian(Real(0), result.location, result.jacobian);
				if (result.jacobian.rows() != dimension || result.jacobian.cols() != dimension) {
					FixedPointAnalysisDetail::Fail(result, AlgorithmStatus::InvalidInput, "Fixed-point Jacobian output has the wrong dimensions");
					return result;
				}
				if (!FixedPointAnalysisDetail::IsFinite(result.jacobian)) {
					FixedPointAnalysisDetail::Fail(result, AlgorithmStatus::NumericalInstability, "Nonfinite fixed-point Jacobian output");
					return result;
				}
				if (converged) {
					ClassifyFixedPoint(result);
					return result;
				}
				try {
					LUSolver<Real> solver(result.jacobian);
					step = solver.Solve(function * static_cast<MML::Real>(-1.0));
				}
				catch (const SingularMatrixError& error) {
					FixedPointAnalysisDetail::Fail(result, AlgorithmStatus::SingularMatrix, error.what());
					return result;
				}
				catch (const MatrixNumericalError& error) {
					FixedPointAnalysisDetail::Fail(result, AlgorithmStatus::NumericalInstability, error.what());
					return result;
				}
				const Vector<Real> nextState = result.location + step;
				if (!FixedPointAnalysisDetail::IsFinite(step) || !FixedPointAnalysisDetail::IsFinite(nextState)) {
					FixedPointAnalysisDetail::Fail(result, AlgorithmStatus::NumericalInstability, "Nonfinite Newton step or updated state");
					return result;
				}
				result.location = nextState;
				++result.iterations;
				result.convergenceResidual = NoEstimate;
				result.jacobian.Fill(NoEstimate);
			}
		}

	public:
		static std::vector<FixedPoint<Real>> FindMultiple(IDynamicalSystem& system,
			const std::vector<Vector<Real>>& initialGuesses,
			Real tolerance = Precision::FixedPointTolerance,
			Real uniqueTolerance = Precision::FixedPointUniquenessTolerance)
		{
			if (!system.isAutonomous() || system.getDim() <= 0 || !std::isfinite(tolerance) || tolerance <= Real(0)
				|| !std::isfinite(uniqueTolerance) || uniqueTolerance <= Real(0))
				throw ArgumentError("FindMultiple requires an autonomous positive-dimensional system and finite positive tolerances");
			std::vector<FixedPoint<Real>> results;
			for (const auto& guess : initialGuesses) {
				auto fixedPoint = Find(system, guess, tolerance);
				if (fixedPoint.status != AlgorithmStatus::Success || !std::isfinite(fixedPoint.convergenceResidual)
					|| fixedPoint.convergenceResidual >= tolerance || !FixedPointAnalysisDetail::IsFinite(fixedPoint.location))
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

		/// Classify the supplied Jacobian only; this does not verify a root at location.
		static void ClassifyFixedPoint(FixedPoint<Real>& fixedPoint)
		{
			fixedPoint.type = FixedPointType::Unknown;
			fixedPoint.isStable = false;
			fixedPoint.eigenvalues.clear();
			const int dimension = fixedPoint.jacobian.rows();
			if (dimension <= 0 || fixedPoint.jacobian.cols() != dimension) {
				FixedPointAnalysisDetail::Fail(fixedPoint, AlgorithmStatus::InvalidInput, "Classification requires a nonempty square Jacobian");
				return;
			}
			if (!FixedPointAnalysisDetail::IsFinite(fixedPoint.jacobian)) {
				FixedPointAnalysisDetail::Fail(fixedPoint, AlgorithmStatus::NumericalInstability, "Classification requires a finite Jacobian");
				return;
			}
			const auto eigenResult = EigenSolver::Solve(fixedPoint.jacobian);
			if (!eigenResult.converged || !IsUsable(eigenResult.status) || !std::isfinite(eigenResult.maxResidual)
				|| eigenResult.eigenvalues.size() != static_cast<std::size_t>(dimension)
				|| !std::all_of(eigenResult.eigenvalues.begin(), eigenResult.eigenvalues.end(), [](const auto& value) {
					return std::isfinite(value.real) && std::isfinite(value.imag);
				})) {
				const auto status = IsUsable(eigenResult.status) ? AlgorithmStatus::NumericalInstability : eigenResult.status;
				FixedPointAnalysisDetail::Fail(fixedPoint, status, "Fixed-point Jacobian eigensolve is unusable: "
					+ (eigenResult.error_message.empty() ? to_string(status) : eigenResult.error_message));
				return;
			}
			for (const auto& eigenvalue : eigenResult.eigenvalues)
				fixedPoint.eigenvalues.emplace_back(eigenvalue.real, eigenvalue.imag);

			int positiveReal = 0;
			int negativeReal = 0;
			bool hasComplex = false;
			bool allOscillatory = true;
			for (const auto& eigenvalue : fixedPoint.eigenvalues) {
				const Real realPart = eigenvalue.real();
				const Real imaginaryPart = eigenvalue.imag();
				if (std::abs(imaginaryPart) > Precision::DefaultToleranceStrict)
					hasComplex = true;
				else
					allOscillatory = false;
				if (realPart > Precision::DefaultToleranceStrict)
					++positiveReal;
				else if (realPart < -Precision::DefaultToleranceStrict)
					++negativeReal;
			}

			fixedPoint.status = AlgorithmStatus::Success;
			fixedPoint.error_message.clear();
			const int neutralReal = dimension - positiveReal - negativeReal;
			if (neutralReal > 0) {
				// A linear center label is not a nonlinear stability or periodic-orbit proof.
				if (neutralReal == dimension && allOscillatory)
					fixedPoint.type = FixedPointType::Center;
				return;
			}
			if (positiveReal > 0 && negativeReal > 0)
				fixedPoint.type = hasComplex ? FixedPointType::SaddleFocus : FixedPointType::Saddle;
			else if (negativeReal == dimension) {
				fixedPoint.type = hasComplex ? FixedPointType::StableFocus : FixedPointType::StableNode;
				fixedPoint.isStable = true;
			}
			else if (positiveReal == dimension)
				fixedPoint.type = hasComplex ? FixedPointType::UnstableFocus : FixedPointType::UnstableNode;
		}
	};

	template<typename Type = Real>
	struct FindFixedPointResult : public EvaluationResultBase
	{
		FixedPoint<Type> fixed_point;
	};

	inline FindFixedPointResult<Real> FindFixedPointDetailed(IDynamicalSystem& system,
		const Vector<Real>& initialGuess, Real tolerance = Precision::DefaultToleranceStrict,
		int maxIterations = Precision::FixedPointMaxIterations, const DynSysConfig& config = {})
	{
		return DynSysDetail::ExecuteDynSysDetailed<FindFixedPointResult<Real>>(
			"FixedPointFinder", config, [&](FindFixedPointResult<Real>& result) {
				result.fixed_point = FixedPointFinder::FindWithWork(system, initialGuess,
					tolerance, maxIterations, result.function_evaluations);
				result.status = result.fixed_point.status;
				result.error_message = result.fixed_point.error_message;
			});
	}
}

#endif // MML_FIXED_POINT_ANALYSIS_H
