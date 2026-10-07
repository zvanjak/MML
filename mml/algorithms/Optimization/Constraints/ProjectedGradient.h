///////////////////////////////////////////////////////////////////////////////////////////
///                         MinimalMathLibrary (MML)                                  ///
///                                                                                   ///
///  File:        Optimization/Constraints/ProjectedGradient.h                        ///
///  Description: Projected-gradient minimization for bound constraints               ///
///                                                                                   ///
///  Copyright:   (c) 2024-2026 Zvonimir Vanjak                                       ///
///  License:     MIT License (see LICENSE.md)                                         ///
///                                                                                   ///
///////////////////////////////////////////////////////////////////////////////////////////
#if !defined MML_PROJECTED_GRADIENT_H
#define MML_PROJECTED_GRADIENT_H

#include <mml/algorithms/Optimization/Constraints/BoundConstraints.h>
#include <mml/algorithms/Optimization/Multidim/LineSearch.h>

#include <chrono>
#include <cmath>
#include <string>

namespace MML::Optimization
{
	struct ProjectedGradientConfig
	{
		Real gradient_tolerance = PrecisionValues<Real>::OptimizationGradientTolerance;
		Real step_tolerance = PrecisionValues<Real>::OptimizationTolerance;
		Real initial_step = REAL(1.0);
		Real min_step = REAL(1e-14);
		Real contraction = REAL(0.5);
		Real armijo = REAL(1e-4);
		int max_iterations = 1000;
		bool verbose = false;
	};

	template<int N>
	struct ProjectedGradientResult
	{
		VectorN<Real, N> xmin;
		Real fmin = REAL(0.0);
		Real projected_gradient_norm = REAL(0.0);
		Real step_norm = REAL(0.0);
		int iterations = 0;
		int function_evaluations = 0;
		int gradient_evaluations = 0;
		bool converged = false;
		AlgorithmStatus status = AlgorithmStatus::Success;
		std::string error_message;
		double elapsed_time_ms = 0.0;
	};

	namespace Detail
	{
		template<int N>
		Real Dot(const VectorN<Real, N>& left, const VectorN<Real, N>& right)
		{
			Real result = REAL(0.0);
			for (int index = 0; index < N; ++index)
				result += left[index] * right[index];
			return result;
		}

		template<int N>
		VectorN<Real, N> Difference(const VectorN<Real, N>& left, const VectorN<Real, N>& right)
		{
			VectorN<Real, N> result;
			for (int index = 0; index < N; ++index)
				result[index] = left[index] - right[index];
			return result;
		}

		template<int N>
		VectorN<Real, N> GradientStepPoint(const VectorN<Real, N>& point, const VectorN<Real, N>& gradient, Real step)
		{
			VectorN<Real, N> result;
			for (int index = 0; index < N; ++index)
				result[index] = point[index] - step * gradient[index];
			return result;
		}
	}

	class ProjectedGradient
	{
		ProjectedGradientConfig _config;

		void validateConfig() const
		{
			ValidateMultidimTolerance(_config.gradient_tolerance, "ProjectedGradient gradient_tolerance");
			ValidateMultidimTolerance(_config.step_tolerance, "ProjectedGradient step_tolerance");
			ValidateMultidimTolerance(_config.initial_step, "ProjectedGradient initial_step");
			ValidateMultidimTolerance(_config.min_step, "ProjectedGradient min_step");
			if (_config.contraction <= REAL(0.0) || _config.contraction >= REAL(1.0))
				throw MultidimOptimizationInputError("ProjectedGradient contraction must be in (0, 1)");
			if (_config.armijo <= REAL(0.0) || _config.armijo >= REAL(1.0))
				throw MultidimOptimizationInputError("ProjectedGradient armijo must be in (0, 1)");
			if (_config.max_iterations <= 0)
				throw MultidimOptimizationInputError("ProjectedGradient max_iterations must be positive");
		}

	public:
		ProjectedGradient() = default;
		explicit ProjectedGradient(const ProjectedGradientConfig& config)
			: _config(config)
		{
			validateConfig();
		}

		const ProjectedGradientConfig& config() const noexcept { return _config; }

		template<int N>
		ProjectedGradientResult<N> Minimize(const IDifferentiableScalarFunction<N>& func,
			const VectorN<Real, N>& start,
			const BoundConstraints& bounds) const
		{
			validateConfig();
			bounds.Project(start); // dimension validation
			ValidateVectorFinite<N>(start, "ProjectedGradient starting point");

			auto startTime = std::chrono::steady_clock::now();
			ProjectedGradientResult<N> result;
			VectorN<Real, N> x = bounds.Project(start);
			VectorN<Real, N> gradient;
			Real f = func(x);
			ValidateMultidimFunctionValue(f, "ProjectedGradient initial evaluation");
			result.function_evaluations = 1;

			for (int iteration = 0; iteration < _config.max_iterations; ++iteration) {
				func.Gradient(x, gradient);
				result.gradient_evaluations++;
				ValidateVectorFinite<N>(gradient, "ProjectedGradient gradient");

				VectorN<Real, N> unitStep = bounds.Project(Detail::GradientStepPoint(x, gradient, REAL(1.0)));
				VectorN<Real, N> projectedGradient = Detail::Difference(x, unitStep);
				result.projected_gradient_norm = projectedGradient.NormL2();
				if (result.projected_gradient_norm <= _config.gradient_tolerance) {
					result.converged = true;
					result.status = AlgorithmStatus::Success;
					result.iterations = iteration;
					break;
				}

				Real step = _config.initial_step;
				bool accepted = false;
				while (step >= _config.min_step) {
					VectorN<Real, N> trial = bounds.Project(Detail::GradientStepPoint(x, gradient, step));
					VectorN<Real, N> direction = Detail::Difference(trial, x);
					Real stepNorm = direction.NormL2();
					if (stepNorm <= _config.step_tolerance) {
						result.step_norm = stepNorm;
						break;
					}

					Real directionalDerivative = Detail::Dot(gradient, direction);
					Real trialValue = func(trial);
					result.function_evaluations++;
					if (!std::isfinite(trialValue)) {
						step *= _config.contraction;
						continue;
					}

					if (directionalDerivative >= REAL(0.0)) {
						step *= _config.contraction;
						continue;
					}

					if (trialValue <= f + _config.armijo * directionalDerivative) {
						x = trial;
						f = trialValue;
						result.step_norm = stepNorm;
						accepted = true;
						break;
					}

					step *= _config.contraction;
				}

				result.iterations = iteration + 1;
				if (!accepted) {
					result.status = AlgorithmStatus::Stalled;
					result.error_message = "ProjectedGradient line search stalled";
					break;
				}
			}

			if (!result.converged && result.status == AlgorithmStatus::Success) {
				result.status = AlgorithmStatus::MaxIterationsExceeded;
				result.error_message = "maximum iterations reached";
			}

			result.xmin = x;
			result.fmin = f;
			auto endTime = std::chrono::steady_clock::now();
			result.elapsed_time_ms = std::chrono::duration<double, std::milli>(endTime - startTime).count();
			return result;
		}
	};

	template<int N>
	ProjectedGradientResult<N> ProjectedGradientMinimize(const IDifferentiableScalarFunction<N>& func,
		const VectorN<Real, N>& start,
		const BoundConstraints& bounds,
		const ProjectedGradientConfig& config = ProjectedGradientConfig())
	{
		ProjectedGradient optimizer(config);
		return optimizer.Minimize(func, start, bounds);
	}
}

#endif // MML_PROJECTED_GRADIENT_H
