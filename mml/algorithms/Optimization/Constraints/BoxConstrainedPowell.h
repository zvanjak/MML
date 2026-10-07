///////////////////////////////////////////////////////////////////////////////////////////
///                         MinimalMathLibrary (MML)                                  ///
///                                                                                   ///
///  File:        Optimization/Constraints/BoxConstrainedPowell.h                     ///
///  Description: Conservative box-constrained Powell wrapper                         ///
///                                                                                   ///
///  Copyright:   (c) 2024-2026 Zvonimir Vanjak                                       ///
///  License:     MIT License (see LICENSE.md)                                         ///
///                                                                                   ///
///////////////////////////////////////////////////////////////////////////////////////////
#if !defined MML_BOX_CONSTRAINED_POWELL_H
#define MML_BOX_CONSTRAINED_POWELL_H

#include <mml/algorithms/Optimization/Constraints/BoundConstraints.h>
#include <mml/algorithms/Optimization/Multidim/Powell.h>

#include <cmath>
#include <string>

namespace MML::Optimization
{
	namespace Detail
	{
		template<int N>
		Matrix<Real> BuildPowellDirectionMatrix(const VectorN<Real, N>& start, const BoundConstraints& bounds)
		{
			Matrix<Real> directions(N, N, REAL(0.0));
			for (int component = 0; component < N; ++component) {
				Real sign = REAL(1.0);
				if (bounds.IsUpperActive(start, component) && !bounds.IsLowerActive(start, component))
					sign = -REAL(1.0);
				directions(component, component) = sign;
			}
			return directions;
		}

		template<int N>
		class PowellProjectedScalarFunction : public IScalarFunction<N>
		{
			const IScalarFunction<N>& _func;
			const BoundConstraints& _bounds;

		public:
			PowellProjectedScalarFunction(const IScalarFunction<N>& func, const BoundConstraints& bounds)
				: _func(func), _bounds(bounds) { }

			Real operator()(const VectorN<Real, N>& point) const override
			{
				return _func(_bounds.Project(point));
			}
		};
	}

	class BoxConstrainedPowell
	{
		MultidimOptimizationConfig _config;

		void validateConfig() const
		{
			ValidateMultidimTolerance(_config.tolerance, "BoxConstrainedPowell tolerance");
			if (_config.max_iterations <= 0)
				throw MultidimOptimizationInputError("BoxConstrainedPowell max_iterations must be positive");
		}

	public:
		BoxConstrainedPowell() = default;
		explicit BoxConstrainedPowell(const MultidimOptimizationConfig& config)
			: _config(config)
		{
			validateConfig();
		}

		const MultidimOptimizationConfig& config() const noexcept { return _config; }

		template<int N>
		MultidimMinimizationResult Minimize(const IScalarFunction<N>& func,
			const VectorN<Real, N>& start,
			const BoundConstraints& bounds) const
		{
			validateConfig();
			bounds.Project(start); // dimension validation
			ValidateVectorFinite<N>(start, "BoxConstrainedPowell starting point");
			AlgorithmTimer timer;

			VectorN<Real, N> projectedStart = bounds.Project(start);
			Detail::PowellProjectedScalarFunction<N> projectedFunc(func, bounds);
			Powell optimizer(_config.tolerance, _config.max_iterations);
			Matrix<Real> directions = Detail::BuildPowellDirectionMatrix(projectedStart, bounds);
			MultidimMinimizationResult result = optimizer.Minimize(projectedFunc, projectedStart, directions);

			VectorN<Real, N> rawXmin;
			for (int component = 0; component < N; ++component)
				rawXmin[component] = result.xmin[component];
			VectorN<Real, N> projectedXmin = bounds.Project(rawXmin);

			result.xmin = Vector<Real>(N);
			for (int component = 0; component < N; ++component)
				result.xmin[component] = projectedXmin[component];
			result.fmin = func(projectedXmin);
			ValidateMultidimFunctionValue(result.fmin, "BoxConstrainedPowell final evaluation");
			result.algorithm_name = "BoxConstrainedPowell";
			result.elapsed_time_ms = timer.elapsed_ms();
			result.function_evaluations = result.iterations;
			if (!result.converged) {
				result.status = AlgorithmStatus::MaxIterationsExceeded;
				result.error_message = "Box-constrained Powell did not converge within " + std::to_string(_config.max_iterations) + " iterations";
			}
			return result;
		}
	};

	template<int N>
	MultidimMinimizationResult BoxConstrainedPowellMinimize(const IScalarFunction<N>& func,
		const VectorN<Real, N>& start,
		const BoundConstraints& bounds,
		const MultidimOptimizationConfig& config = MultidimOptimizationConfig())
	{
		BoxConstrainedPowell optimizer(config);
		return optimizer.Minimize(func, start, bounds);
	}
}

#endif // MML_BOX_CONSTRAINED_POWELL_H
