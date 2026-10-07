///////////////////////////////////////////////////////////////////////////////////////////
///                         MinimalMathLibrary (MML)                                  ///
///                                                                                   ///
///  File:        Optimization/Constraints/BoxConstrainedNelderMead.h                 ///
///  Description: Conservative box-constrained Nelder-Mead wrapper                    ///
///                                                                                   ///
///  Copyright:   (c) 2024-2026 Zvonimir Vanjak                                       ///
///  License:     MIT License (see LICENSE.md)                                         ///
///                                                                                   ///
///////////////////////////////////////////////////////////////////////////////////////////
#if !defined MML_BOX_CONSTRAINED_NELDER_MEAD_H
#define MML_BOX_CONSTRAINED_NELDER_MEAD_H

#include <mml/algorithms/Optimization/Constraints/BoundConstraints.h>
#include <mml/algorithms/Optimization/Multidim/NelderMead.h>

#include <cmath>
#include <string>

namespace MML::Optimization
{
	namespace Detail
	{
		template<int N>
		class ProjectedScalarFunction : public IScalarFunction<N>
		{
			const IScalarFunction<N>& _func;
			const BoundConstraints& _bounds;

		public:
			ProjectedScalarFunction(const IScalarFunction<N>& func, const BoundConstraints& bounds)
				: _func(func), _bounds(bounds) { }

			Real operator()(const VectorN<Real, N>& point) const override
			{
				return _func(_bounds.Project(point));
			}
		};

		template<int N>
		Matrix<Real> BuildBoundedSimplex(const VectorN<Real, N>& start, const BoundConstraints& bounds, const Vector<Real>& deltas)
		{
			if (deltas.size() != N)
				throw MultidimOptimizationError("BoxConstrainedNelderMead deltas vector dimension mismatch");

			ValidateVectorFinite<N>(start, "BoxConstrainedNelderMead starting point");
			VectorN<Real, N> base = bounds.Project(start);
			Matrix<Real> simplex(N + 1, N);

			for (int component = 0; component < N; ++component)
				simplex(0, component) = base[component];

			for (int vertexIndex = 1; vertexIndex <= N; ++vertexIndex) {
				int component = vertexIndex - 1;
				Real delta = deltas[component];
				if (!std::isfinite(delta) || delta == REAL(0.0))
					throw MultidimOptimizationInputError("BoxConstrainedNelderMead: deltas[" + std::to_string(component) + "] must be finite and non-zero");

				VectorN<Real, N> vertex = base;
				vertex[component] += delta;
				VectorN<Real, N> projected = bounds.Project(vertex);

				if (projected[component] == base[component] && !bounds.IsFixed(component)) {
					vertex = base;
					vertex[component] -= delta;
					projected = bounds.Project(vertex);
				}

				for (int column = 0; column < N; ++column)
					simplex(vertexIndex, column) = projected[column];
			}

			return simplex;
		}
	}

	class BoxConstrainedNelderMead
	{
		MultidimOptimizationConfig _config;

		void validateConfig() const
		{
			ValidateMultidimTolerance(_config.tolerance, "BoxConstrainedNelderMead tolerance");
			if (_config.max_iterations <= 0)
				throw MultidimOptimizationInputError("BoxConstrainedNelderMead max_iterations must be positive");
			if (!std::isfinite(_config.initial_delta) || _config.initial_delta == REAL(0.0))
				throw MultidimOptimizationInputError("BoxConstrainedNelderMead initial_delta must be finite and non-zero");
		}

	public:
		BoxConstrainedNelderMead() = default;
		explicit BoxConstrainedNelderMead(const MultidimOptimizationConfig& config)
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
			Vector<Real> deltas(N, _config.initial_delta);
			return Minimize(func, start, bounds, deltas);
		}

		template<int N>
		MultidimMinimizationResult Minimize(const IScalarFunction<N>& func,
			const VectorN<Real, N>& start,
			const BoundConstraints& bounds,
			const Vector<Real>& deltas) const
		{
			validateConfig();
			bounds.Project(start); // dimension validation
			AlgorithmTimer timer;

			Matrix<Real> simplex = Detail::BuildBoundedSimplex(start, bounds, deltas);
			Detail::ProjectedScalarFunction<N> projectedFunc(func, bounds);
			NelderMead optimizer(_config.tolerance, _config.max_iterations);
			MultidimMinimizationResult result = optimizer.Minimize(projectedFunc, simplex);

			VectorN<Real, N> rawXmin;
			for (int component = 0; component < N; ++component)
				rawXmin[component] = result.xmin[component];
			VectorN<Real, N> projectedXmin = bounds.Project(rawXmin);

			result.xmin = Vector<Real>(N);
			for (int component = 0; component < N; ++component)
				result.xmin[component] = projectedXmin[component];
			result.fmin = func(projectedXmin);
			ValidateMultidimFunctionValue(result.fmin, "BoxConstrainedNelderMead final evaluation");
			result.algorithm_name = "BoxConstrainedNelderMead";
			result.elapsed_time_ms = timer.elapsed_ms();
			result.function_evaluations = optimizer.getNumFuncEvals() + 1;
			if (!result.converged) {
				result.status = AlgorithmStatus::MaxIterationsExceeded;
				result.error_message = "Box-constrained Nelder-Mead did not converge within " + std::to_string(_config.max_iterations) + " iterations";
			}
			return result;
		}
	};

	template<int N>
	MultidimMinimizationResult BoxConstrainedNelderMeadMinimize(const IScalarFunction<N>& func,
		const VectorN<Real, N>& start,
		const BoundConstraints& bounds,
		const MultidimOptimizationConfig& config = MultidimOptimizationConfig())
	{
		BoxConstrainedNelderMead optimizer(config);
		return optimizer.Minimize(func, start, bounds);
	}

	template<int N>
	MultidimMinimizationResult BoxConstrainedNelderMeadMinimize(const IScalarFunction<N>& func,
		const VectorN<Real, N>& start,
		const BoundConstraints& bounds,
		const Vector<Real>& deltas,
		const MultidimOptimizationConfig& config = MultidimOptimizationConfig())
	{
		BoxConstrainedNelderMead optimizer(config);
		return optimizer.Minimize(func, start, bounds, deltas);
	}
}

#endif // MML_BOX_CONSTRAINED_NELDER_MEAD_H
