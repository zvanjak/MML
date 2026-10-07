///////////////////////////////////////////////////////////////////////////////////////////
///                         MinimalMathLibrary (MML)                                  ///
///                                                                                   ///
///  File:        NonlinearSystemSolvers.h                                            ///
///  Description: Dense Newton solvers for nonlinear systems F(x) = 0                 ///
///                                                                                   ///
///  Copyright:   (c) 2024-2026 Zvonimir Vanjak                                       ///
///  License:     MIT License (see LICENSE.md)                                         ///
///                                                                                   ///
///////////////////////////////////////////////////////////////////////////////////////////
#if !defined MML_NONLINEAR_SYSTEM_SOLVERS_H
#define MML_NONLINEAR_SYSTEM_SOLVERS_H

#include <mml/MMLBase.h>
#include <mml/base/AlgorithmTypes.h>
#include <mml/base/Matrix/Matrix.h>
#include <mml/base/Matrix/MatrixNM.h>
#include <mml/base/Vector/Vector.h>
#include <mml/base/Vector/VectorN.h>
#include <mml/core/Derivation/Jacobians.h>
#include <mml/core/LinAlgEqSolvers/LinAlgDirect.h>
#include <mml/interfaces/IFunction.h>

#include <algorithm>
#include <cmath>
#include <functional>
#include <limits>
#include <string>
#include <utility>
#include <vector>

namespace MML::RootFinding {

	struct NonlinearSystemConfig {
		Real residual_tolerance = std::max(REAL(1e-10), REAL(100.0) * std::numeric_limits<Real>::epsilon());
		Real step_tolerance = std::max(REAL(1e-12), REAL(10.0) * std::numeric_limits<Real>::epsilon());
		Real relative_step_tolerance = 0.0;
		int max_iterations = 100;
		Real numerical_jacobian_step = 0.0; // Absolute step; zero uses coordinate-scaled automatic steps
		bool use_backtracking = true;
		Real backtracking_factor = 0.5;
		Real minimum_damping = 1e-6;
		int max_backtracking_steps = 20;
		bool store_trace = false;
		bool check_finite = true;
	};

	template<typename VectorType>
	struct NonlinearSystemIteration {
		int iteration = 0;
		VectorType point{};
		Real residual_norm = 0.0;
		Real step_norm = 0.0;
		Real damping = 1.0;
	};

	template<typename VectorType>
	struct BasicNonlinearSystemResult : public IterativeResultBase {
		VectorType solution{};
		VectorType residual{};
		Real residual_norm = std::numeric_limits<Real>::infinity();
		Real step_norm = std::numeric_limits<Real>::infinity();
		int jacobian_evaluations = 0;
		int linear_solves = 0;
		std::vector<NonlinearSystemIteration<VectorType>> trace;

		[[nodiscard]] bool IsSuccess() const noexcept {
			return converged && status == AlgorithmStatus::Success;
		}

		explicit operator bool() const noexcept { return IsSuccess(); }
	};

	using NonlinearSystemResult = BasicNonlinearSystemResult<Vector<Real>>;
	template<int N>
	using NonlinearSystemResultN = BasicNonlinearSystemResult<VectorN<Real, N>>;

	using DynamicSystemFunction = std::function<Vector<Real>(const Vector<Real>&)>;
	using DynamicJacobianFunction = std::function<Matrix<Real>(const Vector<Real>&)>;

	inline bool IsValidConfig(const NonlinearSystemConfig& config) {
		return std::isfinite(config.residual_tolerance) && config.residual_tolerance > 0.0
			&& std::isfinite(config.step_tolerance) && config.step_tolerance > 0.0
			&& std::isfinite(config.relative_step_tolerance) && config.relative_step_tolerance >= 0.0
			&& config.max_iterations > 0
			&& std::isfinite(config.numerical_jacobian_step) && config.numerical_jacobian_step >= 0.0
			&& std::isfinite(config.backtracking_factor) && config.backtracking_factor > 0.0
			&& config.backtracking_factor < 1.0
			&& std::isfinite(config.minimum_damping) && config.minimum_damping > 0.0
			&& config.minimum_damping <= 1.0
			&& config.max_backtracking_steps >= 0;
	}

	namespace Detail {
		inline bool IsFiniteVector(const Vector<Real>& vector) {
			for (int index = 0; index < vector.size(); ++index)
				if (!std::isfinite(vector[index])) return false;
			return true;
		}

		inline bool IsFiniteMatrix(const Matrix<Real>& matrix) {
			for (int row = 0; row < matrix.rows(); ++row)
				for (int column = 0; column < matrix.cols(); ++column)
					if (!std::isfinite(matrix(row, column))) return false;
			return true;
		}

		inline Real EffectiveStepTolerance(const NonlinearSystemConfig& config,
			const Vector<Real>& point) {
			return config.step_tolerance
				+ config.relative_step_tolerance * std::max(Real(1.0), point.NormL2());
		}

	} // namespace Detail

	inline NonlinearSystemResult SolveNonlinearSystemNewton(
		const DynamicSystemFunction& function,
		const DynamicJacobianFunction& jacobianFunction,
		const Vector<Real>& initialGuess,
		const NonlinearSystemConfig& config = {}) {
		AlgorithmTimer timer;
		NonlinearSystemResult result;
		result.algorithm_name = jacobianFunction ? "NonlinearNewtonAnalyticalJacobian"
			: "NonlinearNewtonNumericalJacobian";
		result.solution = initialGuess;
		if (!function || !IsValidConfig(config) || initialGuess.size() == 0
			|| (config.check_finite && !Detail::IsFiniteVector(initialGuess))) {
			result.status = AlgorithmStatus::InvalidInput;
			result.error_message = "SolveNonlinearSystemNewton: invalid configuration or initial guess";
			result.elapsed_time_ms = timer.elapsed_ms();
			return result;
		}

		try {
			++result.function_evaluations;
			result.residual = function(result.solution);
		} catch (const std::exception& error) {
			result.status = AlgorithmStatus::AlgorithmSpecificFailure;
			result.error_message = error.what();
			result.elapsed_time_ms = timer.elapsed_ms();
			return result;
		}
		if (result.residual.size() != result.solution.size()) {
			result.status = AlgorithmStatus::InvalidInput;
			result.error_message = "SolveNonlinearSystemNewton: function output dimension mismatch";
			result.elapsed_time_ms = timer.elapsed_ms();
			return result;
		}
		if (config.check_finite && !Detail::IsFiniteVector(result.residual)) {
			result.status = AlgorithmStatus::NumericalInstability;
			result.error_message = "SolveNonlinearSystemNewton: non-finite initial residual";
			result.elapsed_time_ms = timer.elapsed_ms();
			return result;
		}
		result.residual_norm = result.residual.NormL2();
		result.step_norm = 0.0;
		if (result.residual_norm <= config.residual_tolerance) {
			result.converged = true;
			result.status = AlgorithmStatus::Success;
			result.achieved_tolerance = result.residual_norm;
			result.elapsed_time_ms = timer.elapsed_ms();
			return result;
		}

		for (int iteration = 0; iteration < config.max_iterations; ++iteration) {
			Matrix<Real> jacobian;
			try {
				++result.jacobian_evaluations;
				jacobian = jacobianFunction
					? jacobianFunction(result.solution)
					: Derivation::calcJacobianDyn(function, result.solution,
						config.numerical_jacobian_step, &result.function_evaluations);
			} catch (const std::exception& error) {
				result.status = AlgorithmStatus::NumericalInstability;
				result.error_message = error.what();
				break;
			}
			if (jacobian.rows() != result.solution.size() || jacobian.cols() != result.solution.size()) {
				result.status = AlgorithmStatus::InvalidInput;
				result.error_message = "SolveNonlinearSystemNewton: Jacobian dimension mismatch";
				break;
			}
			if (config.check_finite && !Detail::IsFiniteMatrix(jacobian)) {
				result.status = AlgorithmStatus::NumericalInstability;
				result.error_message = "SolveNonlinearSystemNewton: non-finite Jacobian";
				break;
			}

			Vector<Real> rightHandSide = -result.residual;
			Vector<Real> step;
			try {
				LUSolver<Real> solver(jacobian);
				step = solver.Solve(rightHandSide);
				++result.linear_solves;
			} catch (const std::exception& error) {
				result.status = AlgorithmStatus::SingularMatrix;
				result.error_message = error.what();
				break;
			}
			if (config.check_finite && !Detail::IsFiniteVector(step)) {
				result.status = AlgorithmStatus::NumericalInstability;
				result.error_message = "SolveNonlinearSystemNewton: non-finite Newton step";
				break;
			}

			Real damping = 1.0;
			Vector<Real> trialPoint;
			Vector<Real> trialResidual;
			Real trialNorm = std::numeric_limits<Real>::infinity();
			bool accepted = false;
			const int attempts = config.use_backtracking ? config.max_backtracking_steps + 1 : 1;
			for (int attempt = 0; attempt < attempts; ++attempt) {
				trialPoint = result.solution + damping * step;
				try {
					++result.function_evaluations;
					trialResidual = function(trialPoint);
				} catch (const std::exception&) {
					trialResidual = Vector<Real>();
				}
				if (trialResidual.size() == result.solution.size()
					&& (!config.check_finite || Detail::IsFiniteVector(trialResidual))) {
					trialNorm = trialResidual.NormL2();
					if (!config.use_backtracking || trialNorm < result.residual_norm) {
						accepted = true;
						break;
					}
				}
				damping *= config.backtracking_factor;
				if (damping < config.minimum_damping) break;
			}
			if (!accepted) {
				result.status = AlgorithmStatus::Stalled;
				result.error_message = "SolveNonlinearSystemNewton: backtracking failed to reduce residual";
				break;
			}

			result.step_norm = damping * step.NormL2();
			result.solution = trialPoint;
			result.residual = trialResidual;
			result.residual_norm = trialNorm;
			result.iterations_used = iteration + 1;
			if (config.store_trace)
				result.trace.push_back({iteration + 1, result.solution, result.residual_norm,
					result.step_norm, damping});

			const Real stepTolerance = Detail::EffectiveStepTolerance(config, result.solution);
			if (result.residual_norm <= config.residual_tolerance) {
				result.converged = true;
				result.status = AlgorithmStatus::Success;
				result.achieved_tolerance = result.residual_norm;
				break;
			}
			if (result.step_norm <= stepTolerance) {
				result.status = AlgorithmStatus::Stalled;
				result.error_message = "SolveNonlinearSystemNewton: step tolerance reached before residual tolerance";
				break;
			}
		}

		if (!result.converged && result.status == AlgorithmStatus::Success) {
			result.status = AlgorithmStatus::MaxIterationsExceeded;
			result.error_message = "SolveNonlinearSystemNewton: maximum iterations exceeded";
		}
		result.achieved_tolerance = result.residual_norm;
		result.elapsed_time_ms = timer.elapsed_ms();
		return result;
	}

	inline NonlinearSystemResult SolveNonlinearSystemNewton(
		const DynamicSystemFunction& function,
		const Vector<Real>& initialGuess,
		const NonlinearSystemConfig& config = {}) {
		return SolveNonlinearSystemNewton(function, DynamicJacobianFunction{}, initialGuess, config);
	}

	namespace Detail {
		template<int N>
		Vector<Real> ToDynamic(const VectorN<Real, N>& vector) {
			Vector<Real> dynamic(N);
			for (int index = 0; index < N; ++index) dynamic[index] = vector[index];
			return dynamic;
		}

		template<int N>
		VectorN<Real, N> ToStatic(const Vector<Real>& vector) {
			VectorN<Real, N> fixed;
				for (int index = 0; index < N; ++index)
					fixed[index] = index < vector.size() ? vector[index] : Real(0.0);
			return fixed;
		}

		template<int N>
		Matrix<Real> ToDynamic(const MatrixNM<Real, N, N>& matrix) {
			Matrix<Real> dynamic(N, N);
			for (int row = 0; row < N; ++row)
				for (int column = 0; column < N; ++column)
					dynamic(row, column) = matrix(row, column);
			return dynamic;
		}

		template<int N>
		NonlinearSystemResultN<N> ToStaticResult(const NonlinearSystemResult& dynamic) {
			NonlinearSystemResultN<N> result;
			result.solution = ToStatic<N>(dynamic.solution);
			result.residual = ToStatic<N>(dynamic.residual);
			result.residual_norm = dynamic.residual_norm;
			result.step_norm = dynamic.step_norm;
			result.jacobian_evaluations = dynamic.jacobian_evaluations;
			result.linear_solves = dynamic.linear_solves;
			result.converged = dynamic.converged;
			result.iterations_used = dynamic.iterations_used;
			result.achieved_tolerance = dynamic.achieved_tolerance;
			result.status = dynamic.status;
			result.error_message = dynamic.error_message;
			result.algorithm_name = dynamic.algorithm_name;
			result.elapsed_time_ms = dynamic.elapsed_time_ms;
			result.function_evaluations = dynamic.function_evaluations;
			result.trace.reserve(dynamic.trace.size());
			for (const auto& entry : dynamic.trace)
				result.trace.push_back({entry.iteration, ToStatic<N>(entry.point),
					entry.residual_norm, entry.step_norm, entry.damping});
			return result;
		}
	} // namespace Detail

	template<int N, typename Function, typename Jacobian>
	NonlinearSystemResultN<N> SolveNonlinearSystemNewton(Function&& function,
		Jacobian&& jacobian, const VectorN<Real, N>& initialGuess,
		const NonlinearSystemConfig& config = {}) {
		DynamicSystemFunction dynamicFunction = [function = std::forward<Function>(function)](
			const Vector<Real>& point) mutable {
			VectorN<Real, N> staticPoint = Detail::ToStatic<N>(point);
			return Detail::ToDynamic<N>(std::invoke(function, staticPoint));
		};
		DynamicJacobianFunction dynamicJacobian = [jacobian = std::forward<Jacobian>(jacobian)](
			const Vector<Real>& point) mutable {
			VectorN<Real, N> staticPoint = Detail::ToStatic<N>(point);
			return Detail::ToDynamic<N>(std::invoke(jacobian, staticPoint));
		};
		return Detail::ToStaticResult<N>(SolveNonlinearSystemNewton(
			dynamicFunction, dynamicJacobian, Detail::ToDynamic<N>(initialGuess), config));
	}

	template<int N, typename Function>
	NonlinearSystemResultN<N> SolveNonlinearSystemNewton(Function&& function,
		const VectorN<Real, N>& initialGuess, const NonlinearSystemConfig& config = {}) {
		DynamicSystemFunction dynamicFunction = [function = std::forward<Function>(function)](
			const Vector<Real>& point) mutable {
			VectorN<Real, N> staticPoint = Detail::ToStatic<N>(point);
			return Detail::ToDynamic<N>(std::invoke(function, staticPoint));
		};
		return Detail::ToStaticResult<N>(SolveNonlinearSystemNewton(
			dynamicFunction, Detail::ToDynamic<N>(initialGuess), config));
	}

	template<int N>
	NonlinearSystemResultN<N> SolveNonlinearSystemNewton(const IVectorFunction<N>& function,
		const VectorN<Real, N>& initialGuess, const NonlinearSystemConfig& config = {}) {
		return SolveNonlinearSystemNewton<N>(
			[&function](const VectorN<Real, N>& point) { return function(point); },
			initialGuess, config);
	}

	template<int N, typename Jacobian>
	NonlinearSystemResultN<N> SolveNonlinearSystemNewton(const IVectorFunction<N>& function,
		Jacobian&& jacobian, const VectorN<Real, N>& initialGuess,
		const NonlinearSystemConfig& config = {}) {
		return SolveNonlinearSystemNewton<N>(
			[&function](const VectorN<Real, N>& point) { return function(point); },
			std::forward<Jacobian>(jacobian), initialGuess, config);
	}

} // namespace MML::RootFinding

#endif // MML_NONLINEAR_SYSTEM_SOLVERS_H
