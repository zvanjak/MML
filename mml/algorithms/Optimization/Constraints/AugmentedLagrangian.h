///////////////////////////////////////////////////////////////////////////////////////////
///                         MinimalMathLibrary (MML)                                  ///
///                                                                                   ///
///  File:        Optimization/Constraints/AugmentedLagrangian.h                      ///
///  Description: Augmented Lagrangian Method (ALM) for general constrained           ///
///               optimization - converts constrained problems into a sequence of     ///
///               Nelder-Mead subproblems, optionally restricted to a box               ///
///                                                                                   ///
///  Copyright:   (c) 2024-2026 Zvonimir Vanjak                                       ///
///  License:     MIT License (see LICENSE.md)                                         ///
///                                                                                   ///
///////////////////////////////////////////////////////////////////////////////////////////
//
// Problem form:
//   minimize   f(x)
//   subject to g_i(x) <= 0, i = 1,...,m  (inequality)
//              h_j(x) = 0,  j = 1,...,p  (equality)
//
// Augmented Lagrangian function:
//   L_A(x, lambda, mu, rho) = f(x)
//       + Sum_i (lambda_i * max(g_i(x), -lambda_i/rho) + rho/2 * max(g_i(x), -lambda_i/rho)^2)
//       + Sum_j (mu_j * h_j(x) + rho/2 * h_j(x)^2)
//
// References:
// - Nocedal & Wright (2006). "Numerical Optimization", Chapter 17
// - Conn, Gould & Toint (1991). "A Globally Convergent Augmented Lagrangian Algorithm"
//
// History: repatriated from MML-Packages optimization (2026-10-03) - it is built on
// VectorN/IScalarFunction and MML's NelderMead, and completes the Constraints/ tier
// (bounds -> projected gradient -> box-constrained -> general nonlinear constraints).
///////////////////////////////////////////////////////////////////////////////////////////
#if !defined MML_AUGMENTED_LAGRANGIAN_H
#define MML_AUGMENTED_LAGRANGIAN_H

#include <mml/MMLBase.h>

#include <mml/base/AlgorithmTypes.h>
#include <mml/base/Function.h>
#include <mml/base/Vector/VectorN.h>
#include <mml/interfaces/IFunction.h>

#include <mml/algorithms/Optimization/Constraints/BoxConstrainedNelderMead.h>

#include <algorithm>
#include <cmath>
#include <functional>
#include <limits>
#include <optional>
#include <ostream>
#include <string>
#include <type_traits>
#include <utility>
#include <vector>

namespace MML::Optimization
{
	/// Configuration for the Augmented Lagrangian Method.
	struct ALMConfig {
		// Outer loop
		int max_outer_iterations = 100;         ///< Maximum outer (multiplier-update) iterations
		Real feasibility_tolerance = Precision::LooseTolerance;      ///< Constraint violation tolerance
		Real stationarity_tolerance = Precision::LooseTolerance;     ///< Optimality tolerance

		// Penalty parameter
		Real initial_penalty = 10.0;            ///< Initial penalty parameter rho_0
		Real penalty_multiplier = 10.0;         ///< Penalty increase factor
		Real max_penalty = Precision::PenaltyParameterMax;                ///< Maximum penalty parameter

		// Multiplier updates
		Real multiplier_lower_bound = -Precision::LagrangeMultiplierBound;    ///< Lower bound on Lagrange multipliers
		Real multiplier_upper_bound = Precision::LagrangeMultiplierBound;     ///< Upper bound on Lagrange multipliers

		// Inner solver
		int max_inner_iterations = 1000;        ///< Nelder-Mead evaluation/iteration budget (legacy field name; not a strict call cap)
		Real inner_tolerance = Precision::ModerateTolerance;            ///< Tolerance for the inner solver
		Real inner_simplex_delta = 1.0;         ///< Initial simplex size for the inner solver

		/// Optional progress log (contract policy section 7: no std::cout in library code)
		std::ostream* verbose_stream = nullptr;
	};

	/// Result of an Augmented Lagrangian optimization.
	/// Extends the standard iterative-result contract; `iterations_used` counts outer
	/// iterations, `achieved_tolerance` reports the final max constraint violation.
	template<int N>
		requires (N > 0)
	struct [[nodiscard]] ALMResult : public IterativeResultBase {
		VectorN<Real, N> x;                     ///< Optimal solution
		Real objective_value = NoEstimate;      ///< Objective function value at the returned point
		Real constraint_violation = NoEstimate; ///< Maximum nonlinear and bound violation
		Real step_norm = NoEstimate;            ///< Max-norm change in the last outer iteration
		Real stationarity_measure = NoEstimate; ///< Relative outer-step proxy after a converged inner solve; not a KKT residual
		std::optional<AlgorithmStatus> last_inner_status; ///< Absent until an inner solve is attempted
		int outer_iterations = 0;               ///< Outer iterations performed
		int inner_iterations = 0;               ///< Total inner solver iterations

		// Lagrange multiplier estimates
		std::vector<Real> inequality_multipliers;   ///< lambda for g(x) <= 0
		std::vector<Real> equality_multipliers;     ///< mu for h(x) = 0
	};

	/// @brief Augmented Lagrangian Method for general nonlinear constraints.
	/// Constraints are callables Real(const VectorN<Real,N>&); the objective may be an
	/// IScalarFunction<N> or any fixed-N callable. Inner subproblems use Nelder-Mead.
	template<int N>
		requires (N > 0)
	class AugmentedLagrangian {
	public:
		using Vec = VectorN<Real, N>;
		using Result = ALMResult<N>;
		using Config = ALMConfig;

	private:
		Config _config;

		int _numInequality = 0;
		int _numEquality = 0;
		std::vector<std::function<Real(const Vec&)>> _inequalityConstraints;
		std::vector<std::function<Real(const Vec&)>> _equalityConstraints;

		std::optional<BoundConstraints> _bounds;

		std::vector<Real> _lambda;   // Inequality multipliers
		std::vector<Real> _mu;       // Equality multipliers
		Real _rho;                   // Penalty parameter

		static void ValidateConfig(const Config& config) {
			if (config.max_outer_iterations < 0 || config.max_inner_iterations <= 0)
				throw MultidimOptimizationInputError("ALM requires a nonnegative outer budget and a positive inner budget");
			ValidateMultidimTolerance(config.feasibility_tolerance, "ALM feasibility_tolerance");
			ValidateMultidimTolerance(config.stationarity_tolerance, "ALM stationarity_tolerance");
			ValidateMultidimTolerance(config.inner_tolerance, "ALM inner_tolerance");
			if (!std::isfinite(config.initial_penalty) || config.initial_penalty <= 0 ||
				!std::isfinite(config.max_penalty) || config.max_penalty < config.initial_penalty ||
				!std::isfinite(config.penalty_multiplier) || config.penalty_multiplier < 1)
				throw MultidimOptimizationInputError("ALM requires finite positive initial_penalty <= max_penalty and penalty_multiplier >= 1");
			if (!std::isfinite(config.multiplier_lower_bound) || !std::isfinite(config.multiplier_upper_bound) ||
				config.multiplier_lower_bound > config.multiplier_upper_bound ||
				config.multiplier_lower_bound > 0 || config.multiplier_upper_bound < 0)
				throw MultidimOptimizationInputError("ALM multiplier bounds must be finite, ordered and contain zero");
			if (!std::isfinite(config.inner_simplex_delta) || config.inner_simplex_delta == 0)
				throw MultidimOptimizationInputError("ALM inner_simplex_delta must be finite and nonzero");
		}

	public:
		AugmentedLagrangian() : AugmentedLagrangian(Config()) {}

		explicit AugmentedLagrangian(const Config& config)
			: _config(config), _rho(config.initial_penalty) {
			ValidateConfig(config);
		}

		void SetConfig(const Config& config) {
			ValidateConfig(config);
			_config = config;
			_rho = config.initial_penalty;
		}

		const Config& GetConfig() const { return _config; }

		/// @brief Set variable bounds
		void SetBounds(const Vec& lower, const Vec& upper) {
			BoundConstraints bounds(lower, upper);
			_bounds = std::move(bounds);
		}

		/// @brief Set uniform bounds for all variables
		void SetBounds(Real lower, Real upper) {
			Vec lo, hi;
			for (int i = 0; i < N; ++i) {
				lo[i] = lower;
				hi[i] = upper;
			}
			SetBounds(lo, hi);
		}

		/// @brief Add inequality constraint g(x) <= 0
		void AddInequalityConstraint(std::function<Real(const Vec&)> g) {
			if (!g)
				throw MultidimOptimizationInputError("ALM inequality constraint callable cannot be empty");
			_inequalityConstraints.push_back(std::move(g));
			_numInequality++;
			_lambda.push_back(0.0);
		}

		/// @brief Add equality constraint h(x) = 0
		void AddEqualityConstraint(std::function<Real(const Vec&)> h) {
			if (!h)
				throw MultidimOptimizationInputError("ALM equality constraint callable cannot be empty");
			_equalityConstraints.push_back(std::move(h));
			_numEquality++;
			_mu.push_back(0.0);
		}

		/// @brief Clear all constraints
		void ClearConstraints() {
			_inequalityConstraints.clear();
			_equalityConstraints.clear();
			_lambda.clear();
			_mu.clear();
			_numInequality = 0;
			_numEquality = 0;
		}

		/// @brief Optimize using any objective callable Real(const Vec&)
		template<typename ObjectiveFunc>
			requires (std::is_invocable_r_v<Real, ObjectiveFunc&, const Vec&> ||
				std::is_invocable_r_v<Real, ObjectiveFunc&, Vector<Real>&>)
		Result Optimize(ObjectiveFunc&& objective, const Vec& x0) {
			ValidateConfig(_config);
			ValidateVectorFinite<N>(x0, "ALM starting point");
			AlgorithmTimer timer;

			Result result;
			result.algorithm_name = "AugmentedLagrangian";
			result.status = AlgorithmStatus::MaxIterationsExceeded;

			_lambda.assign(_numInequality, 0.0);
			_mu.assign(_numEquality, 0.0);
			_rho = _config.initial_penalty;

			Vec xk = _bounds ? _bounds->Project(x0) : x0;
			result.constraint_violation = ComputeMaxViolation(xk);
			auto countedObjective = [&](const Vec& point) {
				ValidateVectorFinite<N>(point, "ALM objective point");
				if (result.function_evaluations == std::numeric_limits<int>::max())
					throw MultidimOptimizationError("ALM objective-call counter exceeded its representable range");
				++result.function_evaluations;
				const Real value = CallObjective(objective, point);
				ValidateMultidimFunctionValue(value, "ALM objective");
				return value;
			};

			for (int outer = 0; outer < _config.max_outer_iterations; ++outer) {
				result.outer_iterations = outer + 1;

				auto augLagrangianFunc = [this, &countedObjective](const Vec& x) -> Real {
					return EvaluateAugmentedLagrangian(x, countedObjective);
				};
				ScalarFunctionFromStdFunc<N> augLagrangianWrapper(augLagrangianFunc);

				MultidimOptimizationResult innerResult;
				if (_bounds) {
					MultidimOptimizationConfig innerConfig(_config.inner_tolerance, _config.max_inner_iterations);
					innerConfig.initial_delta = _config.inner_simplex_delta;
					innerResult = BoxConstrainedNelderMead(innerConfig).Minimize(augLagrangianWrapper, xk, *_bounds);
				} else {
					NelderMead innerSolver(_config.inner_tolerance, _config.max_inner_iterations);
					innerResult = innerSolver.Minimize(augLagrangianWrapper, xk, _config.inner_simplex_delta);
				}
				result.last_inner_status = innerResult.status;
				result.stationarity_measure = NoEstimate;

				const Vec xPrev = xk;
				for (int i = 0; i < N; ++i)
					xk[i] = innerResult.xmin[i];
				ValidateVectorFinite<N>(xk, "ALM inner minimizer");
				if (innerResult.iterations > std::numeric_limits<int>::max() - result.inner_iterations)
					throw MultidimOptimizationError("ALM inner-iteration counter exceeded its representable range");
				result.inner_iterations += innerResult.iterations;

				Real maxViolation = ComputeMaxViolation(xk);
				result.constraint_violation = maxViolation;

				Real stepNorm = static_cast<MML::Real>(0.0);
				Real xScale = static_cast<MML::Real>(1.0);
				for (int i = 0; i < N; ++i) {
					stepNorm = std::max(stepNorm, std::abs(xk[i] - xPrev[i]));
					xScale = std::max(xScale, std::abs(xk[i]));
				}
				result.step_norm = stepNorm;
				ValidateMultidimFunctionValue(stepNorm, "ALM outer step norm");

				if (!innerResult.converged || innerResult.status != AlgorithmStatus::Success) {
					result.status = IsUsable(innerResult.status) ? AlgorithmStatus::AlgorithmSpecificFailure : innerResult.status;
					result.error_message = "ALM inner solve failed at outer iteration " + std::to_string(outer + 1) +
						" (" + MML::to_string(innerResult.status) + "): " + innerResult.error_message;
					break;
				}
				result.stationarity_measure = stepNorm / xScale;
				const bool stationary = outer > 0 && result.stationarity_measure <= _config.stationarity_tolerance;

				if (_config.verbose_stream)
					(*_config.verbose_stream) << "ALM outer iter " << (outer + 1)
						<< ": violation=" << maxViolation << ", step=" << stepNorm << ", rho=" << _rho << "\n";

				if (maxViolation <= _config.feasibility_tolerance && stationary) {
					result.converged = true;
					result.status = AlgorithmStatus::Success;
					break;
				}

				UpdateMultipliers(xk);

				if (outer > 0 && maxViolation > _config.feasibility_tolerance)
					_rho = _rho >= _config.max_penalty / _config.penalty_multiplier ?
						_config.max_penalty : _rho * _config.penalty_multiplier;
			}

			result.x = xk;
			result.objective_value = countedObjective(xk);
			result.iterations_used = result.outer_iterations;
			result.achieved_tolerance = result.constraint_violation;
			result.inequality_multipliers = _lambda;
			result.equality_multipliers = _mu;
			result.elapsed_time_ms = timer.elapsed_ms();

			if (!result.converged && result.error_message.empty()) {
				result.error_message = "Failed to reach feasibility and stationarity tolerances after "
					+ std::to_string(result.outer_iterations) + " outer iterations (violation: "
					+ std::to_string(result.constraint_violation) + ", last step: "
					+ std::to_string(result.step_norm) + ")";
			}

			return result;
		}

		/// @brief Optimize an IScalarFunction<N> objective
		Result Optimize(IScalarFunction<N>& objective, const Vec& x0) {
			auto wrapper = [&objective](const Vec& x) { return objective(x); };
			return Optimize(wrapper, x0);
		}

	private:
		Real ConstraintValue(const std::function<Real(const Vec&)>& constraint, const Vec& x,
			const char* kind, int index) const {
			const Real value = constraint(x);
			if (!std::isfinite(value))
				throw MultidimOptimizationInputError(std::string("ALM ") + kind + " constraint " +
					std::to_string(index) + " returned a non-finite value");
			return value;
		}

		/// @brief Call objective (fixed-N callables or dynamic Vector<Real> objectives)
		template<typename ObjectiveFunc>
		Real CallObjective(ObjectiveFunc& obj, const Vec& x) const {
			if constexpr (std::is_invocable_r_v<Real, ObjectiveFunc, const Vec&>) {
				return obj(x);
			} else {
				Vector<Real> xd(N);
				for (int i = 0; i < N; ++i) xd[i] = x[i];
				return obj(xd);
			}
		}

		/// @brief Evaluate the augmented Lagrangian at x
		template<typename ObjectiveFunc>
		Real EvaluateAugmentedLagrangian(const Vec& x, ObjectiveFunc& objective) const {
			Real LA = CallObjective(objective, x);

			// Inequality: lambda_i * max(g_i, -lambda_i/rho) + (rho/2) * max(g_i, -lambda_i/rho)^2
			for (int i = 0; i < _numInequality; ++i) {
				Real g = ConstraintValue(_inequalityConstraints[i], x, "inequality", i);
				Real threshold = -_lambda[i] / _rho;
				Real shifted = std::max(g, threshold);
				LA += _lambda[i] * shifted + (_rho / MML::Real(2.0)) * shifted * shifted;
			}

			// Equality: mu_j * h_j + (rho/2) * h_j^2
			for (int j = 0; j < _numEquality; ++j) {
				Real h = ConstraintValue(_equalityConstraints[j], x, "equality", j);
				LA += _mu[j] * h + (_rho / MML::Real(2.0)) * h * h;
			}

			ValidateMultidimFunctionValue(LA, "ALM augmented Lagrangian");
			return LA;
		}

		/// @brief Maximum constraint violation at x
		Real ComputeMaxViolation(const Vec& x) const {
			ValidateVectorFinite<N>(x, "ALM feasibility point");
			Real maxViol = _bounds ? _bounds->MaxViolation(x) : Real(0);

			for (int i = 0; i < _numInequality; ++i) {
				Real g = ConstraintValue(_inequalityConstraints[i], x, "inequality", i);
				maxViol = std::max(maxViol, std::max(Real(0), g));
			}
			for (int j = 0; j < _numEquality; ++j) {
				Real h = ConstraintValue(_equalityConstraints[j], x, "equality", j);
				maxViol = std::max(maxViol, std::abs(h));
			}

			return maxViol;
		}

		/// @brief First-order Lagrange multiplier updates
		void UpdateMultipliers(const Vec& x) {
			// Inequality: lambda = clamp(max(0, lambda + rho * g(x)))
			for (int i = 0; i < _numInequality; ++i) {
				Real g = ConstraintValue(_inequalityConstraints[i], x, "inequality", i);
				const Real update = _lambda[i] + _rho * g;
				ValidateMultidimFunctionValue(update, "ALM inequality multiplier update");
				_lambda[i] = std::max(Real(0), update);
				_lambda[i] = std::clamp(_lambda[i], _config.multiplier_lower_bound, _config.multiplier_upper_bound);
			}
			// Equality: mu = clamp(mu + rho * h(x))
			for (int j = 0; j < _numEquality; ++j) {
				Real h = ConstraintValue(_equalityConstraints[j], x, "equality", j);
				_mu[j] = _mu[j] + _rho * h;
				ValidateMultidimFunctionValue(_mu[j], "ALM equality multiplier update");
				_mu[j] = std::clamp(_mu[j], _config.multiplier_lower_bound, _config.multiplier_upper_bound);
			}
		}
	};

	/// @brief Convenience: solve a constrained problem with the Augmented Lagrangian Method
	template<int N, typename ObjectiveFunc>
		requires (N > 0 && (std::is_invocable_r_v<Real, ObjectiveFunc&, const VectorN<Real, N>&> ||
			std::is_invocable_r_v<Real, ObjectiveFunc&, Vector<Real>&>))
	ALMResult<N> SolveConstrained(ObjectiveFunc&& objective, const VectorN<Real, N>& x0,
		const std::vector<std::function<Real(const VectorN<Real, N>&)>>& inequalityConstraints = {},
		const std::vector<std::function<Real(const VectorN<Real, N>&)>>& equalityConstraints = {},
		const ALMConfig& config = ALMConfig()) {
		AugmentedLagrangian<N> alm(config);

		for (const auto& g : inequalityConstraints)
			alm.AddInequalityConstraint(g);
		for (const auto& h : equalityConstraints)
			alm.AddEqualityConstraint(h);

		return alm.Optimize(std::forward<ObjectiveFunc>(objective), x0);
	}

} // namespace MML::Optimization
#endif // MML_AUGMENTED_LAGRANGIAN_H
