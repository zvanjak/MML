///////////////////////////////////////////////////////////////////////////////////////////
///                         MinimalMathLibrary (MML)                                  ///
///                                                                                   ///
///  File:        Optimization/Constraints/AugmentedLagrangian.h                      ///
///  Description: Augmented Lagrangian Method (ALM) for general constrained           ///
///               optimization - converts constrained problems into a sequence of     ///
///               unconstrained subproblems solved by Nelder-Mead                     ///
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

#include <mml/algorithms/Optimization/Multidim/NelderMead.h>

#include <algorithm>
#include <cmath>
#include <functional>
#include <limits>
#include <ostream>
#include <string>
#include <vector>

namespace MML::Optimization
{
	/// Configuration for the Augmented Lagrangian Method.
	struct ALMConfig {
		// Outer loop
		int max_outer_iterations = 100;         ///< Maximum outer (multiplier-update) iterations
		Real feasibility_tolerance = 1e-6;      ///< Constraint violation tolerance
		Real stationarity_tolerance = 1e-6;     ///< Optimality tolerance

		// Penalty parameter
		Real initial_penalty = 10.0;            ///< Initial penalty parameter rho_0
		Real penalty_multiplier = 10.0;         ///< Penalty increase factor
		Real max_penalty = 1e12;                ///< Maximum penalty parameter

		// Multiplier updates
		Real multiplier_lower_bound = -1e10;    ///< Lower bound on Lagrange multipliers
		Real multiplier_upper_bound = 1e10;     ///< Upper bound on Lagrange multipliers

		// Inner (unconstrained) solver
		int max_inner_iterations = 1000;        ///< Max iterations for the Nelder-Mead subproblem
		Real inner_tolerance = 1e-8;            ///< Tolerance for the inner solver
		Real inner_simplex_delta = 1.0;         ///< Initial simplex size for the inner solver

		/// Optional progress log (contract policy section 7: no std::cout in library code)
		std::ostream* verbose_stream = nullptr;
	};

	/// Result of an Augmented Lagrangian optimization.
	/// Extends the standard iterative-result contract; `iterations_used` counts outer
	/// iterations, `achieved_tolerance` reports the final max constraint violation.
	template<int N>
	struct ALMResult : public IterativeResultBase {
		VectorN<Real, N> x;                     ///< Optimal solution
		Real objective_value = 0;               ///< Objective function value f(x*)
		Real constraint_violation = 0;          ///< Final maximum constraint violation
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

		Vec _lowerBounds;
		Vec _upperBounds;
		bool _hasBounds = false;

		std::vector<Real> _lambda;   // Inequality multipliers
		std::vector<Real> _mu;       // Equality multipliers
		Real _rho;                   // Penalty parameter

	public:
		AugmentedLagrangian() : _rho(10.0) { InitializeBounds(); }

		explicit AugmentedLagrangian(const Config& config)
			: _config(config), _rho(config.initial_penalty) {
			InitializeBounds();
		}

		void SetConfig(const Config& config) {
			_config = config;
			_rho = config.initial_penalty;
		}

		const Config& GetConfig() const { return _config; }

		/// @brief Set variable bounds
		void SetBounds(const Vec& lower, const Vec& upper) {
			_lowerBounds = lower;
			_upperBounds = upper;
			_hasBounds = true;
		}

		/// @brief Set uniform bounds for all variables
		void SetBounds(Real lower, Real upper) {
			for (int i = 0; i < N; ++i) {
				_lowerBounds[i] = lower;
				_upperBounds[i] = upper;
			}
			_hasBounds = true;
		}

		/// @brief Add inequality constraint g(x) <= 0
		void AddInequalityConstraint(std::function<Real(const Vec&)> g) {
			_inequalityConstraints.push_back(std::move(g));
			_numInequality++;
			_lambda.push_back(0.0);
		}

		/// @brief Add equality constraint h(x) = 0
		void AddEqualityConstraint(std::function<Real(const Vec&)> h) {
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
		Result Optimize(ObjectiveFunc& objective, const Vec& x0) {
			AlgorithmTimer timer;

			Result result;
			result.algorithm_name = "AugmentedLagrangian";
			result.x = x0;

			_lambda.assign(_numInequality, 0.0);
			_mu.assign(_numEquality, 0.0);
			_rho = _config.initial_penalty;

			int totalFuncEvals = 0;
			Vec xk = x0;

			for (int outer = 0; outer < _config.max_outer_iterations; ++outer) {
				result.outer_iterations = outer + 1;

				// Inner unconstrained subproblem: minimize the augmented Lagrangian
				auto augLagrangianFunc = [this, &objective](const Vec& x) -> Real {
					return EvaluateAugmentedLagrangian(x, objective);
				};
				ScalarFunctionFromStdFunc<N> augLagrangianWrapper(augLagrangianFunc);

				NelderMead innerSolver(_config.inner_tolerance, _config.max_inner_iterations);
				auto innerResult = innerSolver.Minimize(augLagrangianWrapper, xk, _config.inner_simplex_delta);

				for (int i = 0; i < N; ++i)
					xk[i] = innerResult.xmin[i];
				result.inner_iterations += innerResult.iterations;
				totalFuncEvals += innerResult.iterations;

				Real maxViolation = ComputeMaxViolation(xk);
				result.constraint_violation = maxViolation;

				if (_config.verbose_stream)
					(*_config.verbose_stream) << "ALM outer iter " << (outer + 1)
						<< ": violation=" << maxViolation << ", rho=" << _rho << "\n";

				if (maxViolation <= _config.feasibility_tolerance) {
					result.converged = true;
					break;
				}

				UpdateMultipliers(xk);

				if (outer > 0)
					_rho = std::min(_rho * _config.penalty_multiplier, _config.max_penalty);
			}

			result.x = xk;
			result.objective_value = CallObjective(objective, xk);
			result.iterations_used = result.outer_iterations;
			result.achieved_tolerance = result.constraint_violation;
			result.function_evaluations = totalFuncEvals;
			result.inequality_multipliers = _lambda;
			result.equality_multipliers = _mu;
			result.elapsed_time_ms = timer.elapsed_ms();

			if (result.converged) {
				result.status = AlgorithmStatus::Success;
			} else {
				result.status = AlgorithmStatus::MaxIterationsExceeded;
				result.error_message = "Failed to reach feasibility tolerance after "
					+ std::to_string(result.outer_iterations) + " outer iterations (violation: "
					+ std::to_string(result.constraint_violation) + ")";
			}

			return result;
		}

		/// @brief Optimize an IScalarFunction<N> objective
		Result Optimize(IScalarFunction<N>& objective, const Vec& x0) {
			auto wrapper = [&objective](const Vec& x) { return objective(x); };
			return Optimize(wrapper, x0);
		}

	private:
		void InitializeBounds() {
			for (int i = 0; i < N; ++i) {
				_lowerBounds[i] = -std::numeric_limits<Real>::infinity();
				_upperBounds[i] = std::numeric_limits<Real>::infinity();
			}
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
				Real g = _inequalityConstraints[i](x);
				Real threshold = -_lambda[i] / _rho;
				Real shifted = std::max(g, threshold);
				LA += _lambda[i] * shifted + (_rho / 2.0) * shifted * shifted;
			}

			// Equality: mu_j * h_j + (rho/2) * h_j^2
			for (int j = 0; j < _numEquality; ++j) {
				Real h = _equalityConstraints[j](x);
				LA += _mu[j] * h + (_rho / 2.0) * h * h;
			}

			return LA;
		}

		/// @brief Maximum constraint violation at x
		Real ComputeMaxViolation(const Vec& x) const {
			Real maxViol = 0;

			for (int i = 0; i < _numInequality; ++i) {
				Real g = _inequalityConstraints[i](x);
				maxViol = std::max(maxViol, std::max(Real(0), g));
			}
			for (int j = 0; j < _numEquality; ++j) {
				Real h = _equalityConstraints[j](x);
				maxViol = std::max(maxViol, std::abs(h));
			}

			return maxViol;
		}

		/// @brief First-order Lagrange multiplier updates
		void UpdateMultipliers(const Vec& x) {
			// Inequality: lambda = clamp(max(0, lambda + rho * g(x)))
			for (int i = 0; i < _numInequality; ++i) {
				Real g = _inequalityConstraints[i](x);
				_lambda[i] = std::max(Real(0), _lambda[i] + _rho * g);
				_lambda[i] = std::clamp(_lambda[i], _config.multiplier_lower_bound, _config.multiplier_upper_bound);
			}
			// Equality: mu = clamp(mu + rho * h(x))
			for (int j = 0; j < _numEquality; ++j) {
				Real h = _equalityConstraints[j](x);
				_mu[j] = _mu[j] + _rho * h;
				_mu[j] = std::clamp(_mu[j], _config.multiplier_lower_bound, _config.multiplier_upper_bound);
			}
		}
	};

	/// @brief Convenience: solve a constrained problem with the Augmented Lagrangian Method
	template<int N, typename ObjectiveFunc>
	ALMResult<N> SolveConstrained(ObjectiveFunc& objective, const VectorN<Real, N>& x0,
		const std::vector<std::function<Real(const VectorN<Real, N>&)>>& inequalityConstraints = {},
		const std::vector<std::function<Real(const VectorN<Real, N>&)>>& equalityConstraints = {},
		const ALMConfig& config = ALMConfig()) {
		AugmentedLagrangian<N> alm(config);

		for (const auto& g : inequalityConstraints)
			alm.AddInequalityConstraint(g);
		for (const auto& h : equalityConstraints)
			alm.AddEqualityConstraint(h);

		return alm.Optimize(objective, x0);
	}

} // namespace MML::Optimization
#endif // MML_AUGMENTED_LAGRANGIAN_H
