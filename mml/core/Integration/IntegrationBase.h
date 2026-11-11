#ifndef MML_INTEGRATION_BASE_H
#define MML_INTEGRATION_BASE_H

#include <mml/MMLBase.h>
#include <mml/base/AlgorithmTypes.h>
#include <mml/MMLNumericValidation.h>

#include <algorithm>
#include <array>
#include <atomic>
#include <cmath>
#include <functional>
#include <initializer_list>
#include <limits>
#include <numeric>
#include <ostream>
#include <queue>
#include <string>
#include <type_traits>
#include <vector>

namespace MML
{
	/// Configuration for integration Detailed APIs
	struct IntegrationConfig : public EvaluationConfigBase {
		// Inherits: estimate_error, check_finite, exception_policy
	};

	/******************************************************************************/
	/*****             The one integration result type                         *****/
	/******************************************************************************/

	/// Result of a numerical integration (RESULT_POLICY: one type for every integration layer).
	///
	/// - `converged`        true only if the method's convergence criterion was checked and met
	/// - `error_estimated`  false for fixed rules (Gauss10, Gauss-Legendre, midpoint sums, ...)
	/// - `error_estimate`   a computed bound, or NoEstimate when none exists (HasEstimate())
	/// - `status`           Success (criterion met) / Completed (fixed rule) / MaxIterationsExceeded ...
	/// - `IsSuccess()`      the value is usable (Success or Completed)
	///
	/// Plain integration functions return the value and throw ConvergenceError when the result is
	/// not usable; the `...Detailed` functions return this struct. Conversion to the value is
	/// explicit (`static_cast<TValue>(r)`) or via `.value` so that discarding diagnostics is
	/// visible at the call site.
	template<typename TValue>
	struct [[nodiscard]] BasicIntegrationResult : public EvaluationResultBase {
		TValue value{};                       ///< Computed integral value
		Real error_estimate = NoEstimate;     ///< Estimated absolute error (NoEstimate until computed)
		int iterations = 0;                   ///< Iterations / refinement levels performed (1 for fixed rules)
		bool converged = false;               ///< Criterion checked and met
		bool error_estimated = false;         ///< An error estimate exists (false for fixed rules)

		BasicIntegrationResult() = default;

		/// Field-wise constructor used by the integrator cores. The status is derived:
		/// converged -> Success; no estimate -> Completed; otherwise MaxIterationsExceeded.
		BasicIntegrationResult(TValue val, Real err, int iter, bool conv, bool errEstimated = true)
			: value(val), error_estimate(errEstimated ? err : NoEstimate), iterations(iter)
			, converged(conv), error_estimated(errEstimated)
		{
			if (conv)                 status = AlgorithmStatus::Success;
			else if (!errEstimated)   status = AlgorithmStatus::Completed;
			else                      status = AlgorithmStatus::MaxIterationsExceeded;
		}

		/// Result of a fixed quadrature rule: nothing was checked, there is no estimate (RESULT_POLICY §2).
		static BasicIntegrationResult FixedRule(TValue val, int function_evals = 0)
		{
			BasicIntegrationResult r(val, NoEstimate, 1, /*converged*/ false, /*error_estimated*/ false);
			r.function_evaluations = function_evals;
			return r;
		}

		/// Mark an existing result as produced by a fixed rule.
		void MarkFixedRule(int function_evals)
		{
			error_estimate = NoEstimate;
			error_estimated = false;
			converged = false;
			iterations = 1;
			status = AlgorithmStatus::Completed;
			function_evaluations = function_evals;
		}

		/// True if the value may be used: the criterion was met, or the method has none.
		/// A checked-and-missed criterion (MaxIterationsExceeded, ...) is not usable.
		[[nodiscard]] bool IsUsable() const noexcept { return IsUsable_(status) && (status == AlgorithmStatus::Completed || converged || !error_estimated); }

		/// Returns the value, or throws ConvergenceError when the result is not usable.
		/// This is what the plain (non-Detailed) integration functions return.
		TValue ValueOrThrow(const char* context) const
		{
			if (IsUsable())
				return value;
			if (!error_message.empty())
				throw ConvergenceError(std::string(context) + ": " + error_message, iterations, error_estimate);
			throw ConvergenceError(std::string(context) + " did not converge after " + std::to_string(iterations) +
				" refinements (error estimate " + std::to_string(error_estimate) + ")", iterations, error_estimate);
		}

		/// Explicit conversion to the integral value.
		/// @warning Discards error_estimate, iterations, converged and status - check them first.
		explicit operator TValue() const { return value; }

		/// Prints the value; a result whose convergence criterion was checked and not met is marked.
		friend std::ostream& operator<<(std::ostream& out, const BasicIntegrationResult& r)
		{
			out << r.value;
			if (r.status != AlgorithmStatus::Completed && r.error_estimated && !r.converged)
				out << " (not converged)";
			return out;
		}

	private:
		static bool IsUsable_(AlgorithmStatus s) noexcept { return MML::IsUsable(s); }
	};

	/// Real-valued integration result (1D, multi-dimensional, path, surface, Monte Carlo base)
	using IntegrationResult = BasicIntegrationResult<Real>;

	/// @deprecated Same type as IntegrationResult since the result unification (wewi.3.3).
	using IntegrationDetailedResult = IntegrationResult;
	/******************************************************************************/
	/*****             Integration Detail Execution Helper                    *****/
	/******************************************************************************/

	namespace IntegrationDetail
	{
		struct Sum {
			Real sum = 0, correction = 0;
			void Add(Real value) {
				const Real next = sum + value;
				correction += std::abs(sum) >= std::abs(value) ? (sum - next) + value : (value - next) + sum;
				sum = next;
			}
			Real Value() const { return sum + correction; }
		};

		struct ComplexSum {
			Sum real, imaginary;
			void Add(Complex value) { real.Add(value.real()); imaginary.Add(value.imag()); }
			Complex Value() const { return {real.Value(), imaginary.Value()}; }
		};

		inline Real Product(std::initializer_list<Real> factors, Real denominator = 1) {
			int exponent = 0;
			Real fraction = 1 / std::frexp(denominator, &exponent);
			exponent = -exponent;
			for (Real factor : factors) {
				int shift = 0;
				fraction *= std::frexp(factor, &shift);
				exponent += shift;
			}
			return std::ldexp(fraction, exponent);
		}

		inline Real HalfWidth(Real a, Real b) {
			ValidateBounds(a, b, "Integration");
			return std::signbit(a) != std::signbit(b) ? b / 2 - a / 2 : (b - a) / 2;
		}

		inline void RequireFinite(Real value, const char* context) {
			if (!std::isfinite(value))
				throw NumericalMethodError(std::string(context) + ": non-finite integral or error estimate");
		}

		inline void RequireFinite(Complex value, const char* context) {
			RequireFinite(value.real(), context);
			RequireFinite(value.imag(), context);
		}

		template<class Function, class... Args>
		Real Sample(Function& function, int& evaluations, bool checkFinite, Args&&... args) {
			if constexpr (requires { static_cast<bool>(function); })
				if (!function) throw NotImplementedError("Integration requires a callable");
			if (evaluations == std::numeric_limits<int>::max())
				throw ResourceLimitError("Integration evaluation counter exhausted");
			++evaluations;
			const Real value = static_cast<Real>(std::invoke(function, std::forward<Args>(args)...));
			if (checkFinite) RequireFinite(value, "Integration integrand");
			return value;
		}

		inline void CountEvaluation(std::atomic<int>& evaluations) {
			int used = evaluations.load(std::memory_order_relaxed);
			do {
				if (used == std::numeric_limits<int>::max())
					throw ResourceLimitError("Parallel integration evaluation counter exhausted");
			} while (!evaluations.compare_exchange_weak(used, used + 1, std::memory_order_relaxed));
		}

		inline void ValidateAdaptive(Real absolute, Real relative, int depth, int evaluations) {
			if (!std::isfinite(absolute) || !std::isfinite(relative) || absolute < 0 || relative < 0 ||
				(absolute == 0 && relative == 0) || depth < 0 || evaluations <= 0)
				throw ArgumentError("Adaptive integration requires nonnegative finite tolerances, a positive target, nonnegative depth and positive evaluation budget");
		}

		inline IntegrationResult CombineIntegrals(const IntegrationResult& left, const IntegrationResult& right, Real tolerance,
			bool checkTolerance = true) {
			if (left.function_evaluations > std::numeric_limits<int>::max() - right.function_evaluations ||
				left.iterations > std::numeric_limits<int>::max() - right.iterations)
				throw ResourceLimitError("Combined integration counters exceed supported range");
			IntegrationResult result(left.value + right.value, left.error_estimate + right.error_estimate,
				left.iterations + right.iterations, left.converged && right.converged,
				left.error_estimated && right.error_estimated);
			result.function_evaluations = left.function_evaluations + right.function_evaluations;
			RequireFinite(result.value, "Combined integration");
			if (result.error_estimated) RequireFinite(result.error_estimate, "Combined integration error");
			for (const auto* part : {&left, &right})
				if (!part->IsUsable()) {
					result.status = part->status;
					result.converged = false;
					result.error_message = part->error_message;
					return result;
				}
			const Real floor = Precision::IntegrationRoundoffFloorScale * std::abs(left.value) +
				Precision::IntegrationRoundoffFloorScale * std::abs(right.value);
			if (checkTolerance && result.converged && result.error_estimate > tolerance * std::abs(result.value) + floor) {
				result.converged = false;
				result.status = AlgorithmStatus::MaxIterationsExceeded;
				result.error_message = "Combined integration uncertainty exceeds the global tolerance";
			}
			return result;
		}

		template<size_t Dimension>
		std::vector<std::array<std::array<Real, 2>, Dimension>> SplitBox(
			const std::array<std::array<Real, 2>, Dimension>& box) {
			using Box = std::array<std::array<Real, 2>, Dimension>;
			std::array<Real, Dimension> midpoints;
			for (size_t i = 0; i < Dimension; ++i) {
				midpoints[i] = std::midpoint(box[i][0], box[i][1]);
				if (midpoints[i] == box[i][0] || midpoints[i] == box[i][1]) return {};
			}
			std::vector<Box> children(size_t(1) << Dimension);
			for (size_t mask = 0; mask < children.size(); ++mask)
				for (size_t axis = 0; axis < Dimension; ++axis)
					children[mask][axis] = (mask & (size_t(1) << axis))
						? std::array<Real, 2>{midpoints[axis], box[axis][1]}
						: std::array<Real, 2>{box[axis][0], midpoints[axis]};
			return children;
		}

		template<class Result, class Domain, class Evaluate, class Split>
		void RefineCells(Result& result, const std::vector<Domain>& domains, Evaluate&& evaluate, Split&& split,
			Real absolute, Real relative, int maxDepth, int maxEvaluations, int ruleCost) {
			ValidateAdaptive(absolute, relative, maxDepth, maxEvaluations);
			using Value = decltype(result.value);
			using ValueSum = std::conditional_t<std::is_same_v<Value, Complex>, ComplexSum, Sum>;
			struct Cell {
				Domain domain;
				BasicIntegrationResult<Value> estimate;
				int depth;
				bool operator<(const Cell& other) const { return estimate.error_estimate < other.estimate.error_estimate; }
			};
			std::priority_queue<Cell> pending;
			ValueSum values, terminalValues;
			Sum errors, terminalErrors;
			int deepest = 0, subdivisions = 0, stalled = 0;
			result.value = NoEstimate;
			result.error_estimate = NoEstimate;
			result.error_estimated = false;
			auto fail = [&](AlgorithmStatus status, const char* message) {
				result.status = status;
				result.converged = false;
				result.error_message = message;
			};
			auto publish = [&]() {
				result.value = values.Value();
				result.error_estimate = errors.Value();
				RequireFinite(result.value, "Adaptive integration");
				RequireFinite(result.error_estimate, "Adaptive integration");
				result.error_estimated = true;
				result.iterations = deepest;
				if constexpr (requires { result.cells_subdivided; result.max_depth_reached; }) {
					result.cells_subdivided = subdivisions;
					result.max_depth_reached = deepest;
				}
			};
			if (domains.size() > static_cast<size_t>(maxEvaluations / ruleCost)) {
				fail(AlgorithmStatus::ResourceLimitExceeded, "Evaluation budget cannot cover the initial integration domain");
				return;
			}
			for (const auto& domain : domains) {
				auto estimate = evaluate(domain);
				RequireFinite(estimate.value, "Adaptive cell");
				RequireFinite(estimate.error_estimate, "Adaptive cell");
				values.Add(estimate.value);
				errors.Add(estimate.error_estimate);
				if (maxDepth > 0) pending.push({domain, estimate, 0});
				else {
					terminalValues.Add(estimate.value);
					terminalErrors.Add(estimate.error_estimate);
				}
			}
			publish();
			for (;;) {
				Real target = std::max(absolute, relative * std::abs(result.value));
				if (result.error_estimate <= target) {
					// Re-sum the retained frontier before approval; removed parent errors must not cancel leaf uncertainty.
					auto frontier = pending;
					values = terminalValues;
					errors = terminalErrors;
					while (!frontier.empty()) {
						values.Add(frontier.top().estimate.value);
						errors.Add(frontier.top().estimate.error_estimate);
						frontier.pop();
					}
					publish();
					target = std::max(absolute, relative * std::abs(result.value));
					if (result.error_estimate <= target) {
						result.status = AlgorithmStatus::Success;
						result.converged = true;
						result.error_message.clear();
						return;
					}
				}
				if (pending.empty()) {
					fail(AlgorithmStatus::MaxIterationsExceeded, "Adaptive integration depth exhausted before the global tolerance was met");
					return;
				}
				const auto parent = pending.top();
				const auto children = split(parent.domain);
				if (children.empty()) {
					fail(AlgorithmStatus::ToleranceUnachievable, "Adaptive integration cannot subdivide floating-point coordinates further");
					return;
				}
				if (children.size() > static_cast<size_t>((maxEvaluations - result.function_evaluations) / ruleCost)) {
					fail(AlgorithmStatus::ResourceLimitExceeded, "Evaluation budget cannot complete another whole-domain refinement");
					return;
				}
				std::vector<Cell> replacements;
				for (const auto& child : children) {
					auto estimate = evaluate(child);
					RequireFinite(estimate.value, "Adaptive cell");
					RequireFinite(estimate.error_estimate, "Adaptive cell");
					replacements.push_back({child, estimate, parent.depth + 1});
				}
				const Real previousError = result.error_estimate;
				pending.pop();
				values.Add(-parent.estimate.value);
				errors.Add(-parent.estimate.error_estimate);
				for (const auto& child : replacements) {
					values.Add(child.estimate.value);
					errors.Add(child.estimate.error_estimate);
					deepest = std::max(deepest, child.depth);
					if (child.depth < maxDepth) pending.push(child);
					else {
						terminalValues.Add(child.estimate.value);
						terminalErrors.Add(child.estimate.error_estimate);
					}
				}
				++subdivisions;
				publish();
				stalled = result.error_estimate >= previousError ? stalled + 1 : 0;
				if (stalled >= Precision::IntegrationStagnationLimit &&
					result.error_estimate > std::max(absolute, relative * std::abs(result.value))) {
					fail(AlgorithmStatus::ToleranceUnachievable, "Adaptive integration error estimate stopped improving before the global tolerance was met");
					return;
				}
			}
		}

		/// Copy the integration payload and diagnostics of a core result into a Detailed result
		/// while keeping the wrapper's algorithm_name / elapsed_time_ms, and derive the status
		/// per RESULT_POLICY: converged -> Success, no estimate -> Completed, checked and missed ->
		/// MaxIterationsExceeded with an actionable message.
		template<typename TValue>
		inline void Adopt(BasicIntegrationResult<TValue>& out,
		                  const BasicIntegrationResult<TValue>& r,
		                  const char* algorithm_name)
		{
			out.value = r.value;
			out.iterations = r.iterations;
			out.converged = r.converged;
			out.error_estimated = r.error_estimated;
			out.error_estimate = r.error_estimated ? r.error_estimate : NoEstimate;
			out.error_message = r.error_message;
			if (r.function_evaluations > 0)
				out.function_evaluations = r.function_evaluations;
			if (!IsUsable(r.status)) {
				out.status = r.status;
				out.error_message = r.error_message;
				return;
			}
			if (r.status == AlgorithmStatus::Completed) {
				out.status = AlgorithmStatus::Completed;
			} else if (r.converged) {
				out.status = AlgorithmStatus::Success;
			} else if (!r.error_estimated) {
				out.status = AlgorithmStatus::Completed;
			} else {
				out.status = AlgorithmStatus::MaxIterationsExceeded;
				out.error_message = std::string(algorithm_name) + " did not converge after " +
					std::to_string(r.iterations) + " refinements (error estimate " +
					std::to_string(r.error_estimate) + ")";
			}
		}

		/// Execute an integration Detailed operation with timing and exception handling.
		///
		/// @tparam ResultType The detailed result type (e.g. IntegrationDetailedResult)
		/// @tparam ComputeFn Lambda that fills in result fields
		/// @param algorithm_name Name of the algorithm for diagnostics
		/// @param config Configuration controlling exception policy
		/// @param compute Lambda receiving ResultType& to populate
		template<typename ResultType, typename ComputeFn>
		ResultType ExecuteIntegrationDetailed(const char* algorithm_name,
		                                     const EvaluationConfigBase& config,
		                                     ComputeFn&& compute)
		{
			AlgorithmTimer timer;
			ResultType result = MakeEvaluationSuccessResult<ResultType>(algorithm_name);
			result.value = NoEstimate;
			auto execute = [&]() {
				compute(result);
				if (MML::IsUsable(result.status)) {
					if constexpr (requires { result.value.real(); result.value.imag(); }) {
						RequireFinite(result.value.real(), algorithm_name);
						RequireFinite(result.value.imag(), algorithm_name);
					} else RequireFinite(result.value, algorithm_name);
					if (result.error_estimated) RequireFinite(result.error_estimate, algorithm_name);
				}
			};
			auto finish = [&]() {
				if (!config.estimate_error) {
					result.error_estimated = false;
					result.error_estimate = NoEstimate;
				}
				result.elapsed_time_ms = timer.elapsed_ms();
				return result;
			};
			auto fail = [&](AlgorithmStatus status, const std::exception& error) {
				result.status = status;
				result.converged = false;
				result.error_message = error.what();
			};
			if (config.exception_policy == EvaluationExceptionPolicy::Propagate) {
				execute();
				return finish();
			}
			try {
				execute();
			}
			catch (const NumericalMethodError& ex) {
				fail(AlgorithmStatus::NumericalInstability, ex);
			}
			catch (const ResourceLimitError& ex) {
				fail(AlgorithmStatus::ResourceLimitExceeded, ex);
			}
			catch (const NumericInputError& ex) {
				fail(AlgorithmStatus::InvalidInput, ex);
			}
			catch (const std::invalid_argument& ex) {
				fail(AlgorithmStatus::InvalidInput, ex);
			}
			catch (const std::exception& ex) {
				fail(AlgorithmStatus::AlgorithmSpecificFailure, ex);
			}
			return finish();
		}
	} // namespace IntegrationDetail
}

#endif  // MML_INTEGRATION_BASE_H
