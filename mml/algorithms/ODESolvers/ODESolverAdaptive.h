///////////////////////////////////////////////////////////////////////////////////////////
///                         MinimalMathLibrary (MML)                                  ///
///                                                                                   ///
///  File:        ODESolverAdaptive.h                                             ///
///  Description: Production-ready adaptive ODE integration with dense output        ///
///               Features: FSAL optimization, Hermite interpolation, diagnostics    ///
///                                                                                   ///
///  Copyright:   (c) 2024-2026 Zvonimir Vanjak                                       ///
///  License:     MIT License (see LICENSE.md)                   ///
///                                                                                   ///
///////////////////////////////////////////////////////////////////////////////////////////
#if !defined MML_ODE_SOLVER_ADAPTIVE_H
#define MML_ODE_SOLVER_ADAPTIVE_H

#include <mml/MMLBase.h>
#include <algorithm>
#include <mml/base/AlgorithmTypes.h>
#include <mml/interfaces/IODESystem.h>
#include <mml/base/ODESystemSolution.h>
#include <mml/algorithms/ODESolvers/ODESteppers.h>

#include <optional>
#include <type_traits>

namespace MML {

	/******************************************************************************/
	/*****            ODE Integrator Configuration                            *****/
	/******************************************************************************/

	/// Configuration for ODE integrators.
	///
	/// Provides user control over step size bounds, tolerances, and integration limits.
	/// All fields have sensible defaults suitable for most problems.
	///
	/// @example
	/// ODEIntegratorConfig config;
	/// config.tolerance = 1e-12;        // Higher precision
	/// config.max_steps = 100000;       // Allow more steps for stiff problems
	/// config.initial_step_size = 0.0;  // Auto-estimate (recommended)
	/// auto result = integrator.integrate(x0, t0, tEnd, config);
	struct ODEIntegratorConfig {
		/// Initial step size (default: 0.0 = auto-estimate)
		/// Set to 0 to use Hairer's automatic step size estimation
		Real initial_step_size = 0.0;

		/// Minimum allowed step size (default: type-adaptive)
		/// Integration fails if step size falls below this
		Real min_step_size = Precision::ODEMinStepSize;

		/// Maximum allowed step size (default: infinity)
		/// Useful to limit step growth in smooth regions
		Real max_step_size = std::numeric_limits<Real>::infinity();

		/// Error tolerance for step acceptance (default: type-adaptive)
		/// Used for both absolute and relative error control
		Real tolerance = Precision::ODEDefaultTolerance;

		/// Maximum number of attempted steps, accepted plus rejected (default: 100000)
		/// Integration fails if this is exceeded
		int max_steps = 100000;

		/// Output interval for dense output (default: 0.0 = output at integration steps)
		/// If > 0, solutions are interpolated at fixed intervals
		Real output_interval = 0.0;

		/// Enable verbose output for debugging (default: false)
		bool verbose = false;

		/// Factory method for high-precision configuration
		static ODEIntegratorConfig HighPrecision() {
			ODEIntegratorConfig config;
			config.tolerance = Precision::ODEDefaultTolerance / Real(100);
			config.max_steps = 500000;
			return config;
		}

		/// Factory method for fast (lower precision) configuration
		static ODEIntegratorConfig Fast() {
			ODEIntegratorConfig config;
			config.tolerance = Precision::DefaultTolerance;
			config.max_steps = 10000;
			return config;
		}

		/// Factory method for stiff problems
		static ODEIntegratorConfig Stiff() {
			ODEIntegratorConfig config;
			config.tolerance = Precision::DefaultToleranceStrict;
			config.max_steps = 1000000;
			config.min_step_size = Precision::DivisionSafetyThreshold;
			return config;
		}
	};

	//===================================================================================
	//                           SolutionStatistics
	//===================================================================================
	/// @brief Diagnostics for adaptive integration with timing information
	struct SolutionStatistics {
		// === Step Counts ===
		int acceptedSteps = 0;	///< Number of accepted steps
		int rejectedSteps = 0;	///< Number of rejected steps
		int totalFuncEvals = 0; ///< Total function evaluations

		// === Step Size Statistics ===
		Real minStepSize = 0;		///< Smallest step size used
		Real maxStepSize = 0;		///< Largest step size used

		// === Diagnostics ===
		std::string algorithm_name = "DormandPrince5";  ///< Stepper algorithm used
		AlgorithmStatus status = AlgorithmStatus::Success;
		std::string error_message;  ///< Error description (empty on success)
		double elapsed_time_ms = 0.0;  ///< Total integration time
		int attemptedSteps = 0; ///< Trials started, including an incomplete trial if a callback fails

		void reset() {
			acceptedSteps = rejectedSteps = totalFuncEvals = attemptedSteps = 0;
			minStepSize = std::numeric_limits<Real>::max();
			maxStepSize = 0;
			status = AlgorithmStatus::Success;
			error_message.clear();
			elapsed_time_ms = 0.0;
		}

		void recordStep(const StepResult& result) {
			if (result.attemptedSteps < 1 || result.funcEvals < 0)
				throw ODESolverError("Invalid adaptive ODE work count");
			if (result.attemptedSteps > std::numeric_limits<int>::max() - attemptedSteps ||
				result.funcEvals > std::numeric_limits<int>::max() - totalFuncEvals)
				throw ResourceLimitError("Adaptive ODE work counter exhausted");
			if (result.accepted) {
				acceptedSteps++;
				minStepSize = std::min(minStepSize, std::abs(result.hDone));
				maxStepSize = std::max(maxStepSize, std::abs(result.hDone));
			}
			attemptedSteps += result.attemptedSteps;
			rejectedSteps += result.rejectedSteps();
			totalFuncEvals += result.funcEvals;
		}

		Real acceptanceRate() const {
			int total = acceptedSteps + rejectedSteps;
			return total > 0 ? static_cast<Real>(acceptedSteps) / total : MML::Real(1.0);
		}

		int totalSteps() const {
			return acceptedSteps + rejectedSteps;
		}
	};

	struct ODEAdaptiveConfig : public EvaluationConfigBase {
		/// Embedded integrator configuration (step sizes, tolerance, etc.)
		ODEIntegratorConfig integrator;
	};
	struct ODEAdaptiveResult;
	// Function-template defaults must appear on the first declaration.
	template<typename Stepper = DormandPrince5_Stepper>
	ODEAdaptiveResult ODEAdaptiveIntegrateDetailed(const IODESystem& sys,
		const Vector<Real>& x0, Real t0, Real tEnd, const ODEAdaptiveConfig& config = {});

	//===================================================================================
	//                           ODEAdaptiveIntegrator
	//===================================================================================
	/// @brief Adaptive ODE integrator with dense output support
	///
	/// This integrator uses an adaptive stepper to efficiently solve ODEs while
	/// providing solution output at user-specified fixed intervals (dense output).
	///
	/// Key features:
	/// - Adaptive step size control for efficiency and accuracy
	/// - Dense output: solution at arbitrary times without extra function evals
	/// - FSAL optimization when supported by stepper
	/// - Detailed statistics for performance analysis
	template<typename Stepper = DormandPrince5_Stepper>
	class ODEAdaptiveIntegrator {
		template<typename S>
		friend ODEAdaptiveResult ODEAdaptiveIntegrateDetailed(const IODESystem&,
			const Vector<Real>&, Real, Real, const ODEAdaptiveConfig&);

	protected:
		const IODESystem& _sys;
		Stepper _stepper;
		SolutionStatistics _stats;

	public:
		/// @brief Construct integrator for given ODE system
		explicit ODEAdaptiveIntegrator(const IODESystem& sys)
				: _sys(sys)
				, _stepper(sys) {}

		/// @brief Get integration statistics from last solve
		const SolutionStatistics& getStatistics() const { return _stats; }

		/// @brief Estimate initial step size
		/// Uses Hairer's formula based on local error estimation
		Real estimateInitialStep(Real t0, const Vector<Real>& x0, Real tEnd, Real eps) {
			_stats.reset();
			return EstimateInitialStepBounded(t0, x0, tEnd, eps, std::numeric_limits<Real>::infinity());
		}

	protected:
		struct StepWorkRecorder {
			SolutionStatistics& stats;
			Stepper& stepper;
			const StepResult& result;
			ODESystemSolution& solution;
			int previousEvaluations;
			int previousAttempts;
			int attemptLimit;
			bool returned = false;
			bool recorded = false;
			StepWorkRecorder(SolutionStatistics& statistics, Stepper& method, const StepResult& outcome,
				ODESystemSolution& trajectory, int maxSteps = 0)
				: stats(statistics), stepper(method), result(outcome), solution(trajectory),
				  previousEvaluations(statistics.totalFuncEvals), previousAttempts(statistics.attemptedSteps),
				  attemptLimit((maxSteps > 0 ? maxSteps : std::numeric_limits<int>::max()) - previousAttempts) {
				if (attemptLimit <= 0)
					throw ResourceLimitError("Adaptive ODE step attempt counter exhausted");
				if constexpr (std::is_base_of_v<IAdaptiveStepper, Stepper>) {
					stepper._lastFunctionEvaluations = stepper._lastAttemptedSteps = stepper._lastRejectedSteps = 0;
					stepper._functionEvaluationLimit = std::numeric_limits<int>::max() - previousEvaluations;
					stepper._stepAttemptLimit = attemptLimit;
				}
			}
			void recordReturn() {
				if (result.funcEvals < 0)
					throw ODESolverError("Negative adaptive ODE RHS evaluation count");
				if (result.funcEvals > std::numeric_limits<int>::max() - previousEvaluations)
					throw ResourceLimitError("Adaptive ODE RHS evaluation counter exhausted");
				if (result.attemptedSteps < 1)
					throw ODESolverError("Invalid adaptive ODE step attempt count");
				if (result.attemptedSteps > attemptLimit)
					throw ResourceLimitError("Adaptive ODE step attempt budget exceeded");
				returned = true;
			}
			void recordStep() {
				stats.recordStep(result);
				recorded = true;
				if (result.accepted) solution.incrementSuccessfulSteps();
				for (int i = 0; i < result.rejectedSteps(); ++i)
					solution.incrementRejectedSteps();
			}
			~StepWorkRecorder() {
				if (returned) {
					stats.totalFuncEvals = previousEvaluations + result.funcEvals;
					stats.attemptedSteps = previousAttempts + result.attemptedSteps;
				} else if constexpr (std::is_base_of_v<IAdaptiveStepper, Stepper>) {
					stats.totalFuncEvals = previousEvaluations + stepper.getLastFunctionEvaluations();
					stats.attemptedSteps = previousAttempts + stepper.getLastAttemptedSteps();
				}
				if constexpr (std::is_base_of_v<IAdaptiveStepper, Stepper>) {
					if (!recorded) {
						stats.rejectedSteps += stepper.getLastRejectedSteps();
						for (int i = 0; i < stepper.getLastRejectedSteps(); ++i)
							solution.incrementRejectedSteps();
					}
					stepper._functionEvaluationLimit = std::numeric_limits<int>::max();
					stepper._stepAttemptLimit = std::numeric_limits<int>::max();
				}
			}
		};

		Real EstimateInitialStepBounded(Real t0, const Vector<Real>& x0, Real tEnd, Real eps, Real maxStepSize) {
			ValidateIntegrationInputs(x0, t0, tEnd, eps, 0, 0, maxStepSize);
			int n = _sys.getDim();
			Vector<Real> f0(n), f1(n), x1(n);

			const ODEStepperDetail::CountingODESystem system(_sys, _stats.totalFuncEvals);
			ODEStepperDetail::EvaluateRHS(system, t0, x0, f0);

			// Estimate scale of problem
			Real d0 = 0, d1 = 0;
			for (int i = 0; i < n; i++) {
				Real scale = std::abs(x0[i]) + Precision::InitialStepScale;
				d0 = std::max(d0, std::abs(x0[i]) / scale);
				d1 = std::max(d1, std::abs(f0[i]) / scale);
			}
			ODEStepperDetail::RequireFinite(d1, "Non-finite ODE initial derivative scale");

			// Initial guess
			Real h0 = (d0 < Real(1e-5) || d1 < Real(1e-5)) ? Real(1e-6) : Real(0.01) * d0 / d1;   // tol-ok: Hairer-Norsett-Wanner initial-step heuristic (I.4, algorithm 4.14)
			h0 = std::min({h0, std::abs(tEnd - t0), maxStepSize});

			// Take one explicit Euler step
			for (int i = 0; i < n; i++)
				x1[i] = x0[i] + std::copysign(h0, tEnd - t0) * f0[i];
			ODEStepperDetail::EvaluateRHS(system, t0 + std::copysign(h0, tEnd - t0), x1, f1);

			// Estimate second derivative
			Real d2 = 0;
			for (int i = 0; i < n; i++) {
				Real scale = std::abs(x0[i]) + Precision::InitialStepScale;
				d2 = std::max(d2, std::abs(f1[i] - f0[i]) / scale / h0);
			}
			ODEStepperDetail::RequireFinite(d2, "Non-finite ODE initial curvature estimate");

			// Optimal step for 5th order method
			Real h1;
			if (std::max(d1, d2) <= Real(1e-15)) {   // tol-ok: Hairer initial-step heuristic
				h1 = std::max(Real(1e-6), h0 * Real(1e-3));   // tol-ok: Hairer initial-step heuristic
			} else {
				h1 = std::pow(Real(0.01) / std::max(d1, d2), Real(1.0) / Real(5.0));
			}

			return std::min({100 * h0, h1, std::abs(tEnd - t0), maxStepSize});
		}

	public:
		/// @brief Integrate ODE system with dense output at fixed intervals
		/// @param x0 Initial conditions
		/// @param t0 Initial time
		/// @param tEnd Final time
		/// @param outputInterval Fixed output spacing, or zero for accepted-step output
		/// @param eps Error tolerance (default 1e-10)
		/// @param h0 Initial step size (0 = auto-estimate)
		/// @param minStepSize Minimum allowed step size (default 1e-15)
		/// @return Solution with points at regular intervals
		ODESystemSolution integrate(const Vector<Real>& x0, Real t0, Real tEnd, Real outputInterval, Real eps = Precision::ODEDefaultTolerance, Real h0 = 0, Real minStepSize = Precision::ODEMinStepSize, int maxSteps = 0) {
			return IntegrateWithStoragePolicy(x0, t0, tEnd, outputInterval, eps, h0, minStepSize, maxSteps, false);
		}

	private:
		void ValidateIntegrationInputs(const Vector<Real>& x0, Real t0, Real tEnd, Real eps,
			Real initialStep, Real minStepSize, Real maxStepSize) const {
			ODEStepperDetail::RequireFinite(t0, "Non-finite ODE initial time");
			ODEStepperDetail::RequireFinite(tEnd, "Non-finite ODE final time");
			if (!std::isfinite(tEnd - t0) || tEnd == t0)
				throw ODESolverError("integrate: requires distinct finite endpoints");
			ODEStepperDetail::RequireVector(x0, _sys.getDim(), "Invalid ODE initial state");
			if (!std::isfinite(eps) || eps <= 0 || !std::isfinite(initialStep) || initialStep < 0 ||
				!std::isfinite(minStepSize) || minStepSize < 0 || std::isnan(maxStepSize) ||
				maxStepSize <= 0 || minStepSize > maxStepSize)
				throw ODESolverError("Invalid adaptive ODE tolerance or step bounds");
		}

		void ValidateStepResult(const StepResult& result, const Vector<Real>& state, const Vector<Real>& derivative,
			Real attemptedStep, Real direction) const {
			ODEStepperDetail::RequireFinite(result.errMax, "Non-finite adaptive ODE error ratio");
			ODEStepperDetail::RequireFinite(result.hDone, "Non-finite adaptive ODE completed step");
			ODEStepperDetail::RequireFinite(result.hNext, "Non-finite adaptive ODE next step");
			if (result.errMax < 0 || direction * result.hNext <= 0)
				throw ODESolverError("Invalid adaptive ODE step proposal");
			if (result.accepted) {
				if (result.errMax > 1 || direction * result.hDone <= 0 ||
					std::abs(result.hDone) > std::abs(attemptedStep))
					throw ODESolverError("Invalid adaptive ODE step acceptance");
				ODEStepperDetail::RequireVector(state, _sys.getDim(), "Invalid adaptive ODE accepted state");
				ODEStepperDetail::RequireVector(derivative, _sys.getDim(), "Invalid adaptive ODE accepted derivative");
			}
		}

		static Real ValidateTimeAdvance(Real time, Real step, Real limit, Real direction) {
			const Real next = time + step;
			if (!std::isfinite(next) || direction * (next - time) <= 0 || direction * (next - limit) > 0)
				throw ODESolverError("Adaptive step cannot advance within the requested floating-point time interval");
			return next;
		}

		ODESystemSolution IntegrateWithStoragePolicy(const Vector<Real>& x0, Real t0, Real tEnd,
			Real outputInterval, Real eps, Real h0, Real minStepSize, int maxSteps, bool retainPartial,
			Real maxStepSize = std::numeric_limits<Real>::infinity(),
			std::optional<ODESystemSolution>* failureSolution = nullptr) {
			_stats.reset();
			ValidateIntegrationInputs(x0, t0, tEnd, eps, h0, minStepSize, maxStepSize);
			if (!std::isfinite(outputInterval) || outputInterval < 0)
				throw ODESolverError("outputInterval must be finite and non-negative");

			int n = _sys.getDim();
			const Real direction = tEnd > t0 ? Real(1) : Real(-1);
			ODESystemSolution sol(t0, tEnd, n, 0);
			_stepper.resetFSAL();
			auto storePoint = [&](Real time, const Vector<Real>& state) {
				ODEStepperDetail::RequireFinite(time, "Non-finite adaptive ODE output time");
				ODEStepperDetail::RequireVector(state, n, "Invalid adaptive ODE output state");
				if (sol.size() > 0 && direction * (time - sol.getTValue(sol.size() - 1)) <= 0)
					throw ODESolverError("Adaptive ODE output times must advance strictly");
				sol.fillValues(sol.size(), time, state);
			};

			Vector<Real> x = x0;
			Vector<Real> dxdt(n);
			Real t = t0;
			storePoint(t, x);

			auto retainEndpoint = [&]() {
				if (sol.getTValue(sol.size() - 1) != t) {
					AlgorithmStatus storageStatus = AlgorithmStatus::Success;
					std::string storageMessage;
					if (detail::PrepareTrajectoryPoint(sol, sol.size(), storageStatus, storageMessage))
						sol.fillValues(sol.size(), t, x);
					else
						_stats.error_message += "; last accepted endpoint could not be stored: " + storageMessage;
				}
			};
			auto captureFailure = [&](AlgorithmStatus status, const char* message) {
				_stats.status = status;
				_stats.error_message = message;
				if (failureSolution) {
					retainEndpoint();
					failureSolution->emplace(std::move(sol));
				}
			};

			try {
				const bool stepOutput = outputInterval == 0;
				int outputIndex = 1;
				auto nextOutput = [&]() {
					// Anchor each sample to the origin instead of accumulating rounding drift.
					const Real next = outputInterval >= std::abs(tEnd - t0) / Real(outputIndex)
						? tEnd : std::fma(direction * outputInterval, Real(outputIndex), t0);
					const Real bounded = direction * (next - tEnd) >= 0 ? tEnd : next;
					if (!std::isfinite(bounded) || direction * (bounded - sol.getTValue(sol.size() - 1)) <= 0)
						throw ODESolverError("Output interval cannot advance floating-point time");
					return bounded;
				};
				Real tNextOutput = stepOutput ? tEnd : nextOutput();
				const ODEStepperDetail::CountingODESystem system(_sys, _stats.totalFuncEvals);
				ODEStepperDetail::EvaluateRHS(system, t, x, dxdt);
				Real h = h0 > 0 ? h0 : EstimateInitialStepBounded(t0, x0, tEnd, eps, maxStepSize);

				while (direction * (tEnd - t) > 0) {
					if (maxSteps > 0 && _stats.attemptedSteps >= maxSteps) {
						_stats.status = AlgorithmStatus::MaxIterationsExceeded;
						_stats.error_message = "Maximum steps (" + std::to_string(maxSteps) + ") exceeded";
						break;
					}
					h = direction * std::min({std::abs(h), maxStepSize, std::abs(tEnd - t)});
					ValidateTimeAdvance(t, h, tEnd, direction);

					// Keep the last accepted endpoint intact even if a trial mutates then fails.
					Vector<Real> trial = x, trialDerivative = dxdt;
					StepResult result;
					StepWorkRecorder work(_stats, _stepper, result, sol, maxSteps);
					result = _stepper.doStep(t, trial, trialDerivative, h, eps);
					work.recordReturn();
					ValidateStepResult(result, trial, trialDerivative, h, direction);
					if (result.accepted) {
						const Real tNew = ValidateTimeAdvance(t, result.hDone, tEnd, direction);
						x = std::move(trial);
						dxdt = std::move(trialDerivative);
						t = tNew;
					}
					work.recordStep();
					h = result.hNext;

					if (result.accepted) {
						if (stepOutput) storePoint(t, x);
						else {
							while (direction * (tNextOutput - t) <= 0) {
								storePoint(tNextOutput, tNextOutput == t ? x : _stepper.interpolate(tNextOutput));
								if (tNextOutput == tEnd) break;
								++outputIndex;
								tNextOutput = nextOutput();
							}
						}
					}

					if (std::abs(h) < minStepSize && direction * (tEnd - t) > 0)
						throw ODESolverError("Step size too small in ODEAdaptiveIntegrator");
				}

				if (sol.getTValue(sol.size() - 1) != t) storePoint(t, x);
				sol.setFinalSize(sol.size() - 1);
			}
			catch (const ResourceLimitError& error) {
				_stats.status = AlgorithmStatus::ResourceLimitExceeded;
				_stats.error_message = error.what();
				if (!retainPartial) throw;
				retainEndpoint();
			}
			catch (const ODESolverError& error) {
				captureFailure(AlgorithmStatus::NumericalInstability, error.what());
				throw;
			}
			catch (const std::exception& error) {
				captureFailure(AlgorithmStatus::AlgorithmSpecificFailure, error.what());
				throw;
			}
			return sol;
		}

	public:
		/// @brief Integrate and return solution at specified times
		/// @param x0 Initial conditions
		/// @param times Strictly monotone output times (ascending or descending)
		/// @param eps Error tolerance
		/// @return Solution at specified times
		ODESystemSolution integrateAt(const Vector<Real>& x0, const Vector<Real>& times, Real eps = Precision::ODEDefaultTolerance, Real h0 = 0, Real minStepSize = Precision::ODEMinStepSize) {
			_stats.reset();
			if (times.size() < 2) {
				throw ODESolverError("Need at least 2 output times");
			}

			int n = _sys.getDim();
			int numPoints = static_cast<int>(times.size());
			Real t0 = times[0];
			Real tEnd = times[numPoints - 1];
			ValidateIntegrationInputs(x0, t0, tEnd, eps, h0, minStepSize, std::numeric_limits<Real>::infinity());
			const Real direction = tEnd > t0 ? Real(1) : Real(-1);
			for (int i = 1; i < numPoints; ++i)
				if (!std::isfinite(times[i - 1]) || !std::isfinite(times[i]) ||
					direction * (times[i] - times[i - 1]) <= 0)
					throw ODESolverError("Output times must be finite and strictly monotone");

			ODESystemSolution sol(t0, tEnd, n, 0);
			_stepper.resetFSAL();

			Vector<Real> x = x0;
			Vector<Real> dxdt(n);
			Real t = t0;

			const ODEStepperDetail::CountingODESystem system(_sys, _stats.totalFuncEvals);
			ODEStepperDetail::EvaluateRHS(system, t, x, dxdt);

			Real h = (h0 > 0) ? h0 : EstimateInitialStepBounded(t0, x0, tEnd, eps, std::numeric_limits<Real>::infinity());
			h *= direction;

			// Store initial point
			int outputIdx = 0;
			sol.fillValues(outputIdx++, t0, x0);

			// Integration loop
			while (outputIdx < numPoints) {
				Real tTarget = times[outputIdx];

				h = direction * std::min(std::abs(h), std::abs(tTarget - t));
				ValidateTimeAdvance(t, h, tTarget, direction);

				Vector<Real> trial = x, trialDerivative = dxdt;
				StepResult result;
				StepWorkRecorder work(_stats, _stepper, result, sol);
				result = _stepper.doStep(t, trial, trialDerivative, h, eps);
				work.recordReturn();
				ValidateStepResult(result, trial, trialDerivative, h, direction);
				if (result.accepted) {
					Real tNew = ValidateTimeAdvance(t, result.hDone, tTarget, direction);
					x = std::move(trial);
					dxdt = std::move(trialDerivative);
					work.recordStep();

					// Output at any times within this step
					while (outputIdx < numPoints && direction * (times[outputIdx] - tNew) <= 0) {
						Real tOut = times[outputIdx];
						if (tOut == tNew) {
							sol.fillValues(outputIdx++, tOut, x);
						} else {
							Vector<Real> xInterp = _stepper.interpolate(tOut);
							ODEStepperDetail::RequireVector(xInterp, n, "Invalid adaptive ODE output state");
							sol.fillValues(outputIdx++, tOut, xInterp);
						}
					}

					t = tNew;
					h = result.hNext;
				} else {
					work.recordStep();
					h = result.hNext;
				}

				if (std::abs(h) < minStepSize && direction * (tEnd - t) > 0) {
					throw ODESolverError("Step size too small in ODEAdaptiveIntegrator");
				}
			}

			sol.setFinalSize(sol.size() - 1);
			return sol;
		}

		/// @brief Integrate ODE system using configuration object
		/// 
		/// @param x0 Initial conditions
		/// @param t0 Initial time
		/// @param tEnd Final time
		/// @param config Integration configuration
		/// @return Solution with timing and diagnostic information in statistics
		ODESystemSolution integrate(const Vector<Real>& x0, Real t0, Real tEnd, const ODEIntegratorConfig& config) {
			AlgorithmTimer timer;  // Starts automatically

			ODESystemSolution sol = IntegrateWithStoragePolicy(x0, t0, tEnd, config.output_interval,
				config.tolerance, config.initial_step_size, config.min_step_size, config.max_steps, true, config.max_step_size);

			// Populate diagnostic fields in statistics
			_stats.elapsed_time_ms = timer.elapsed_ms();

			return sol;
		}

};

	// Convenient type aliases
	using DormandPrince5Integrator = ODEAdaptiveIntegrator<DormandPrince5_Stepper>;
	using CashKarpIntegrator = ODEAdaptiveIntegrator<CashKarp_Stepper>;
	using DormandPrince8Integrator = ODEAdaptiveIntegrator<DormandPrince8_Stepper>;
	using BulirschStoerIntegrator = ODEAdaptiveIntegrator<BulirschStoer_Stepper>;
	using BulirschStoerRationalIntegrator = ODEAdaptiveIntegrator<BulirschStoerRational_Stepper>;
	using BulirschStoerBulirschSeqIntegrator = ODEAdaptiveIntegrator<BulirschStoerRational_Stepper>; // Alias with clearer name


	///////////////////////////////////////////////////////////////////////////////////////////
	// ODEAdaptiveResult - Result type for adaptive ODE detailed APIs
	///////////////////////////////////////////////////////////////////////////////////////////
	struct ODEAdaptiveResult : public EvaluationResultBase {
		std::optional<ODESystemSolution> solution;

		int accepted_steps = 0;
		int rejected_steps = 0;
		int total_func_evals = 0;
		Real min_step_size = 0.0;
		Real max_step_size = 0.0;
		int attempted_steps = 0;
	};

	///////////////////////////////////////////////////////////////////////////////////////////
	// ODEAdaptiveDetail - Internal helpers for Detailed API execution
	///////////////////////////////////////////////////////////////////////////////////////////
	namespace ODEAdaptiveDetail
	{
		template<typename ResultType, typename ComputeFn>
		ResultType ExecuteODEAdaptiveDetailed(const char* algorithm_name,
		                                      const EvaluationConfigBase& config,
		                                      ComputeFn&& compute)
		{
			AlgorithmTimer timer;
			ResultType result = MakeEvaluationSuccessResult<ResultType>(algorithm_name);
			auto execute = [&]() {
				compute(result);
			};

			if (config.exception_policy == EvaluationExceptionPolicy::Propagate)
				execute();
			else {
				try {
					execute();
				}
				catch (const ResourceLimitError& ex) {
					result.status = AlgorithmStatus::ResourceLimitExceeded;
					result.error_message = ex.what();
				}
				catch (const ODESolverError& ex) {
					result.status = AlgorithmStatus::NumericalInstability;
					result.error_message = ex.what();
				}
				catch (const std::exception& ex) {
					result.status = AlgorithmStatus::AlgorithmSpecificFailure;
					result.error_message = ex.what();
				}
			}
			result.elapsed_time_ms = timer.elapsed_ms();
			return result;
		}
	} // namespace ODEAdaptiveDetail

	///////////////////////////////////////////////////////////////////////////////////////////
	// Detailed free function API
	///////////////////////////////////////////////////////////////////////////////////////////

	/// Adaptive ODE integration with full instrumentation.
	/// Returns ODEAdaptiveResult with solution, step statistics, timing, and AlgorithmStatus.
	template<typename Stepper>
	ODEAdaptiveResult ODEAdaptiveIntegrateDetailed(
		const IODESystem& sys,
		const Vector<Real>& x0, Real t0, Real tEnd,
		const ODEAdaptiveConfig& config)
	{
		std::optional<ODEAdaptiveIntegrator<Stepper>> integrator;
		auto result = ODEAdaptiveDetail::ExecuteODEAdaptiveDetailed<ODEAdaptiveResult>(
			"ODEAdaptiveIntegrator", config,
			[&](ODEAdaptiveResult& output) {
				integrator.emplace(sys);
				const auto& options = config.integrator;
				output.solution = integrator->IntegrateWithStoragePolicy(x0, t0, tEnd, options.output_interval,
					options.tolerance, options.initial_step_size, options.min_step_size, options.max_steps,
					true, options.max_step_size, &output.solution);
				output.status = integrator->getStatistics().status;
				output.error_message = integrator->getStatistics().error_message;
			});
		if (!integrator) return result;
		const auto& stats = integrator->getStatistics();
		if (!stats.error_message.empty()) result.error_message = stats.error_message;
		result.accepted_steps = stats.acceptedSteps;
		result.rejected_steps = stats.rejectedSteps;
		result.attempted_steps = stats.attemptedSteps;
		result.total_func_evals = stats.totalFuncEvals;
		result.min_step_size = stats.acceptedSteps > 0 ? stats.minStepSize : Real(0);
		result.max_step_size = stats.maxStepSize;
		result.function_evaluations = stats.totalFuncEvals;
		return result;
	}

} // namespace MML

#endif // MML_ODE_SOLVER_ADAPTIVE_H
