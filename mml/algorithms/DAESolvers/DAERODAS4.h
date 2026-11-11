///////////////////////////////////////////////////////////////////////////////////////////
///                         MinimalMathLibrary (MML)                                  ///
///  RODAS4: six-stage, L-stable Rosenbrock 4(3) for semi-explicit index-1 DAEs        ///
///  Copyright: (c) 2024-2026 Zvonimir Vanjak. MIT License (see LICENSE.md).           ///
///////////////////////////////////////////////////////////////////////////////////////////
#if !defined MML_DAE_RODAS4_H
#define MML_DAE_RODAS4_H

#include "DAESolverBase.h"
#include <algorithm>
#include <array>

namespace MML {

	namespace DAERODASDetail {
		// Hairer-Wanner RODAS, ROCOE METH=1 (not RODAS4P):
		// https://www.unige.ch/~hairer/prog/stiff/rodas.f
		// (M/(h*gamma)-J) k_i = F_i + M/h * sum(C_ij*k_j) + h*D_i*F_t.
		inline constexpr Real Gamma = static_cast<MML::Real>(0.25);
		inline constexpr Real Times[6] = {0, static_cast<MML::Real>(0.386), static_cast<MML::Real>(0.21), static_cast<MML::Real>(0.63), 1, 1};
		inline constexpr Real D[6] = {static_cast<MML::Real>(0.25), static_cast<MML::Real>(-0.1043), static_cast<MML::Real>(0.1035), static_cast<MML::Real>(-0.03620000000000023), 0, 0};
		inline constexpr Real A[6][6] = {
			{},
			{static_cast<MML::Real>(1.544)},
			{static_cast<MML::Real>(0.9466785280815826), static_cast<MML::Real>(0.2557011698983284)},
			{static_cast<MML::Real>(3.314825187068521), static_cast<MML::Real>(2.896124015972201), static_cast<MML::Real>(0.9986419139977817)},
			{static_cast<MML::Real>(1.221224509226641), static_cast<MML::Real>(6.019134481288629), static_cast<MML::Real>(12.53708332932087), static_cast<MML::Real>(-0.6878860361058950)},
			{static_cast<MML::Real>(1.221224509226641), static_cast<MML::Real>(6.019134481288629), static_cast<MML::Real>(12.53708332932087), static_cast<MML::Real>(-0.6878860361058950), 1}};
		inline constexpr Real C[6][6] = {
			{},
			{static_cast<MML::Real>(-5.6688)},
			{static_cast<MML::Real>(-2.430093356833875), static_cast<MML::Real>(-0.2063599157091915)},
			{static_cast<MML::Real>(-0.1073529058151375), static_cast<MML::Real>(-9.594562251023355), static_cast<MML::Real>(-20.47028614809616)},
			{static_cast<MML::Real>(7.496443313967647), static_cast<MML::Real>(-10.24680431464352), static_cast<MML::Real>(-33.99990352819905), static_cast<MML::Real>(11.70890893206160)},
			{static_cast<MML::Real>(8.083246795921522), static_cast<MML::Real>(-7.981132988064893), static_cast<MML::Real>(-31.52159432874371), static_cast<MML::Real>(16.31930543123136), static_cast<MML::Real>(-6.058818238834054)}};

		inline void Evaluate(IODESystemDAEWithJacobian& system, Real t, const Vector<Real>& x, const Vector<Real>& y, Vector<Real>& values) {
			Vector<Real> f(x.size()), g(y.size());
			system.diffEqs(t, x, y, f);
			system.algConstraints(t, x, y, g);
			for (int i = 0; i < x.size(); ++i)
				values[i] = f[i];
			for (int i = 0; i < y.size(); ++i)
				values[x.size() + i] = g[i];
		}

		inline Vector<Real> TimeDerivative(IODESystemDAEWithJacobian& system, Real t, Real end, const Vector<Real>& x, const Vector<Real>& y,
																			 const Vector<Real>& base) {
			const Real h = end - t;
			const Real delta = std::copysign(std::min(Precision::CbrtEps * std::max(Real(1), std::abs(t)), std::abs(h) / 2), h);
			const Real first = t + delta;
			const Real second = h > 0 ? std::min(t + 2 * delta, end) : std::max(t + 2 * delta, end);
			Vector<Real> f1(base.size()), f2(base.size()), derivative(base.size());
			// One-sided, fixed-state differences stay inside the requested time interval.
			// With fewer than three representable clocks, use the endpoint secant instead.
			if (first == t || second == first) {
				Evaluate(system, end, x, y, f1);
				for (int i = 0; i < base.size(); ++i)
					derivative[i] = (f1[i] - base[i]) / h;
			} else {
				Evaluate(system, first, x, y, f1);
				Evaluate(system, second, x, y, f2);
				const Real a = first - t, b = second - t;
				for (int i = 0; i < base.size(); ++i)
					derivative[i] = ((b / a) * (f1[i] - base[i]) - (a / b) * (f2[i] - base[i])) / (b - a);
			}
			return derivative;
		}
	} // namespace DAERODASDetail

	/// Adaptive Hairer-Wanner RODAS4 for x'=f(t,x,y), 0=g(t,x,y).
	/// Uses one frozen Jacobian/LU and six linear solves per attempted step, no Newton
	/// iterations. The sixth increment estimates the embedded third-order error.
	/// Both differential and algebraic errors control the positive step magnitude.
	/// Non-autonomous partial time derivatives are computed by one-sided differences.
	/// The legacy SolveDAERODAS alias still selects the old linearly-implicit trapezoid.
	inline DAESolverResult SolveDAERODAS4(IODESystemDAEWithJacobian& system, Real t0, const Vector<Real>& x0, const Vector<Real>& y0,
																				Real t_end, const DAESolverConfig& config = DAESolverConfig()) {
		const int diffDim = system.getDiffDim(), algDim = system.getAlgDim();
		if (!std::isfinite(t0) || !std::isfinite(t_end) || !std::isfinite(t_end - t0) || t_end == t0 || !std::isfinite(config.step_size) ||
				config.step_size <= 0 || !std::isfinite(config.min_step_size) || config.min_step_size <= 0 ||
				!std::isfinite(config.max_step_size) || config.max_step_size < config.min_step_size || !std::isfinite(config.abs_tolerance) ||
				config.abs_tolerance <= 0 || !std::isfinite(config.rel_tolerance) || config.rel_tolerance < 0 ||
				!std::isfinite(config.constraint_tol) || config.constraint_tol <= 0 || !std::isfinite(config.step_safety) ||
				config.step_safety <= 0 || config.step_safety >= 1 || config.max_steps <= 0 || diffDim <= 0 || algDim < 0 || x0.size() != diffDim ||
				y0.size() != algDim)
			throw ArgumentError("SolveDAERODAS4: requires distinct finite endpoints, matching dimensions, "
													"positive steps/absolute and constraint tolerances/budget, nonnegative relative tolerance, "
													"ordered step bounds and 0 < step_safety < 1");

		AlgorithmTimer timer;
		DAESolverResult result(t0, t_end, diffDim, algDim, 0);
		result.algorithm_name = "DAERODAS4";
		result.solution.fillValues(0, t0, x0, y0);
		if (!ValidateDAEInitialState(system, t0, x0, y0, config, result)) {
			result.elapsed_time_ms = timer.elapsed_ms();
			return result;
		}

		using namespace DAERODASDetail;
		const int dimension = diffDim + algDim;
		const Real direction = t_end > t0 ? Real(1) : Real(-1);
		Real t = t0;
		Real h = direction * std::clamp(config.step_size, config.min_step_size, config.max_step_size);
		Vector<Real> x = x0, y = y0;
		Vector<Real> base(dimension), values(dimension), rhs(dimension), constraints(algDim);
		Vector<Real> trialX(diffDim), trialY(algDim);
		Matrix<Real> fx(diffDim, diffDim), fy(diffDim, algDim), gx(algDim, diffDim), gy(algDim, algDim);
		Matrix<Real> matrix(dimension, dimension);
		std::array<Vector<Real>, 6> stages;
		auto rejectAttempt = [&]() {
			++result.rejected_steps;
			result.solution.incrementRejectedSteps();
		};

		while (direction * (t_end - t) > 0 && result.accepted_steps + result.rejected_steps < config.max_steps) {
			if (!PrepareDAETrajectoryPoint(result, result.accepted_steps + 1)) break;
			h = direction * std::min(std::abs(h), std::abs(t_end - t));
			const Real nextTime = t + h;
			if (!ValidateDAETimeStep(t, nextTime, result, "RODAS4"))
				break;

			Evaluate(system, t, x, y, base);
			const Vector<Real> timeDerivative = TimeDerivative(system, t, nextTime, x, y, base);
			system.allJacobians(t, x, y, fx, fy, gx, gy);
			++result.jacobian_evaluations;
			AssembleDAEAugmentedJacobian(fx, fy, gx, gy, h * Gamma, matrix);
			bool finite = IsFiniteDAEVector(base) && IsFiniteDAEVector(timeDerivative);
			Detail::ImplicitNewtonMatrix iterationMatrix;
			try {
				iterationMatrix.Factor(matrix);
				++result.lu_factorizations;
			} catch (const SingularMatrixError& error) {
				result.status = AlgorithmStatus::SingularMatrix;
				result.failure_reason = DAEFailureReason::SingularJacobian;
				result.error_message = std::string("RODAS4 stage matrix: ") + error.what();
			} catch (const MatrixNumericalError& error) {
				result.status = AlgorithmStatus::NumericalInstability;
				result.failure_reason = DAEFailureReason::NonFiniteState;
				result.error_message = std::string("RODAS4 stage matrix: ") + error.what();
			}
			if (result.status != AlgorithmStatus::Success) {
				rejectAttempt();
				break;
			}
			for (int stage = 0; stage < 6 && finite; ++stage) {
				trialX = x;
				trialY = y;
				for (int previous = 0; previous < stage; ++previous) {
					for (int i = 0; i < diffDim; ++i)
						trialX[i] += A[stage][previous] * stages[previous][i];
					for (int i = 0; i < algDim; ++i)
						trialY[i] += A[stage][previous] * stages[previous][diffDim + i];
				}
				if (!IsFiniteDAEVector(trialX) || !IsFiniteDAEVector(trialY)) {
					finite = false;
					break;
				}
				if (stage == 0)
					values = base;
				else
					Evaluate(system, t + Times[stage] * h, trialX, trialY, values);
				// Row-scale M/(h*gamma)-J: differential rows by h*gamma,
				// algebraic rows by -1. M=diag(I,0), so C never multiplies algebraic rows.
				for (int i = 0; i < diffDim; ++i) {
					rhs[i] = h * Gamma * (values[i] + h * D[stage] * timeDerivative[i]);
					for (int previous = 0; previous < stage; ++previous)
						rhs[i] += Gamma * C[stage][previous] * stages[previous][i];
				}
				for (int i = diffDim; i < dimension; ++i)
					rhs[i] = -(values[i] + h * D[stage] * timeDerivative[i]);
				finite = IsFiniteDAEVector(values) && IsFiniteDAEVector(rhs);
				if (!finite)
					break;
				try {
					stages[stage] = iterationMatrix.Solve(rhs);
					++result.linear_solves;
				} catch (const SingularMatrixError& error) {
					result.status = AlgorithmStatus::SingularMatrix;
					result.failure_reason = DAEFailureReason::SingularJacobian;
					result.error_message = std::string("RODAS4 stage solve: ") + error.what();
					break;
				} catch (const MatrixNumericalError& error) {
					result.status = AlgorithmStatus::NumericalInstability;
					result.failure_reason = DAEFailureReason::NonFiniteState;
					result.error_message = std::string("RODAS4 stage solve: ") + error.what();
					break;
				}
				finite = IsFiniteDAEVector(stages[stage]);
			}
			if (result.status != AlgorithmStatus::Success) {
				rejectAttempt();
				break;
			}
			if (!finite) {
				result.status = AlgorithmStatus::NumericalInstability;
				result.failure_reason = DAEFailureReason::NonFiniteState;
				result.error_message = "RODAS4: non-finite state, right-hand side or stage";
				rejectAttempt();
				break;
			}
			for (int i = 0; i < diffDim; ++i)
				trialX[i] += stages[5][i];
			for (int i = 0; i < algDim; ++i)
				trialY[i] += stages[5][diffDim + i];
			system.algConstraints(nextTime, trialX, trialY, constraints);
			if (!IsFiniteDAEVector(trialX) || !IsFiniteDAEVector(trialY) || !IsFiniteDAEVector(constraints)) {
				result.status = AlgorithmStatus::NumericalInstability;
				result.failure_reason = DAEFailureReason::NonFiniteState;
				result.error_message = "RODAS4: non-finite candidate state or constraint";
				rejectAttempt();
				break;
			}

			Real errorRatio = 0;
			for (int i = 0; i < dimension; ++i) {
				const Real oldValue = i < diffDim ? x[i] : y[i - diffDim];
				const Real newValue = i < diffDim ? trialX[i] : trialY[i - diffDim];
				const Real scale = config.abs_tolerance + config.rel_tolerance * std::max(std::abs(oldValue), std::abs(newValue));
				if (!std::isfinite(scale)) {
					errorRatio = std::numeric_limits<Real>::infinity();
					break;
				}
				errorRatio = std::max(errorRatio, std::abs(stages[5][i]) / scale);
			}
			result.last_error_norm = errorRatio;
			result.max_error_norm = std::max(result.max_error_norm, errorRatio);
			const Real constraintNorm = constraints.NormL2();
			if (!std::isfinite(errorRatio) || !std::isfinite(constraintNorm)) {
				result.status = AlgorithmStatus::NumericalInstability;
				result.failure_reason = DAEFailureReason::NonFiniteState;
				result.error_message = "RODAS4: non-finite error or constraint norm";
				rejectAttempt();
				break;
			}
			const bool constraintFailed = constraintNorm > config.constraint_tol;
			const bool accepted = errorRatio <= 1 && !constraintFailed;
			Real factor = errorRatio > 0 ? std::clamp(config.step_safety / std::pow(errorRatio, static_cast<MML::Real>(0.25)), Precision::RODASMinStepFactor,
																								Precision::RODASMaxStepFactor)
																	 : Precision::RODASMaxStepFactor;
			if (!accepted) {
				rejectAttempt();
				if (std::abs(h) <= config.min_step_size) {
					result.status = AlgorithmStatus::NumericalInstability;
					result.failure_reason = constraintFailed ? DAEFailureReason::ConstraintViolation : DAEFailureReason::StepSizeUnderflow;
					result.error_message = constraintFailed ? "RODAS4: constraint tolerance unattainable at minimum step"
																									: "RODAS4: embedded error tolerance unattainable at minimum step";
					break;
				}
				factor = std::min(factor, Precision::RODASRejectMaxFactor);
			} else {
				x = trialX;
				y = trialY;
				t = nextTime;
				result.final_constraint_norm = constraintNorm;
				result.max_constraint_violation = std::max(result.max_constraint_violation, constraintNorm);
				result.solution.fillValues(++result.accepted_steps, t, x, y);
				result.solution.incrementSuccessfulSteps();
			}
			h = direction * std::clamp(std::abs(h) * factor, config.min_step_size, config.max_step_size);
		}
		result.total_steps = result.accepted_steps;
		if (direction * (t_end - t) > 0 && result.status == AlgorithmStatus::Success) {
			result.status = AlgorithmStatus::MaxIterationsExceeded;
			result.failure_reason = DAEFailureReason::MaxStepsExceeded;
			result.error_message = "RODAS4: attempted-step budget exhausted at t=" + std::to_string(t);
		}
		result.elapsed_time_ms = timer.elapsed_ms();
		return result;
	}
} // namespace MML

#endif // MML_DAE_RODAS4_H
