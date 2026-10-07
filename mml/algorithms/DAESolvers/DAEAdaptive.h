///////////////////////////////////////////////////////////////////////////////////////////
///                         MinimalMathLibrary (MML)                                  ///
///  Adaptive Backward Euler for index-1 DAEs using step doubling                    ///
///////////////////////////////////////////////////////////////////////////////////////////
#if !defined MML_DAE_ADAPTIVE_H
#define MML_DAE_ADAPTIVE_H

#include "DAESolverBase.h"

namespace MML {

	inline DAESolverResult SolveDAEBackwardEulerAdaptive(IODESystemDAEWithJacobian& system,
		Real t0, const Vector<Real>& x0, const Vector<Real>& y0, Real t_end,
		const DAESolverConfig& config = DAESolverConfig()) {
		if (t_end <= t0)
			throw ArgumentError("SolveDAEBackwardEulerAdaptive: t_end must be greater than t0 (reverse-time integration is not supported)");
		if (config.step_size <= 0)
			throw ArgumentError("SolveDAEBackwardEulerAdaptive: config.step_size must be positive");

		AlgorithmTimer timer;
		const int diffDim = system.getDiffDim();
		const int algDim = system.getAlgDim();
		DAESolverResult result(t0, t_end, diffDim, algDim, 128);
		result.algorithm_name = "DAEBackwardEulerAdaptive";
		Vector<Real> x = x0, y = y0, constraints(algDim);
		Real t = t0;
		Real h = std::clamp(config.step_size, config.min_step_size, config.max_step_size);
		result.solution.fillValues(0, t, x, y);
		system.algConstraints(t, x, y, constraints);
		result.final_constraint_norm = constraints.NormL2();
		result.max_constraint_violation = result.final_constraint_norm;
		if (!IsFiniteDAEVector(x) || !IsFiniteDAEVector(y) || !IsFiniteDAEVector(constraints)) {
			result.status = AlgorithmStatus::NumericalInstability;
			result.failure_reason = DAEFailureReason::NonFiniteState;
			result.error_message = "Non-finite initial DAE state";
			result.solution.setFinalSize(0);
			result.elapsed_time_ms = timer.elapsed_ms();
			return result;
		}
		if (result.final_constraint_norm > config.constraint_tol) {
			result.status = AlgorithmStatus::InvalidInput;
			result.failure_reason = DAEFailureReason::InconsistentInitialConditions;
			result.error_message = "Initial algebraic constraints are inconsistent";
			result.solution.setFinalSize(0);
			result.elapsed_time_ms = timer.elapsed_ms();
			return result;
		}

		auto takeStep = [&](Real startTime, const Vector<Real>& startX, const Vector<Real>& startY,
			Real stepSize, Vector<Real>& endX, Vector<Real>& endY) {
			Vector<Real> derivative(diffDim), residualValues(diffDim + algDim);
			Matrix<Real> df_dx(diffDim, diffDim), df_dy(diffDim, algDim);
			Matrix<Real> dg_dx(algDim, diffDim), dg_dy(algDim, algDim);
			system.diffEqs(startTime, startX, startY, derivative);
			endX = startX + derivative * stepSize;
			endY = startY;
			const Real endTime = startTime + stepSize;
			auto residual = [&](const Vector<Real>& trialX, const Vector<Real>& trialY, Vector<Real>& values) {
				system.diffEqs(endTime, trialX, trialY, derivative);
				system.algConstraints(endTime, trialX, trialY, constraints);
				for (int i = 0; i < diffDim; ++i) values[i] = trialX[i] - startX[i] - stepSize * derivative[i];
				for (int i = 0; i < algDim; ++i) values[diffDim + i] = constraints[i];
			};
			auto jacobian = [&](const Vector<Real>& trialX, const Vector<Real>& trialY, Matrix<Real>& matrix) {
				system.allJacobians(endTime, trialX, trialY, df_dx, df_dy, dg_dx, dg_dy);
				AssembleDAEAugmentedJacobian(df_dx, df_dy, dg_dx, dg_dy, stepSize, matrix);
			};
			return SolveDAENewton(endX, endY, diffDim, algDim, config, residual, jacobian);
		};

		int savedStep = 1;
		while (t < t_end && result.accepted_steps + result.rejected_steps < config.max_steps) {
			h = std::min(h, t_end - t);
			Vector<Real> fullX, fullY, halfX, halfY, acceptedX, acceptedY;
			DAENewtonResult full = takeStep(t, x, y, h, fullX, fullY);
			DAENewtonResult firstHalf = takeStep(t, x, y, h * REAL(0.5), halfX, halfY);
			DAENewtonResult secondHalf;
			if (firstHalf.converged)
				secondHalf = takeStep(t + h * REAL(0.5), halfX, halfY, h * REAL(0.5), acceptedX, acceptedY);
			AccumulateDAENewtonDiagnostics(result, full);
			AccumulateDAENewtonDiagnostics(result, firstHalf);
			if (firstHalf.converged) AccumulateDAENewtonDiagnostics(result, secondHalf);

			if (!full.converged || !firstHalf.converged || !secondHalf.converged) {
				result.rejected_steps++;
				result.solution.incrementRejectedSteps();
				if (h <= config.min_step_size) {
					const DAENewtonResult& failure = !full.converged ? full : (!firstHalf.converged ? firstHalf : secondHalf);
					result.failure_reason = failure.failure_reason;
					result.status = failure.failure_reason == DAEFailureReason::SingularJacobian
						? AlgorithmStatus::SingularMatrix : AlgorithmStatus::NumericalInstability;
					result.error_message = failure.message;
					break;
				}
				h = std::max(config.min_step_size, h * REAL(0.5));
				continue;
			}

			Real errorRatio = REAL(0.0);
			for (int i = 0; i < diffDim; ++i) {
				Real scale = config.abs_tolerance + config.rel_tolerance * std::max(std::abs(acceptedX[i]), std::abs(x[i]));
				errorRatio = std::max(errorRatio, std::abs(acceptedX[i] - fullX[i]) / scale);
			}
			for (int i = 0; i < algDim; ++i) {
				Real scale = config.abs_tolerance + config.rel_tolerance * std::max(std::abs(acceptedY[i]), std::abs(y[i]));
				errorRatio = std::max(errorRatio, std::abs(acceptedY[i] - fullY[i]) / scale);
			}

			if (errorRatio > REAL(1.0) && h > config.min_step_size) {
				result.rejected_steps++;
				result.solution.incrementRejectedSteps();
				Real factor = std::max(REAL(0.2), config.step_safety / std::sqrt(errorRatio));
				h = std::max(config.min_step_size, h * factor);
				continue;
			}

			x = acceptedX;
			y = acceptedY;
			t += h;
			system.algConstraints(t, x, y, constraints);
			if (!ValidateAcceptedDAEState(x, y, constraints, config, result, "Adaptive Backward Euler step"))
				break;
			result.solution.fillValues(savedStep++, t, x, y);
			result.solution.incrementSuccessfulSteps();
			result.accepted_steps++;
			Real factor = errorRatio > REAL(0.0)
				? std::clamp(config.step_safety / std::sqrt(errorRatio), REAL(0.5), REAL(2.0)) : REAL(2.0);
			h = std::clamp(h * factor, config.min_step_size, config.max_step_size);
		}

		result.solution.setFinalSize(savedStep - 1);
		result.total_steps = result.accepted_steps;
		if (t < t_end && result.status == AlgorithmStatus::Success) {
			result.status = AlgorithmStatus::MaxIterationsExceeded;
			result.failure_reason = DAEFailureReason::MaxStepsExceeded;
			result.error_message = "Maximum adaptive DAE step count reached";
		}
		result.elapsed_time_ms = timer.elapsed_ms();
		return result;
	}
}

#endif