///////////////////////////////////////////////////////////////////////////////////////////
///                         MinimalMathLibrary (MML)                                  ///
///                                                                                   ///
///  File:        DAERadauIIA.h                                                       ///
///  Description: Radau IIA method for semi-explicit index-1 DAE systems              ///
///                                                                                   ///
///  3-stage implicit Runge-Kutta method with Radau nodes.                            ///
///  Gold standard for stiff DAEs - used in SUNDIALS IDA.                             ///
///                                                                                   ///
///  Properties:                                                                      ///
///  - Order 5 accuracy                                                               ///
///  - L-stable (strong damping of stiff components)                                  ///
///  - Excellent for index-1 and some index-2 problems                                ///
///                                                                                   ///
///  Reference: Hairer & Wanner, "Solving ODEs II", Section IV.8                      ///
///                                                                                   ///
///  Copyright:   (c) 2024-2026 Zvonimir Vanjak                                       ///
///  License:     MIT License (see LICENSE.md)                                         ///
///////////////////////////////////////////////////////////////////////////////////////////
#if !defined MML_DAE_RADAU_IIA_H
#define MML_DAE_RADAU_IIA_H

#include "DAESolverBase.h"
#include <array>
#include <cmath>

namespace MML {

	/// @brief Radau IIA - 3-stage implicit Runge-Kutta method for index-1 DAEs
	/// 
	/// Complexity: O((3N)³) = O(27N³) per step where N = diffDim + algDim.
	///            All 3 stages are coupled, requiring Newton on a 3N-dimensional system.
	///            Most expensive DAE solver per step, but highest order (5th).
	///
	/// Radau IIA is the gold standard for stiff DAEs:
	/// - L-stable (strong damping of stiff components)
	/// - Order 5 accuracy
	/// - Excellent for index-1 and some index-2 problems
	/// - Used in SUNDIALS IDA and other production solvers
	///
	/// The 3-stage Radau IIA has coefficients:
	///   c1 = (4 - sqrt(6))/10 ≈ 0.1550510257
	///   c2 = (4 + sqrt(6))/10 ≈ 0.6449489743
	///   c3 = 1
	///
	/// All stages are implicit and coupled, requiring Newton iteration
	/// on a 3*(diffDim+algDim) dimensional system.
	///
	/// Reference: Hairer & Wanner, "Solving ODEs II", Section IV.8
	///
	/// @param system DAE system with Jacobian
	/// @param t0 Initial time
	/// @param x0 Initial differential state
	/// @param y0 Initial algebraic state (must satisfy constraints)
	/// @param t_end Final time
	/// @param config Solver configuration
	/// @return DAESolverResult with solution and diagnostics
	inline DAESolverResult SolveDAERadauIIA(IODESystemDAEWithJacobian& system,
	                                         Real t0, const Vector<Real>& x0, const Vector<Real>& y0,
	                                         Real t_end, const DAESolverConfig& config = DAESolverConfig())
	{
		if (t_end <= t0)
			throw ArgumentError("SolveDAERadauIIA: t_end must be greater than t0 (reverse-time integration is not supported)");
		if (config.step_size <= 0)
			throw ArgumentError("SolveDAERadauIIA: config.step_size must be positive");

		AlgorithmTimer timer;

		int diffDim = system.getDiffDim();
		int algDim = system.getAlgDim();
		int totalDim = diffDim + algDim;
		int num_steps = static_cast<int>((t_end - t0) / config.step_size) + 1;

		DAESolverResult result(t0, t_end, diffDim, algDim, num_steps);
		result.algorithm_name = "DAERadauIIA";
		result.solution.fillValues(0, t0, x0, y0);
		if (!ValidateDAEInitialState(system, t0, x0, y0, config, result)) {
			result.solution.setFinalSize(0);
			result.elapsed_time_ms = timer.elapsed_ms();
			return result;
		}

		// Radau IIA coefficients (3-stage, order 5)
		const Real sqrt6 = std::sqrt(Real(6));
		
		// Abscissae (nodes)
		const Real c1 = (4 - sqrt6) / 10;  // ≈ 0.1550510257
		const Real c2 = (4 + sqrt6) / 10;  // ≈ 0.6449489743
		const Real c3 = 1;

		// Butcher tableau A matrix (for Radau IIA order 5)
		const Real a11 = (88 - 7*sqrt6) / 360;
		const Real a12 = (296 - 169*sqrt6) / 1800;
		const Real a13 = (-2 + 3*sqrt6) / 225;
		
		const Real a21 = (296 + 169*sqrt6) / 1800;
		const Real a22 = (88 + 7*sqrt6) / 360;
		const Real a23 = (-2 - 3*sqrt6) / 225;
		
		const Real a31 = (16 - sqrt6) / 36;
		const Real a32 = (16 + sqrt6) / 36;
		const Real a33 = 1.0 / 9;

		// Weights (b = last row of A for Radau methods)
		// b1 = a31, b2 = a32, b3 = a33

		// Current state
		Vector<Real> x = x0;
		Vector<Real> y = y0;
		Real t = t0;
		Real h = config.step_size;

		// Stage values: K_i for differential, L_i for algebraic
		// K_i = f(t + c_i*h, X_i, Y_i), L_i = y-values at stage
		Vector<Real> X1(diffDim), X2(diffDim), X3(diffDim);  // Stage x-values
		Vector<Real> Y1(algDim), Y2(algDim), Y3(algDim);     // Stage y-values
		Vector<Real> K1(diffDim), K2(diffDim), K3(diffDim);  // f evaluations
		Vector<Real> G1(algDim), G2(algDim), G3(algDim);     // g evaluations

		std::array<Matrix<Real>, 3> dfDx = {
			Matrix<Real>(diffDim, diffDim), Matrix<Real>(diffDim, diffDim), Matrix<Real>(diffDim, diffDim)};
		std::array<Matrix<Real>, 3> dfDy = {
			Matrix<Real>(diffDim, algDim), Matrix<Real>(diffDim, algDim), Matrix<Real>(diffDim, algDim)};
		std::array<Matrix<Real>, 3> dgDx = {
			Matrix<Real>(algDim, diffDim), Matrix<Real>(algDim, diffDim), Matrix<Real>(algDim, diffDim)};
		std::array<Matrix<Real>, 3> dgDy = {
			Matrix<Real>(algDim, algDim), Matrix<Real>(algDim, algDim), Matrix<Real>(algDim, algDim)};

		int step = 1;
		while (t + h/2 <= t_end && step < config.max_steps) {

			// Clamp step size to reach t_end exactly
			Real h_actual = std::min(h, t_end - t);
			if (h_actual < DAETimeTolerance(t_end, h)) break;

			Vector<Real> stageX(3 * diffDim), stageY(3 * algDim);
			for (int stage = 0; stage < 3; ++stage) {
				for (int i = 0; i < diffDim; ++i) stageX[stage * diffDim + i] = x[i];
				for (int i = 0; i < algDim; ++i) stageY[stage * algDim + i] = y[i];
			}

			auto unpackStages = [&](const Vector<Real>& valuesX, const Vector<Real>& valuesY) {
				for (int i = 0; i < diffDim; ++i) {
					X1[i] = valuesX[i]; X2[i] = valuesX[diffDim + i]; X3[i] = valuesX[2 * diffDim + i];
				}
				for (int i = 0; i < algDim; ++i) {
					Y1[i] = valuesY[i]; Y2[i] = valuesY[algDim + i]; Y3[i] = valuesY[2 * algDim + i];
				}
			};
			auto evaluateResidual = [&](const Vector<Real>& valuesX, const Vector<Real>& valuesY, Vector<Real>& residual) {
				unpackStages(valuesX, valuesY);
				system.diffEqs(t + c1 * h_actual, X1, Y1, K1);
				system.diffEqs(t + c2 * h_actual, X2, Y2, K2);
				system.diffEqs(t + c3 * h_actual, X3, Y3, K3);
				system.algConstraints(t + c1 * h_actual, X1, Y1, G1);
				system.algConstraints(t + c2 * h_actual, X2, Y2, G2);
				system.algConstraints(t + c3 * h_actual, X3, Y3, G3);
				for (int i = 0; i < diffDim; ++i) {
					residual[i] = X1[i] - x[i] - h_actual * (a11 * K1[i] + a12 * K2[i] + a13 * K3[i]);
					residual[diffDim + i] = X2[i] - x[i] - h_actual * (a21 * K1[i] + a22 * K2[i] + a23 * K3[i]);
					residual[2 * diffDim + i] = X3[i] - x[i] - h_actual * (a31 * K1[i] + a32 * K2[i] + a33 * K3[i]);
				}
				const int algebraicOffset = 3 * diffDim;
				for (int i = 0; i < algDim; ++i) {
					residual[algebraicOffset + i] = G1[i];
					residual[algebraicOffset + algDim + i] = G2[i];
					residual[algebraicOffset + 2 * algDim + i] = G3[i];
				}
			};
			auto evaluateJacobian = [&](const Vector<Real>& valuesX, const Vector<Real>& valuesY, Matrix<Real>& jacobian) {
				unpackStages(valuesX, valuesY);
				const std::array<Real, 3> stageTimes = {t + c1 * h_actual, t + c2 * h_actual, t + c3 * h_actual};
				const std::array<Vector<Real>*, 3> stageXs = {&X1, &X2, &X3};
				const std::array<Vector<Real>*, 3> stageYs = {&Y1, &Y2, &Y3};
				for (int stage = 0; stage < 3; ++stage)
					system.allJacobians(stageTimes[stage], *stageXs[stage], *stageYs[stage],
						dfDx[stage], dfDy[stage], dgDx[stage], dgDy[stage]);

				for (int row = 0; row < 3 * totalDim; ++row)
					for (int col = 0; col < 3 * totalDim; ++col) jacobian(row, col) = REAL(0.0);
				const Real coefficients[3][3] = {{a11, a12, a13}, {a21, a22, a23}, {a31, a32, a33}};
				const int algebraicOffset = 3 * diffDim;
				for (int si = 0; si < 3; ++si) {
					for (int sj = 0; sj < 3; ++sj) {
						for (int i = 0; i < diffDim; ++i) {
							for (int j = 0; j < diffDim; ++j)
								jacobian(si * diffDim + i, sj * diffDim + j) =
									(si == sj && i == j ? REAL(1.0) : REAL(0.0))
									- h_actual * coefficients[si][sj] * dfDx[sj](i, j);
							for (int j = 0; j < algDim; ++j)
								jacobian(si * diffDim + i, algebraicOffset + sj * algDim + j) =
									-h_actual * coefficients[si][sj] * dfDy[sj](i, j);
						}
					}
					for (int i = 0; i < algDim; ++i) {
						for (int j = 0; j < diffDim; ++j)
							jacobian(algebraicOffset + si * algDim + i, si * diffDim + j) = dgDx[si](i, j);
						for (int j = 0; j < algDim; ++j)
							jacobian(algebraicOffset + si * algDim + i, algebraicOffset + si * algDim + j) = dgDy[si](i, j);
					}
				}
			};

			DAENewtonResult newton = SolveDAENewton(stageX, stageY, 3 * diffDim, 3 * algDim,
				config, evaluateResidual, evaluateJacobian);
			AccumulateDAENewtonDiagnostics(result, newton);
			if (!newton.converged) {
				SetDAENewtonFailure(result, newton, "Radau IIA stage Newton failure");
				break;
			}
			unpackStages(stageX, stageY);

			// Update solution using stage 3 (for Radau IIA, c3 = 1, so X3, Y3 is the solution at t+h_actual)
			x = X3;
			y = Y3;
			t += h_actual;

			// Track constraint violation
			Vector<Real> g_check(algDim);
			system.algConstraints(t, x, y, g_check);
			if (!ValidateAcceptedDAEState(x, y, g_check, config, result, "Radau IIA step"))
				break;

			result.solution.fillValues(step, t, x, y);
			result.solution.incrementSuccessfulSteps();
			result.accepted_steps++;
			++step;
		}

		// Finalize
		result.solution.setFinalSize(step - 1);
		result.total_steps = step - 1;
		if (t_end - t > DAETimeTolerance(t_end, config.step_size) && result.status == AlgorithmStatus::Success) {
			result.status = AlgorithmStatus::MaxIterationsExceeded;
			result.failure_reason = DAEFailureReason::MaxStepsExceeded;
			result.error_message = "Maximum Radau IIA step count reached before t_end";
		}
		result.elapsed_time_ms = timer.elapsed_ms();

		return result;
	}

} // namespace MML

#endif // MML_DAE_RADAU_IIA_H
