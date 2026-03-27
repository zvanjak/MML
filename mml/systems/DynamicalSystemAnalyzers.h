///////////////////////////////////////////////////////////////////////////////////////////
///                         MinimalMathLibrary (MML)                                  ///
///                                                                                   ///
///  File:        DynamicalSystemAnalyzers.h                                          ///
///  Description: Analysis tools for dynamical systems                                ///
///               Fixed points, Lyapunov exponents, bifurcations, phase space         ///
///                                                                                   ///
///  Copyright:   (c) 2024-2026 Zvonimir Vanjak                                       ///
///  License:     MIT License (see LICENSE.md)                                        ///
///                                                                                   ///
///////////////////////////////////////////////////////////////////////////////////////////
#if !defined MML_DYNAMICAL_SYSTEM_ANALYZERS_H
#define MML_DYNAMICAL_SYSTEM_ANALYZERS_H

#include "MMLBase.h"
#include "base/Vector/Vector.h"
#include "base/Matrix/Matrix.h"
#include "interfaces/IDynamicalSystem.h"
#include "systems/DynamicalSystemTypes.h"
#include "algorithms/EigenSystemSolvers.h"
#include "core/LinAlgEqSolvers.h"
#include "core/AlgorithmTypes.h"

#include <vector>
#include <cmath>
#include <algorithm>

namespace MML::Systems 
{
	//=============================================================================
	// FIXED POINT FINDER
	//=============================================================================

	/// @brief Find and classify fixed points of a dynamical system
	class FixedPointFinder {
	public:
		/// @brief Find a fixed point near initial guess using Newton's method
		/// @param sys The dynamical system
		/// @param initialGuess Starting point for Newton iteration
		/// @param tol Convergence tolerance
		/// @param maxIter Maximum iterations
		/// @return Fixed point result (check convergenceResidual for success)
		static FixedPoint<Real> Find(IDynamicalSystem& sys, const Vector<Real>& initialGuess, Real tol = Precision::DefaultToleranceStrict, int maxIter = 50) {
			int n = sys.getDim();
			FixedPoint<Real> result;
			result.location = initialGuess;
			result.jacobian.Resize(n, n);

			Vector<Real> x = initialGuess;
			Vector<Real> f(n), dx(n);
			Matrix<Real> J(n, n);

			for (int iter = 0; iter < maxIter; ++iter) {
				// Evaluate f(x) = derivs(x) - should be zero at fixed point
				sys.derivs(0.0, x, f);

				Real residual = f.NormL2();
				if (residual < tol) {
					result.location = x;
					result.convergenceResidual = residual;
					result.iterations = iter;

					// Compute Jacobian and classify
					sys.jacobian(0.0, x, result.jacobian);
					ClassifyFixedPoint(result);
					return result;
				}

				// Newton step: J * dx = -f
				sys.jacobian(0.0, x, J);

				// Solve using LU
				try {
					LUSolver<Real> solver(J);
					Vector<Real> negF = f * (-1.0);
					dx = solver.Solve(negF);
				} catch (...) {
					// Jacobian singular - can't continue
					result.type = FixedPointType::Unknown;
					result.convergenceResidual = residual;
					result.iterations = iter;
					return result;
				}

				// Update
				x = x + dx;
			}

			// Didn't converge
			sys.derivs(0.0, x, f);
			result.location = x;
			result.convergenceResidual = f.NormL2();
			result.iterations = maxIter;
			result.type = FixedPointType::Unknown;
			return result;
		}

		/// @brief Find multiple fixed points from several initial guesses
		static std::vector<FixedPoint<Real>> FindMultiple(IDynamicalSystem& sys, const std::vector<Vector<Real>>& initialGuesses,
																											Real tol = 1e-10, Real uniqueTol = 1e-6) {
			std::vector<FixedPoint<Real>> results;

			for (const auto& guess : initialGuesses) {
				auto fp = Find(sys, guess, tol);

				if (fp.convergenceResidual < tol) {
					// Check if this is a new fixed point
					bool isNew = true;
					for (const auto& existing : results) {
						if ((fp.location - existing.location).NormL2() < uniqueTol) {
							isNew = false;
							break;
						}
					}
					if (isNew)
						results.push_back(fp);
				}
			}

			return results;
		}

		/// @brief Classify fixed point based on eigenvalues
		static void ClassifyFixedPoint(FixedPoint<Real>& fp) {
			int n = fp.jacobian.rows();

			// Compute eigenvalues
			auto eigenResult = EigenSolver::Solve(fp.jacobian);

			fp.eigenvalues.clear();
			for (const auto& ev : eigenResult.eigenvalues)
				fp.eigenvalues.push_back(std::complex<Real>(ev.real, ev.imag));

			// Classify based on eigenvalues
			int numPositiveReal = 0;
			int numNegativeReal = 0;
			bool hasComplex = false;
			Real maxRealPart = -1e30;
			Real minRealPart = 1e30;

			for (const auto& ev : fp.eigenvalues) {
				Real re = ev.real();
				Real im = ev.imag();

				maxRealPart = std::max(maxRealPart, re);
				minRealPart = std::min(minRealPart, re);

				if (std::abs(im) > Precision::DefaultToleranceStrict)
					hasComplex = true;

				if (re > Precision::DefaultToleranceStrict)
					numPositiveReal++;
				else if (re < -Precision::DefaultToleranceStrict)
					numNegativeReal++;
			}

			// Determine stability
			fp.isStable = (maxRealPart < -Precision::DefaultToleranceStrict);

			// Classify type
			if (hasComplex) {
				if (std::abs(maxRealPart) < Precision::DefaultToleranceStrict && std::abs(minRealPart) < Precision::DefaultToleranceStrict)
					fp.type = FixedPointType::Center;
				else if (numPositiveReal > 0 && numNegativeReal > 0)
					fp.type = FixedPointType::SaddleFocus;
				else if (maxRealPart < -Precision::DefaultToleranceStrict)
					fp.type = FixedPointType::StableFocus;
				else
					fp.type = FixedPointType::UnstableFocus;
			} else {
				if (numPositiveReal > 0 && numNegativeReal > 0)
					fp.type = FixedPointType::Saddle;
				else if (numPositiveReal == 0 && numNegativeReal == n)
					fp.type = FixedPointType::StableNode;
				else if (numPositiveReal == n && numNegativeReal == 0)
					fp.type = FixedPointType::UnstableNode;
				else
					fp.type = FixedPointType::Unknown;
			}
		}
	};

	//=============================================================================
	// LYAPUNOV ANALYZER
	//=============================================================================

	/// @brief Compute Lyapunov exponents via variational equations
	///
	/// Uses the Benettin algorithm:
	/// 1. Integrate trajectory and n perturbation vectors simultaneously
	/// 2. Periodically orthonormalize using Gram-Schmidt
	/// 3. Accumulate stretching factors
	class LyapunovAnalyzer {
	public:
		/// @brief Compute all Lyapunov exponents
		/// @param sys The dynamical system
		/// @param x0 Initial condition
		/// @param tTotal Total integration time
		/// @param dtOrthonormalize Time between orthonormalizations
		/// @param h Integration step size
		/// @return LyapunovResult with exponents and related quantities
		static LyapunovResult<Real> Compute(IDynamicalSystem& sys, const Vector<Real>& x0, Real tTotal, Real dtOrthonormalize = 1.0,
																				Real h = 0.01) {
			int n = sys.getDim();
			LyapunovResult<Real> result;
			result.exponents.Resize(n);

			// State vector and perturbation frame
			Vector<Real> x = x0;
			Matrix<Real> Q = Matrix<Real>::Identity(n); // Orthonormal frame
			Vector<Real> lyapunovSums(n, 0.0);

			Real t = 0.0;
			int numOrth = 0;

			// Create RK4 stepper for trajectory
			while (t < tTotal) {
				// Integrate for dtOrthonormalize
				Real tEnd = std::min(t + dtOrthonormalize, tTotal);
				IntegrateWithVariational(sys, x, Q, t, tEnd, h);

				// Gram-Schmidt orthonormalization and accumulate stretching
				GramSchmidtQR(Q, lyapunovSums);

				numOrth++;
				t = tEnd;
			}

			// Compute final exponents
			for (int i = 0; i < n; ++i)
				result.exponents[i] = lyapunovSums[i] / tTotal;

			// Sort descending
			std::vector<Real> sorted(n);
			for (int i = 0; i < n; ++i)
				sorted[i] = result.exponents[i];
			std::sort(sorted.begin(), sorted.end(), std::greater<Real>());
			for (int i = 0; i < n; ++i)
				result.exponents[i] = sorted[i];

			result.maxExponent = result.exponents[0];
			result.sum = 0;
			for (int i = 0; i < n; ++i)
				result.sum += result.exponents[i];

			result.isChaotic = (result.maxExponent > Precision::DefaultTolerance);
			result.kaplanYorkeDimension = ComputeKaplanYorkeDimension(result.exponents);
			result.numOrthonormalizations = numOrth;
			result.totalTime = tTotal;

			return result;
		}

	private:
		/// @brief Integrate trajectory and variational equations
		static void IntegrateWithVariational(IDynamicalSystem& sys, Vector<Real>& x, Matrix<Real>& Q, Real t0, Real t1, Real h) {
			int n = sys.getDim();
			Vector<Real> k1(n), k2(n), k3(n), k4(n);
			Vector<Real> xTemp(n), xOld(n), xMid(n), xEnd(n);
			Matrix<Real> J(n, n);
			Matrix<Real> dQ1(n, n), dQ2(n, n), dQ3(n, n), dQ4(n, n);
			Matrix<Real> QTemp(n, n);

			Real t = t0;
			while (t < t1 - Precision::NumericalZeroThreshold) {
				Real dt = std::min(h, t1 - t);

				// Save x(t) before trajectory update
				for (int i = 0; i < n; ++i)
					xOld[i] = x[i];

				// RK4 for trajectory
				sys.derivs(t, x, k1);
				for (int i = 0; i < n; ++i)
					xTemp[i] = x[i] + 0.5 * dt * k1[i];
				// Save midpoint for Jacobian evaluation
				for (int i = 0; i < n; ++i)
					xMid[i] = xTemp[i];
				sys.derivs(t + 0.5 * dt, xTemp, k2);
				for (int i = 0; i < n; ++i)
					xTemp[i] = x[i] + 0.5 * dt * k2[i];
				sys.derivs(t + 0.5 * dt, xTemp, k3);
				for (int i = 0; i < n; ++i)
					xTemp[i] = x[i] + dt * k3[i];
				// Save endpoint for Jacobian evaluation
				for (int i = 0; i < n; ++i)
					xEnd[i] = xTemp[i];
				sys.derivs(t + dt, xTemp, k4);

				for (int i = 0; i < n; ++i)
					x[i] += dt * (k1[i] + 2 * k2[i] + 2 * k3[i] + k4[i]) / 6.0;

				// RK4 for variational equations: dQ/dt = J(x(t)) * Q
				// Jacobian evaluated at trajectory intermediate points (not post-update x)
				sys.jacobian(t, xOld, J);
				MultiplyJQ(J, Q, dQ1);

				for (int i = 0; i < n; ++i)
					for (int j = 0; j < n; ++j)
						QTemp(i, j) = Q(i, j) + 0.5 * dt * dQ1(i, j);
				sys.jacobian(t + 0.5 * dt, xMid, J);
				MultiplyJQ(J, QTemp, dQ2);

				for (int i = 0; i < n; ++i)
					for (int j = 0; j < n; ++j)
						QTemp(i, j) = Q(i, j) + 0.5 * dt * dQ2(i, j);
				MultiplyJQ(J, QTemp, dQ3);

				for (int i = 0; i < n; ++i)
					for (int j = 0; j < n; ++j)
						QTemp(i, j) = Q(i, j) + dt * dQ3(i, j);
				sys.jacobian(t + dt, xEnd, J);
				MultiplyJQ(J, QTemp, dQ4);

				for (int i = 0; i < n; ++i)
					for (int j = 0; j < n; ++j)
						Q(i, j) += dt * (dQ1(i, j) + 2 * dQ2(i, j) + 2 * dQ3(i, j) + dQ4(i, j)) / 6.0;

				t += dt;
			}
		}

		/// @brief Multiply J * Q (matrix-matrix)
		static void MultiplyJQ(const Matrix<Real>& J, const Matrix<Real>& Q, Matrix<Real>& result) {
			int n = J.rows();
			for (int i = 0; i < n; ++i) {
				for (int j = 0; j < n; ++j) {
					Real sum = 0;
					for (int k = 0; k < n; ++k)
						sum += J(i, k) * Q(k, j);
					result(i, j) = sum;
				}
			}
		}

		/// @brief Gram-Schmidt with accumulation of stretching factors
		static void GramSchmidtQR(Matrix<Real>& Q, Vector<Real>& lyapunovSums) {
			int n = Q.rows();

			for (int j = 0; j < n; ++j) {
				// Get column j
				Vector<Real> v(n);
				for (int i = 0; i < n; ++i)
					v[i] = Q(i, j);

				// Orthogonalize against previous columns
				for (int k = 0; k < j; ++k) {
					Real dot = 0;
					for (int i = 0; i < n; ++i)
						dot += v[i] * Q(i, k);
					for (int i = 0; i < n; ++i)
						v[i] -= dot * Q(i, k);
				}

				// Compute norm (stretching factor)
				Real norm = 0;
				for (int i = 0; i < n; ++i)
					norm += v[i] * v[i];
				norm = std::sqrt(norm);

				// Accumulate log of stretching
				if (norm > Precision::DivisionSafetyThreshold)
					lyapunovSums[j] += std::log(norm);

				// Normalize
				for (int i = 0; i < n; ++i)
					Q(i, j) = v[i] / norm;
			}
		}

		/// @brief Compute Kaplan-Yorke dimension
		///
		/// D_KY = j + sum(λ_1..λ_j) / |λ_{j+1}|
		/// where j is largest int such that sum(λ_1..λ_j) >= 0
		static Real ComputeKaplanYorkeDimension(const Vector<Real>& exponents) {
			int n = exponents.size();
			if (n == 0)
				return 0;

			Real cumSum = 0;
			int j = 0;

			for (int i = 0; i < n; ++i) {
				cumSum += exponents[i];
				if (cumSum >= 0)
					j = i + 1;
				else
					break;
			}

			if (j == n)
				return static_cast<Real>(n); // All positive
			if (j == 0)
				return 0; // All negative

			// Recompute partial sum up to j
			Real partialSum = 0;
			for (int i = 0; i < j; ++i)
				partialSum += exponents[i];

			return j + partialSum / std::abs(exponents[j]);
		}
	};

	//=============================================================================
	// BIFURCATION ANALYZER
	//=============================================================================

	/// @brief Generate bifurcation diagrams via parameter sweeps
	class BifurcationAnalyzer {
	public:
		/// @brief Sweep parameter and record attractor values
		/// @param sys The dynamical system (will be modified)
		/// @param paramIndex Which parameter to sweep
		/// @param paramMin Starting value
		/// @param paramMax Ending value
		/// @param numSteps Number of parameter values
		/// @param x0 Initial condition
		/// @param component Which state component to record
		/// @param tTransient Time to integrate before recording
		/// @param tRecord Time to record
		/// @param h Integration step size
		static BifurcationDiagram<Real> Sweep(IDynamicalSystem& sys, int paramIndex, Real paramMin, Real paramMax, int numSteps,
																					Vector<Real> x0, int component, Real tTransient = 100.0, Real tRecord = 50.0, Real h = 0.01) {
			BifurcationDiagram<Real> diagram;
			diagram.parameterName = sys.getParamName(paramIndex);
			Real originalParam = sys.getParam(paramIndex);			Real dParam = (paramMax - paramMin) / (numSteps - 1);

			for (int step = 0; step < numSteps; ++step) {
				Real param = paramMin + step * dParam;
				sys.setParam(paramIndex, param);

				// Integrate through transient
				Vector<Real> x = x0;
				IntegrateTo(sys, x, tTransient, h);

				// Record local maxima in component
				std::vector<Real> maxima = FindLocalMaxima(sys, x, component, tRecord, h);

				diagram.parameterValues.push_back(param);
				diagram.attractorValues.push_back(maxima);

				// Use final state as initial for next parameter (continuation)
				x0 = x;
			}

			sys.setParam(paramIndex, originalParam);
			return diagram;
		}

	private:
		/// @brief Simple RK4 integration
		static void IntegrateTo(IDynamicalSystem& sys, Vector<Real>& x, Real tTotal, Real h) {
			int n = sys.getDim();
			Vector<Real> k1(n), k2(n), k3(n), k4(n), xTemp(n);

			Real t = 0;
			while (t < tTotal) {
				Real dt = std::min(h, tTotal - t);

				sys.derivs(t, x, k1);
				for (int i = 0; i < n; ++i)
					xTemp[i] = x[i] + 0.5 * dt * k1[i];
				sys.derivs(t + 0.5 * dt, xTemp, k2);
				for (int i = 0; i < n; ++i)
					xTemp[i] = x[i] + 0.5 * dt * k2[i];
				sys.derivs(t + 0.5 * dt, xTemp, k3);
				for (int i = 0; i < n; ++i)
					xTemp[i] = x[i] + dt * k3[i];
				sys.derivs(t + dt, xTemp, k4);

				for (int i = 0; i < n; ++i)
					x[i] += dt * (k1[i] + 2 * k2[i] + 2 * k3[i] + k4[i]) / 6.0;

				t += dt;
			}
		}

		/// @brief Find local maxima of a component during integration
		static std::vector<Real> FindLocalMaxima(IDynamicalSystem& sys, Vector<Real>& x, int component, Real tTotal, Real h) {
			std::vector<Real> maxima;
			int n = sys.getDim();
			Vector<Real> k1(n), k2(n), k3(n), k4(n), xTemp(n);

			Real prevVal = x[component];
			Real prevPrevVal = prevVal;

			Real t = 0;
			while (t < tTotal) {
				Real dt = std::min(h, tTotal - t);

				sys.derivs(t, x, k1);
				for (int i = 0; i < n; ++i)
					xTemp[i] = x[i] + 0.5 * dt * k1[i];
				sys.derivs(t + 0.5 * dt, xTemp, k2);
				for (int i = 0; i < n; ++i)
					xTemp[i] = x[i] + 0.5 * dt * k2[i];
				sys.derivs(t + 0.5 * dt, xTemp, k3);
				for (int i = 0; i < n; ++i)
					xTemp[i] = x[i] + dt * k3[i];
				sys.derivs(t + dt, xTemp, k4);

				for (int i = 0; i < n; ++i)
					x[i] += dt * (k1[i] + 2 * k2[i] + 2 * k3[i] + k4[i]) / 6.0;

				Real val = x[component];

				// Check for local maximum
				if (prevVal > prevPrevVal && prevVal > val)
					maxima.push_back(prevVal);

				prevPrevVal = prevVal;
				prevVal = val;
				t += dt;
			}

			return maxima;
		}
	};

	//=============================================================================
	// PHASE SPACE ANALYZER
	//=============================================================================

	/// @brief Phase space analysis tools (Poincaré sections, trajectories)
	class PhaseSpaceAnalyzer {
	public:
		/// @brief Compute Poincaré section intersections
		static std::vector<Vector<Real>> ComputePoincareSection(IDynamicalSystem& sys, const Vector<Real>& x0,
																														const PoincareSection<Real>& section, int numIntersections, Real h = 0.01) {
			std::vector<Vector<Real>> intersections;
			int n = sys.getDim();

			Vector<Real> x = x0;
			Vector<Real> xPrev = x;
			Vector<Real> k1(n), k2(n), k3(n), k4(n), xTemp(n);

			Real t = 0;
			int count = 0;

			// Skip initial transient
			Real tTransient = 100.0;
			while (t < tTransient) {
				xPrev = x;
				RK4Step(sys, x, t, h, k1, k2, k3, k4, xTemp);
				t += h;
			}

			// Collect intersections
			while (count < numIntersections) {
				xPrev = x;
				Real valPrev = xPrev[section.variable] - section.value;

				RK4Step(sys, x, t, h, k1, k2, k3, k4, xTemp);
				t += h;

				Real val = x[section.variable] - section.value;

				// Check for crossing
				if (valPrev * val < 0) {
					// Check direction
					bool correctDirection =
						(section.direction == 0) || (section.direction > 0 && val > valPrev) || (section.direction < 0 && val < valPrev);

					if (correctDirection) {
						// Linear interpolation to find crossing point
						Real alpha = -valPrev / (val - valPrev);
						Vector<Real> xCross(n);
						for (int i = 0; i < n; ++i)
							xCross[i] = xPrev[i] + alpha * (x[i] - xPrev[i]);

						intersections.push_back(xCross);
						count++;
					}
				}
			}

			return intersections;
		}

		/// @brief Integrate trajectory and return points
		static std::vector<Vector<Real>> IntegrateTrajectory(IDynamicalSystem& sys, const Vector<Real>& x0, Real tTotal, Real dtOutput,
																												 Real h = 0.01) {
			std::vector<Vector<Real>> trajectory;
			int n = sys.getDim();

			Vector<Real> x = x0;
			Vector<Real> k1(n), k2(n), k3(n), k4(n), xTemp(n);

			Real t = 0;
			Real nextOutput = 0;

			trajectory.push_back(x);

			while (t < tTotal) {
				RK4Step(sys, x, t, h, k1, k2, k3, k4, xTemp);
				t += h;

				if (t >= nextOutput) {
					trajectory.push_back(x);
					nextOutput += dtOutput;
				}
			}

			return trajectory;
		}

	private:
		static void RK4Step(IDynamicalSystem& sys, Vector<Real>& x, Real t, Real h, Vector<Real>& k1, Vector<Real>& k2, Vector<Real>& k3,
												Vector<Real>& k4, Vector<Real>& xTemp) {
			int n = sys.getDim();

			sys.derivs(t, x, k1);
			for (int i = 0; i < n; ++i)
				xTemp[i] = x[i] + 0.5 * h * k1[i];
			sys.derivs(t + 0.5 * h, xTemp, k2);
			for (int i = 0; i < n; ++i)
				xTemp[i] = x[i] + 0.5 * h * k2[i];
			sys.derivs(t + 0.5 * h, xTemp, k3);
			for (int i = 0; i < n; ++i)
				xTemp[i] = x[i] + h * k3[i];
			sys.derivs(t + h, xTemp, k4);

			for (int i = 0; i < n; ++i)
				x[i] += h * (k1[i] + 2 * k2[i] + 2 * k3[i] + k4[i]) / 6.0;
		}
	};

	///////////////////////////////////////////////////////////////////////////////////////////
	// DynSysConfig - Configuration for dynamical system detailed APIs
	///////////////////////////////////////////////////////////////////////////////////////////
	struct DynSysConfig : public EvaluationConfigBase {
		// Inherits: estimate_error, check_finite, exception_policy
	};

	///////////////////////////////////////////////////////////////////////////////////////////
	// Result types for dynamical system detailed APIs
	///////////////////////////////////////////////////////////////////////////////////////////

	/// Result of fixed point finding with full instrumentation
	template<typename Type = Real>
	struct FindFixedPointResult : public EvaluationResultBase {
		/// The found fixed point (location, eigenvalues, stability, etc.)
		FixedPoint<Type> fixed_point;
	};

	/// Result of Lyapunov exponent computation with full instrumentation
	template<typename Type = Real>
	struct LyapunovAnalysisResult : public EvaluationResultBase {
		/// Lyapunov exponents and related quantities
		LyapunovResult<Type> lyapunov;
	};

	/// Result of bifurcation analysis with full instrumentation
	template<typename Type = Real>
	struct BifurcationAnalysisResult : public EvaluationResultBase {
		/// Bifurcation diagram data
		BifurcationDiagram<Type> diagram;
	};

	/// Result of Poincare section computation with full instrumentation
	template<typename Type = Real>
	struct PoincareSectionResult : public EvaluationResultBase {
		/// Section intersection points
		std::vector<Vector<Type>> intersections;
	};

	/// Result of trajectory integration with full instrumentation
	template<typename Type = Real>
	struct TrajectoryResult : public EvaluationResultBase {
		/// Trajectory points
		std::vector<Vector<Type>> points;
	};

	///////////////////////////////////////////////////////////////////////////////////////////
	// DynSysDetail - Internal helpers for Detailed API execution
	///////////////////////////////////////////////////////////////////////////////////////////
	namespace DynSysDetail
	{
		/// Execute a dynamical system Detailed operation with timing and exception handling.
		template<typename ResultType, typename ComputeFn>
		ResultType ExecuteDynSysDetailed(const char* algorithm_name,
		                                 const DynSysConfig& config,
		                                 ComputeFn&& compute)
		{
			auto execute = [&]() {
				AlgorithmTimer timer;
				ResultType result = MakeEvaluationSuccessResult<ResultType>(algorithm_name);
				compute(result);
				result.elapsed_time_ms = timer.elapsed_ms();
				return result;
			};

			if (config.exception_policy == EvaluationExceptionPolicy::Propagate)
				return execute();

			try {
				return execute();
			}
			catch (const std::exception& ex) {
				return MakeEvaluationFailureResult<ResultType>(
					AlgorithmStatus::AlgorithmSpecificFailure, ex.what(), algorithm_name);
			}
		}
	} // namespace DynSysDetail

	///////////////////////////    DETAILED API    ///////////////////////////

	/// Find a fixed point with full instrumentation.
	/// Maps convergence failure to AlgorithmStatus::MaxIterationsExceeded.
	inline FindFixedPointResult<Real> FindFixedPointDetailed(
		IDynamicalSystem& sys,
		const Vector<Real>& initialGuess,
		Real tol = Precision::DefaultToleranceStrict,
		int maxIter = 50,
		const DynSysConfig& config = {})
	{
		return DynSysDetail::ExecuteDynSysDetailed<FindFixedPointResult<Real>>(
			"FixedPointFinder", config,
			[&](FindFixedPointResult<Real>& result) {
				result.fixed_point = FixedPointFinder::Find(sys, initialGuess, tol, maxIter);
				result.function_evaluations = result.fixed_point.iterations;

				if (result.fixed_point.convergenceResidual >= tol) {
					result.status = AlgorithmStatus::MaxIterationsExceeded;
					result.error_message = "Newton iteration did not converge (residual="
						+ std::to_string(result.fixed_point.convergenceResidual) + ")";
				}
			});
	}

	/// Compute all Lyapunov exponents with full instrumentation.
	inline LyapunovAnalysisResult<Real> ComputeLyapunovDetailed(
		IDynamicalSystem& sys,
		const Vector<Real>& x0,
		Real tTotal,
		Real dtOrthonormalize = 1.0,
		Real h = 0.01,
		const DynSysConfig& config = {})
	{
		return DynSysDetail::ExecuteDynSysDetailed<LyapunovAnalysisResult<Real>>(
			"LyapunovAnalyzer", config,
			[&](LyapunovAnalysisResult<Real>& result) {
				result.lyapunov = LyapunovAnalyzer::Compute(sys, x0, tTotal, dtOrthonormalize, h);
				result.function_evaluations = result.lyapunov.numOrthonormalizations;
			});
	}

	/// Bifurcation parameter sweep with full instrumentation.
	inline BifurcationAnalysisResult<Real> SweepBifurcationDetailed(
		IDynamicalSystem& sys,
		int paramIndex,
		Real paramMin,
		Real paramMax,
		int numSteps,
		Vector<Real> x0,
		int component,
		Real tTransient = 100.0,
		Real tRecord = 50.0,
		Real h = 0.01,
		const DynSysConfig& config = {})
	{
		return DynSysDetail::ExecuteDynSysDetailed<BifurcationAnalysisResult<Real>>(
			"BifurcationAnalyzer", config,
			[&](BifurcationAnalysisResult<Real>& result) {
				result.diagram = BifurcationAnalyzer::Sweep(sys, paramIndex, paramMin, paramMax,
					numSteps, x0, component, tTransient, tRecord, h);
				result.function_evaluations = numSteps;
			});
	}

	/// Compute Poincare section with full instrumentation.
	inline PoincareSectionResult<Real> ComputePoincareSectionDetailed(
		IDynamicalSystem& sys,
		const Vector<Real>& x0,
		const PoincareSection<Real>& section,
		int numIntersections,
		Real h = 0.01,
		const DynSysConfig& config = {})
	{
		return DynSysDetail::ExecuteDynSysDetailed<PoincareSectionResult<Real>>(
			"PoincareSection", config,
			[&](PoincareSectionResult<Real>& result) {
				result.intersections = PhaseSpaceAnalyzer::ComputePoincareSection(
					sys, x0, section, numIntersections, h);
				result.function_evaluations = numIntersections;
			});
	}

	/// Integrate trajectory with full instrumentation.
	inline TrajectoryResult<Real> IntegrateTrajectoryDetailed(
		IDynamicalSystem& sys,
		const Vector<Real>& x0,
		Real tTotal,
		Real dtOutput,
		Real h = 0.01,
		const DynSysConfig& config = {})
	{
		return DynSysDetail::ExecuteDynSysDetailed<TrajectoryResult<Real>>(
			"TrajectoryIntegration", config,
			[&](TrajectoryResult<Real>& result) {
				result.points = PhaseSpaceAnalyzer::IntegrateTrajectory(sys, x0, tTotal, dtOutput, h);
				result.function_evaluations = static_cast<int>(result.points.size());
			});
	}

} // namespace MML::Systems
#endif // MML_DYNAMICAL_SYSTEM_ANALYZERS_H
