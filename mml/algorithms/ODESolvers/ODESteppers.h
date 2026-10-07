///////////////////////////////////////////////////////////////////////////////////////////
///                         MinimalMathLibrary (MML)                                  ///
///                                                                                   ///
///  File:        ODESteppers.h                                                 ///
///  Description: Adaptive ODE stepper implementations for use with integrators      ///
///               Includes: DormandPrince5, CashKarp, DormandPrince8, BulirschStoer   ///
///                                                                                   ///
///  Copyright:   (c) 2024-2026 Zvonimir Vanjak                                       ///
///  License:     MIT License (see LICENSE.md)                   ///
///                                                                                   ///
///////////////////////////////////////////////////////////////////////////////////////////
#if !defined MML_ODE_STEPPERS_H
#define MML_ODE_STEPPERS_H

#include <mml/MMLBase.h>
#include <mml/base/AlgorithmTypes.h>
#include <mml/interfaces/IODESystem.h>
#include "ODERKCoefficients.h"
#include "ODEStepperInfrastructure.h"

namespace MML {

	//===================================================================================
	//                              StepResult
	//===================================================================================
	/// @brief Result of an adaptive step attempt with diagnostic information
	struct StepResult {
		// === Step Outcome ===
		bool accepted = false; ///< Was the step accepted?
		Real hDone = 0.0;	   ///< Actual step size taken (or attempted if rejected)
		Real hNext = 0.0;	   ///< Suggested next step size

		// === Error Information ===
		Real errMax = 0.0;	   ///< Maximum error ratio (for diagnostics)
		int funcEvals = 0;	   ///< Number of function evaluations used

		// === Diagnostics (for single-step analysis) ===
		AlgorithmStatus status = AlgorithmStatus::Success;  ///< Step status
		std::string error_message;  ///< Error description (usually empty for steps)
	};

	//===================================================================================
	//                           IAdaptiveStepper
	//===================================================================================
	/// @brief Interface for adaptive steppers with dense output capability
	///
	/// An adaptive stepper takes a trial step and returns whether it was accepted.
	/// Key features:
	/// - Error estimation and step size control
	/// - Dense output (interpolation within accepted steps)
	/// - FSAL (First Same As Last) optimization support
	class IAdaptiveStepper {
	public:
		virtual ~IAdaptiveStepper() = default;

		/// @brief Attempt a step from current state
		/// @param t Current time
		/// @param x Current state (modified in place on accepted step)
		/// @param dxdt Current derivative (modified in place on accepted step)
		/// @param htry Trial step size
		/// @param eps Error tolerance
		/// @return StepResult with acceptance status and step info
		virtual StepResult doStep(Real t, Vector<Real>& x, Vector<Real>& dxdt, Real htry, Real eps) = 0;

		/// @brief Interpolate solution at time t within the last accepted step
		/// @param t Time to interpolate at (must be in [t_old, t_old + hDone])
		/// @return Interpolated state vector
		/// @pre Must be called after a successful doStep()
		virtual Vector<Real> interpolate(Real t) const = 0;

		/// @brief Check if stepper supports FSAL optimization
		/// @return true if the final derivative can be reused as initial derivative
		virtual bool isFSAL() const { return false; }

		/// @brief Get the final derivative from last step (for FSAL)
		/// @return Derivative at end of last step
		/// @pre Must be called after a successful doStep() when isFSAL() is true
		virtual const Vector<Real>& getFinalDeriv() const = 0;

		/// @brief Get the number of stages in this method
		virtual int stageCount() const = 0;

		/// @brief Get the order of the method
		virtual int order() const = 0;

		/// @brief Reset FSAL state (call when restarting integration)
		virtual void resetFSAL() = 0;
	};

	//===================================================================================
	//                         DormandPrince5_Stepper
	//===================================================================================
	/// @brief Dormand-Prince 5(4) adaptive stepper with FSAL and dense output
	///
	/// This is the workhorse method used by MATLAB's ode45 and SciPy's RK45.
	/// Features:
	/// - 5th order accurate solution with 4th order error estimate
	/// - FSAL: First Same As Last optimization (6 evals per step instead of 7)
	/// - 4th order dense output using Shampine's stage-based quartic extension
	/// - PI step size control for smooth adaptation
	class DormandPrince5_Stepper : public IAdaptiveStepper {
	private:
		const IODESystem& _sys;
		int _n; // System dimension
		using DP = RKCoeff::DormandPrince5;
		using StageEvaluator = ExplicitRKStageEvaluator<DP>;

		std::vector<Vector<Real>> _stages;
		Vector<Real> _workspace;
		Vector<Real> _xNew;
		Vector<Real> _error;
		StepSizeController _controller;
		DormandPrince5Interpolator _interpolator;

		// FSAL state
		bool _haveFSAL;

	public:
		explicit DormandPrince5_Stepper(const IODESystem& sys)
			: _sys(sys)
			, _n(sys.getDim())
			, _haveFSAL(false) {
			_stages.resize(DP::stages, Vector<Real>(_n));
			_workspace.Resize(_n);
			_xNew.Resize(_n);
			_error.Resize(_n);
		}

		StepResult doStep(Real t, Vector<Real>& x, Vector<Real>& dxdt, Real htry, Real eps) override {
			StepResult result;
			result.accepted = false;
			result.funcEvals = 0;

			Real h = htry;
			// Main step loop (retry with smaller h if rejected)
			while (true) {
				// FSAL reuses the previous accepted step's final stage as this step's k1.
				const Vector<Real>& initialDerivative = _haveFSAL ? _stages.back() : dxdt;
				if (!_haveFSAL)
					result.funcEvals++;
				StageEvaluator::evaluate(_sys, t, x, h, initialDerivative, _stages, _workspace);
				StageEvaluator::combineSolution(x, h, _stages, _xNew);
				StageEvaluator::combineError(h, _stages, _error);
				result.funcEvals += DP::stages - 1;

				// Error estimation
				Real errMax = 0.0;
				for (int i = 0; i < _n; i++) {
					Real scale = std::abs(x[i]) + std::abs(h * _stages[0][i]) + Precision::DivisionSafetyThreshold;
					errMax = std::max(errMax, std::abs(_error[i]) / scale);
				}
				errMax /= eps;
				result.errMax = errMax;

				if (errMax <= 1.0) {
					// Step accepted
					result.accepted = true;
					result.hDone = h;

					_interpolator.setStep(t, h, x, _xNew, _stages);
					x = _xNew;
					dxdt = _stages.back();
					_haveFSAL = true;
					result.hNext = _controller.acceptedStep(h, errMax);
					break;
				}

				h = _controller.rejectedStep(h, errMax);
				_haveFSAL = false;

				if (std::abs(h) < Constants::Eps) {
					throw ODESolverError("Step size underflow in DormandPrince5_Stepper");
				}
			}

			return result;
		}

		Vector<Real> interpolate(Real t) const override {
			return _interpolator.interpolate(t);
		}

		bool isFSAL() const override { return true; }

		const Vector<Real>& getFinalDeriv() const override { return _stages.back(); }

		int stageCount() const override { return 7; }

		int order() const override { return 5; }

		void resetFSAL() override {
			_haveFSAL = false;
			_controller.reset();
			_interpolator.reset();
		}
	};

	//===================================================================================
	//                           CashKarp_Stepper
	//===================================================================================
	/// @brief Cash-Karp 5(4) adaptive Runge-Kutta stepper
	///
	/// Classic adaptive RK method from Numerical Recipes.
	/// Features:
	/// - 5th order accurate solution with 4th order error estimate
	/// - Well-suited for general non-stiff ODEs
	/// - Efficient 6-stage method
	/// - Smooth step size adaptation
	///
	/// Note: Does not support FSAL optimization (simpler but less efficient than DP5).
	class CashKarp_Stepper : public IAdaptiveStepper {
	private:
		const IODESystem& _sys;
		int _n;
		using CK = RKCoeff::CashKarp5;
		using StageEvaluator = ExplicitRKStageEvaluator<CK>;

		std::vector<Vector<Real>> _stages;
		Vector<Real> _workspace;
		Vector<Real> _xNew;
		Vector<Real> _error;
		Vector<Real> _dxNew;
		StepSizeController _controller;
		HermiteInterpolator _interpolator;

	public:
		explicit CashKarp_Stepper(const IODESystem& sys)
			: _sys(sys)
			, _n(sys.getDim()) {
			_stages.resize(CK::stages, Vector<Real>(_n));
			_workspace.Resize(_n);
			_xNew.Resize(_n);
			_error.Resize(_n);
			_dxNew.Resize(_n);
		}

		StepResult doStep(Real t, Vector<Real>& x, Vector<Real>& dxdt, Real htry, Real eps) override {
			StepResult result;
			result.accepted = false;

			Real h = htry;

			while (true) {
				StageEvaluator::evaluate(_sys, t, x, h, dxdt, _stages, _workspace);
				StageEvaluator::combineSolution(x, h, _stages, _xNew);
				StageEvaluator::combineError(h, _stages, _error);
				result.funcEvals += CK::stages;

				// Error estimate
				Real errMax = 0.0;
				for (int i = 0; i < _n; i++) {
					Real scale = std::abs(x[i]) + std::abs(h * dxdt[i]) + Precision::DivisionSafetyThreshold;
					errMax = std::max(errMax, std::abs(_error[i]) / scale);
				}
				errMax /= eps;
				result.errMax = errMax;

				if (errMax <= 1.0) {
					result.accepted = true;
					result.hDone = h;

					// Cash-Karp is not FSAL; evaluate the endpoint derivative for Hermite output.
					_sys.derivs(t + h, _xNew, _dxNew);
					result.funcEvals++;
					_interpolator.setStep(t, h, x, _xNew, _stages[0], _dxNew);
					x = _xNew;
					dxdt = _dxNew;
					result.hNext = _controller.acceptedStep(h, errMax);

					break;
				}

				h = _controller.rejectedStep(h, errMax);

				if (std::abs(h) < Constants::Eps)
					throw ODESolverError("Step size underflow in CashKarp_Stepper");
			}

			return result;
		}

		Vector<Real> interpolate(Real t) const override {
			return _interpolator.interpolate(t);
		}

		bool isFSAL() const override { return false; }

		const Vector<Real>& getFinalDeriv() const override {
			return _dxNew;
		}

		int stageCount() const override { return 6; }

		int order() const override { return 5; }

		void resetFSAL() override {
			_controller.reset();
			_interpolator.reset();
		}
	};

	//===================================================================================
	//                         DormandPrince8_Stepper
	//===================================================================================
	/// @brief DOP853 8(5,3) high-order adaptive stepper
	///
	/// High-order method for problems requiring very high accuracy.
	/// Features:
	/// - 8th order solution with blended 5th and 3rd order error estimates
	/// - 12 propagation evaluations and one endpoint derivative per attempt
	/// - 7th order dense output using three interpolation-only stages
	/// - Excellent for smooth problems, astronomical trajectories
	/// - Higher cost per step but much larger accurate step sizes
	class DormandPrince8_Stepper : public IAdaptiveStepper {
	private:
		const IODESystem& _sys;
		int _n;
		using DP8 = RKCoeff::DormandPrince8;
		using StageEvaluator = ExplicitRKStageEvaluator<DP8>;

		std::vector<Vector<Real>> _k;
		Vector<Real> _workspace;
		Vector<Real> _xNew;
		StepSizeController _controller;
		DormandPrince8Interpolator _interpolator;

		static StepSizeControllerConfig controllerConfig() {
			// Retain DP8's conservative eighth-order exponent and historical growth cap.
			StepSizeControllerConfig config;
			config.alpha = Real(1.0 / 8.0);
			config.maxFactor = Real(6.0);
			return config;
		}

	public:
		explicit DormandPrince8_Stepper(const IODESystem& sys)
			: _sys(sys)
			, _n(sys.getDim())
			, _controller(controllerConfig()) {
			_k.resize(DP8::extended_stages);
			for (int i = 0; i < DP8::extended_stages; ++i)
				_k[i].Resize(_n);
			_workspace.Resize(_n);
			_xNew.Resize(_n);
		}

		StepResult doStep(Real t, Vector<Real>& x, Vector<Real>& dxdt, Real htry, Real eps) override {
			StepResult result;
			result.accepted = false;
			result.funcEvals = 0;

			Real h = htry;

			while (true) {
				StageEvaluator::evaluate(_sys, t, x, h, dxdt, _k, _workspace);
				StageEvaluator::combineSolution(x, h, _k, _xNew);
				_sys.derivs(t + h, _xNew, _k[12]);
				result.funcEvals += DP8::stages;

				Real error5NormSquared = 0;
				Real error3NormSquared = 0;
				for (int i = 0; i < _n; i++) {
					Real scale = std::abs(x[i]) + std::abs(h * _k[0][i]) + Precision::DivisionSafetyThreshold;
					Real error5 = 0;
					Real error3 = 0;
					for (int stage = 0; stage <= DP8::stages; ++stage) {
						error5 += DP8::error5Weight(stage) * _k[stage][i];
						error3 += DP8::error3Weight(stage) * _k[stage][i];
					}
					const Real normalizedError5 = error5 / scale;
					const Real normalizedError3 = error3 / scale;
					error5NormSquared += normalizedError5 * normalizedError5;
					error3NormSquared += normalizedError3 * normalizedError3;
				}
				const Real denominator = error5NormSquared + Real(0.01) * error3NormSquared;
				const Real errMax = denominator == 0 ? 0 :
					std::abs(h) * error5NormSquared / std::sqrt(denominator * _n) / eps;
				result.errMax = errMax;

				if (errMax <= 1.0) {
					result.accepted = true;
					result.hDone = h;

					result.funcEvals += _interpolator.setStep(_sys, t, h, x, _xNew, _k, _workspace);
					x = _xNew;
					dxdt = _k[12];
					result.hNext = _controller.acceptedStep(h, errMax);
					break;
				}

				h = _controller.rejectedStep(h, errMax);

				if (std::abs(h) < Constants::Eps)
					throw ODESolverError("Step size underflow in DormandPrince8_Stepper");
			}

			return result;
		}

		Vector<Real> interpolate(Real t) const override {
			return _interpolator.interpolate(t);
		}

		bool isFSAL() const override { return false; }

		const Vector<Real>& getFinalDeriv() const override { return _k[12]; }

		int stageCount() const override { return DP8::stages; }

		int order() const override { return 8; }

		void resetFSAL() override {
			_controller.reset();
			_interpolator.reset();
		}
	};

	//===================================================================================
	//                         BulirschStoer_Stepper
	//===================================================================================
	/// @brief Bulirsch-Stoer adaptive stepper with polynomial extrapolation
	///
	/// High-order adaptive method based on modified midpoint rule and extrapolation.
	/// Excellent for smooth problems requiring high accuracy.
	///
	/// Method:
	/// 1. Modified midpoint method with n substeps gives O(h²) error
	/// 2. Polynomial extrapolation to h→0 using Neville's algorithm
	/// 3. Adaptive convergence: keep adding sequence elements until error acceptable
	///
	/// References:
	/// - Hairer, Norsett, Wanner: "Solving Ordinary Differential Equations I"
	/// - Press et al.: "Numerical Recipes", Chapter 17
	class BulirschStoer_Stepper : public IAdaptiveStepper {
	private:
		const IODESystem& _sys;
		int _n; ///< System dimension

		// Substep sequence (even numbers for modified midpoint)
		static constexpr int KMAXX = 8; // Maximum columns in extrapolation
		int _nseq[KMAXX + 1] = {2, 4, 6, 8, 10, 12, 14, 16, 18};

		// Working arrays
		mutable Vector<Real> _xOld; ///< State at t
		mutable Vector<Real> _xNew; ///< State at t+h
		mutable Vector<Real> _err;	///< Error estimate
		mutable Vector<Real> _dydxOld;	 ///< Derivative at start of last accepted step
		mutable Vector<Real> _dydxFinal; ///< Derivative at end of last accepted step
		StepSizeController _controller;
		HermiteInterpolator _interpolator;

		// Extrapolation tableau - stores results at each level
		mutable std::vector<Vector<Real>> _d; // Differences for Neville
		mutable std::vector<Real> _xCoords;	  // x-coordinates for extrapolation

		/// @brief Modified midpoint method (Gragg's method)
		/// Computes y(x+H) using nstep substeps of size h = H/nstep
		void mmid(const Vector<Real>& y, const Vector<Real>& dydx, Real xs, Real htot, int nstep, Vector<Real>& yout) const {
			Real h = htot / nstep;
			Real h2 = 2.0 * h;

			Vector<Real> ym = y;
			Vector<Real> yn = y + dydx * h; // First step

			Real x = xs + h;
			Vector<Real> dyn(_n);
			_sys.derivs(x, yn, dyn);

			// General step
			for (int n = 1; n < nstep; n++) {
				Vector<Real> swap = ym + dyn * h2;
				ym = yn;
				yn = swap;
				x += h;
				_sys.derivs(x, yn, dyn);
			}

			// Last step - smoothing
			yout = (ym + yn + dyn * h) * 0.5;
		}

		/// @brief Polynomial extrapolation using Neville's algorithm
		/// Extrapolates sequence of estimates to step size h→0
		void pzextr(int iest, Real xest, const Vector<Real>& yest, Vector<Real>& yz, Vector<Real>& dy) const {
			// xest = (h/nseq[iest])^2 is the "x-coordinate" for extrapolation
			// yest is the estimate from modified midpoint with nseq[iest] steps
			// yz returns the extrapolated value, dy returns the error estimate

			_xCoords[iest] = xest;

			if (iest == 0) {
				// First estimate - just copy
				for (int j = 0; j < _n; j++) {
					yz[j] = yest[j];
					dy[j] = yest[j];
					_d[0][j] = yest[j];
				}
			} else {
				// Use Neville's algorithm
				Vector<Real> c = yest;

				for (int k = 0; k < iest; k++) {
					Real delta = 1.0 / (_xCoords[iest - k - 1] - xest);
					Real f1 = xest * delta;
					Real f2 = _xCoords[iest - k - 1] * delta;

					for (int j = 0; j < _n; j++) {
						Real q = _d[k][j];
						_d[k][j] = dy[j];
						delta = c[j] - q;
						dy[j] = f1 * delta;
						c[j] = f2 * delta;
					}
				}

				for (int j = 0; j < _n; j++) {
					yz[j] += dy[j];
					_d[iest][j] = dy[j];
				}
			}
		}

	public:
		explicit BulirschStoer_Stepper(const IODESystem& sys)
			: _sys(sys)
			, _n(sys.getDim())
			, _xOld(_n)
			, _xNew(_n)
			, _err(_n)
			, _dydxOld(_n)
			, _dydxFinal(_n) {
			// Initialize extrapolation tableau
			_d.resize(KMAXX + 1);
			for (int i = 0; i <= KMAXX; ++i) {
				_d[i] = Vector<Real>(_n);
			}
			_xCoords.resize(KMAXX + 1);
		}

		StepResult doStep(Real t, Vector<Real>& x, Vector<Real>& dxdt, Real htry, Real eps) override {
			StepResult result;
			result.accepted = false;
			result.funcEvals = 1; // We already have dxdt

			_xOld = x;
			_dydxOld = dxdt; // Store initial derivative for Hermite interpolation
			Real h = htry;

			Vector<Real> ysav = x;
			Vector<Real> yseq(_n), yest(_n), yerr(_n);

			// Try step with current h, building up extrapolation tableau
			for (int k = 0; k <= KMAXX; k++) {
				// Modified midpoint with _nseq[k] substeps
				mmid(ysav, dxdt, t, h, _nseq[k], yseq);
				result.funcEvals += _nseq[k];

				// The "x-coordinate" for extrapolation is (h/n)^2
				Real xest = (h / _nseq[k]) * (h / _nseq[k]);

				// Extrapolate
				pzextr(k, xest, yseq, yest, yerr);

				if (k > 0) { // Need at least 2 points to estimate error
					// Compute scaled error
					Real errMax = 0.0;
					for (int i = 0; i < _n; ++i) {
						Real scale = eps * (std::abs(ysav[i]) + std::abs(yest[i]) + Precision::DivisionSafetyThreshold);
						Real errRatio = std::abs(yerr[i]) / scale;
						errMax = std::max<Real>(errMax, errRatio);
					}

					result.errMax = errMax;

					if (errMax < 1.0) {
						// Converged!
						result.accepted = true;
						_xNew = yest;
						_err = yerr;

						_sys.derivs(t + h, yest, dxdt);
						_dydxFinal = dxdt;
						result.funcEvals++;
						result.hDone = h;
						_interpolator.setStep(t, h, x, yest, _dydxOld, _dydxFinal);
						x = yest;
						// Extrapolation order grows with the converged tableau column.
						result.hNext = _controller.acceptedStep(h, errMax, Real(1.0) / (2 * k + 1));

						return result;
					}
				}
			}

			// Failed to converge - reduce step size
			result.accepted = false;
			result.hDone = h;
			result.errMax = Real(999.0);
			result.hNext = _controller.rejectedStep(h, result.errMax, Real(1.0) / (2 * KMAXX + 1));

			return result;
		}

		Vector<Real> interpolate(Real t) const override {
			return _interpolator.interpolate(t);
		}

		bool isFSAL() const override { return false; }

		const Vector<Real>& getFinalDeriv() const override {
			return _dydxFinal;
		}

		int stageCount() const override { return KMAXX; }

		int order() const override { return 2 * KMAXX; }

		void resetFSAL() override {
			_controller.reset();
			_interpolator.reset();
		}
	};

	//===================================================================================
	//                      BulirschStoerRational_Stepper
	//===================================================================================
	/// @brief Bulirsch-Stoer stepper with Bulirsch's original step sequence
	///
	/// Variant using Bulirsch's recommended sequence {2,4,6,8,12,16,24,32,48}
	/// instead of the simple even sequence {2,4,6,8,10,12,14,16,18}.
	/// Uses polynomial extrapolation (rational extrapolation is numerically unstable).
	///
	/// The Bulirsch sequence can be more efficient for certain problems as it
	/// allows larger jumps in the extrapolation tableau.
	///
	/// References:
	/// - Stoer & Bulirsch: "Introduction to Numerical Analysis"
	/// - Press et al.: "Numerical Recipes", Chapter 17
	class BulirschStoerRational_Stepper : public IAdaptiveStepper {
	private:
		const IODESystem& _sys;
		int _n; ///< System dimension

		// Bulirsch's original step sequence (more aggressive growth)
		static constexpr int KMAXX = 8; // Maximum columns in extrapolation
		int _nseq[KMAXX + 1] = {2, 4, 6, 8, 12, 16, 24, 32, 48};

		// Working arrays
		mutable Vector<Real> _xOld; ///< State at t
		mutable Vector<Real> _xNew; ///< State at t+h
		mutable Vector<Real> _err;	///< Error estimate
		mutable Vector<Real> _dydxOld;	 ///< Derivative at start of last accepted step
		mutable Vector<Real> _dydxFinal; ///< Derivative at end of last accepted step
		StepSizeController _controller;
		HermiteInterpolator _interpolator;

		// Extrapolation tableau for rational extrapolation
		mutable std::vector<Vector<Real>> _d; // Differences for rational extrapolation
		mutable std::vector<Real> _xCoords;	  // x-coordinates for extrapolation

		/// @brief Modified midpoint method (Gragg's method)
		/// Computes y(x+H) using nstep substeps of size h = H/nstep
		void mmid(const Vector<Real>& y, const Vector<Real>& dydx, Real xs, Real htot, int nstep, Vector<Real>& yout) const {
			Real h = htot / nstep;
			Real h2 = 2.0 * h;

			Vector<Real> ym = y;
			Vector<Real> yn = y + dydx * h; // First step

			Real x = xs + h;
			Vector<Real> dyn(_n);
			_sys.derivs(x, yn, dyn);

			// General step
			for (int n = 1; n < nstep; n++) {
				Vector<Real> swap = ym + dyn * h2;
				ym = yn;
				yn = swap;
				x += h;
				_sys.derivs(x, yn, dyn);
			}

			// Last step - smoothing
			yout = (ym + yn + dyn * h) * 0.5;
		}

		/// @brief Polynomial extrapolation (same as pzextr but for rational stepper)
		/// Note: True rational extrapolation is numerically unstable in many cases.
		/// This uses polynomial extrapolation which is more robust.
		/// The difference from BulirschStoer_Stepper is the step sequence used.
		void rzextr(int iest, Real xest, const Vector<Real>& yest, Vector<Real>& yz, Vector<Real>& dy) const {
			// Store x-coordinate for this estimate
			_xCoords[iest] = xest;

			if (iest == 0) {
				// First point - just copy
				for (int j = 0; j < _n; j++) {
					yz[j] = yest[j];
					dy[j] = yest[j];
					_d[0][j] = yest[j];
				}
				return;
			}

			// Use Neville's polynomial extrapolation (same as pzextr)
			// This is more stable than rational extrapolation
			Vector<Real> c = yest;

			for (int k = 0; k < iest; k++) {
				Real delta = 1.0 / (_xCoords[iest - k - 1] - xest);
				Real f1 = xest * delta;
				Real f2 = _xCoords[iest - k - 1] * delta;

				for (int j = 0; j < _n; j++) {
					Real q = _d[k][j];
					_d[k][j] = dy[j];
					Real diff = c[j] - q;
					dy[j] = f1 * diff;
					c[j] = f2 * diff;
				}
			}

			for (int j = 0; j < _n; j++) {
				yz[j] += dy[j];
				_d[iest][j] = dy[j];
			}
		}

	public:
		explicit BulirschStoerRational_Stepper(const IODESystem& sys)
			: _sys(sys)
			, _n(sys.getDim())
			, _xOld(_n)
			, _xNew(_n)
			, _err(_n)
			, _dydxOld(_n)
			, _dydxFinal(_n) {
			// Initialize extrapolation tableau
			_d.resize(KMAXX + 1);
			for (int i = 0; i <= KMAXX; ++i) {
				_d[i] = Vector<Real>(_n);
			}
			_xCoords.resize(KMAXX + 1);
		}

		StepResult doStep(Real t, Vector<Real>& x, Vector<Real>& dxdt, Real htry, Real eps) override {
			StepResult result;
			result.accepted = false;
			result.funcEvals = 1; // We already have dxdt

			_xOld = x;
			_dydxOld = dxdt; // Store initial derivative for Hermite interpolation
			Real h = htry;

			Vector<Real> ysav = x;
			Vector<Real> yseq(_n), yest(_n), yerr(_n);

			// Try step with current h, building up extrapolation tableau
			for (int k = 0; k <= KMAXX; k++) {
				// Modified midpoint with _nseq[k] substeps
				mmid(ysav, dxdt, t, h, _nseq[k], yseq);
				result.funcEvals += _nseq[k];

				// The "x-coordinate" for extrapolation is (h/n)^2
				Real xest = (h / _nseq[k]) * (h / _nseq[k]);

				// Rational extrapolation
				rzextr(k, xest, yseq, yest, yerr);

				if (k > 0) { // Need at least 2 points to estimate error
					// Compute scaled error
					Real errMax = 0.0;
					for (int i = 0; i < _n; ++i) {
						Real scale = eps * (std::abs(ysav[i]) + std::abs(yest[i]) + Precision::DivisionSafetyThreshold);
						Real errRatio = std::abs(yerr[i]) / scale;
						errMax = std::max<Real>(errMax, errRatio);
					}

					result.errMax = errMax;

					if (errMax < 1.0) {
						// Converged!
						result.accepted = true;
						_xNew = yest;
						_err = yerr;

						_sys.derivs(t + h, yest, dxdt);
						_dydxFinal = dxdt;
						result.funcEvals++;
						result.hDone = h;
						_interpolator.setStep(t, h, x, yest, _dydxOld, _dydxFinal);
						x = yest;
						// Extrapolation order grows with the converged tableau column.
						result.hNext = _controller.acceptedStep(h, errMax, Real(1.0) / (2 * k + 1));

						return result;
					}
				}
			}

			// Failed to converge - reduce step size
			result.accepted = false;
			result.hDone = h;
			result.errMax = Real(999.0);
			result.hNext = _controller.rejectedStep(h, result.errMax, Real(1.0) / (2 * KMAXX + 1));

			return result;
		}

		Vector<Real> interpolate(Real t) const override {
			return _interpolator.interpolate(t);
		}

		bool isFSAL() const override { return false; }

		const Vector<Real>& getFinalDeriv() const override {
			return _dydxFinal;
		}

		int stageCount() const override { return KMAXX; }

		int order() const override { return 2 * KMAXX; }

		void resetFSAL() override {
			_controller.reset();
			_interpolator.reset();
		}
	};

} // namespace MML

#endif // MML_ODE_STEPPERS_H
