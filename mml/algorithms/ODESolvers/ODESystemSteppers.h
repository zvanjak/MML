///////////////////////////////////////////////////////////////////////////////////////////
///                         MinimalMathLibrary (MML)                                  ///
///                                                                                   ///
///  File:        ODESystemSteppers.h                                                 ///
///  Description: Adaptive ODE stepper implementations for use with integrators      ///
///               Includes: DormandPrince5, CashKarp, DormandPrince8, BulirschStoer   ///
///                                                                                   ///
///  Copyright:   (c) 2024-2026 Zvonimir Vanjak                                       ///
///  License:     MIT License (see LICENSE.md)                   ///
///                                                                                   ///
///////////////////////////////////////////////////////////////////////////////////////////
#if !defined MML_ODE_SYSTEM_STEPPERS_H
#define MML_ODE_SYSTEM_STEPPERS_H

#include "mml/MMLBase.h"
#include "mml/core/AlgorithmTypes.h"
#include "mml/interfaces/IODESystem.h"
#include "ODERKCoefficients.h"

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
	/// - 4th order dense output using Hermite interpolation
	/// - PI step size control for smooth adaptation
	class DormandPrince5_Stepper : public IAdaptiveStepper {
	private:
		const IODESystem& _sys;
		int _n; // System dimension

		// Stage vectors (stored for dense output)
		Vector<Real> _k1, _k2, _k3, _k4, _k5, _k6, _k7;
		Vector<Real> _xtemp; // Temporary for stage computation

		// Dense output data
		Real _tOld;			// Start time of last step
		Real _hDone;		// Size of last accepted step
		Vector<Real> _xOld; // State at start of last step
		bool _stepReady;	// Is dense output data valid?

		// FSAL state
		bool _haveFSAL; // Do we have a valid k7 to reuse?

		// Butcher tableau coefficients from centralized header
		using DP = RKCoeff::DormandPrince5;
		// Time nodes (c vector)
		static constexpr Real a2 = DP::c2;
		static constexpr Real a3 = DP::c3;
		static constexpr Real a4 = DP::c4;
		static constexpr Real a5 = DP::c5;
		static constexpr Real a6 = DP::c6;
		static constexpr Real a7 = DP::c7;
		// Stage coefficients (a matrix)
		static constexpr Real b21 = DP::a21;
		static constexpr Real b31 = DP::a31, b32 = DP::a32;
		static constexpr Real b41 = DP::a41, b42 = DP::a42, b43 = DP::a43;
		static constexpr Real b51 = DP::a51, b52 = DP::a52, b53 = DP::a53, b54 = DP::a54;
		static constexpr Real b61 = DP::a61, b62 = DP::a62, b63 = DP::a63, b64 = DP::a64, b65 = DP::a65;
		static constexpr Real b71 = DP::b1, b73 = DP::b3, b74 = DP::b4, b75 = DP::b5, b76 = DP::b6;
		// Error coefficients
		static constexpr Real e1 = DP::e1;
		static constexpr Real e3 = DP::e3;
		static constexpr Real e4 = DP::e4;
		static constexpr Real e5 = DP::e5;
		static constexpr Real e6 = DP::e6;
		static constexpr Real e7 = DP::e7;

		// Step size control
		static constexpr Real SAFETY = 0.9;
		static constexpr Real MIN_FACTOR = 0.2;
		static constexpr Real MAX_FACTOR = 10.0;
		static constexpr Real BETA = 0.04;	// For PI controller
		static constexpr Real ALPHA = 0.17; // For PI controller

		Real _errOld = 1.0; // Error from previous step (for PI controller)

	public:
		explicit DormandPrince5_Stepper(const IODESystem& sys)
			: _sys(sys)
			, _n(sys.getDim())
			, _tOld(0)
			, _hDone(0)
			, _stepReady(false)
			, _haveFSAL(false) {
			_k1.Resize(_n);
			_k2.Resize(_n);
			_k3.Resize(_n);
			_k4.Resize(_n);
			_k5.Resize(_n);
			_k6.Resize(_n);
			_k7.Resize(_n);
			_xtemp.Resize(_n);
			_xOld.Resize(_n);
		}

		StepResult doStep(Real t, Vector<Real>& x, Vector<Real>& dxdt, Real htry, Real eps) override {
			StepResult result;
			result.accepted = false;
			result.funcEvals = 0;

			Real h = htry;
			Vector<Real> xNew(_n);

			// Main step loop (retry with smaller h if rejected)
			while (true) {
				// Use FSAL if available
				if (_haveFSAL) {
					_k1 = _k7; // Reuse last stage
				} else {
					_k1 = dxdt;
					result.funcEvals++;
				}

				// Stage 2
				for (int i = 0; i < _n; i++) {
					_xtemp[i] = x[i] + h * b21 * _k1[i];
				}
				_sys.derivs(t + a2 * h, _xtemp, _k2);

				// Stage 3
				for (int i = 0; i < _n; i++) {
					_xtemp[i] = x[i] + h * (b31 * _k1[i] + b32 * _k2[i]);
				}
				_sys.derivs(t + a3 * h, _xtemp, _k3);

				// Stage 4
				for (int i = 0; i < _n; i++) {
					_xtemp[i] = x[i] + h * (b41 * _k1[i] + b42 * _k2[i] + b43 * _k3[i]);
				}
				_sys.derivs(t + a4 * h, _xtemp, _k4);

				// Stage 5
				for (int i = 0; i < _n; i++) {
					_xtemp[i] = x[i] + h * (b51 * _k1[i] + b52 * _k2[i] + b53 * _k3[i] + b54 * _k4[i]);
				}
				_sys.derivs(t + a5 * h, _xtemp, _k5);

				// Stage 6
				for (int i = 0; i < _n; i++) {
					_xtemp[i] = x[i] + h * (b61 * _k1[i] + b62 * _k2[i] + b63 * _k3[i] + b64 * _k4[i] + b65 * _k5[i]);
				}
				_sys.derivs(t + a6 * h, _xtemp, _k6);

				// Stage 7 (5th order solution, also FSAL)
				for (int i = 0; i < _n; i++) {
					xNew[i] = x[i] + h * (b71 * _k1[i] + b73 * _k3[i] + b74 * _k4[i] + b75 * _k5[i] + b76 * _k6[i]);
				}
				_sys.derivs(t + h, xNew, _k7);

				result.funcEvals += 6;

				// Error estimation
				Real errMax = 0.0;
				for (int i = 0; i < _n; i++) {
					Real err = h * (e1 * _k1[i] + e3 * _k3[i] + e4 * _k4[i] + e5 * _k5[i] + e6 * _k6[i] + e7 * _k7[i]);
					Real scale = std::abs(x[i]) + std::abs(h * _k1[i]) + 1e-30; // Avoid division by zero
					errMax = std::max(errMax, std::abs(err) / scale);
				}
				errMax /= eps;
				result.errMax = errMax;

				if (errMax <= 1.0) {
					// Step accepted
					result.accepted = true;
					result.hDone = h;

					// Store for dense output
					_tOld = t;
					_hDone = h;
					_xOld = x;
					_stepReady = true;

					// Update state
					x = xNew;
					dxdt = _k7;
					_haveFSAL = true;

					// PI controller for next step
					Real factor = SAFETY * std::pow(errMax, -ALPHA) * std::pow(_errOld, BETA);
					factor = std::max(MIN_FACTOR, std::min(MAX_FACTOR, factor));
					result.hNext = h * factor;

					_errOld = std::max<Real>(errMax, Real(1e-4)); // Prevent extreme growth
					break;
				}

				// Step rejected - reduce h and try again
				Real factor = SAFETY * std::pow(errMax, -ALPHA);
				factor = std::max<Real>(MIN_FACTOR, factor);
				h *= factor;
				_haveFSAL = false;

				if (std::abs(h) < Constants::Eps) {
					throw ODESolverError("Step size underflow in DormandPrince5_Stepper");
				}
			}

			return result;
		}

		Vector<Real> interpolate(Real t) const override {
			if (!_stepReady) {
				throw ODESolverError("No valid step data for interpolation");
			}

			// 4th order Hermite interpolation
			Real theta = (t - _tOld) / _hDone;
			Real theta2 = theta * theta;
			Real theta3 = theta2 * theta;

			// Hermite basis polynomials
			Real h00 = 2 * theta3 - 3 * theta2 + 1;		  // 1 - 3t² + 2t³
			Real h10 = theta3 - 2 * theta2 + theta;		  // t - 2t² + t³
			Real h01 = -2 * theta3 + 3 * theta2;		  // 3t² - 2t³
			Real h11 = theta3 - theta2;					  // t³ - t²

			Vector<Real> result(_n);
			for (int i = 0; i < _n; i++) {
				Real xEnd = _xOld[i] + _hDone * (b71 * _k1[i] + b73 * _k3[i] + b74 * _k4[i] + b75 * _k5[i] + b76 * _k6[i]);
				Real f0 = _hDone * _k1[i];
				Real f1 = _hDone * _k7[i];
				result[i] = h00 * _xOld[i] + h10 * f0 + h01 * xEnd + h11 * f1;
			}

			return result;
		}

		bool isFSAL() const override { return true; }

		const Vector<Real>& getFinalDeriv() const override { return _k7; }

		int stageCount() const override { return 7; }

		int order() const override { return 5; }

		void resetFSAL() override {
			_haveFSAL = false;
			_stepReady = false;
			_errOld = 1.0;
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

		// Stage vectors
		Vector<Real> _k1, _k2, _k3, _k4, _k5, _k6;
		Vector<Real> _xtemp;
		Vector<Real> _xNew;

		// Dense output data
		Real _tOld;
		Real _hDone;
		Vector<Real> _xOld;
		bool _stepReady;

		// Butcher tableau from centralized header
		using CK = RKCoeff::CashKarp5;
		// Time nodes (c vector)
		static constexpr Real a2 = CK::c2;
		static constexpr Real a3 = CK::c3;
		static constexpr Real a4 = CK::c4;
		static constexpr Real a5 = CK::c5;
		static constexpr Real a6 = CK::c6;
		// Stage coefficients (a matrix)
		static constexpr Real b21 = CK::a21;
		static constexpr Real b31 = CK::a31, b32 = CK::a32;
		static constexpr Real b41 = CK::a41, b42 = CK::a42, b43 = CK::a43;
		static constexpr Real b51 = CK::a51, b52 = CK::a52, b53 = CK::a53, b54 = CK::a54;
		static constexpr Real b61 = CK::a61, b62 = CK::a62, b63 = CK::a63, b64 = CK::a64, b65 = CK::a65;
		// 5th order coefficients (b weights)
		static constexpr Real c1 = CK::b1;
		static constexpr Real c3 = CK::b3;
		static constexpr Real c4 = CK::b4;
		static constexpr Real c6 = CK::b6;
		// Error coefficients
		static constexpr Real dc1 = CK::e1;
		static constexpr Real dc3 = CK::e3;
		static constexpr Real dc4 = CK::e4;
		static constexpr Real dc5 = CK::e5;
		static constexpr Real dc6 = CK::e6;

		static constexpr Real SAFETY = 0.9;
		static constexpr Real PGROW = -0.2;
		static constexpr Real PSHRINK = -0.25;
		static constexpr Real ERRCON = 1.89e-4; // (5/SAFETY)^(1/PGROW)

	public:
		explicit CashKarp_Stepper(const IODESystem& sys)
			: _sys(sys)
			, _n(sys.getDim())
			, _tOld(0)
			, _hDone(0)
			, _stepReady(false) {
			_k1.Resize(_n);
			_k2.Resize(_n);
			_k3.Resize(_n);
			_k4.Resize(_n);
			_k5.Resize(_n);
			_k6.Resize(_n);
			_xtemp.Resize(_n);
			_xNew.Resize(_n);
			_xOld.Resize(_n);
		}

		StepResult doStep(Real t, Vector<Real>& x, Vector<Real>& dxdt, Real htry, Real eps) override {
			StepResult result;
			result.accepted = false;

			Real h = htry;

			while (true) {
				// Compute all stages
				_k1 = dxdt;

				for (int i = 0; i < _n; i++)
					_xtemp[i] = x[i] + h * b21 * _k1[i];
				_sys.derivs(t + a2 * h, _xtemp, _k2);

				for (int i = 0; i < _n; i++)
					_xtemp[i] = x[i] + h * (b31 * _k1[i] + b32 * _k2[i]);
				_sys.derivs(t + a3 * h, _xtemp, _k3);

				for (int i = 0; i < _n; i++)
					_xtemp[i] = x[i] + h * (b41 * _k1[i] + b42 * _k2[i] + b43 * _k3[i]);
				_sys.derivs(t + a4 * h, _xtemp, _k4);

				for (int i = 0; i < _n; i++)
					_xtemp[i] = x[i] + h * (b51 * _k1[i] + b52 * _k2[i] + b53 * _k3[i] + b54 * _k4[i]);
				_sys.derivs(t + a5 * h, _xtemp, _k5);

				for (int i = 0; i < _n; i++)
					_xtemp[i] = x[i] + h * (b61 * _k1[i] + b62 * _k2[i] + b63 * _k3[i] + b64 * _k4[i] + b65 * _k5[i]);
				_sys.derivs(t + a6 * h, _xtemp, _k6);

				result.funcEvals = 6;

				// 5th order solution
				for (int i = 0; i < _n; i++)
					_xNew[i] = x[i] + h * (c1 * _k1[i] + c3 * _k3[i] + c4 * _k4[i] + c6 * _k6[i]);

				// Error estimate
				Real errMax = 0.0;
				for (int i = 0; i < _n; i++) {
					Real err = h * (dc1 * _k1[i] + dc3 * _k3[i] + dc4 * _k4[i] + dc5 * _k5[i] + dc6 * _k6[i]);
					Real scale = std::abs(x[i]) + std::abs(h * dxdt[i]) + 1e-30;
					errMax = std::max(errMax, std::abs(err) / scale);
				}
				errMax /= eps;
				result.errMax = errMax;

				if (errMax <= 1.0) {
					result.accepted = true;
					result.hDone = h;

					// Store for dense output
					_tOld = t;
					_hDone = h;
					_xOld = x;
					_stepReady = true;

					// Update state
					x = _xNew;
					_sys.derivs(t + h, x, dxdt);
					result.funcEvals++;

					// Next step size
					if (errMax > ERRCON)
						result.hNext = SAFETY * h * std::pow(errMax, PGROW);
					else
						result.hNext = Real(5.0) * h;

					break;
				}

				// Reduce step size
				Real htemp = SAFETY * h * std::pow(errMax, PSHRINK);
				h = (h >= 0) ? std::max<Real>(htemp, Real(0.1) * h) : std::min<Real>(htemp, Real(0.1) * h);

				if (std::abs(h) < Constants::Eps)
					throw ODESolverError("Step size underflow in CashKarp_Stepper");
			}

			return result;
		}

		Vector<Real> interpolate(Real t) const override {
			if (!_stepReady)
				throw ODESolverError("No valid step data for interpolation");

			// Simple Hermite interpolation
			Real theta = (t - _tOld) / _hDone;
			Real theta2 = theta * theta;
			Real theta3 = theta2 * theta;

			Real h00 = 2 * theta3 - 3 * theta2 + 1;
			Real h10 = theta3 - 2 * theta2 + theta;
			Real h01 = -2 * theta3 + 3 * theta2;
			Real h11 = theta3 - theta2;

			Vector<Real> result(_n);
			for (int i = 0; i < _n; i++) {
				Real f0 = _hDone * _k1[i];
				Real f1 = _hDone * (c1 * _k1[i] + c3 * _k3[i] + c4 * _k4[i] + c6 * _k6[i]); // Approximate
				result[i] = h00 * _xOld[i] + h10 * f0 + h01 * _xNew[i] + h11 * f1;
			}

			return result;
		}

		bool isFSAL() const override { return false; }

		const Vector<Real>& getFinalDeriv() const override {
			static Vector<Real> dummy;
			return dummy;
		}

		int stageCount() const override { return 6; }

		int order() const override { return 5; }

		void resetFSAL() override {
			_stepReady = false;
		}
	};

	//===================================================================================
	//                         DormandPrince8_Stepper
	//===================================================================================
	/// @brief Dormand-Prince 8(7) high-order adaptive stepper
	///
	/// High-order method for problems requiring very high accuracy.
	/// Features:
	/// - 8th order accurate solution with 7th order error estimate
	/// - 13 stages (FSAL optimization applies)
	/// - Excellent for smooth problems, astronomical trajectories
	/// - Higher cost per step but much larger accurate step sizes
	class DormandPrince8_Stepper : public IAdaptiveStepper {
	private:
		const IODESystem& _sys;
		int _n;

		static constexpr int NSTAGES = 13;
		std::vector<Vector<Real>> _k; // Stage vectors
		Vector<Real> _xNew;

		// Dense output data
		Real _tOld;
		Real _hDone;
		Vector<Real> _xOld;
		bool _stepReady;
		bool _haveFSAL;
		Real _errOld;

		// Dormand-Prince 8(7) coefficients from centralized header
		using DP8 = RKCoeff::DormandPrince8;
		// Time nodes
		static constexpr auto& c = DP8::c;
		// Stage coefficients (a matrix, row-by-row)
		static constexpr auto& a2 = DP8::a2;
		static constexpr auto& a3 = DP8::a3;
		static constexpr auto& a4 = DP8::a4;
		static constexpr auto& a5 = DP8::a5;
		static constexpr auto& a6 = DP8::a6;
		static constexpr auto& a7 = DP8::a7;
		static constexpr auto& a8 = DP8::a8;
		static constexpr auto& a9 = DP8::a9;
		static constexpr auto& a10 = DP8::a10;
		static constexpr auto& a11 = DP8::a11;
		static constexpr auto& a12 = DP8::a12;
		static constexpr auto& a13 = DP8::a13;
		// Solution weights
		static constexpr auto& b8 = DP8::b8;
		static constexpr auto& b7 = DP8::b7;

		static constexpr Real SAFETY = 0.9;
		static constexpr Real MIN_FACTOR = 0.2;
		static constexpr Real MAX_FACTOR = 6.0;
		static constexpr Real ALPHA = 1.0 / 8.0;

	public:
		explicit DormandPrince8_Stepper(const IODESystem& sys)
			: _sys(sys)
			, _n(sys.getDim())
			, _tOld(0)
			, _hDone(0)
			, _stepReady(false)
			, _haveFSAL(false)
			, _errOld(1.0) {
			_k.resize(NSTAGES);
			for (int i = 0; i < NSTAGES; ++i)
				_k[i].Resize(_n);
			_xNew.Resize(_n);
			_xOld.Resize(_n);
		}

		StepResult doStep(Real t, Vector<Real>& x, Vector<Real>& dxdt, Real htry, Real eps) override {
			StepResult result;
			result.accepted = false;
			result.funcEvals = 0;

			Real h = htry;

			while (true) {
				// Use FSAL if available
				if (_haveFSAL) {
					_k[0] = _k[12];
				} else {
					_k[0] = dxdt;
				}

				// Stage 2
				for (int i = 0; i < _n; i++)
					_xNew[i] = x[i] + h * a2[0] * _k[0][i];
				_sys.derivs(t + c[1] * h, _xNew, _k[1]);

				// Stage 3
				for (int i = 0; i < _n; i++)
					_xNew[i] = x[i] + h * (a3[0] * _k[0][i] + a3[1] * _k[1][i]);
				_sys.derivs(t + c[2] * h, _xNew, _k[2]);

				// Stage 4
				for (int i = 0; i < _n; i++)
					_xNew[i] = x[i] + h * (a4[0] * _k[0][i] + a4[2] * _k[2][i]);
				_sys.derivs(t + c[3] * h, _xNew, _k[3]);

				// Stage 5
				for (int i = 0; i < _n; i++)
					_xNew[i] = x[i] + h * (a5[0] * _k[0][i] + a5[2] * _k[2][i] + a5[3] * _k[3][i]);
				_sys.derivs(t + c[4] * h, _xNew, _k[4]);

				// Stage 6
				for (int i = 0; i < _n; i++)
					_xNew[i] = x[i] + h * (a6[0] * _k[0][i] + a6[3] * _k[3][i] + a6[4] * _k[4][i]);
				_sys.derivs(t + c[5] * h, _xNew, _k[5]);

				// Stage 7
				for (int i = 0; i < _n; i++)
					_xNew[i] = x[i] + h * (a7[0] * _k[0][i] + a7[3] * _k[3][i] + a7[4] * _k[4][i] + a7[5] * _k[5][i]);
				_sys.derivs(t + c[6] * h, _xNew, _k[6]);

				// Stage 8
				for (int i = 0; i < _n; i++)
					_xNew[i] = x[i] + h * (a8[0] * _k[0][i] + a8[3] * _k[3][i] + a8[4] * _k[4][i] + a8[5] * _k[5][i] + a8[6] * _k[6][i]);
				_sys.derivs(t + c[7] * h, _xNew, _k[7]);

				// Stage 9
				for (int i = 0; i < _n; i++)
					_xNew[i] = x[i] + h * (a9[0] * _k[0][i] + a9[3] * _k[3][i] + a9[4] * _k[4][i] + a9[5] * _k[5][i] + a9[6] * _k[6][i] + a9[7] * _k[7][i]);
				_sys.derivs(t + c[8] * h, _xNew, _k[8]);

				// Stage 10
				for (int i = 0; i < _n; i++)
					_xNew[i] = x[i] + h * (a10[0] * _k[0][i] + a10[3] * _k[3][i] + a10[4] * _k[4][i] + a10[5] * _k[5][i] + a10[6] * _k[6][i] + a10[7] * _k[7][i] + a10[8] * _k[8][i]);
				_sys.derivs(t + c[9] * h, _xNew, _k[9]);

				// Stage 11
				for (int i = 0; i < _n; i++)
					_xNew[i] = x[i] + h * (a11[0] * _k[0][i] + a11[3] * _k[3][i] + a11[4] * _k[4][i] + a11[5] * _k[5][i] + a11[6] * _k[6][i] + a11[7] * _k[7][i] + a11[8] * _k[8][i] + a11[9] * _k[9][i]);
				_sys.derivs(t + c[10] * h, _xNew, _k[10]);

				// Stage 12
				for (int i = 0; i < _n; i++)
					_xNew[i] = x[i] + h * (a12[0] * _k[0][i] + a12[3] * _k[3][i] + a12[4] * _k[4][i] + a12[5] * _k[5][i] + a12[6] * _k[6][i] + a12[7] * _k[7][i] + a12[8] * _k[8][i] + a12[9] * _k[9][i] + a12[10] * _k[10][i]);
				_sys.derivs(t + c[11] * h, _xNew, _k[11]);

				// Stage 13 (final, also FSAL)
				for (int i = 0; i < _n; i++)
					_xNew[i] = x[i] + h * (a13[0] * _k[0][i] + a13[3] * _k[3][i] + a13[4] * _k[4][i] + a13[5] * _k[5][i] + a13[6] * _k[6][i] + a13[7] * _k[7][i] + a13[8] * _k[8][i] + a13[9] * _k[9][i] + a13[10] * _k[10][i]);
				_sys.derivs(t + c[12] * h, _xNew, _k[12]);

				result.funcEvals = 13;

				// 8th order solution
				for (int i = 0; i < _n; i++) {
					_xNew[i] = x[i];
					for (int j = 0; j < NSTAGES; j++)
						_xNew[i] += h * b8[j] * _k[j][i];
				}

				// Error estimate (8th - 7th order)
				Real errMax = 0.0;
				for (int i = 0; i < _n; i++) {
					Real err = 0.0;
					for (int j = 0; j < NSTAGES; j++)
						err += (b8[j] - b7[j]) * _k[j][i];
					err *= h;
					Real scale = std::abs(x[i]) + std::abs(h * _k[0][i]) + 1e-30;
					errMax = std::max(errMax, std::abs(err) / scale);
				}
				errMax /= eps;
				result.errMax = errMax;

				if (errMax <= 1.0) {
					result.accepted = true;
					result.hDone = h;

					// Store for dense output
					_tOld = t;
					_hDone = h;
					_xOld = x;
					_stepReady = true;

					// Update state
					x = _xNew;
					dxdt = _k[12];
					_haveFSAL = true;

					// Step size control for next step
					Real factor = SAFETY * std::pow(errMax, -ALPHA);
					factor = std::max<Real>(MIN_FACTOR, std::min<Real>(MAX_FACTOR, factor));
					result.hNext = h * factor;

					_errOld = std::max<Real>(errMax, Real(1e-4));
					break;
				}

				// Reduce step size
				Real factor = SAFETY * std::pow(errMax, -ALPHA);
				factor = std::max<Real>(MIN_FACTOR, factor);
				h *= factor;
				_haveFSAL = false;

				if (std::abs(h) < Constants::Eps)
					throw ODESolverError("Step size underflow in DormandPrince8_Stepper");
			}

			return result;
		}

		Vector<Real> interpolate(Real t) const override {
			if (!_stepReady)
				throw ODESolverError("No valid step data for interpolation");

			Real theta = (t - _tOld) / _hDone;
			Real theta2 = theta * theta;
			Real theta3 = theta2 * theta;

			// Hermite interpolation
			Real h00 = 2 * theta3 - 3 * theta2 + 1;
			Real h10 = theta3 - 2 * theta2 + theta;
			Real h01 = -2 * theta3 + 3 * theta2;
			Real h11 = theta3 - theta2;

			Vector<Real> result(_n);

			for (int i = 0; i < _n; i++) {
				Real hf0 = _hDone * _k[0][i];
				Real hf1 = _hDone * _k[12][i];
				result[i] = h00 * _xOld[i] + h10 * hf0 + h01 * _xNew[i] + h11 * hf1;
			}

			return result;
		}

		bool isFSAL() const override { return true; }

		const Vector<Real>& getFinalDeriv() const override { return _k[12]; }

		int stageCount() const override { return NSTAGES; }

		int order() const override { return 8; }

		void resetFSAL() override {
			_haveFSAL = false;
			_stepReady = false;
			_errOld = 1.0;
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
		mutable Real _tOld;			///< Time at start of step
		mutable Real _hDone;		///< Actual step size taken

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
			, _tOld(0)
			, _hDone(0)
			, _xOld(_n)
			, _xNew(_n)
			, _err(_n) {
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

			_tOld = t;
			_xOld = x;
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
						Real scale = eps * (std::abs(ysav[i]) + std::abs(yest[i]) + Real(1e-30));
						Real errRatio = std::abs(yerr[i]) / scale;
						errMax = std::max<Real>(errMax, errRatio);
					}

					result.errMax = errMax;

					if (errMax < 1.0) {
						// Converged!
						result.accepted = true;
						_xNew = yest;
						_err = yerr;
						_hDone = h;

						x = yest;
						_sys.derivs(t + h, x, dxdt);
						result.funcEvals++;
						result.hDone = h;

						// Step size control - be conservative
						Real factor = Real(0.9) * std::pow(errMax, Real(-1.0) / (2 * k + 1));
						factor = std::max<Real>(Real(0.1), std::min<Real>(factor, Real(4.0)));
						result.hNext = h * factor;

						return result;
					}
				}
			}

			// Failed to converge - reduce step size
			result.accepted = false;
			result.hDone = h;
			result.hNext = h * Real(0.5);
			result.errMax = Real(999.0);

			return result;
		}

		Vector<Real> interpolate(Real t) const override {
			// Simple linear interpolation
			Real theta = (t - _tOld) / _hDone;
			return _xOld * (1.0 - theta) + _xNew * theta;
		}

		bool isFSAL() const override { return false; }

		const Vector<Real>& getFinalDeriv() const override {
			static Vector<Real> dummy;
			return dummy;
		}

		int stageCount() const override { return KMAXX; }

		int order() const override { return 2 * KMAXX; }

		void resetFSAL() override {}
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
		mutable Real _tOld;			///< Time at start of step
		mutable Real _hDone;		///< Actual step size taken

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
			, _tOld(0)
			, _hDone(0)
			, _xOld(_n)
			, _xNew(_n)
			, _err(_n) {
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

			_tOld = t;
			_xOld = x;
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
						Real scale = eps * (std::abs(ysav[i]) + std::abs(yest[i]) + Real(1e-30));
						Real errRatio = std::abs(yerr[i]) / scale;
						errMax = std::max<Real>(errMax, errRatio);
					}

					result.errMax = errMax;

					if (errMax < 1.0) {
						// Converged!
						result.accepted = true;
						_xNew = yest;
						_err = yerr;
						_hDone = h;

						x = yest;
						_sys.derivs(t + h, x, dxdt);
						result.funcEvals++;
						result.hDone = h;

						// Step size control - be conservative
						Real factor = Real(0.9) * std::pow(errMax, Real(-1.0) / (2 * k + 1));
						factor = std::max<Real>(Real(0.1), std::min<Real>(factor, Real(4.0)));
						result.hNext = h * factor;

						return result;
					}
				}
			}

			// Failed to converge - reduce step size
			result.accepted = false;
			result.hDone = h;
			result.hNext = h * Real(0.5);
			result.errMax = Real(999.0);

			return result;
		}

		Vector<Real> interpolate(Real t) const override {
			// Simple linear interpolation
			Real theta = (t - _tOld) / _hDone;
			return _xOld * (1.0 - theta) + _xNew * theta;
		}

		bool isFSAL() const override { return false; }

		const Vector<Real>& getFinalDeriv() const override {
			static Vector<Real> dummy;
			return dummy;
		}

		int stageCount() const override { return KMAXX; }

		int order() const override { return 2 * KMAXX; }

		void resetFSAL() override {}
	};

} // namespace MML

#endif // MML_ODE_SYSTEM_STEPPERS_H
