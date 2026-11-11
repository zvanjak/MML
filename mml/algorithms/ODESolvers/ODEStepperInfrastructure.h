///////////////////////////////////////////////////////////////////////////////////////////
///                         MinimalMathLibrary (MML)                                  ///
///                                                                                   ///
///  File:        ODEStepperInfrastructure.h                                         ///
///  Description: Shared adaptive ODE stepper control and interpolation machinery    ///
///                                                                                   ///
///  Copyright:   (c) 2024-2026 Zvonimir Vanjak                                      ///
///  License:     MIT License (see LICENSE.md)                                        ///
///////////////////////////////////////////////////////////////////////////////////////////
#if !defined MML_ODE_STEPPER_INFRASTRUCTURE_H
#define MML_ODE_STEPPER_INFRASTRUCTURE_H

#include <mml/MMLBase.h>
#include <mml/interfaces/IODESystem.h>
#include <mml/algorithms/ODESolvers/ODERKCoefficients.h>

#include <algorithm>
#include <array>
#include <cmath>
#include <vector>

namespace MML {

	namespace ODEStepperDetail {
		class CountingODESystem final : public IODESystem {
			const IODESystem& _system;
			int& _evaluations;
			int _limit;
		public:
			CountingODESystem(const IODESystem& system, int& evaluations, int limit = std::numeric_limits<int>::max())
				: _system(system), _evaluations(evaluations), _limit(limit) {}
			int getDim() const override { return _system.getDim(); }
			void derivs(Real time, const Vector<Real>& state, Vector<Real>& derivative) const override {
				if (_evaluations >= _limit)
					throw ResourceLimitError("Adaptive ODE RHS evaluation counter exhausted");
				++_evaluations;
				_system.derivs(time, state, derivative);
			}
		};

		inline void RequireFinite(Real value, const char* message) {
			if (!std::isfinite(value))
				throw ODESolverError(message);
		}

		inline void RequireVector(const Vector<Real>& values, int dimension, const char* message) {
			if (dimension <= 0 || values.size() != dimension)
				throw ODESolverError(message);
			for (int i = 0; i < dimension; ++i)
				RequireFinite(values[i], message);
		}

		inline void EvaluateRHS(const IODESystem& system, Real time, const Vector<Real>& state, Vector<Real>& derivative) {
			RequireFinite(time, "Non-finite ODE evaluation time");
			RequireVector(state, system.getDim(), "Invalid ODE evaluation state");
			system.derivs(time, state, derivative);
			RequireVector(derivative, system.getDim(), "Invalid ODE right-hand side");
		}

		inline void ValidateStep(int dimension, Real time, const Vector<Real>& state,
			const Vector<Real>& derivative, Real step, Real tolerance) {
			RequireFinite(time, "Non-finite ODE step time");
			RequireFinite(step, "Non-finite ODE step size");
			RequireFinite(time + step, "Non-finite ODE step endpoint");
			RequireFinite(tolerance, "Non-finite ODE tolerance");
			if (step == 0 || tolerance <= 0)
				throw ODESolverError("ODE step size must be nonzero and tolerance must be positive");
			RequireVector(state, dimension, "Invalid ODE step state");
			RequireVector(derivative, dimension, "Invalid ODE initial derivative");
		}

		inline Real RKErrorScale(Real state, Real derivative, Real step, Real candidate) {
			RequireFinite(candidate, "Non-finite ODE candidate for error scaling");
			const Real history = std::abs(state) + std::abs(step * derivative) + Precision::DivisionSafetyThreshold;
			RequireFinite(history, "Non-finite ODE error scale");
			// A stationary zero can grow at higher order than the Euler predictor.
			return std::max(history, std::abs(candidate));
		}

		inline Real ErrorRatio(Real error, Real scale) {
			RequireFinite(error, "Non-finite ODE error estimate");
			RequireFinite(scale, "Non-finite ODE error scale");
			if (scale <= 0)
				throw ODESolverError("ODE error scale must be positive");
			const Real ratio = std::abs(error) / scale;
			RequireFinite(ratio, "Non-finite ODE normalized error");
			return ratio;
		}
	}

	/// Tunable limits and exponents for adaptive step-size control.
	/// Errors passed to the controller are normalized, so error <= 1 means acceptable.
	struct StepSizeControllerConfig {
		Real safety = Real(0.9);       ///< Conservative multiplier applied to every proposal.
		Real alpha = Real(0.17);       ///< Proportional exponent for the current error.
		Real beta = Real(0.04);        ///< Integral exponent for the previous accepted error.
		Real minFactor = Real(0.2);    ///< Maximum shrinkage allowed in one proposal.
		Real maxFactor = Real(10.0);   ///< Maximum growth allowed in one proposal.
		Real errorFloor = Precision::StepControllerErrorFloor;  ///< Prevents zero error from producing an infinite factor.
	};

	/// Stateful PI controller shared by adaptive ODE steppers.
	///
	/// Accepted steps use safety * error^(-alpha) * previousError^beta. Rejected
	/// steps omit the history term because a rejected estimate must not enter the
	/// accepted-error history. The sign of h is preserved for backward integration.
	class StepSizeController {
	private:
		StepSizeControllerConfig _config;
		Real _previousError = Real(1.0);

		static void validate(const StepSizeControllerConfig& config) {
			for (Real value : {config.safety, config.alpha, config.beta, config.minFactor,
				config.maxFactor, config.errorFloor})
				ODEStepperDetail::RequireFinite(value, "Non-finite step-size controller configuration");
			if (config.safety <= 0 || config.alpha <= 0 || config.beta < 0 ||
				config.minFactor <= 0 || config.maxFactor < config.minFactor || config.errorFloor <= 0) {
				throw ODESolverError("Invalid step-size controller configuration");
			}
		}

		Real boundedFactor(Real factor) const {
			return std::max(_config.minFactor, std::min(_config.maxFactor, factor));
		}

	public:
		explicit StepSizeController(const StepSizeControllerConfig& config = {})
			: _config(config) {
			validate(_config);
		}

		Real acceptedStep(Real h, Real error) {
			return acceptedStep(h, error, _config.alpha);
		}

		Real acceptedStep(Real h, Real error, Real alpha) {
			ODEStepperDetail::RequireFinite(h, "Non-finite ODE controller step");
			ODEStepperDetail::RequireFinite(error, "Non-finite ODE controller error");
			ODEStepperDetail::RequireFinite(alpha, "Non-finite ODE controller exponent");
			if (alpha <= 0)
				throw ODESolverError("Step-size controller exponent must be positive");
			const Real boundedError = std::max(std::abs(error), _config.errorFloor);
			const Real factor = boundedFactor(_config.safety * std::pow(boundedError, -alpha) *
													 std::pow(_previousError, _config.beta));
			const Real next = h * factor;
			ODEStepperDetail::RequireFinite(next, "Non-finite ODE step proposal");
			if (next == 0)
				throw ODESolverError("ODE step proposal underflow");
			_previousError = boundedError;
			return next;
		}

		Real rejectedStep(Real h, Real error) const {
			return rejectedStep(h, error, _config.alpha);
		}

		Real rejectedStep(Real h, Real error, Real alpha) const {
			ODEStepperDetail::RequireFinite(h, "Non-finite ODE controller step");
			ODEStepperDetail::RequireFinite(error, "Non-finite ODE controller error");
			ODEStepperDetail::RequireFinite(alpha, "Non-finite ODE controller exponent");
			if (alpha <= 0)
				throw ODESolverError("Step-size controller exponent must be positive");
			const Real boundedError = std::max(std::abs(error), _config.errorFloor);
			const Real factor = boundedFactor(_config.safety * std::pow(boundedError, -alpha));
			const Real next = h * factor;
			ODEStepperDetail::RequireFinite(next, "Non-finite ODE step proposal");
			if (next == 0)
				throw ODESolverError("ODE step proposal underflow");
			return next;
		}

		void reset() { _previousError = Real(1.0); }

		Real previousError() const { return _previousError; }
		const StepSizeControllerConfig& config() const { return _config; }
	};

	/// Cubic Hermite dense output over the last accepted ODE step.
	///
	/// Endpoint states and derivatives are copied so interpolation remains valid
	/// after the stepper updates its working state. Queries outside the interval
	/// clamp to the nearest endpoint, including when stepSize is negative.
	class HermiteInterpolator {
	private:
		Real _tStart = 0;
		Real _stepSize = 0;
		Vector<Real> _xStart;
		Vector<Real> _xEnd;
		Vector<Real> _dxStart;
		Vector<Real> _dxEnd;
		bool _ready = false;

	public:
		void setStep(Real tStart, Real stepSize, const Vector<Real>& xStart, const Vector<Real>& xEnd,
					 const Vector<Real>& dxStart, const Vector<Real>& dxEnd) {
			_ready = false;
			if (stepSize == 0 || !std::isfinite(stepSize))
				throw ODESolverError("Hermite interpolation requires a finite, non-zero step size");
			if (xStart.size() != xEnd.size() || xStart.size() != dxStart.size() || xStart.size() != dxEnd.size())
				throw ODESolverError("Hermite interpolation data dimensions do not match");
			ODEStepperDetail::ValidateStep(xStart.size(), tStart, xStart, dxStart, stepSize, Real(1));
			ODEStepperDetail::RequireVector(xEnd, xStart.size(), "Invalid Hermite endpoint state");
			ODEStepperDetail::RequireVector(dxEnd, xStart.size(), "Invalid Hermite endpoint derivative");
			for (int i = 0; i < xStart.size(); ++i) {
				ODEStepperDetail::RequireFinite(stepSize * dxStart[i], "Non-finite Hermite start slope");
				ODEStepperDetail::RequireFinite(stepSize * dxEnd[i], "Non-finite Hermite endpoint slope");
			}

			_tStart = tStart;
			_stepSize = stepSize;
			_xStart = xStart;
			_xEnd = xEnd;
			_dxStart = dxStart;
			_dxEnd = dxEnd;
			_ready = true;
		}

		Vector<Real> interpolate(Real t) const {
			if (!_ready)
				throw ODESolverError("No valid step data for interpolation");

			ODEStepperDetail::RequireFinite(t, "Non-finite Hermite interpolation time");
			// Normalizing by signed stepSize gives the same [0,1] interval in both directions.
			const Real theta = (t - _tStart) / _stepSize;
			if (theta <= 0)
				return _xStart;
			if (theta >= 1)
				return _xEnd;

			const Real theta2 = theta * theta;
			const Real theta3 = theta2 * theta;
			// Standard Hermite basis for values (h00, h01) and derivatives (h10, h11).
			const Real h00 = 2 * theta3 - 3 * theta2 + 1;
			const Real h10 = theta3 - 2 * theta2 + theta;
			const Real h01 = -2 * theta3 + 3 * theta2;
			const Real h11 = theta3 - theta2;

			Vector<Real> result(_xStart.size());
			for (int i = 0; i < _xStart.size(); ++i) {
				result[i] = h00 * _xStart[i] + h10 * _stepSize * _dxStart[i] +
							h01 * _xEnd[i] + h11 * _stepSize * _dxEnd[i];
			}
			ODEStepperDetail::RequireVector(result, _xStart.size(), "Non-finite Hermite interpolated state");
			return result;
		}

		void reset() { _ready = false; }
		bool ready() const { return _ready; }
	};

	/// Shampine's quartic continuous extension for Dormand-Prince 5(4).
	///
	/// This is fourth-order dense output with O(h^5) local interpolation error.
	/// The four polynomial coefficient vectors are formed from the seven RK
	/// stages already available after an accepted step, so dense output adds no
	/// derivative evaluations. Endpoint states are retained explicitly to make
	/// clamping exact and symmetric for positive and negative steps.
	class DormandPrince5Interpolator {
	private:
		Real _tStart = 0;
		Real _stepSize = 0;
		Vector<Real> _xStart;
		Vector<Real> _xEnd;
		std::array<Vector<Real>, 4> _coefficients;
		bool _ready = false;

	public:
		void setStep(Real tStart, Real stepSize, const Vector<Real>& xStart, const Vector<Real>& xEnd,
					 const std::vector<Vector<Real>>& stages) {
			_ready = false;
			if (stepSize == 0 || !std::isfinite(stepSize))
				throw ODESolverError("Dormand-Prince interpolation requires a finite, non-zero step size");
			if (stages.size() != 7 || xStart.size() != xEnd.size())
				throw ODESolverError("Dormand-Prince interpolation data dimensions do not match");
			for (const auto& stage : stages) {
				if (stage.size() != xStart.size())
					throw ODESolverError("Dormand-Prince interpolation data dimensions do not match");
				ODEStepperDetail::RequireVector(stage, xStart.size(), "Invalid Dormand-Prince interpolation stage");
			}
			ODEStepperDetail::ValidateStep(xStart.size(), tStart, xStart, stages.front(), stepSize, Real(1));
			ODEStepperDetail::RequireVector(xEnd, xStart.size(), "Invalid Dormand-Prince endpoint state");

			// Rows correspond to k1...k7 and columns to theta...theta^4.
			static constexpr Real weights[7][4] = {
				{MML::Real(1.0), -MML::Real(8048581381.0) / MML::Real(2820520608.0), MML::Real(8663915743.0) / MML::Real(2820520608.0), -MML::Real(12715105075.0) / MML::Real(11282082432.0)},
				{0.0, 0.0, 0.0, 0.0},
				{MML::Real(0.0), MML::Real(131558114200.0) / MML::Real(32700410799.0), -MML::Real(68118460800.0) / MML::Real(10900136933.0), MML::Real(87487479700.0) / MML::Real(32700410799.0)},
				{MML::Real(0.0), -MML::Real(1754552775.0) / MML::Real(470086768.0), MML::Real(14199869525.0) / MML::Real(1410260304.0), -MML::Real(10690763975.0) / MML::Real(1880347072.0)},
				{MML::Real(0.0), MML::Real(127303824393.0) / MML::Real(49829197408.0), -MML::Real(318862633887.0) / MML::Real(49829197408.0), MML::Real(701980252875.0) / MML::Real(199316789632.0)},
				{MML::Real(0.0), -MML::Real(282668133.0) / MML::Real(205662961.0), MML::Real(2019193451.0) / MML::Real(616988883.0), -MML::Real(1453857185.0) / MML::Real(822651844.0)},
				{MML::Real(0.0), MML::Real(40617522.0) / MML::Real(29380423.0), -MML::Real(110615467.0) / MML::Real(29380423.0), MML::Real(69997945.0) / MML::Real(29380423.0)}
			};

			_tStart = tStart;
			_stepSize = stepSize;
			_xStart = xStart;
			_xEnd = xEnd;
			for (int power = 0; power < 4; ++power) {
				_coefficients[power].resize(xStart.size());
				for (int component = 0; component < xStart.size(); ++component) {
					Real coefficient = 0;
					for (int stage = 0; stage < 7; ++stage)
						coefficient += weights[stage][power] * stages[stage][component];
					_coefficients[power][component] = coefficient;
					ODEStepperDetail::RequireFinite(coefficient, "Non-finite Dormand-Prince dense coefficient");
					ODEStepperDetail::RequireFinite(stepSize * coefficient, "Non-finite Dormand-Prince scaled dense coefficient");
				}
			}
			_ready = true;
		}

		Vector<Real> interpolate(Real t) const {
			if (!_ready)
				throw ODESolverError("No valid step data for interpolation");

			const Real theta = (t - _tStart) / _stepSize;
			if (theta <= 0)
				return _xStart;
			if (theta >= 1)
				return _xEnd;

			Vector<Real> result(_xStart);
			for (int component = 0; component < result.size(); ++component) {
				// Horner evaluation of h * (q1*theta + ... + q4*theta^4).
				result[component] += _stepSize * theta *
					(_coefficients[0][component] + theta *
					 (_coefficients[1][component] + theta *
					  (_coefficients[2][component] + theta * _coefficients[3][component])));
			}
			ODEStepperDetail::RequireVector(result, _xStart.size(), "Non-finite Dormand-Prince interpolated state");
			return result;
		}

		void reset() { _ready = false; }
		bool ready() const { return _ready; }
	};

	/// Hairer's seventh-order continuous extension for DOP853.
	///
	/// The accepted-step endpoint derivative and three interpolation-only stages
	/// complete the 16 derivative vectors used to form the seven coefficients.
	/// Once cached, interpolation does not evaluate the ODE system.
	class DormandPrince8Interpolator {
	private:
		using DP8 = RKCoeff::DormandPrince8;

		Real _tStart = 0;
		Real _stepSize = 0;
		Vector<Real> _xStart;
		Vector<Real> _xEnd;
		std::array<Vector<Real>, 7> _coefficients;
		bool _ready = false;

	public:
		int setStep(const IODESystem& system, Real tStart, Real stepSize,
					const Vector<Real>& xStart, const Vector<Real>& xEnd,
					std::vector<Vector<Real>>& stages, Vector<Real>& workspace) {
			_ready = false;
			if (stepSize == 0 || !std::isfinite(stepSize))
				throw ODESolverError("DOP853 interpolation requires a finite, non-zero step size");
			if (stages.size() != DP8::extended_stages || xStart.size() != xEnd.size())
				throw ODESolverError("DOP853 interpolation data dimensions do not match");

			const int dimension = system.getDim();
			for (auto& stage : stages) {
				if (stage.size() != dimension)
					stage.resize(dimension);
			}
			if (workspace.size() != dimension)
				workspace.resize(dimension);

			ODEStepperDetail::ValidateStep(dimension, tStart, xStart, stages.front(), stepSize, Real(1));
			ODEStepperDetail::RequireVector(xEnd, dimension, "Invalid DOP853 endpoint state");
			for (int stage = 0; stage <= DP8::stages; ++stage)
				ODEStepperDetail::RequireVector(stages[stage], dimension, "Invalid DOP853 interpolation stage");

			// Stages 13 through 15 are used only by the continuous extension.
			for (int stage = DP8::stages + 1; stage < DP8::extended_stages; ++stage) {
				for (int component = 0; component < dimension; ++component) {
					Real increment = 0;
					for (int previousStage = 0; previousStage < stage; ++previousStage)
						increment += DP8::a[stage][previousStage] * stages[previousStage][component];
					workspace[component] = xStart[component] + stepSize * increment;
				}
				ODEStepperDetail::EvaluateRHS(system, tStart + DP8::c[stage] * stepSize, workspace, stages[stage]);
			}

			_tStart = tStart;
			_stepSize = stepSize;
			_xStart = xStart;
			_xEnd = xEnd;
			for (auto& coefficient : _coefficients)
				coefficient.resize(dimension);

			for (int component = 0; component < dimension; ++component) {
				const Real delta = xEnd[component] - xStart[component];
				_coefficients[0][component] = delta;
				_coefficients[1][component] = stepSize * stages[0][component] - delta;
				_coefficients[2][component] = 2 * delta - stepSize * (stages[12][component] + stages[0][component]);
				for (int coefficient = 0; coefficient < 4; ++coefficient) {
					Real value = 0;
					for (int stage = 0; stage < DP8::extended_stages; ++stage)
						value += DP8::denseWeight(coefficient, stage) * stages[stage][component];
					_coefficients[coefficient + 3][component] = stepSize * value;
				}
			}
			for (const auto& coefficient : _coefficients)
				ODEStepperDetail::RequireVector(coefficient, dimension, "Non-finite DOP853 dense coefficient");
			_ready = true;
			return 3;
		}

		Vector<Real> interpolate(Real t) const {
			if (!_ready)
				throw ODESolverError("No valid step data for interpolation");

			const Real theta = (t - _tStart) / _stepSize;
			if (theta <= 0)
				return _xStart;
			if (theta >= 1)
				return _xEnd;

			Vector<Real> result(_xStart.size());
			for (int component = 0; component < result.size(); ++component) {
				Real value = 0;
				for (int iteration = 0; iteration < 7; ++iteration) {
					value += _coefficients[6 - iteration][component];
					value *= iteration % 2 == 0 ? theta : 1 - theta;
				}
				result[component] = _xStart[component] + value;
			}
			ODEStepperDetail::RequireVector(result, _xStart.size(), "Non-finite DOP853 interpolated state");
			return result;
		}

		void reset() { _ready = false; }
		bool ready() const { return _ready; }
	};

	/// Tableau-driven stage engine for explicit embedded Runge-Kutta methods.
	///
	/// Tableau must expose stages, node(), stageCoefficient(), solutionWeight(),
	/// and errorWeight(). The caller supplies k1 because adaptive steppers may
	/// reuse it through FSAL; evaluate() therefore performs stages - 1 RHS calls.
	template<typename Tableau>
	class ExplicitRKStageEvaluator {
	public:
		static int evaluate(const IODESystem& system, Real t, const Vector<Real>& x, Real h,
							const Vector<Real>& initialDerivative, std::vector<Vector<Real>>& stages,
							Vector<Real>& workspace) {
			const int dimension = system.getDim();
			ODEStepperDetail::RequireFinite(t, "Non-finite Runge-Kutta start time");
			ODEStepperDetail::RequireFinite(h, "Non-finite Runge-Kutta step");
			ODEStepperDetail::RequireVector(x, dimension, "Invalid Runge-Kutta start state");
			ODEStepperDetail::RequireVector(initialDerivative, dimension, "Invalid Runge-Kutta initial derivative");
			if (stages.size() < Tableau::stages)
				stages.resize(Tableau::stages);
			for (auto& stage : stages) {
				if (stage.size() != dimension)
					stage.resize(dimension);
			}
			if (workspace.size() != dimension)
				workspace.resize(dimension);

			stages[0] = initialDerivative;
			// Explicit tableaux are strictly lower triangular: stage i uses only k[0..i-1].
			for (int stage = 1; stage < Tableau::stages; ++stage) {
				for (int component = 0; component < dimension; ++component) {
					Real increment = 0;
					for (int previousStage = 0; previousStage < stage; ++previousStage)
						increment += Tableau::stageCoefficient(stage, previousStage) * stages[previousStage][component];
					workspace[component] = x[component] + h * increment;
				}
				ODEStepperDetail::EvaluateRHS(system, t + Tableau::node(stage) * h, workspace, stages[stage]);
			}
			return Tableau::stages - 1;
		}

		static void combineSolution(const Vector<Real>& x, Real h, const std::vector<Vector<Real>>& stages,
								Vector<Real>& result) {
			if (result.size() != x.size())
				result.resize(x.size());
			for (int component = 0; component < x.size(); ++component) {
				Real increment = 0;
				for (int stage = 0; stage < Tableau::stages; ++stage)
					increment += Tableau::solutionWeight(stage) * stages[stage][component];
				result[component] = x[component] + h * increment;
			}
			ODEStepperDetail::RequireVector(result, x.size(), "Non-finite Runge-Kutta candidate");
		}

		static void combineError(Real h, const std::vector<Vector<Real>>& stages, Vector<Real>& error) {
			const int dimension = stages[0].size();
			if (error.size() != dimension)
				error.resize(dimension);
			for (int component = 0; component < dimension; ++component) {
				Real estimate = 0;
				// errorWeight is the difference between the high- and low-order weights.
				for (int stage = 0; stage < Tableau::stages; ++stage)
					estimate += Tableau::errorWeight(stage) * stages[stage][component];
				error[component] = h * estimate;
			}
			ODEStepperDetail::RequireVector(error, dimension, "Non-finite Runge-Kutta error estimate");
		}
	};

} // namespace MML

#endif // MML_ODE_STEPPER_INFRASTRUCTURE_H