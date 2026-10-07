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

#include <algorithm>
#include <array>
#include <cmath>
#include <vector>

namespace MML {

	/// Tunable limits and exponents for adaptive step-size control.
	/// Errors passed to the controller are normalized, so error <= 1 means acceptable.
	struct StepSizeControllerConfig {
		Real safety = Real(0.9);       ///< Conservative multiplier applied to every proposal.
		Real alpha = Real(0.17);       ///< Proportional exponent for the current error.
		Real beta = Real(0.04);        ///< Integral exponent for the previous accepted error.
		Real minFactor = Real(0.2);    ///< Maximum shrinkage allowed in one proposal.
		Real maxFactor = Real(10.0);   ///< Maximum growth allowed in one proposal.
		Real errorFloor = Real(1e-4);  ///< Prevents zero error from producing an infinite factor.
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
			if (alpha <= 0)
				throw ODESolverError("Step-size controller exponent must be positive");
			const Real boundedError = std::max(std::abs(error), _config.errorFloor);
			const Real factor = boundedFactor(_config.safety * std::pow(boundedError, -alpha) *
													 std::pow(_previousError, _config.beta));
			_previousError = boundedError;
			return h * factor;
		}

		Real rejectedStep(Real h, Real error) const {
			return rejectedStep(h, error, _config.alpha);
		}

		Real rejectedStep(Real h, Real error, Real alpha) const {
			if (alpha <= 0)
				throw ODESolverError("Step-size controller exponent must be positive");
			const Real boundedError = std::max(std::abs(error), _config.errorFloor);
			const Real factor = boundedFactor(_config.safety * std::pow(boundedError, -alpha));
			return h * factor;
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
			if (stepSize == 0 || !std::isfinite(stepSize))
				throw ODESolverError("Hermite interpolation requires a finite, non-zero step size");
			if (xStart.size() != xEnd.size() || xStart.size() != dxStart.size() || xStart.size() != dxEnd.size())
				throw ODESolverError("Hermite interpolation data dimensions do not match");

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
			if (stepSize == 0 || !std::isfinite(stepSize))
				throw ODESolverError("Dormand-Prince interpolation requires a finite, non-zero step size");
			if (stages.size() != 7 || xStart.size() != xEnd.size())
				throw ODESolverError("Dormand-Prince interpolation data dimensions do not match");
			for (const auto& stage : stages) {
				if (stage.size() != xStart.size())
					throw ODESolverError("Dormand-Prince interpolation data dimensions do not match");
			}

			// Rows correspond to k1...k7 and columns to theta...theta^4.
			static constexpr Real weights[7][4] = {
				{1.0, -8048581381.0 / 2820520608.0, 8663915743.0 / 2820520608.0, -12715105075.0 / 11282082432.0},
				{0.0, 0.0, 0.0, 0.0},
				{0.0, 131558114200.0 / 32700410799.0, -68118460800.0 / 10900136933.0, 87487479700.0 / 32700410799.0},
				{0.0, -1754552775.0 / 470086768.0, 14199869525.0 / 1410260304.0, -10690763975.0 / 1880347072.0},
				{0.0, 127303824393.0 / 49829197408.0, -318862633887.0 / 49829197408.0, 701980252875.0 / 199316789632.0},
				{0.0, -282668133.0 / 205662961.0, 2019193451.0 / 616988883.0, -1453857185.0 / 822651844.0},
				{0.0, 40617522.0 / 29380423.0, -110615467.0 / 29380423.0, 69997945.0 / 29380423.0}
			};

			_tStart = tStart;
			_stepSize = stepSize;
			_xStart = xStart;
			_xEnd = xEnd;
			for (int power = 0; power < 4; ++power) {
				_coefficients[power].Resize(xStart.size());
				for (int component = 0; component < xStart.size(); ++component) {
					Real coefficient = 0;
					for (int stage = 0; stage < 7; ++stage)
						coefficient += weights[stage][power] * stages[stage][component];
					_coefficients[power][component] = coefficient;
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
			if (stepSize == 0 || !std::isfinite(stepSize))
				throw ODESolverError("DOP853 interpolation requires a finite, non-zero step size");
			if (stages.size() != DP8::extended_stages || xStart.size() != xEnd.size())
				throw ODESolverError("DOP853 interpolation data dimensions do not match");

			const int dimension = system.getDim();
			for (auto& stage : stages) {
				if (stage.size() != dimension)
					stage.Resize(dimension);
			}
			if (workspace.size() != dimension)
				workspace.Resize(dimension);

			// Stages 13 through 15 are used only by the continuous extension.
			for (int stage = DP8::stages + 1; stage < DP8::extended_stages; ++stage) {
				for (int component = 0; component < dimension; ++component) {
					Real increment = 0;
					for (int previousStage = 0; previousStage < stage; ++previousStage)
						increment += DP8::a[stage][previousStage] * stages[previousStage][component];
					workspace[component] = xStart[component] + stepSize * increment;
				}
				system.derivs(tStart + DP8::c[stage] * stepSize, workspace, stages[stage]);
			}

			_tStart = tStart;
			_stepSize = stepSize;
			_xStart = xStart;
			_xEnd = xEnd;
			for (auto& coefficient : _coefficients)
				coefficient.Resize(dimension);

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
			if (stages.size() < Tableau::stages)
				stages.resize(Tableau::stages);
			for (auto& stage : stages) {
				if (stage.size() != dimension)
					stage.Resize(dimension);
			}
			if (workspace.size() != dimension)
				workspace.Resize(dimension);

			stages[0] = initialDerivative;
			// Explicit tableaux are strictly lower triangular: stage i uses only k[0..i-1].
			for (int stage = 1; stage < Tableau::stages; ++stage) {
				for (int component = 0; component < dimension; ++component) {
					Real increment = 0;
					for (int previousStage = 0; previousStage < stage; ++previousStage)
						increment += Tableau::stageCoefficient(stage, previousStage) * stages[previousStage][component];
					workspace[component] = x[component] + h * increment;
				}
				system.derivs(t + Tableau::node(stage) * h, workspace, stages[stage]);
			}
			return Tableau::stages - 1;
		}

		static void combineSolution(const Vector<Real>& x, Real h, const std::vector<Vector<Real>>& stages,
								Vector<Real>& result) {
			if (result.size() != x.size())
				result.Resize(x.size());
			for (int component = 0; component < x.size(); ++component) {
				Real increment = 0;
				for (int stage = 0; stage < Tableau::stages; ++stage)
					increment += Tableau::solutionWeight(stage) * stages[stage][component];
				result[component] = x[component] + h * increment;
			}
		}

		static void combineError(Real h, const std::vector<Vector<Real>>& stages, Vector<Real>& error) {
			const int dimension = stages[0].size();
			if (error.size() != dimension)
				error.Resize(dimension);
			for (int component = 0; component < dimension; ++component) {
				Real estimate = 0;
				// errorWeight is the difference between the high- and low-order weights.
				for (int stage = 0; stage < Tableau::stages; ++stage)
					estimate += Tableau::errorWeight(stage) * stages[stage][component];
				error[component] = h * estimate;
			}
		}
	};

} // namespace MML

#endif // MML_ODE_STEPPER_INFRASTRUCTURE_H