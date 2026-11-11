///////////////////////////////////////////////////////////////////////////////////////////
///                         MinimalMathLibrary (MML)                                  ///
///                                                                                   ///
///  File:        LyapunovAnalysis.h                                                 ///
///  Description: Lyapunov exponents and attractor-dimension analysis                ///
///                                                                                   ///
///  Copyright:   (c) 2024-2026 Zvonimir Vanjak                                       ///
///  License:     MIT License (see LICENSE.md)                                         ///
///                                                                                   ///
///////////////////////////////////////////////////////////////////////////////////////////
#if !defined MML_LYAPUNOV_ANALYSIS_H
#define MML_LYAPUNOV_ANALYSIS_H

#include <mml/systems/DynamicalSystem/DynamicalAnalysisCommon.h>
#include <mml/systems/DynamicalSystem/DynamicalSystemTypes.h>
#include <mml/systems/DynamicalSystem/LyapunovQR.h>
#include <mml/algorithms/ODESolvers/ODESolverFixedStep.h>
#include <mml/algorithms/ODESolvers/ODEStepCalculators.h>
#include <mml/base/Matrix/Matrix.h>
#include <mml/interfaces/IODESystem.h>

#include <algorithm>
#include <cmath>
#include <functional>
#include <limits>
#include <type_traits>
#include <utility>
#include <vector>

namespace MML::Systems
{
	namespace DynamicalAnalysisDetail
	{
		class TangentFrameResolutionError : public NumericalMethodError
		{
		public:
			using NumericalMethodError::NumericalMethodError;
		};

		class VariationalSystem final : public IODESystem
		{
			IDynamicalSystem& _system;
			int _dimension;
			int& _functionEvaluations;

		public:
			VariationalSystem(IDynamicalSystem& system, int& functionEvaluations)
				: _system(system), _dimension(system.getDim()), _functionEvaluations(functionEvaluations) {}

			int getDim() const override { return _dimension + _dimension * _dimension; }

			void derivs(Real t, const Vector<Real>& state, Vector<Real>& derivative) const override
			{
				Vector<Real> trajectory(_dimension);
				for (int row = 0; row < _dimension; ++row)
					trajectory[row] = state[row];
				DynSysDetail::CheckState(trajectory, _dimension, "Continuous Lyapunov trajectory");

				if (_functionEvaluations == std::numeric_limits<int>::max())
					throw ResourceLimitError("Lyapunov: RHS evaluation counter exhausted");
				++_functionEvaluations;
				Vector<Real> trajectoryDerivative(_dimension, std::numeric_limits<Real>::quiet_NaN());
				_system.derivs(t, trajectory, trajectoryDerivative);
				DynSysDetail::CheckState(trajectoryDerivative, _dimension, "Continuous Lyapunov RHS output");
				for (int row = 0; row < _dimension; ++row)
					derivative[row] = trajectoryDerivative[row];

				Matrix<Real> jacobian(_dimension, _dimension, std::numeric_limits<Real>::quiet_NaN());
				_system.jacobian(t, trajectory, jacobian);
				DynSysDetail::CheckJacobian(jacobian, _dimension);
				for (int row = 0; row < _dimension; ++row) {
					for (int column = 0; column < _dimension; ++column) {
						Real value = 0;
						for (int inner = 0; inner < _dimension; ++inner)
							value += jacobian(row, inner) * state[_dimension + inner * _dimension + column];
						if (!std::isfinite(value))
							throw TangentFrameResolutionError("Lyapunov: non-finite variational derivative; shorten the QR block");
						derivative[_dimension + row * _dimension + column] = value;
					}
				}
			}
		};
	class LyapunovComputation
	{
	public:
		static LyapunovResult<Real> Compute(IDynamicalSystem& system, const Vector<Real>& initialState,
			Real totalTime, Real orthonormalizationInterval, Real initialStep, int& functionEvaluations)
		{
			const int dimension = system.getDim();
			DynSysDetail::CheckVariationalDimension(dimension);
			if (!std::isfinite(totalTime) || totalTime <= 0 ||
				!std::isfinite(orthonormalizationInterval) || orthonormalizationInterval <= 0 ||
				!std::isfinite(initialStep) || initialStep <= 0)
				throw ArgumentError("Lyapunov: total time, QR interval and initial step must be finite and positive");
			DynSysDetail::CheckState<ArgumentError>(initialState, dimension, "Continuous Lyapunov initial state");
			LyapunovResult<Real> result;
			result.exponents.resize(dimension);
			Vector<Real> state = initialState;
			Matrix<Real> frame = Matrix<Real>::Identity(dimension);
			Vector<Real> sums(dimension, static_cast<MML::Real>(0.0));
			Real time = static_cast<MML::Real>(0.0);
			int orthonormalizations = 0;
			Real blockInterval = std::min(orthonormalizationInterval, totalTime);
			// Keep accepted tangent stretches above the IVP absolute-error floor.
			// This is a resolution safeguard, not a Lyapunov-exponent error estimate.
			const Real logStretchLimit = -std::log(std::sqrt(IntegrationTolerance()));

			while (time < totalTime) {
				int refinements = 0;
				while (true) {
					const Real endTime = time + std::min(blockInterval, totalTime - time);
					if (!(endTime > time))
						throw NumericalMethodError("Lyapunov: QR block cannot advance floating-point time");
					Vector<Real> trialState = state;
					Matrix<Real> trialFrame = frame;
					bool resolved = true;
					try {
						IntegrateWithVariational(system, trialState, trialFrame, time, endTime, initialStep, functionEvaluations);
					}
					catch (const TangentFrameResolutionError&) {
						resolved = false;
					}
					Vector<Real> logarithms;
					if (resolved) {
						try {
							logarithms = LyapunovDetail::Orthonormalize(trialFrame);
						}
						catch (const SingularMatrixError&) {
							resolved = false;
						}
					}
					if (resolved)
						for (Real logarithm : logarithms)
							if (std::abs(logarithm) > logStretchLimit)
								resolved = false;
					if (resolved) {
						if (orthonormalizations == std::numeric_limits<int>::max())
							throw ResourceLimitError("Lyapunov: QR counter exhausted");
						for (int column = 0; column < dimension; ++column)
							sums[column] += logarithms[column];
						state = std::move(trialState);
						frame = std::move(trialFrame);
						++orthonormalizations;
						time = endTime;
						break;
					}
					if (refinements == Precision::LyapunovMaxQRRefinements)
						throw NumericalMethodError("Lyapunov: tangent frame unresolved after QR refinement budget; shorten the requested interval");
					if (result.qr_refinements == std::numeric_limits<int>::max())
						throw ResourceLimitError("Lyapunov: QR refinement counter exhausted");
					++refinements;
					++result.qr_refinements;
					blockInterval = (endTime - time) / 2;
				}
			}

			std::vector<Real> sorted(dimension);
			for (int index = 0; index < dimension; ++index)
				sorted[index] = sums[index] / totalTime;
			std::sort(sorted.begin(), sorted.end(), std::greater<Real>());
			for (int index = 0; index < dimension; ++index)
				result.exponents[index] = sorted[index];

			result.maxExponent = result.exponents[0];
			result.sum = static_cast<MML::Real>(0.0);
			for (int index = 0; index < dimension; ++index)
				result.sum += result.exponents[index];
			result.isChaotic = result.maxExponent > Precision::LyapunovChaosThreshold;
			result.kaplanYorkeDimension = ComputeKaplanYorkeDimension(result.exponents);
			result.numOrthonormalizations = orthonormalizations;
			result.totalTime = totalTime;
			result.function_evaluations = functionEvaluations;
			if (!std::isfinite(result.maxExponent) || !std::isfinite(result.sum) ||
				!std::isfinite(result.kaplanYorkeDimension))
				throw NumericalMethodError("Lyapunov: non-finite spectrum or derived quantity");
			return result;
		}

	private:
		static Real IntegrationTolerance()
		{
			return std::max(Precision::ODEDefaultTolerance,
				Precision::LyapunovToleranceEpsilonFactor * std::numeric_limits<Real>::epsilon());
		}

		static void IntegrateWithVariational(IDynamicalSystem& system, Vector<Real>& state,
			Matrix<Real>& frame, Real startTime, Real endTime, Real initialStep, int& functionEvaluations)
		{
			if (endTime < startTime)
				throw ArgumentError("IntegrateWithVariational: t1 must be >= t0 (reverse-time integration is not supported)");
			if (initialStep <= 0)
				throw ArgumentError("IntegrateWithVariational: step size h must be positive");

			const int dimension = system.getDim();
			Vector<Real> augmented(dimension + dimension * dimension);
			for (int row = 0; row < dimension; ++row) {
				augmented[row] = state[row];
				for (int column = 0; column < dimension; ++column)
					augmented[dimension + row * dimension + column] = frame(row, column);
			}

			VariationalSystem variational(system, functionEvaluations);
			const auto solution = [&]() {
				if constexpr (std::is_same_v<Real, float>) {
					const Real requestedSteps = std::ceil((endTime - startTime) / initialStep);
					if (!std::isfinite(requestedSteps) ||
						static_cast<long double>(requestedSteps) > std::numeric_limits<int>::max())
						throw ResourceLimitError("Lyapunov: fixed-step count exceeds the integer range");
					const int numSteps = std::max(1, static_cast<int>(requestedSteps));
					ODESystemFixedStepSolver solver(variational, StepCalculators::RK4_Basic);
					return solver.integrate(augmented, startTime, endTime, numSteps);
				}
				else {
					ODEAdaptiveIntegrator<> integrator(variational);
					return integrator.integrate(augmented, startTime, endTime, endTime - startTime,
						IntegrationTolerance(), std::min(initialStep, endTime - startTime));
				}
			}();
			const Vector<Real> finalState = solution.getXValuesAtEnd();
			for (int row = 0; row < dimension; ++row) {
				state[row] = finalState[row];
				for (int column = 0; column < dimension; ++column)
					frame(row, column) = finalState[dimension + row * dimension + column];
			}
			DynSysDetail::CheckState(state, dimension, "Continuous Lyapunov final trajectory");
		}

		static Real ComputeKaplanYorkeDimension(const Vector<Real>& exponents)
		{
			const int dimension = exponents.size();
			if (dimension == 0)
				return static_cast<MML::Real>(0.0);
			Real cumulative = static_cast<MML::Real>(0.0);
			int index = 0;
			for (int current = 0; current < dimension; ++current) {
				cumulative += exponents[current];
				if (cumulative >= 0)
					index = current + 1;
				else
					break;
			}
			if (index == dimension)
				return static_cast<Real>(dimension);
			if (index == 0)
				return static_cast<MML::Real>(0.0);
			Real partialSum = static_cast<MML::Real>(0.0);
			for (int current = 0; current < index; ++current)
				partialSum += exponents[current];
			return index + partialSum / std::abs(exponents[index]);
		}
	};
	}

	class LyapunovAnalyzer
	{
	public:
		/// Finite-time spectrum from the full variational flow; the requested QR interval
		/// is a maximum. Throws if bounded subdivision cannot resolve a finite full frame.
		static LyapunovResult<Real> Compute(IDynamicalSystem& system, const Vector<Real>& initialState,
			Real totalTime, Real orthonormalizationInterval = Real(1), Real initialStep = Real(0.01))
		{
			int functionEvaluations = 0;
			return DynamicalAnalysisDetail::LyapunovComputation::Compute(system, initialState, totalTime,
				orthonormalizationInterval, initialStep, functionEvaluations);
		}
	};

	template<typename Type = Real>
	struct LyapunovAnalysisResult : public EvaluationResultBase
	{
		LyapunovResult<Type> lyapunov;
	};

	inline LyapunovAnalysisResult<Real> ComputeLyapunovDetailed(IDynamicalSystem& system,
		const Vector<Real>& initialState, Real totalTime, Real orthonormalizationInterval = static_cast<MML::Real>(1.0),
		Real initialStep = static_cast<MML::Real>(0.01), const DynSysConfig& config = {})
	{
		return DynSysDetail::ExecuteDynSysDetailed<LyapunovAnalysisResult<Real>>(
			"LyapunovAnalyzer", config, [&](LyapunovAnalysisResult<Real>& result) {
				result.lyapunov = DynamicalAnalysisDetail::LyapunovComputation::Compute(system, initialState, totalTime,
					orthonormalizationInterval, initialStep, result.function_evaluations);
			});
	}
}

#endif // MML_LYAPUNOV_ANALYSIS_H
