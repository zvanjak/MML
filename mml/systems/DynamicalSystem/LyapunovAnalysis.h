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
#include <mml/algorithms/ODESolvers/ODESolverFixedStep.h>
#include <mml/algorithms/ODESolvers/ODEStepCalculators.h>
#include <mml/base/Matrix/Matrix.h>
#include <mml/interfaces/IODESystem.h>

#include <algorithm>
#include <cmath>
#include <functional>
#include <limits>
#include <type_traits>
#include <vector>

namespace MML::Systems
{
	namespace DynamicalAnalysisDetail
	{
		class VariationalSystem final : public IODESystem
		{
			IDynamicalSystem& _system;
			int _dimension;

		public:
			explicit VariationalSystem(IDynamicalSystem& system)
				: _system(system), _dimension(system.getDim()) {}

			int getDim() const override { return _dimension + _dimension * _dimension; }

			void derivs(Real t, const Vector<Real>& state, Vector<Real>& derivative) const override
			{
				Vector<Real> trajectory(_dimension);
				for (int row = 0; row < _dimension; ++row)
					trajectory[row] = state[row];

				Vector<Real> trajectoryDerivative(_dimension);
				_system.derivs(t, trajectory, trajectoryDerivative);
				for (int row = 0; row < _dimension; ++row)
					derivative[row] = trajectoryDerivative[row];

				Matrix<Real> jacobian(_dimension, _dimension);
				_system.jacobian(t, trajectory, jacobian);
				for (int row = 0; row < _dimension; ++row) {
					for (int column = 0; column < _dimension; ++column) {
						Real value = 0;
						for (int inner = 0; inner < _dimension; ++inner)
							value += jacobian(row, inner) * state[_dimension + inner * _dimension + column];
						derivative[_dimension + row * _dimension + column] = value;
					}
				}
			}
		};
	}

	class LyapunovAnalyzer
	{
	public:
		static LyapunovResult<Real> Compute(IDynamicalSystem& system, const Vector<Real>& initialState,
			Real totalTime, Real orthonormalizationInterval = REAL(1.0), Real initialStep = REAL(0.01))
		{
			const int dimension = system.getDim();
			LyapunovResult<Real> result;
			result.exponents.Resize(dimension);
			Vector<Real> state = initialState;
			Matrix<Real> frame = Matrix<Real>::Identity(dimension);
			Vector<Real> sums(dimension, REAL(0.0));
			Real time = REAL(0.0);
			int orthonormalizations = 0;

			while (time < totalTime) {
				const Real endTime = std::min(time + orthonormalizationInterval, totalTime);
				IntegrateWithVariational(system, state, frame, time, endTime, initialStep);
				GramSchmidtQR(frame, sums);
				++orthonormalizations;
				time = endTime;
			}

			std::vector<Real> sorted(dimension);
			for (int index = 0; index < dimension; ++index)
				sorted[index] = sums[index] / totalTime;
			std::sort(sorted.begin(), sorted.end(), std::greater<Real>());
			for (int index = 0; index < dimension; ++index)
				result.exponents[index] = sorted[index];

			result.maxExponent = result.exponents[0];
			result.sum = REAL(0.0);
			for (int index = 0; index < dimension; ++index)
				result.sum += result.exponents[index];
			result.isChaotic = result.maxExponent > Precision::LyapunovChaosThreshold;
			result.kaplanYorkeDimension = ComputeKaplanYorkeDimension(result.exponents);
			result.numOrthonormalizations = orthonormalizations;
			result.totalTime = totalTime;
			return result;
		}

	private:
		static void IntegrateWithVariational(IDynamicalSystem& system, Vector<Real>& state,
			Matrix<Real>& frame, Real startTime, Real endTime, Real initialStep)
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

			DynamicalAnalysisDetail::VariationalSystem variational(system);
			const auto solution = [&]() {
				if constexpr (std::is_same_v<Real, float>) {
					const int numSteps = std::max(1, static_cast<int>(std::ceil((endTime - startTime) / initialStep)));
					ODESystemFixedStepSolver solver(variational, StepCalculators::RK4_Basic);
					return solver.integrate(augmented, startTime, endTime, numSteps);
				}
				else {
					ODEAdaptiveIntegrator<> integrator(variational);
					const Real tolerance = std::max(Precision::ODEDefaultTolerance,
						Real(1000000) * std::numeric_limits<Real>::epsilon());
					return integrator.integrate(augmented, startTime, endTime, endTime - startTime,
						tolerance, std::min(initialStep, endTime - startTime));
				}
			}();
			const Vector<Real> finalState = solution.getXValuesAtEnd();
			for (int row = 0; row < dimension; ++row) {
				state[row] = finalState[row];
				for (int column = 0; column < dimension; ++column)
					frame(row, column) = finalState[dimension + row * dimension + column];
			}
		}

		static void GramSchmidtQR(Matrix<Real>& frame, Vector<Real>& sums)
		{
			const int dimension = frame.rows();
			for (int column = 0; column < dimension; ++column) {
				Vector<Real> vector(dimension);
				for (int row = 0; row < dimension; ++row)
					vector[row] = frame(row, column);

				for (int previous = 0; previous < column; ++previous) {
					Real dot = REAL(0.0);
					for (int row = 0; row < dimension; ++row)
						dot += vector[row] * frame(row, previous);
					for (int row = 0; row < dimension; ++row)
						vector[row] -= dot * frame(row, previous);
				}

				Real norm = REAL(0.0);
				for (int row = 0; row < dimension; ++row)
					norm += vector[row] * vector[row];
				norm = std::sqrt(norm);
				if (norm > Precision::DivisionSafetyThreshold) {
					sums[column] += std::log(norm);
					for (int row = 0; row < dimension; ++row)
						frame(row, column) = vector[row] / norm;
				}
				else {
					for (int row = 0; row < dimension; ++row)
						frame(row, column) = row == column ? REAL(1.0) : REAL(0.0);
				}
			}
		}

		static Real ComputeKaplanYorkeDimension(const Vector<Real>& exponents)
		{
			const int dimension = exponents.size();
			if (dimension == 0)
				return REAL(0.0);
			Real cumulative = REAL(0.0);
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
				return REAL(0.0);
			Real partialSum = REAL(0.0);
			for (int current = 0; current < index; ++current)
				partialSum += exponents[current];
			return index + partialSum / std::abs(exponents[index]);
		}
	};

	template<typename Type = Real>
	struct LyapunovAnalysisResult : public EvaluationResultBase
	{
		LyapunovResult<Type> lyapunov;
	};

	inline LyapunovAnalysisResult<Real> ComputeLyapunovDetailed(IDynamicalSystem& system,
		const Vector<Real>& initialState, Real totalTime, Real orthonormalizationInterval = REAL(1.0),
		Real initialStep = REAL(0.01), const DynSysConfig& config = {})
	{
		return DynSysDetail::ExecuteDynSysDetailed<LyapunovAnalysisResult<Real>>(
			"LyapunovAnalyzer", config, [&](LyapunovAnalysisResult<Real>& result) {
				result.lyapunov = LyapunovAnalyzer::Compute(system, initialState, totalTime,
					orthonormalizationInterval, initialStep);
				result.function_evaluations = result.lyapunov.numOrthonormalizations;
			});
	}
}

#endif // MML_LYAPUNOV_ANALYSIS_H