///////////////////////////////////////////////////////////////////////////////////////////
///                         MinimalMathLibrary (MML)                                  ///
///                                                                                   ///
///  File:        Geodesic.h                                                          ///
///  Description: Geodesic equation ODE adapter and fixed-step integration helpers    ///
///                                                                                   ///
///  Copyright:   (c) 2024-2026 Zvonimir Vanjak                                       ///
///  License:     MIT License (see LICENSE.md)                                         ///
///                                                                                   ///
///////////////////////////////////////////////////////////////////////////////////////////
#if !defined MML_GEODESIC_H
#define MML_GEODESIC_H

#include <mml/base/ODESystem.h>
#include <mml/base/ODESystemSolution.h>

#include <mml/algorithms/ODESolvers/ODESolverFixedStep.h>
#include <mml/algorithms/ODESolvers/ODEStepCalculators.h>

#include <mml/core/DifferentialGeometry/InducedMetric.h>
#include <mml/core/MetricTensor.h>

namespace MML
{
	template<int N>
	Vector<Real> MakeGeodesicState(const VectorN<Real, N>& position, const VectorN<Real, N>& velocity)
	{
		Vector<Real> state(2 * N);
		for (int i = 0; i < N; i++)
		{
			state[i] = position[i];
			state[N + i] = velocity[i];
		}
		return state;
	}

	template<int N>
	VectorN<Real, N> GeodesicPositionFromState(const Vector<Real>& state)
	{
		VectorN<Real, N> position;
		for (int i = 0; i < N; i++)
			position[i] = state[i];
		return position;
	}

	template<int N>
	VectorN<Real, N> GeodesicVelocityFromState(const Vector<Real>& state)
	{
		VectorN<Real, N> velocity;
		for (int i = 0; i < N; i++)
			velocity[i] = state[N + i];
		return velocity;
	}

	/// @brief First-order ODE system for the geodesic equation.
	/// @details State layout is (q^0, ..., q^(N-1), v^0, ..., v^(N-1)), with
	///          dq^i/dlambda = v^i and dv^i/dlambda = -Gamma^i_jk v^j v^k.
	template<int N>
	class GeodesicEquationSystem : public IODESystem
	{
		const MetricTensorField<N>& _metric;

	public:
		GeodesicEquationSystem(const MetricTensorField<N>& metric) : _metric(metric) { }

		int getDim() const override { return 2 * N; }

		void derivs(const Real lambda, const Vector<Real>& state, Vector<Real>& dstate) const override
		{
			(void)lambda;
			VectorN<Real, N> position = GeodesicPositionFromState<N>(state);

			for (int i = 0; i < N; i++)
				dstate[i] = state[N + i];

			for (int i = 0; i < N; i++)
			{
				Real acceleration = REAL(0.0);
				for (int j = 0; j < N; j++)
					for (int k = 0; k < N; k++)
						acceleration -= _metric.GetChristoffelSymbolSecondKind(i, j, k, position) * state[N + j] * state[N + k];

				dstate[N + i] = acceleration;
			}
		}

		std::string getVarName(int ind) const override
		{
			if (ind < 0 || ind >= 2 * N)
				return IODESystem::getVarName(ind);
			return ind < N ? "q" + std::to_string(ind) : "v" + std::to_string(ind - N);
		}
	};

	template<int N>
	ODESystemSolution IntegrateGeodesicFixedStep(const MetricTensorField<N>& metric,
		const VectorN<Real, N>& initialPosition,
		const VectorN<Real, N>& initialVelocity,
		Real lambdaStart,
		Real lambdaEnd,
		int numSteps,
		const IODESystemStepCalculator& stepCalculator = StepCalculators::RK4_Basic)
	{
		GeodesicEquationSystem<N> system(metric);
		ODESystemFixedStepSolver solver(system, stepCalculator);
		return solver.integrate(MakeGeodesicState<N>(initialPosition, initialVelocity), lambdaStart, lambdaEnd, numSteps);
	}

	/// @brief ODE system for parallel transport of a contravariant vector along a curve.
	/// @details State layout is (V^0, ..., V^(N-1)), with dV^i/dlambda = -Gamma^i_jk qdot^j V^k.
	template<int N>
	class ParallelTransportEquationSystem : public IODESystem
	{
		const MetricTensorField<N>& _metric;
		const IParametricCurve<N>& _curve;
		Real _curveDerivativeStep;

	public:
		ParallelTransportEquationSystem(const MetricTensorField<N>& metric,
			const IParametricCurve<N>& curve,
			Real curveDerivativeStep = PrecisionValues<Real>::DerivativeStepSize)
			: _metric(metric), _curve(curve), _curveDerivativeStep(curveDerivativeStep) { }

		int getDim() const override { return N; }

		void derivs(const Real lambda, const Vector<Real>& vector, Vector<Real>& dvector) const override
		{
			VectorN<Real, N> position = _curve(lambda);
			VectorN<Real, N> tangent = Derivation::NDer4(_curve, lambda, _curveDerivativeStep);

			for (int i = 0; i < N; i++)
			{
				Real derivative = REAL(0.0);
				for (int j = 0; j < N; j++)
					for (int k = 0; k < N; k++)
						derivative -= _metric.GetChristoffelSymbolSecondKind(i, j, k, position) * tangent[j] * vector[k];

				dvector[i] = derivative;
			}
		}

		std::string getVarName(int ind) const override
		{
			if (ind < 0 || ind >= N)
				return IODESystem::getVarName(ind);
			return "V" + std::to_string(ind);
		}
	};

	template<int N>
	Vector<Real> MakeParallelTransportState(const VectorN<Real, N>& vector)
	{
		Vector<Real> state(N);
		for (int i = 0; i < N; i++)
			state[i] = vector[i];
		return state;
	}

	template<int N>
	VectorN<Real, N> ParallelTransportVectorFromState(const Vector<Real>& state)
	{
		VectorN<Real, N> vector;
		for (int i = 0; i < N; i++)
			vector[i] = state[i];
		return vector;
	}

	template<int N>
	ODESystemSolution IntegrateParallelTransportFixedStep(const MetricTensorField<N>& metric,
		const IParametricCurve<N>& curve,
		const VectorN<Real, N>& initialVector,
		Real lambdaStart,
		Real lambdaEnd,
		int numSteps,
		const IODESystemStepCalculator& stepCalculator = StepCalculators::RK4_Basic,
		Real curveDerivativeStep = PrecisionValues<Real>::DerivativeStepSize)
	{
		ParallelTransportEquationSystem<N> system(metric, curve, curveDerivativeStep);
		ODESystemFixedStepSolver solver(system, stepCalculator);
		return solver.integrate(MakeParallelTransportState<N>(initialVector), lambdaStart, lambdaEnd, numSteps);
	}
}

namespace MML::DifferentialGeometry
{
	inline ODESystemSolution IntegrateSurfaceGeodesicFixedStep(const Surfaces::ISurfaceCartesian& surface,
		const VectorN<Real, 2>& initialPosition,
		const VectorN<Real, 2>& initialVelocity,
		Real lambdaStart,
		Real lambdaEnd,
		int numSteps,
		const IODESystemStepCalculator& stepCalculator = StepCalculators::RK4_Basic)
	{
		InducedMetric2D metric(surface);
		return IntegrateGeodesicFixedStep<2>(metric, initialPosition, initialVelocity, lambdaStart, lambdaEnd, numSteps, stepCalculator);
	}
}

#endif // MML_GEODESIC_H