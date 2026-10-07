///////////////////////////////////////////////////////////////////////////////////////////
///                         MinimalMathLibrary (MML)                                  ///
///                                                                                   ///
///  File:        DifferentialGeometry/Frames.h                                       ///
///  Description: Frenet and Darboux frame helpers                                    ///
///                                                                                   ///
///  Copyright:   (c) 2024-2026 Zvonimir Vanjak                                       ///
///  License:     MIT License (see LICENSE.md)                                         ///
///                                                                                   ///
///////////////////////////////////////////////////////////////////////////////////////////
#if !defined MML_CORE_DIFFERENTIAL_GEOMETRY_FRAMES_H
#define MML_CORE_DIFFERENTIAL_GEOMETRY_FRAMES_H

#include <mml/MMLBase.h>

#include <mml/core/Derivation.h>
#include <mml/core/Surfaces.h>

namespace MML::DifferentialGeometry
{
	template<int N>
	struct FrenetFrame
	{
		VectorN<Real, N> tangent;
		VectorN<Real, N> normal;
		Real speed = REAL(0.0);
		Real curvature = REAL(0.0);
	};

	template<int N>
	FrenetFrame<N> ComputeFrenetFrame(const IParametricCurve<N>& curve, Real t, Real derivativeStep = REAL(0.0))
	{
		VectorN<Real, N> first = derivativeStep > REAL(0.0)
			? Derivation::NDer4(curve, t, derivativeStep)
			: Derivation::NDer4(curve, t);
		VectorN<Real, N> second = derivativeStep > REAL(0.0)
			? Derivation::NSecDer4(curve, t, derivativeStep)
			: Derivation::NSecDer4(curve, t);

		FrenetFrame<N> frame;
		frame.speed = first.NormL2();
		if (frame.speed <= PrecisionValues<Real>::NumericalZeroThreshold)
			throw DomainError("ComputeFrenetFrame - curve speed is zero");

		frame.tangent = first / frame.speed;
		Real tangentialAcceleration = Utils::ScalarProduct(second, frame.tangent);
		VectorN<Real, N> normalComponent = second - tangentialAcceleration * frame.tangent;
		Real normalMagnitude = normalComponent.NormL2();
		frame.curvature = normalMagnitude / (frame.speed * frame.speed);
		if (normalMagnitude > PrecisionValues<Real>::NumericalZeroThreshold)
			frame.normal = normalComponent / normalMagnitude;

		return frame;
	}

	struct DarbouxFrame
	{
		VectorN<Real, 3> tangent;
		VectorN<Real, 3> surfaceNormal;
		VectorN<Real, 3> tangentNormal;
		VectorN<Real, 3> curvatureVector;
		Real curvature = REAL(0.0);
		Real geodesicCurvature = REAL(0.0);
		Real normalCurvature = REAL(0.0);
	};

	namespace Detail
	{
		inline VectorN<Real, 3> Cross3D(const VectorN<Real, 3>& left, const VectorN<Real, 3>& right)
		{
			return {
				left[1] * right[2] - left[2] * right[1],
				left[2] * right[0] - left[0] * right[2],
				left[0] * right[1] - left[1] * right[0]
			};
		}

		class SurfaceParameterCurve3D : public IParametricCurve<3>
		{
			const Surfaces::ISurfaceCartesian& _surface;
			const IParametricCurve<2>& _parameterCurve;

		public:
			SurfaceParameterCurve3D(const Surfaces::ISurfaceCartesian& surface, const IParametricCurve<2>& parameterCurve)
				: _surface(surface), _parameterCurve(parameterCurve)
			{
			}

			Real getMinT() const override { return _parameterCurve.getMinT(); }
			Real getMaxT() const override { return _parameterCurve.getMaxT(); }

			VectorN<Real, 3> operator()(Real t) const override
			{
				VectorN<Real, 2> params = _parameterCurve(t);
				return _surface(params[0], params[1]);
			}
		};
	}

	inline DarbouxFrame ComputeDarbouxFrame(const Surfaces::ISurfaceCartesian& surface,
		const IParametricCurve<2>& parameterCurve,
		Real t,
		Real derivativeStep = REAL(0.0))
	{
		Detail::SurfaceParameterCurve3D ambientCurve(surface, parameterCurve);
		VectorN<Real, 2> params = parameterCurve(t);

		VectorN<Real, 3> first = derivativeStep > REAL(0.0)
			? Derivation::NDer4(ambientCurve, t, derivativeStep)
			: Derivation::NDer4(ambientCurve, t);
		VectorN<Real, 3> second = derivativeStep > REAL(0.0)
			? Derivation::NSecDer4(ambientCurve, t, derivativeStep)
			: Derivation::NSecDer4(ambientCurve, t);

		DarbouxFrame frame;
		Real speed = first.NormL2();
		if (speed <= PrecisionValues<Real>::NumericalZeroThreshold)
			throw DomainError("ComputeDarbouxFrame - surface curve speed is zero");

		frame.tangent = first / speed;
		frame.surfaceNormal = surface.Normal(params[0], params[1]);
		frame.tangentNormal = Detail::Cross3D(frame.surfaceNormal, frame.tangent);
		Real tangentNormalNorm = frame.tangentNormal.NormL2();
		if (tangentNormalNorm > PrecisionValues<Real>::NumericalZeroThreshold)
			frame.tangentNormal = frame.tangentNormal / tangentNormalNorm;

		Real tangentialAcceleration = Utils::ScalarProduct(second, frame.tangent);
		frame.curvatureVector = (second - tangentialAcceleration * frame.tangent) / (speed * speed);
		frame.curvature = frame.curvatureVector.NormL2();
		frame.normalCurvature = Utils::ScalarProduct(frame.curvatureVector, frame.surfaceNormal);
		frame.geodesicCurvature = Utils::ScalarProduct(frame.curvatureVector, frame.tangentNormal);

		return frame;
	}
} // namespace MML::DifferentialGeometry

#endif // MML_CORE_DIFFERENTIAL_GEOMETRY_FRAMES_H