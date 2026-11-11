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
#include <mml/core/Derivation/DerivativeProviders.h>
#include <mml/core/Surfaces.h>
#include <mml/core/detail/GeometryChecks.h>

namespace MML::DifferentialGeometry
{
	template<int N>
	struct FrenetFrame
	{
		VectorN<Real, N> tangent;
		VectorN<Real, N> normal;
		Real speed = static_cast<MML::Real>(0.0);
		Real curvature = static_cast<MML::Real>(0.0);
		bool normal_defined = false;
	};

	template<int N>
	FrenetFrame<N> ComputeFrenetFrame(const IParametricCurve<N>& curve, Real t, Real derivativeStep = static_cast<MML::Real>(0.0))
	{
		ValidateFinite(t, "ComputeFrenetFrame parameter");
		ValidateNonNegativeTolerance(derivativeStep, "ComputeFrenetFrame derivative step");
		GeometryDetail::RequireRange(curve.getMinT(), curve.getMaxT(), "ComputeFrenetFrame curve domain");
		if (t < curve.getMinT() || t > curve.getMaxT())
			throw DomainError("ComputeFrenetFrame - parameter outside curve domain");
		VectorN<Real, N> first = derivativeStep > static_cast<MML::Real>(0.0)
			? Derivation::NDer4(curve, t, derivativeStep)
			: Differential::CurveDerivative(curve, t, 1);
		VectorN<Real, N> second = derivativeStep > static_cast<MML::Real>(0.0)
			? Derivation::NSecDer4(curve, t, derivativeStep)
			: Differential::CurveDerivative(curve, t, 2);

		FrenetFrame<N> frame;
		frame.speed = GeometryDetail::Norm(first);
		if (frame.speed == 0)
			throw DomainError("ComputeFrenetFrame - curve speed is zero");

		frame.tangent = GeometryDetail::Unit(first);
		const Real acceleration = GeometryDetail::Norm(second);
		if (acceleration == 0)
			return frame;
		const auto direction = second / acceleration;
		const auto normalComponent = direction - Utils::ScalarProduct(direction, frame.tangent) * frame.tangent;
		const Real normalMagnitude = GeometryDetail::Norm(normalComponent);
		if (normalMagnitude > PrecisionValues<Real>::NumericalZeroThreshold) {
			frame.curvature = GeometryDetail::Checked(((normalMagnitude * acceleration) / frame.speed) / frame.speed);
			frame.normal = normalComponent / normalMagnitude;
			frame.normal_defined = true;
		}

		return frame;
	}

	struct DarbouxFrame
	{
		VectorN<Real, 3> tangent;
		VectorN<Real, 3> surfaceNormal;
		VectorN<Real, 3> tangentNormal;
		VectorN<Real, 3> curvatureVector;
		Real curvature = static_cast<MML::Real>(0.0);
		Real geodesicCurvature = static_cast<MML::Real>(0.0);
		Real normalCurvature = static_cast<MML::Real>(0.0);
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
		Real derivativeStep = static_cast<MML::Real>(0.0))
	{
		Detail::SurfaceParameterCurve3D ambientCurve(surface, parameterCurve);
		const auto frenet = ComputeFrenetFrame(ambientCurve, t, derivativeStep);
		VectorN<Real, 2> params = parameterCurve(t);

		DarbouxFrame frame;
		frame.tangent = frenet.tangent;
		frame.surfaceNormal = GeometryDetail::Unit(surface.Normal(params[0], params[1]));
		frame.tangentNormal = GeometryDetail::Unit(Detail::Cross3D(frame.surfaceNormal, frame.tangent));
		frame.curvatureVector = GeometryDetail::Checked(frenet.curvature * frenet.normal);
		frame.curvature = frenet.curvature;
		frame.normalCurvature = GeometryDetail::Checked(Utils::ScalarProduct(frame.curvatureVector, frame.surfaceNormal));
		frame.geodesicCurvature = GeometryDetail::Checked(Utils::ScalarProduct(frame.curvatureVector, frame.tangentNormal));

		return frame;
	}
} // namespace MML::DifferentialGeometry

#endif // MML_CORE_DIFFERENTIAL_GEOMETRY_FRAMES_H
