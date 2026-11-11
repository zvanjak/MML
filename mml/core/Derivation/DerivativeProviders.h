///////////////////////////////////////////////////////////////////////////////////////////
///                         MinimalMathLibrary (MML)                                  ///
///  File:        DerivativeProviders.h                                               ///
///  Description: Optional Real-valued derivative capabilities and dispatch            ///
///  Copyright:   (c) 2024-2026 Zvonimir Vanjak                                         ///
///  License:     MIT License (see LICENSE.md)                                           ///
///////////////////////////////////////////////////////////////////////////////////////////
#if !defined MML_DERIVATIVE_PROVIDERS_H
#define MML_DERIVATIVE_PROVIDERS_H

#include <mml/core/Derivation.h>
#include <mml/base/Matrix/MatrixNM.h>
#include <array>

namespace MML {
	template<int N>
	class ICurveDerivativeProvider {
	public:
		virtual ~ICurveDerivativeProvider() = default;
		virtual VectorN<Real, N> EvaluateDerivative(Real t, int order) const = 0;
	};

	template<int N>
	class IScalarDerivativeProvider {
	public:
		virtual ~IScalarDerivativeProvider() = default;
		virtual Real EvaluatePartial(const VectorN<Real, N>& pos, int axis) const = 0;
		virtual Real EvaluateSecondPartial(const VectorN<Real, N>& pos, int axis1, int axis2) const = 0;
	};

	template<int N>
	class IVectorDerivativeProvider {
	public:
		virtual ~IVectorDerivativeProvider() = default;
		virtual Real EvaluateComponentPartial(const VectorN<Real, N>& pos, int component, int axis) const = 0;
	};

	template<int N>
	struct MetricJet {
		MatrixNM<Real, N, N> value;
		std::array<MatrixNM<Real, N, N>, N> first;
		// second[a][b](i,j) = partial_a partial_b g_ij.
		std::array<std::array<MatrixNM<Real, N, N>, N>, N> second;
	};

	template<int N>
	class IMetricJetProvider {
	public:
		virtual ~IMetricJetProvider() = default;
		virtual MetricJet<N> EvaluateMetricJet(const VectorN<Real, N>& pos) const = 0;
	};
} // namespace MML

namespace MML::Differential {
	template<int N>
	VectorN<Real, N> CurveDerivative(const IParametricCurve<N>& curve, Real t, int order) {
		if (order < 1 || order > 3)
			throw DomainError("CurveDerivative requires order 1, 2 or 3");
		if (const auto* provider = dynamic_cast<const ICurveDerivativeProvider<N>*>(&curve))
			return provider->EvaluateDerivative(t, order);
		if (order == 1)
			return Derivation::DeriveCurve<N>(curve, t, nullptr);
		if (order == 2)
			return Derivation::DeriveCurveSec<N>(curve, t, nullptr);
		return Derivation::DeriveCurveThird<N>(curve, t, nullptr);
	}

	template<int N>
	Real ScalarPartial(const IScalarFunction<N>& f, int axis, const VectorN<Real, N>& pos) {
		if (const auto* provider = dynamic_cast<const IScalarDerivativeProvider<N>*>(&f))
			return provider->EvaluatePartial(pos, axis);
		return Derivation::DerivePartial<N>(f, axis, pos, nullptr);
	}

	template<int N>
	VectorN<Real, N> ScalarGradient(const IScalarFunction<N>& f, const VectorN<Real, N>& pos) {
		if (const auto* provider = dynamic_cast<const IScalarDerivativeProvider<N>*>(&f)) {
			VectorN<Real, N> result;
			for (int i = 0; i < N; ++i)
				result[i] = provider->EvaluatePartial(pos, i);
			return result;
		}
		return Derivation::DerivePartialAll<N>(f, pos, nullptr);
	}

	template<int N>
	Real ScalarSecondPartial(const IScalarFunction<N>& f, int a, int b, const VectorN<Real, N>& pos) {
		if (const auto* provider = dynamic_cast<const IScalarDerivativeProvider<N>*>(&f))
			return provider->EvaluateSecondPartial(pos, a, b);
		return Derivation::DeriveSecPartial<N>(f, a, b, pos, nullptr);
	}

	template<int N>
	Real VectorPartial(const IVectorFunction<N>& f, int component, int axis, const VectorN<Real, N>& pos) {
		if (const auto* provider = dynamic_cast<const IVectorDerivativeProvider<N>*>(&f))
			return provider->EvaluateComponentPartial(pos, component, axis);
		return Derivation::DeriveVecPartial<N>(f, component, axis, pos, nullptr);
	}
} // namespace MML::Differential
#endif
