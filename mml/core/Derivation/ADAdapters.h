///////////////////////////////////////////////////////////////////////////////////////////
///                         MinimalMathLibrary (MML)                                  ///
///  File:        ADAdapters.h                                                        ///
///  Description: Opt-in forward-AD geometry and field adapters                        ///
///  Copyright:   (c) 2024-2026 Zvonimir Vanjak                                         ///
///  License:     MIT License (see LICENSE.md)                                           ///
///////////////////////////////////////////////////////////////////////////////////////////
#if !defined MML_AD_ADAPTERS_H
#define MML_AD_ADAPTERS_H

#include <mml/core/Derivation/ForwardAD.h>
#include <mml/core/Derivation/DerivativeProviders.h>
#include <mml/core/Curves.h>
#include <mml/core/Surfaces.h>
#include <mml/core/MetricTensor.h>
#include <array>
#include <utility>

namespace MML::AD {
	namespace AdapterDetail {
		inline void RequireFinite(Real x) {
			if (!std::isfinite(x))
				throw DomainError("AD adapter requires finite inputs and results");
		}
		template<typename T>
		void RequireFinite(const Dual<T>& x) {
			RequireFinite(x.value);
			RequireFinite(x.deriv);
		}
		template<int N>
		void RequirePoint(const VectorN<Real, N>& pos) {
			for (int i = 0; i < N; ++i)
				RequireFinite(pos[i]);
		}
		template<int N>
		void RequireIndex(int i) {
			if (i < 0 || i >= N)
				throw IndexError("AD adapter index outside dimension");
		}
		template<int N, typename Values>
		void RequireSize(const Values& values) {
			if (static_cast<size_t>(values.size()) != static_cast<size_t>(N))
				throw DomainError("AD expression output dimension mismatch");
		}
		template<int N, typename Values>
		void RequireMatrixSize(const Values& values) {
			RequireSize<N>(values);
			for (int i = 0; i < N; ++i)
				RequireSize<N>(values[i]);
		}
		inline void RequireBounds(Real lower, Real upper) {
			if (std::isnan(lower) || std::isnan(upper) || !(lower < upper))
				throw DomainError("AD adapter requires ordered parameter bounds");
		}
		template<int N, typename S>
		std::array<S, N> Point(const VectorN<Real, N>& pos) {
			RequirePoint(pos);
			std::array<S, N> result;
			for (int i = 0; i < N; ++i)
				result[i] = S(pos[i]);
			return result;
		}
		template<int N, typename F>
		Real Partial(const F& f, const VectorN<Real, N>& pos, int axis) {
			RequireIndex<N>(axis);
			auto x = Point<N, Dual<Real>>(pos);
			x[axis].deriv = 1;
			const Dual<Real> y = f(x);
			RequireFinite(y.value);
			RequireFinite(y.deriv);
			return y.deriv;
		}
		template<int N, typename F>
		Real SecondPartial(const F& f, const VectorN<Real, N>& pos, int a, int b) {
			RequireIndex<N>(a);
			RequireIndex<N>(b);
			using D = Dual<Real>;
			auto x = Point<N, Dual<D>>(pos);
			x[a].deriv.value = 1;
			x[b].value.deriv = 1;
			const Dual<D> y = f(x);
			RequireFinite(y.value.value);
			RequireFinite(y.value.deriv);
			RequireFinite(y.deriv.value);
			RequireFinite(y.deriv.deriv);
			return y.deriv.deriv;
		}
		template<int N, typename F>
		MetricJet<N> EvaluateMetricJet(const F& expression, const VectorN<Real, N>& pos) {
			MetricJet<N> jet;
			const auto values = expression(Point<N, Real>(pos));
			RequireMatrixSize<N>(values);
			for (int i = 0; i < N; ++i)
				for (int j = 0; j < N; ++j) {
					jet.value(i, j) = values[i][j];
					RequireFinite(jet.value(i, j));
					if (jet.value(i, j) != values[j][i])
						throw DomainError("AD metric requires symmetric covariant components");
				}
			using D = Dual<Real>;
			for (int a = 0; a < N; ++a) {
				auto x = Point<N, D>(pos);
				x[a].deriv = 1;
				const auto first = expression(x);
				RequireMatrixSize<N>(first);
				for (int i = 0; i < N; ++i)
					for (int j = 0; j < N; ++j) {
						RequireFinite(first[i][j]);
						jet.first[a](i, j) = first[i][j].deriv;
					}
				for (int b = a; b < N; ++b) {
					auto nested = Point<N, Dual<D>>(pos);
					nested[a].deriv.value = 1;
					nested[b].value.deriv = 1;
					const auto second = expression(nested);
					RequireMatrixSize<N>(second);
					for (int i = 0; i < N; ++i)
						for (int j = 0; j < N; ++j) {
							RequireFinite(second[i][j]);
							jet.second[a][b](i, j) = jet.second[b][a](i, j) = second[i][j].deriv.deriv;
						}
				}
			}
			return jet;
		}
	} // namespace AdapterDetail

	template<int N, typename Expression>
	class ADScalarField final : public IScalarFunction<N>, public IScalarDerivativeProvider<N> {
		Expression _expression;

	public:
		explicit ADScalarField(Expression expression)
				: _expression(std::move(expression)) {}
		Real operator()(const VectorN<Real, N>& pos) const override {
			const Real result = _expression(AdapterDetail::Point<N, Real>(pos));
			AdapterDetail::RequireFinite(result);
			return result;
		}
		Real EvaluatePartial(const VectorN<Real, N>& pos, int axis) const override { return AdapterDetail::Partial<N>(_expression, pos, axis); }
		Real EvaluateSecondPartial(const VectorN<Real, N>& pos, int a, int b) const override {
			return AdapterDetail::SecondPartial<N>(_expression, pos, a, b);
		}
		MatrixNM<Real, N, N> Hessian(const VectorN<Real, N>& pos) const {
			MatrixNM<Real, N, N> result;
			for (int a = 0; a < N; ++a)
				for (int b = a; b < N; ++b)
					result(a, b) = result(b, a) = EvaluateSecondPartial(pos, a, b);
			return result;
		}
	};
	template<int N, typename Expression>
	auto MakeADScalarField(Expression expression) {
		return ADScalarField<N, Expression>(std::move(expression));
	}

	template<int N, typename Expression>
	class ADVectorField final : public IVectorFunction<N>, public IVectorDerivativeProvider<N> {
		Expression _expression;

	public:
		explicit ADVectorField(Expression expression)
				: _expression(std::move(expression)) {}
		VectorN<Real, N> operator()(const VectorN<Real, N>& pos) const override {
			const auto values = _expression(AdapterDetail::Point<N, Real>(pos));
			AdapterDetail::RequireSize<N>(values);
			VectorN<Real, N> result;
			for (int i = 0; i < N; ++i) {
				result[i] = values[i];
				AdapterDetail::RequireFinite(result[i]);
			}
			return result;
		}
		Real EvaluateComponentPartial(const VectorN<Real, N>& pos, int component, int axis) const override {
			AdapterDetail::RequireIndex<N>(component);
			return AdapterDetail::Partial<N>(
				[&](const auto& x) {
					const auto values = _expression(x);
					AdapterDetail::RequireSize<N>(values);
					return values[component];
				},
				pos, axis);
		}
		MatrixNM<Real, N, N> Jacobian(const VectorN<Real, N>& pos) const {
			MatrixNM<Real, N, N> result;
			for (int j = 0; j < N; ++j) {
				auto x = AdapterDetail::Point<N, Dual<Real>>(pos);
				x[j].deriv = 1;
				const auto values = _expression(x);
				AdapterDetail::RequireSize<N>(values);
				for (int i = 0; i < N; ++i) {
					AdapterDetail::RequireFinite(values[i]);
					result(i, j) = values[i].deriv;
				}
			}
			return result;
		}
	};
	template<int N, typename Expression>
	auto MakeADVectorField(Expression expression) {
		return ADVectorField<N, Expression>(std::move(expression));
	}

	template<int N, typename Expression>
	class ADCurve final : public std::conditional_t<N == 3, Curves::ICurveCartesian3D, Curves::ICurveCartesian2D>,
												public ICurveDerivativeProvider<N> {
		static_assert(N == 2 || N == 3, "ADCurve supports Cartesian 2D and 3D curves");
		Expression _expression;
		Real _minT, _maxT;

	public:
		ADCurve(Expression expression, Real minT, Real maxT)
				: _expression(std::move(expression))
				, _minT(minT)
				, _maxT(maxT) {
			AdapterDetail::RequireBounds(minT, maxT);
		}
		Real getMinT() const override { return _minT; }
		Real getMaxT() const override { return _maxT; }
		VectorN<Real, N> operator()(Real t) const override {
			AdapterDetail::RequireFinite(t);
			const auto values = _expression(t);
			AdapterDetail::RequireSize<N>(values);
			VectorN<Real, N> result;
			for (int i = 0; i < N; ++i) {
				result[i] = values[i];
				AdapterDetail::RequireFinite(result[i]);
			}
			return result;
		}
		VectorN<Real, N> EvaluateDerivative(Real t, int order) const override {
			AdapterDetail::RequireFinite(t);
			if (order < 1 || order > 3)
				throw DomainError("AD curve requires order 1, 2 or 3");
			using D = Dual<Real>;
			using D2 = Dual<D>;
			VectorN<Real, N> result;
			if (order == 1) {
				const auto values = _expression(D(t, 1));
				AdapterDetail::RequireSize<N>(values);
				for (int i = 0; i < N; ++i) {
					AdapterDetail::RequireFinite(values[i]);
					result[i] = values[i].deriv;
				}
			} else if (order == 2) {
				const auto values = _expression(D2(D(t, 1), D(1, 0)));
				AdapterDetail::RequireSize<N>(values);
				for (int i = 0; i < N; ++i) {
					AdapterDetail::RequireFinite(values[i]);
					result[i] = values[i].deriv.deriv;
				}
			} else {
				const auto values = _expression(Dual<D2>(D2(D(t, 1), D(1, 0)), D2(D(1, 0), D(0, 0))));
				AdapterDetail::RequireSize<N>(values);
				for (int i = 0; i < N; ++i) {
					AdapterDetail::RequireFinite(values[i]);
					result[i] = values[i].deriv.deriv.deriv;
				}
			}
			for (int i = 0; i < N; ++i)
				AdapterDetail::RequireFinite(result[i]);
			return result;
		}
	};
	template<int N = 3, typename Expression>
	auto MakeADCurve(Expression expression, Real minT = Constants::NegInf, Real maxT = Constants::PosInf) {
		return ADCurve<N, Expression>(std::move(expression), minT, maxT);
	}

	template<int N, typename Expression>
	class ADMetric final : public MetricTensorField<N>, public IMetricJetProvider<N> {
		Expression _expression;
		typename MetricTensorField<N>::SignatureType _signature;

	public:
		explicit ADMetric(Expression expression, typename MetricTensorField<N>::SignatureType signature)
				: _expression(std::move(expression))
				, _signature(signature) {
			if (signature != MetricTensorField<N>::SignatureType::Riemannian && signature != MetricTensorField<N>::SignatureType::Lorentzian)
				throw DomainError("AD metric signature is invalid");
		}
		typename MetricTensorField<N>::SignatureType Signature() const override { return _signature; }
		Real component(int i, int j, const VectorN<Real, N>& pos) const override {
			AdapterDetail::RequireIndex<N>(i);
			AdapterDetail::RequireIndex<N>(j);
			const auto values = _expression(AdapterDetail::Point<N, Real>(pos));
			AdapterDetail::RequireMatrixSize<N>(values);
			AdapterDetail::RequireFinite(values[i][j]);
			if (values[i][j] != values[j][i])
				throw DomainError("AD metric requires symmetric covariant components");
			return values[i][j];
		}
		MetricJet<N> EvaluateMetricJet(const VectorN<Real, N>& pos) const override {
			return AdapterDetail::EvaluateMetricJet<N>(_expression, pos);
		}
	};
	template<int N, typename Expression>
	auto MakeADMetric(Expression expression,
										typename MetricTensorField<N>::SignatureType signature = MetricTensorField<N>::SignatureType::Riemannian) {
		return ADMetric<N, Expression>(std::move(expression), signature);
	}

	template<typename Expression>
	class ADSurface final : public Surfaces::ISurfaceCartesian, public IMetricJetProvider<2> {
		Expression _expression;
		Real _minU, _maxU, _minW, _maxW;
		VectorN<Real, 3> Partial(Real u, Real w, int a, int b = -1) const {
			VectorN<Real, 3> result;
			for (int i = 0; i < 3; ++i) {
				auto component = [&](const auto& x) {
					const auto values = _expression(x[0], x[1]);
					AdapterDetail::RequireSize<3>(values);
					return values[i];
				};
				result[i] = b < 0 ? AdapterDetail::Partial<2>(component, {u, w}, a) : AdapterDetail::SecondPartial<2>(component, {u, w}, a, b);
			}
			return result;
		}

	public:
		ADSurface(Expression expression, Real minU, Real maxU, Real minW, Real maxW)
				: _expression(std::move(expression))
				, _minU(minU)
				, _maxU(maxU)
				, _minW(minW)
				, _maxW(maxW) {
			AdapterDetail::RequireBounds(minU, maxU);
			AdapterDetail::RequireBounds(minW, maxW);
		}
		Real getMinU() const override { return _minU; }
		Real getMaxU() const override { return _maxU; }
		Real getMinW() const override { return _minW; }
		Real getMaxW() const override { return _maxW; }
		using Surfaces::ISurfaceCartesian::operator();
		using Surfaces::ISurfaceCartesian::GetFirstNormalFormCoefficients;
		using Surfaces::ISurfaceCartesian::GetSecondNormalFormCoefficients;
		VectorN<Real, 3> operator()(Real u, Real w) const override {
			AdapterDetail::RequireFinite(u);
			AdapterDetail::RequireFinite(w);
			const auto values = _expression(u, w);
			AdapterDetail::RequireSize<3>(values);
			VectorN<Real, 3> result;
			for (int i = 0; i < 3; ++i) {
				result[i] = values[i];
				AdapterDetail::RequireFinite(result[i]);
			}
			return result;
		}
		void Tangents(Real u, Real w, VectorN<Real, 3>& tU, VectorN<Real, 3>& tW) const override {
			tU = Partial(u, w, 0);
			tW = Partial(u, w, 1);
		}
		VectorN<Real, 3> Normal(Real u, Real w) const override {
			const Vec3Cart cross = VectorProduct(Vec3Cart(Partial(u, w, 0)), Vec3Cart(Partial(u, w, 1)));
			const Real norm = cross.NormL2();
			if (!std::isfinite(norm) || norm <= Precision::NumericalZeroThreshold)
				throw DomainError("AD surface normal is undefined at a degenerate parameterization");
			return cross / norm;
		}
		void GetFirstNormalFormCoefficients(Real u, Real w, Real& E, Real& F, Real& G) const override {
			const auto a = Partial(u, w, 0), b = Partial(u, w, 1);
			E = Utils::ScalarProduct(a, a);
			F = Utils::ScalarProduct(a, b);
			G = Utils::ScalarProduct(b, b);
			AdapterDetail::RequireFinite(E);
			AdapterDetail::RequireFinite(F);
			AdapterDetail::RequireFinite(G);
		}
		void GetSecondNormalFormCoefficients(Real u, Real w, Real& L, Real& M, Real& N) const override {
			const auto n = Normal(u, w);
			L = Utils::ScalarProduct(Partial(u, w, 0, 0), n);
			M = Utils::ScalarProduct(Partial(u, w, 0, 1), n);
			N = Utils::ScalarProduct(Partial(u, w, 1, 1), n);
			AdapterDetail::RequireFinite(L);
			AdapterDetail::RequireFinite(M);
			AdapterDetail::RequireFinite(N);
		}
		MetricJet<2> EvaluateMetricJet(const VectorN<Real, 2>& pos) const override {
			// Differentiate the induced metric expression, including its tangent seeds.
			auto metric = [&](const auto& x) {
				using S = std::decay_t<decltype(x[0])>;
				std::array<std::array<S, 3>, 2> tangents;
				for (int a = 0; a < 2; ++a) {
					const auto y = _expression(Dual<S>(x[0], S(a == 0)), Dual<S>(x[1], S(a == 1)));
					AdapterDetail::RequireSize<3>(y);
					for (int i = 0; i < 3; ++i)
						tangents[a][i] = y[i].deriv;
				}
				std::array<std::array<S, 2>, 2> g{};
				for (int a = 0; a < 2; ++a)
					for (int b = 0; b < 2; ++b)
						for (int i = 0; i < 3; ++i)
							g[a][b] += tangents[a][i] * tangents[b][i];
				return g;
			};
			return AdapterDetail::EvaluateMetricJet<2>(metric, pos);
		}
	};
	template<typename Expression>
	auto MakeADSurface(Expression expression, Real minU = Constants::NegInf, Real maxU = Constants::PosInf, Real minW = Constants::NegInf,
										 Real maxW = Constants::PosInf) {
		return ADSurface<Expression>(std::move(expression), minU, maxU, minW, maxW);
	}
} // namespace MML::AD
#endif
