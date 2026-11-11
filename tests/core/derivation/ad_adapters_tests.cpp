#include <catch2/catch_all.hpp>
#include <mml/core/Derivation/ADAdapters.h>
#include <mml/core/Derivation/ADJacobians.h>
#include <mml/core/DifferentialGeometry/InducedMetric.h>
#include <mml/core/Fields/FieldOperations.h>
#include <mml/core/Integration/PathIntegration.h>

using namespace MML;

namespace {
	template<typename F>
	void CheckThird(F f, Real x, Real value, Real first, Real second, Real third) {
		using D = AD::Dual<Real>;
		using D2 = AD::Dual<D>;
		const auto y = f(AD::Dual<D2>(D2(D(x, 1), D(1, 0)), D2(D(1, 0), D(0, 0))));
		CHECK(y.value.value.value == Catch::Approx(value).margin(1e-12));
		CHECK(y.deriv.value.value == Catch::Approx(first).margin(1e-12));
		CHECK(y.deriv.deriv.value == Catch::Approx(second).margin(1e-12));
		CHECK(y.deriv.deriv.deriv == Catch::Approx(third).margin(1e-12));
	}

	auto SphereExpression(Real radius) {
		return [radius](auto u, auto w) {
			using std::sin;
			using std::cos;
			using S = decltype(u);
			return std::array<S, 3>{radius * sin(u) * cos(w), radius * sin(u) * sin(w), radius * cos(u)};
		};
	}

	struct NumericalMetric : MetricTensorField<2> {
		Real radius;
		explicit NumericalMetric(Real r)
				: radius(r) {}
		Real component(int i, int j, const VectorN<Real, 2>& x) const override {
			if (i != j)
				return 0;
			return i == 0 ? radius * radius : radius * radius * std::sin(x[0]) * std::sin(x[0]);
		}
	};

	struct FailingScalar : IScalarFunction<3>, IScalarDerivativeProvider<3> {
		mutable int calls = 0;
		Real operator()(const VectorN<Real, 3>& x) const override {
			++calls;
			return x[0];
		}
		Real EvaluatePartial(const VectorN<Real, 3>&, int) const override { throw DomainError("provider failure"); }
		Real EvaluateSecondPartial(const VectorN<Real, 3>&, int, int) const override { throw DomainError("provider failure"); }
	};
} // namespace

TEST_CASE("Nested ForwardAD elementary functions support third derivatives", "[ADAdapters][ForwardAD]") {
	const Real x = 0.3;
	CheckThird(
		[](auto t) {
			using std::sin;
			return sin(t);
		},
		x, std::sin(x), std::cos(x), -std::sin(x), -std::cos(x));
	CheckThird(
		[](auto t) {
			using std::cos;
			return cos(t);
		},
		x, std::cos(x), -std::sin(x), -std::cos(x), std::sin(x));
	CheckThird(
		[](auto t) {
			using std::exp;
			return exp(t);
		},
		x, std::exp(x), std::exp(x), std::exp(x), std::exp(x));
	CheckThird(
		[](auto t) {
			using std::log;
			return log(t);
		},
		2, std::log(2.0), 0.5, -0.25, 0.25);
	CheckThird(
		[](auto t) {
			using std::log10;
			return log10(t);
		},
		2, std::log10(2.0), 0.5 / std::log(10.0), -0.25 / std::log(10.0), 0.25 / std::log(10.0));
	CheckThird(
		[](auto t) {
			using std::sqrt;
			return sqrt(t);
		},
		4, 2, 0.25, -1.0 / 32, 3.0 / 256);
	CheckThird(
		[](auto t) {
			using std::sin;
			using std::asin;
			return sin(asin(t));
		},
		x, x, 1, 0, 0);
	CheckThird(
		[](auto t) {
			using std::cos;
			using std::acos;
			return cos(acos(t));
		},
		x, x, 1, 0, 0);
	CheckThird(
		[](auto t) {
			using std::tan;
			using std::atan;
			return tan(atan(t));
		},
		x, x, 1, 0, 0);
	CheckThird(
		[](auto t) {
			using std::sinh;
			return sinh(t);
		},
		x, std::sinh(x), std::cosh(x), std::sinh(x), std::cosh(x));
	CheckThird(
		[](auto t) {
			using std::cosh;
			return cosh(t);
		},
		x, std::cosh(x), std::sinh(x), std::cosh(x), std::sinh(x));
	CheckThird(
		[](auto t) {
			using std::tanh;
			using std::log;
			return tanh(log(t));
		},
		1, 0, 1, -1, 0);
	CheckThird(
		[](auto t) {
			using std::abs;
			return abs(t);
		},
		-2, 2, -1, 0, 0);
	CheckThird(
		[](auto t) {
			using std::floor;
			return floor(t);
		},
		x, 0, 0, 0, 0);
	CheckThird(
		[](auto t) {
			using std::ceil;
			return ceil(t);
		},
		x, 1, 0, 0, 0);
	CheckThird(
		[](auto t) {
			using std::atan2;
			using S = decltype(t);
			return atan2(t, S(1));
		},
		1, Constants::PI / 4, 0.5, -0.5, 0.5);
	CheckThird(
		[](auto t) {
			using std::pow;
			return pow(t, 3);
		},
		0, 0, 0, 0, 6);
	CheckThird(
		[](auto t) {
			using std::pow;
			return pow(t, 0);
		},
		0, 1, 0, 0, 0);
	CheckThird(
		[](auto t) {
			using std::pow;
			return pow(2, t);
		},
		0, 1, std::log(2.0), std::pow(std::log(2.0), 2), std::pow(std::log(2.0), 3));
	CheckThird(
		[](auto t) {
			using std::fma;
			return fma(t, t, t);
		},
		x, x * x + x, 2 * x + 1, 2, 0);
}

TEST_CASE("Nested variable powers preserve mixed exponent seeds", "[ADAdapters][ForwardAD]") {
	auto f = [](const auto& x) {
		using std::pow;
		return pow(x[0], x[1]);
	};
	auto h = AD::CalcHessianForwardAD(f, std::vector<Real>{2, 3});
	CHECK(h(0, 0) == Catch::Approx(12));
	CHECK(h(0, 1) == Catch::Approx(4 * (1 + 3 * std::log(2.0))));
	CHECK(h(1, 1) == Catch::Approx(8 * std::pow(std::log(2.0), 2)));
}

TEST_CASE("AD curves drive base geometry and arc length without stencils", "[ADAdapters][Curves]") {
	for (Real scale : {0.01, 1.0, 100.0}) {
		int realCalls = 0;
		auto expression = [scale, &realCalls](auto t) {
			using S = decltype(t);
			using std::sin;
			using std::cos;
			if constexpr (std::is_same_v<S, Real>)
				++realCalls;
			return std::array<S, 3>{2 * scale * cos(t), 2 * scale * sin(t), scale * t};
		};
		auto curve = AD::MakeADCurve(expression, 0, 2 * Constants::PI);
		const Curves::ICurveCartesian3D& base = curve;
		for (Real t : {0.0, 0.4, 3.0}) {
			CHECK(base.getCurvature(t) == Catch::Approx(2 / (5 * scale)).epsilon(1e-12));
			CHECK(base.getTorsion(t) == Catch::Approx(1 / (5 * scale)).epsilon(1e-12));
			CHECK(base.getTangent(t).NormL2() == Catch::Approx(std::sqrt(5.0) * scale));
			CHECK(base.getNormal(t).NormL2() == Catch::Approx(1));
			CHECK(base.getBinormal(t).NormL2() == Catch::Approx(1));
		}
		CHECK(PathIntegration::ParametricCurveLength(curve, 0, 1) == Catch::Approx(std::sqrt(5.0) * scale));
		CHECK(realCalls == 0);
		Curves::CurveCartesian3D numerical(0, 2 * Constants::PI, [](Real t) { return VectorN<Real, 3>{2 * std::cos(t), 2 * std::sin(t), t}; });
		CHECK(numerical.getCurvature(0.4) == Catch::Approx(0.4).epsilon(1e-5));
		CHECK(numerical.getTorsion(0.4) == Catch::Approx(0.2).epsilon(1e-4));
	}
	auto circle = AD::MakeADCurve<2>([](auto t) {
		using std::sin;
		using std::cos;
		using S = decltype(t);
		return std::array<S, 2>{cos(t), sin(t)};
	});
	CHECK(Curves::getCurvature2D(circle, 0.3) == Catch::Approx(1));
	CHECK(Curves::getSignedCurvature2D(circle, 0.3) == Catch::Approx(1));
	CHECK_THROWS_AS(circle.EvaluateDerivative(0, 4), DomainError);
}

TEST_CASE("AD surfaces support extrinsic and intrinsic curvature through base references", "[ADAdapters][Surfaces][Metric]") {
	for (Real radius : {0.1, 1.0, 10.0}) {
		auto sphere = AD::MakeADSurface(SphereExpression(radius), 0, Constants::PI, 0, 2 * Constants::PI);
		Surfaces::ISurfaceCartesian& base = sphere;
		const VectorN<Real, 2> x{0.8, 0.4};
		CHECK(base.GaussianCurvature(x[0], x[1]) == Catch::Approx(1 / (radius * radius)).epsilon(1e-12));
		CHECK(std::abs(base.MeanCurvature(x[0], x[1])) == Catch::Approx(1 / radius).epsilon(1e-12));
		CHECK(base.Normal(x[0], x[1]).NormL2() == Catch::Approx(1));
		Real k1, k2;
		base.PrincipalCurvatures(x[0], x[1], k1, k2);
		CHECK(k1 * k2 == Catch::Approx(1 / (radius * radius)).epsilon(1e-12));
		DifferentialGeometry::InducedMetric2D metric(base);
		CHECK(metric.GetRicciScalar(x) == Catch::Approx(2 / (radius * radius)).epsilon(1e-11));
		CHECK(DifferentialGeometry::GaussianCurvatureIntrinsic(base, x[0], x[1]) ==
					Catch::Approx(base.GaussianCurvature(x[0], x[1])).epsilon(1e-11));
		CHECK(metric.GetRiemannCurvatureTensor(0, 1, 0, 1, x) == Catch::Approx(std::pow(std::sin(x[0]), 2)).epsilon(1e-12));
		Surfaces::Sphere numerical(radius);
		const Real numericalK = DifferentialGeometry::GaussianCurvatureIntrinsic(numerical, x[0], x[1]);
		const Real exactK = 1 / (radius * radius);
		CHECK(std::isfinite(numericalK));
		CHECK(std::abs(DifferentialGeometry::GaussianCurvatureIntrinsic(base, x[0], x[1]) - exactK) <=
					std::abs(numericalK - exactK) / 1000 + 64 * Precision::Eps * std::abs(exactK));
	}
	auto plane = AD::MakeADSurface([](auto u, auto v) {
		using S = decltype(u);
		return std::array<S, 3>{u, v, S(0)};
	});
	CHECK(plane.GaussianCurvature(0.3, 0.4) == 0);
	CHECK(DifferentialGeometry::GaussianCurvatureIntrinsic(plane, 0.3, 0.4) == 0);
	auto degenerate = AD::MakeADSurface([](auto u, auto v) {
		using S = decltype(u);
		return std::array<S, 3>{u + v, S(0), S(0)};
	});
	CHECK_THROWS_AS(degenerate.GaussianCurvature(0, 0), DomainError);
	CHECK_THROWS_AS(DifferentialGeometry::GaussianCurvatureIntrinsic(degenerate, 0, 0), SingularMatrixError);
}

TEST_CASE("AD metrics compute connection and Riemann without numerical differentiation", "[ADAdapters][Metric]") {
	int realCalls = 0;
	auto metric = AD::MakeADMetric<2>([&](const auto& x) {
		using S = std::decay_t<decltype(x[0])>;
		using std::sin;
		if constexpr (std::is_same_v<S, Real>)
			++realCalls;
		return std::array<std::array<S, 2>, 2>{{{S(4), S(0)}, {S(0), 4 * sin(x[0]) * sin(x[0])}}};
	});
	const MetricTensorField<2>& base = metric;
	const VectorN<Real, 2> x{0.8, 0.4};
	CHECK(base.GetChristoffelSymbolSecondKind(0, 1, 1, x) == Catch::Approx(-std::sin(x[0]) * std::cos(x[0])));
	CHECK(base.GetRicciScalar(x) == Catch::Approx(0.5).epsilon(1e-12));
	NumericalMetric numerical(2);
	CHECK(base.GetRicciScalar(x) == Catch::Approx(numerical.GetRicciScalar(x)).epsilon(1e-5));
	CHECK(base.GetChristoffelSymbolSecondKindNumerical(1, 0, 1, x) ==
				Catch::Approx(base.GetChristoffelSymbolSecondKind(1, 0, 1, x)).epsilon(1e-8));
	CHECK(realCalls > 0);
	CHECK_THROWS_AS(base.GetRiemannCurvatureTensor(0, 0, 2, 0, x), IndexError);
	CHECK_THROWS_AS(base.GetChristoffelSymbolSecondKind(0, 0, 0, {0, 0}), SingularMatrixError);

	auto flat = AD::MakeADMetric<2>([](const auto& x) {
		using S = std::decay_t<decltype(x[0])>;
		return std::array<std::array<S, 2>, 2>{{{S(1), S(0.5)}, {S(0.5), S(2)}}};
	});
	for (int a = 0; a < 2; ++a)
		for (int b = 0; b < 2; ++b)
			for (int c = 0; c < 2; ++c)
				for (int d = 0; d < 2; ++d)
					CHECK(flat.GetRiemannCurvatureTensor(a, b, c, d, x) == 0);
}

TEST_CASE("AD fields feed existing scalar and vector calculus operations", "[ADAdapters][Fields]") {
	auto potential = AD::MakeADScalarField<3>([](const auto& x) {
		using std::sin;
		return sin(x[0]) + x[1] * x[1] + x[0] * x[2];
	});
	const IScalarFunction<3>& scalar = potential;
	VectorN<Real, 3> x{0.3, 2, 4};
	auto gradient = ScalarFieldOperations::GradientCart(scalar, x);
	CHECK(gradient[0] == Catch::Approx(std::cos(x[0]) + x[2]));
	CHECK(gradient[1] == 4);
	CHECK(gradient[2] == x[0]);
	CHECK(ScalarFieldOperations::LaplacianCart(scalar, x) == Catch::Approx(2 - std::sin(x[0])));
	CHECK(potential.Hessian(x)(0, 2) == 1);
	CHECK(potential.Hessian(x)(2, 0) == 1);
	auto field = AD::MakeADVectorField<3>([](const auto& p) {
		using S = std::decay_t<decltype(p[0])>;
		return std::array<S, 3>{p[1], -p[0], p[2] * p[2]};
	});
	const IVectorFunction<3>& vector = field;
	CHECK(VectorFieldOperations::DivCart(vector, x) == 8);
	CHECK(VectorFieldOperations::CurlCart(vector, x)[2] == -2);
	CHECK(field.Jacobian(x)(2, 2) == 8);
	MetricTensorCartesian3D cart;
	CHECK(VectorFieldOperations::Divergence(vector, x, cart) == 8);
	CHECK(cart.CovariantDerivativeContravar(vector, 2, x)[2] == 8);
	CHECK(ScalarFieldOperations::Gradient(scalar, x, cart)[0] == gradient[0]);
}

TEST_CASE("AD failures propagate without silently invoking real-only stencils", "[ADAdapters][Fields]") {
	FailingScalar broken;
	CHECK_THROWS_AS(ScalarFieldOperations::GradientCart(broken, {1, 2, 3}), DomainError);
	CHECK_THROWS_AS(ScalarFieldOperations::LaplacianCart(broken, {1, 2, 3}), DomainError);
	CHECK(broken.calls == 0);
	auto bad = AD::MakeADScalarField<2>([](const auto& x) {
		using std::sqrt;
		return sqrt(x[0]) + x[1];
	});
	CHECK_THROWS_AS(bad.EvaluatePartial({0, 1}, 0), DomainError);
	CHECK_THROWS_AS(bad.EvaluatePartial({1, 1}, 2), IndexError);
	CHECK_THROWS_AS(bad({Constants::PosInf, 0}), DomainError);
	auto nanCurve = AD::MakeADCurve([](auto t) {
		using S = decltype(t);
		using std::sqrt;
		return std::array<S, 3>{sqrt(t), t, t};
	});
	CHECK_THROWS_AS(nanCurve.getCurvature(0), DomainError);
	FieldOperationConfig config;
	config.exception_policy = EvaluationExceptionPolicy::ConvertToStatus;
	CHECK(ScalarFieldOperations::GradientCartDetailed(broken, {1, 2, 3}, config).status == AlgorithmStatus::InvalidInput);
	CHECK(ScalarFieldOperations::LaplacianCartDetailed(broken, {1, 2, 3}, config).status == AlgorithmStatus::InvalidInput);
	CHECK(broken.calls == 0);
}

TEST_CASE("AD metric inverse derivatives work with variable off-diagonal and Lorentzian metrics", "[ADAdapters][Metric]") {
	auto metric = AD::MakeADMetric<2>([](const auto& x) {
		using S = std::decay_t<decltype(x[0])>;
		using std::exp;
		S factor = exp(2 * (x[0] * x[0] + x[1] * x[1]));
		return std::array<std::array<S, 2>, 2>{{{factor, factor / 2}, {factor / 2, 2 * factor}}};
	});
	VectorN<Real, 2> x{0.3, -0.2};
	CHECK(metric.GetRicciScalar(x) == Catch::Approx(-48.0 / 7 * std::exp(-2 * (x[0] * x[0] + x[1] * x[1]))).epsilon(1e-12));
	auto lorentz = AD::MakeADMetric<2>(
		[](const auto& p) {
			using S = std::decay_t<decltype(p[0])>;
			return std::array<std::array<S, 2>, 2>{{{S(-1), S(0)}, {S(0), S(1)}}};
		},
		MetricTensorField<2>::SignatureType::Lorentzian);
	CHECK(lorentz.IsLorentzian());
	CHECK(lorentz.GetRicciScalar(x) == 0);
	auto asymmetric = AD::MakeADMetric<2>([](const auto& p) {
		using S = std::decay_t<decltype(p[0])>;
		return std::array<std::array<S, 2>, 2>{{{S(1), S(0.1)}, {S(0), S(1)}}};
	});
	CHECK_THROWS_AS(asymmetric.GetRicciScalar(x), DomainError);
}

TEST_CASE("AD surface mixed derivatives recover negative Gaussian curvature", "[ADAdapters][Surfaces]") {
	auto saddle = AD::MakeADSurface([](auto u, auto v) {
		using S = decltype(u);
		return std::array<S, 3>{u, v, u * u - v * v};
	});
	const Real u = 0.3, v = -0.2;
	const Real expected = -4 / std::pow(1 + 4 * u * u + 4 * v * v, 2);
	CHECK(saddle.GaussianCurvature(u, v) == Catch::Approx(expected).epsilon(1e-12));
	CHECK(DifferentialGeometry::GaussianCurvatureIntrinsic(saddle, u, v) == Catch::Approx(expected).epsilon(1e-12));
	CHECK_THROWS_AS(AD::MakeADSurface(SphereExpression(1), 1, 0), DomainError);
}

TEST_CASE("AD field coordinate conventions and detailed diagnostics stay consistent", "[ADAdapters][Fields]") {
	VectorN<Real, 3> x{2, 0.8, 0.4};
	auto radial = AD::MakeADScalarField<3>([](const auto& p) { return p[0] * p[0]; });
	MetricTensorSpherical sphere;
	auto physical = ScalarFieldOperations::GradientSpher(radial, x);
	auto coordinate = ScalarFieldOperations::Gradient(radial, x, sphere);
	CHECK(physical[0] == 4);
	CHECK(coordinate[0] == physical[0]);
	CHECK(ScalarFieldOperations::LaplacianSpher(radial, x) == 6);
	auto cylinderPotential = AD::MakeADScalarField<3>([](const auto& p) { return p[0] * p[0] + p[2] * p[2]; });
	CHECK(ScalarFieldOperations::LaplacianCyl(cylinderPotential, x) == 6);
	auto physicalRadial = AD::MakeADVectorField<3>([](const auto& p) {
		using S = std::decay_t<decltype(p[0])>;
		return std::array<S, 3>{p[0], S(0), S(0)};
	});
	CHECK(VectorFieldOperations::DivSpher(physicalRadial, x) == 3);
	CHECK(VectorFieldOperations::DivCyl(physicalRadial, x) == 2);
	auto rotation = AD::MakeADVectorField<3>([](const auto& p) {
		using S = std::decay_t<decltype(p[0])>;
		return std::array<S, 3>{S(0), p[0], S(0)};
	});
	CHECK(VectorFieldOperations::CurlCyl(rotation, x)[2] == 2);
	auto sphericalRotation = AD::MakeADVectorField<3>([](const auto& p) {
		using S = std::decay_t<decltype(p[0])>;
		using std::sin;
		return std::array<S, 3>{S(0), S(0), p[0] * sin(p[1])};
	});
	const auto curl = VectorFieldOperations::CurlSpher(sphericalRotation, x);
	CHECK(curl[0] == Catch::Approx(2 * std::cos(x[1])));
	CHECK(curl[1] == Catch::Approx(-2 * std::sin(x[1])));
	auto detailed = ScalarFieldOperations::GradientCartDetailed(radial, x);
	CHECK(detailed.value[0] == 4);
	CHECK(detailed.function_evaluations == 3);
	CHECK(std::isnan(detailed.error[0]));
	CHECK(std::isnan(ScalarFieldOperations::LaplacianCartDetailed(radial, x).error));
	CHECK(std::isnan(VectorFieldOperations::DivCartDetailed(physicalRadial, x).error));
	CHECK(std::isnan(VectorFieldOperations::CurlCartDetailed(rotation, x).error[0]));
	CHECK(std::isnan(ScalarFieldOperations::GradientSpherDetailed(radial, x).error[0]));
	CHECK(std::isnan(ScalarFieldOperations::GradientCylDetailed(radial, x).error[0]));
	CHECK(std::isnan(ScalarFieldOperations::LaplacianSpherDetailed(radial, x).error));
	CHECK(std::isnan(ScalarFieldOperations::LaplacianCylDetailed(radial, x).error));
	CHECK(std::isnan(VectorFieldOperations::DivSpherDetailed(physicalRadial, x).error));
	CHECK(std::isnan(VectorFieldOperations::DivCylDetailed(physicalRadial, x).error));
	CHECK(std::isnan(VectorFieldOperations::CurlSpherDetailed(sphericalRotation, x).error[0]));
	CHECK(std::isnan(VectorFieldOperations::CurlCylDetailed(rotation, x).error[0]));
	FieldOperationConfig numerical;
	numerical.derivative_order = 4;
	auto explicitNumerical = ScalarFieldOperations::GradientCartDetailed(radial, x, numerical);
	CHECK(explicitNumerical.value[0] == Catch::Approx(4));
	CHECK(std::isfinite(explicitNumerical.error[0]));
	CHECK_THROWS_AS(ScalarFieldOperations::GradientSpher(radial, {0, 0, 0}), DomainError);
}

TEST_CASE("AD determinant divergence uses metric jets and preserves numerical reference", "[ADAdapters][Fields][Metric]") {
	auto metric = AD::MakeADMetric<2>([](const auto& x) {
		using S = std::decay_t<decltype(x[0])>;
		return std::array<std::array<S, 2>, 2>{{{S(1), S(0)}, {S(0), x[0] * x[0]}}};
	});
	auto field = AD::MakeADVectorField<2>([](const auto& x) {
		using S = std::decay_t<decltype(x[0])>;
		return std::array<S, 2>{x[0], S(0)};
	});
	const VectorN<Real, 2> x{2, 0.4};
	auto result = VectorFieldOperations::DivergenceViaDetDetailed(field, x, metric);
	CHECK(result.value == 2);
	CHECK(result.function_evaluations == 3);
	CHECK(std::isnan(result.error));
	CHECK(VectorFieldOperations::Divergence(field, x, metric) == 2);
	FieldOperationConfig config;
	config.derivative_order = 4;
	auto numerical = VectorFieldOperations::DivergenceViaDetDetailed(field, x, metric, config);
	CHECK(numerical.value == Catch::Approx(2));
	CHECK(std::isfinite(numerical.error));
	config.derivative_order = 0;
	config.exception_policy = EvaluationExceptionPolicy::ConvertToStatus;
	CHECK(VectorFieldOperations::DivergenceViaDetDetailed(field, {0, 0}, metric, config).status == AlgorithmStatus::InvalidInput);
}

TEST_CASE("Geometry AD failures never switch to off-point numerical evaluations", "[ADAdapters][Curves][Surfaces][Metric]") {
	int realCalls = 0;
	auto curve = AD::MakeADCurve([&](auto t) {
		using S = decltype(t);
		if constexpr (std::is_same_v<S, Real>)
			++realCalls;
		else
			throw DomainError("AD expression failed");
		return std::array<S, 3>{t, t, t};
	});
	CHECK_THROWS_AS(curve.getCurvature(1), DomainError);
	CHECK(realCalls == 0);
	auto surface = AD::MakeADSurface([&](auto u, auto v) {
		using S = decltype(u);
		if constexpr (std::is_same_v<S, Real>)
			++realCalls;
		else
			throw DomainError("AD expression failed");
		return std::array<S, 3>{u, v, S(0)};
	});
	CHECK_THROWS_AS(surface.GaussianCurvature(1, 1), DomainError);
	CHECK(realCalls == 0);
	CHECK_THROWS_AS(DifferentialGeometry::GaussianCurvatureIntrinsic(surface, 1, 1), DomainError);
	CHECK(realCalls == 0); // Even the induced metric's values require the surface's AD tangents.
	auto metric = AD::MakeADMetric<2>([&](const auto& x) {
		using S = std::decay_t<decltype(x[0])>;
		if constexpr (std::is_same_v<S, Real>)
			++realCalls;
		else
			throw DomainError("AD expression failed");
		return std::array<std::array<S, 2>, 2>{{{S(1), S(0)}, {S(0), S(1)}}};
	});
	CHECK_THROWS_AS(metric.GetRicciScalar({1, 1}), DomainError);
	CHECK(realCalls == 5); // Four component reads for inversion, then one center jet evaluation; no stencil.
}

TEST_CASE("AD adapters reject wrong output shapes and invalid bounds", "[ADAdapters]") {
	auto vector = AD::MakeADVectorField<3>([](const auto& x) {
		using S = std::decay_t<decltype(x[0])>;
		return std::array<S, 2>{x[0], x[1]};
	});
	CHECK_THROWS_AS(vector({1, 2, 3}), DomainError);
	CHECK_THROWS_AS(vector.Jacobian({1, 2, 3}), DomainError);
	CHECK_THROWS_AS(vector.EvaluateComponentPartial({1, 2, 3}, 0, 0), DomainError);
	auto curveExpression = [](auto t) {
		using S = decltype(t);
		return std::array<S, 2>{t, t};
	};
	auto curve = AD::MakeADCurve<3>(curveExpression);
	CHECK_THROWS_AS(curve(1), DomainError);
	CHECK_THROWS_AS(curve.getTorsion(1), DomainError);
	CHECK_THROWS_AS(AD::MakeADCurve<2>(curveExpression, 2, 1), DomainError);
	auto metric = AD::MakeADMetric<2>([](const auto& x) {
		using S = std::decay_t<decltype(x[0])>;
		return std::array<std::array<S, 1>, 1>{{{S(1)}}};
	});
	CHECK_THROWS_AS(metric.component(0, 0, {1, 2}), DomainError);
	CHECK_THROWS_AS(metric.GetRicciScalar({1, 2}), DomainError);
}

TEST_CASE("Modified induced-metric subclasses keep numerical connection fallback", "[ADAdapters][Metric]") {
	auto surface = AD::MakeADSurface(SphereExpression(1));
	struct ConformalInduced : DifferentialGeometry::InducedMetric2D {
		using InducedMetric2D::InducedMetric2D;
		Real component(int i, int j, const VectorN<Real, 2>& pos) const override {
			return std::exp(2 * pos[0]) * InducedMetric2D::component(i, j, pos);
		}

	} modified(surface);
	VectorN<Real, 2> x{0.8, 0.4};
	const Real expected = 2 * std::exp(-2 * x[0]) * (1 - std::cos(x[0]) / std::sin(x[0]));
	CHECK(modified.GetRicciScalar(x) == Catch::Approx(expected).epsilon(1e-5).margin(1e-7));
}

TEST_CASE("Metric connection derivative extension is virtual and propagates errors", "[ADAdapters][Metric]") {
	struct AnalyticalConnection : MetricTensorField<2> {
		mutable int connectionCalls = 0;
		Real component(int i, int j, const VectorN<Real, 2>&) const override { return i == j ? 1 : 0; }
		Real GetChristoffelSymbolSecondKind(int, int, int, const VectorN<Real, 2>&) const override {
			++connectionCalls;
			return 0;
		}

	protected:
		Real DeriveChristoffelSymbolSecondKind(int, int, int, int, const VectorN<Real, 2>&) const override {
			throw DomainError("Analytical connection derivative failure");
		}
	} metric;
	const MetricTensorField<2>& base = metric;
	CHECK_THROWS_AS(base.GetRiemannCurvatureTensor(0, 1, 0, 1, {1, 1}), DomainError);
	CHECK(metric.connectionCalls == 0);
}
