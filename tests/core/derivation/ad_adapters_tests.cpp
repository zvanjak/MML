#include <catch2/catch_all.hpp>
#include "../../TestPrecision.h"
#include <limits>
#ifndef MML_SINGLE_HEADER_INCLUDE
#include <mml/core/Derivation/ADAdapters.h>
#include <mml/core/Derivation/ADJacobians.h>
#include <mml/core/DifferentialGeometry/InducedMetric.h>
#include <mml/core/Fields/FieldOperations.h>
#include <mml/core/Integration/PathIntegration.h>
#endif

using namespace MML;

namespace {
	Real NumericalMetricBudget(Real coordinate) {
		const Real h = Derivation::ScaleStep(Derivation::NDer4_h, coordinate);
		// Four 2D connection-derivative contributions, each with nested stencil L1 norm (3/2)^2.
		// The smooth metric fixtures vary on scale 2, giving a fourth-order (2h)^4 term.
		return 9 * Precision::Eps / (h * h) + POW4(2 * h);
	}

	template<typename F>
	void CheckThird(F f, Real x, Real value, Real first, Real second, Real third) {
		using D = AD::Dual<Real>;
		using D2 = AD::Dual<D>;
		const auto y = f(AD::Dual<D2>(D2(D(x, 1), D(1, 0)), D2(D(1, 0), D(0, 0))));
		CHECK(y.value.value.value == Catch::Approx(value).margin(TOL(1e-12, Precision::EqualityTolerance)));
		CHECK(y.deriv.value.value == Catch::Approx(first).margin(TOL(1e-12, Precision::EqualityTolerance)));
		CHECK(y.deriv.deriv.value == Catch::Approx(second).margin(TOL(1e-12, Precision::EqualityTolerance)));
		CHECK(y.deriv.deriv.deriv == Catch::Approx(third).margin(TOL(1e-12, Precision::EqualityTolerance)));
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
	const Real x = MML::Real(0.3);
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
		2, std::log(MML::Real(2.0)), MML::Real(0.5), -MML::Real(0.25), MML::Real(0.25));
	CheckThird(
		[](auto t) {
			using std::log10;
			return log10(t);
		},
		2, std::log10(MML::Real(2.0)), MML::Real(0.5) / std::log(MML::Real(10.0)), -MML::Real(0.25) / std::log(MML::Real(10.0)), MML::Real(0.25) / std::log(MML::Real(10.0)));
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
		0, 1, std::log(MML::Real(2.0)), std::pow(std::log(MML::Real(2.0)), MML::Real(2)), std::pow(std::log(MML::Real(2.0)), MML::Real(3)));
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
	for (Real scale : {MML::Real(0.01), MML::Real(1.0), MML::Real(100.0)}) {
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
		for (Real t : {MML::Real(0.0), MML::Real(0.4), MML::Real(3.0)}) {
			CHECK(base.getCurvature(t) == Catch::Approx(2 / (5 * scale)).epsilon(TOL(1e-12, Precision::EqualityTolerance)));
			CHECK(base.getTorsion(t) == Catch::Approx(1 / (5 * scale)).epsilon(TOL(1e-12, Precision::EqualityTolerance)));
			CHECK(base.getTangent(t).NormL2() == Catch::Approx(std::sqrt(5.0) * scale));
			CHECK(base.getNormal(t).NormL2() == Catch::Approx(1));
			CHECK(base.getBinormal(t).NormL2() == Catch::Approx(1));
		}
		CHECK(PathIntegration::ParametricCurveLength(curve, 0, 1) == Catch::Approx(std::sqrt(5.0) * scale));
		CHECK(realCalls == 0);
		Curves::CurveCartesian3D numerical(0, 2 * Constants::PI, [](Real t) { return VectorN<Real, 3>{2 * std::cos(t), 2 * std::sin(t), t}; });
		CHECK(numerical.getCurvature(MML::Real(0.4)) == Catch::Approx(MML::Real(0.4)).epsilon(MML::Real(1e-5)));
		const Real h = Derivation::ScaleStep(Derivation::NDer4_h, Real(0.4));
		// The third-derivative stencil has h^-3 roundoff and fourth-order truncation.
		const Real torsionBudget = std::max(Precision::CoarseTolerance, Precision::Eps / (h * h * h) + POW4(h));
		CHECK(numerical.getTorsion(MML::Real(0.4)) == Catch::Approx(MML::Real(0.2)).epsilon(torsionBudget));
	}
	auto circle = AD::MakeADCurve<2>([](auto t) {
		using std::sin;
		using std::cos;
		using S = decltype(t);
		return std::array<S, 2>{cos(t), sin(t)};
	});
	CHECK(Curves::getCurvature2D(circle, MML::Real(0.3)) == Catch::Approx(1));
	CHECK(Curves::getSignedCurvature2D(circle, MML::Real(0.3)) == Catch::Approx(1));
	CHECK_THROWS_AS(circle.EvaluateDerivative(0, 4), DomainError);
}

TEST_CASE("AD surfaces support extrinsic and intrinsic curvature through base references", "[ADAdapters][Surfaces][Metric]") {
	for (Real radius : {MML::Real(0.1), MML::Real(1.0), MML::Real(10.0)}) {
		auto sphere = AD::MakeADSurface(SphereExpression(radius), 0, Constants::PI, 0, 2 * Constants::PI);
		Surfaces::ISurfaceCartesian& base = sphere;
		const VectorN<Real, 2> x{MML::Real(0.8), MML::Real(0.4)};
		CHECK(base.GaussianCurvature(x[0], x[1]) == Catch::Approx(1 / (radius * radius)).epsilon(TOL(1e-12, Precision::EqualityTolerance)));
		CHECK(std::abs(base.MeanCurvature(x[0], x[1])) == Catch::Approx(1 / radius).epsilon(TOL(1e-12, Precision::EqualityTolerance)));
		CHECK(base.Normal(x[0], x[1]).NormL2() == Catch::Approx(1));
		Real k1, k2;
		base.PrincipalCurvatures(x[0], x[1], k1, k2);
		CHECK(k1 * k2 == Catch::Approx(1 / (radius * radius)).epsilon(TOL(1e-12, Precision::EqualityTolerance)));
		DifferentialGeometry::InducedMetric2D metric(base);
		CHECK(metric.GetRicciScalar(x) == Catch::Approx(2 / (radius * radius)).epsilon(TOL(1e-11, Precision::EqualityTolerance)));
		CHECK(DifferentialGeometry::GaussianCurvatureIntrinsic(base, x[0], x[1]) ==
					Catch::Approx(base.GaussianCurvature(x[0], x[1])).epsilon(TOL(1e-11, Precision::EqualityTolerance)));
		CHECK(metric.GetRiemannCurvatureTensor(0, 1, 0, 1, x) == Catch::Approx(std::pow(std::sin(x[0]), 2)).epsilon(TOL(1e-12, Precision::EqualityTolerance)));
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
	CHECK(plane.GaussianCurvature(MML::Real(0.3), MML::Real(0.4)) == 0);
	CHECK(DifferentialGeometry::GaussianCurvatureIntrinsic(plane, MML::Real(0.3), MML::Real(0.4)) == 0);
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
	const VectorN<Real, 2> x{MML::Real(0.8), MML::Real(0.4)};
	CHECK(base.GetChristoffelSymbolSecondKind(0, 1, 1, x) == Catch::Approx(-std::sin(x[0]) * std::cos(x[0])));
	CHECK(base.GetRicciScalar(x) == Catch::Approx(0.5).epsilon(TOL(1e-12, Precision::EqualityTolerance)));
	NumericalMetric numerical(2);
	CHECK(base.GetRicciScalar(x) == Catch::Approx(numerical.GetRicciScalar(x)).epsilon(TOL(1e-5, NumericalMetricBudget(x[0]))));
	Real derivativeError = 0;
	Derivation::NDer4Partial(base, 1, 1, 0, x, &derivativeError);
	REQUIRE(std::isfinite(derivativeError));
	REQUIRE(derivativeError >= 0);
	const Real exactConnection = base.GetChristoffelSymbolSecondKind(1, 0, 1, x);
	const Real connectionBudget = derivativeError / (2 * base.component(1, 1, x))
		+ Precision::EqualityTolerance * std::abs(exactConnection);
	CHECK(base.GetChristoffelSymbolSecondKindNumerical(1, 0, 1, x) ==
				Catch::Approx(exactConnection).epsilon(1e-8).margin(TOL(Real(0), connectionBudget)));
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
	VectorN<Real, 3> x{MML::Real(0.3), 2, 4};
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

TEST_CASE("AD field derivative indices are rejected before evaluating expressions", "[ADAdapters][Fields][regression]") {
	int scalarCalls = 0;
	auto scalar = AD::MakeADScalarField<2>([&](const auto& x) {
		++scalarCalls;
		return x[0] * x[1];
	});
	int vectorCalls = 0;
	auto vector = AD::MakeADVectorField<2>([&](const auto& x) {
		++vectorCalls;
		using S = std::decay_t<decltype(x[0])>;
		return std::array<S, 2>{x[0] * x[1], x[0] + x[1]};
	});
	const IScalarDerivativeProvider<2>& scalarProvider = scalar;
	const IVectorDerivativeProvider<2>& vectorProvider = vector;
	const VectorN<Real, 2> point{1, 2};
	for (int invalid : {std::numeric_limits<int>::min(), -1, 2, std::numeric_limits<int>::max()}) {
		CHECK_THROWS_AS(scalarProvider.EvaluatePartial(point, invalid), IndexError);
		CHECK_THROWS_AS(scalarProvider.EvaluateSecondPartial(point, invalid, 0), IndexError);
		CHECK_THROWS_AS(scalarProvider.EvaluateSecondPartial(point, 0, invalid), IndexError);
		CHECK_THROWS_AS(vectorProvider.EvaluateComponentPartial(point, invalid, 0), IndexError);
		CHECK_THROWS_AS(vectorProvider.EvaluateComponentPartial(point, 0, invalid), IndexError);
	}
	CHECK(scalarCalls == 0);
	CHECK(vectorCalls == 0);
	CHECK(scalarProvider.EvaluatePartial(point, 0) == 2);
	CHECK(vectorProvider.EvaluateComponentPartial(point, 0, 0) == 2);
	CHECK(scalarCalls == 1);
	CHECK(vectorCalls == 1);
}

TEST_CASE("AD metric inverse derivatives work with variable off-diagonal and Lorentzian metrics", "[ADAdapters][Metric]") {
	auto metric = AD::MakeADMetric<2>([](const auto& x) {
		using S = std::decay_t<decltype(x[0])>;
		using std::exp;
		S factor = exp(2 * (x[0] * x[0] + x[1] * x[1]));
		return std::array<std::array<S, 2>, 2>{{{factor, factor / 2}, {factor / 2, 2 * factor}}};
	});
	VectorN<Real, 2> x{MML::Real(0.3), -MML::Real(0.2)};
	CHECK(metric.GetRicciScalar(x) == Catch::Approx(-48.0 / 7 * std::exp(-2 * (x[0] * x[0] + x[1] * x[1]))).epsilon(TOL(1e-12, Precision::EqualityTolerance)));
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
	const Real u = MML::Real(0.3), v = -MML::Real(0.2);
	const Real expected = -4 / std::pow(1 + 4 * u * u + 4 * v * v, MML::Real(2));
	CHECK(saddle.GaussianCurvature(u, v) == Catch::Approx(expected).epsilon(TOL(1e-12, Precision::EqualityTolerance)));
	CHECK(DifferentialGeometry::GaussianCurvatureIntrinsic(saddle, u, v) == Catch::Approx(expected).epsilon(TOL(1e-12, Precision::EqualityTolerance)));
	CHECK_THROWS_AS(AD::MakeADSurface(SphereExpression(1), 1, 0), DomainError);
}

TEST_CASE("AD field coordinate conventions and detailed diagnostics stay consistent", "[ADAdapters][Fields]") {
	VectorN<Real, 3> x{2, MML::Real(0.8), MML::Real(0.4)};
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
	CHECK(detailed.function_evaluations == 0);
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
	const VectorN<Real, 2> x{2, MML::Real(0.4)};
	auto result = VectorFieldOperations::DivergenceViaDetDetailed(field, x, metric);
	CHECK(result.value == 2);
	CHECK(result.function_evaluations == 1);
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
		Real center = 0;
		mutable int offPointCalls = 0;
		Real component(int i, int j, const VectorN<Real, 2>& pos) const override {
			if (pos[0] != center) ++offPointCalls;
			return std::exp(2 * pos[0]) * InducedMetric2D::component(i, j, pos);
		}

	} modified(surface);
	VectorN<Real, 2> x{MML::Real(0.8), MML::Real(0.4)};
	modified.center = x[0];
	// The conformal factor changes Gamma^0_00 from the sphere's zero to one.
	CHECK(modified.GetChristoffelSymbolSecondKind(0, 0, 0, x) ==
		Catch::Approx(1).epsilon(TOL(1e-7, NumericalMetricBudget(x[0]))));
	CHECK(modified.offPointCalls > 0);
	const Real expected = 2 * std::exp(-2 * x[0]) * (1 - std::cos(x[0]) / std::sin(x[0]));
	CHECK(modified.GetRicciScalar(x) == Catch::Approx(expected).epsilon(1e-5).margin(TOL(1e-7, NumericalMetricBudget(x[0]))));
}

TEST_CASE("Metric connection derivative extension is virtual and propagates errors", "[ADAdapters][Metric]") {
	struct AnalyticalConnection : MetricTensorField<2> {
		mutable int connectionCalls = 0;
		Real component(int i, int j, const VectorN<Real, 2>&) const override { return i == j ? Real(1) : Real(0); }
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

TEST_CASE("AD curve derivative seeds preserve exactly representable polynomials", "[ADAdapters][Curve]") {
	int realCalls = 0;
	auto curve = AD::MakeADCurve<3>([&](auto t) {
		using S = decltype(t);
		if constexpr (std::is_same_v<S, Real>) ++realCalls;
		return std::array<S, 3>{t, t * t, t * t * t};
	});
	const Curves::ICurveCartesian3D& base = curve;
	for (Real t : {Real(-2), Real(-1), Real(0), Real(1), Real(2)}) {
		const std::array<VectorN<Real, 3>, 3> expected{
			VectorN<Real, 3>{1, 2 * t, 3 * t * t},
			VectorN<Real, 3>{0, 2, 6 * t},
			VectorN<Real, 3>{0, 0, 6}
		};
		for (int order = 1; order <= 3; ++order) {
			const auto derivative = Differential::CurveDerivative<3>(base, t, order);
			for (int component = 0; component < 3; ++component) {
				CAPTURE(t, order, component);
				CHECK(derivative[component] == expected[order - 1][component]);
			}
		}
	}
	CHECK(realCalls == 0);
}

TEST_CASE("AD scalar mixed derivatives preserve exact polynomial zeros and coefficients", "[ADAdapters][Field]") {
	int realCalls = 0;
	auto field = AD::MakeADScalarField<3>([&](const auto& x) {
		using S = std::decay_t<decltype(x[0])>;
		if constexpr (std::is_same_v<S, Real>) ++realCalls;
		return x[0] * x[0] * x[1] + x[0] * x[2];
	});
	const VectorN<Real, 3> point{2, 3, 4};
	CHECK(field(point) == 20);
	CHECK(realCalls == 1);
	const IScalarDerivativeProvider<3>& provider = field;
	const VectorN<Real, 3> expectedGradient{16, 4, 2};
	const Real expectedHessian[3][3]{{6, 4, 1}, {4, 0, 0}, {1, 0, 0}};
	const auto hessian = field.Hessian(point);
	for (int row = 0; row < 3; ++row) {
		CHECK(provider.EvaluatePartial(point, row) == expectedGradient[row]);
		for (int column = 0; column < 3; ++column) {
			CAPTURE(row, column);
			CHECK(hessian(row, column) == expectedHessian[row][column]);
			CHECK(provider.EvaluateSecondPartial(point, row, column) == expectedHessian[row][column]);
		}
	}
	CHECK(realCalls == 1);
}
