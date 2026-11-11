#include <catch2/catch_test_macros.hpp>
#include <catch2/matchers/catch_matchers_floating_point.hpp>
#include "../../TestPrecision.h"
#ifndef MML_SINGLE_HEADER_INCLUDE
#include <mml/base/Function.h>
#include <mml/base/Geometry/Geometry3DBodies.h>
#include <mml/core/Integration.h>
#include <mml/core/Integration/SurfaceIntegration.h>
#include <mml/core/Integration/PathIntegration.h>
#endif
#include <algorithm>
#include <atomic>
#include <cmath>
#include <limits>
#include <type_traits>

using namespace MML;
using Catch::Matchers::WithinAbs;
using Catch::Matchers::WithinRel;

namespace {
struct ContextRestore {
	AlgorithmContext saved = AlgorithmContext::Get();
	~ContextRestore() { AlgorithmContext::Get() = saved; }
};
Real Lower(Real) { return 0; }
Real Upper(Real) { return 1; }
Real Lower3(Real, Real) { return 0; }
Real Upper3(Real, Real) { return 1; }

struct HugeInterval {
	Real lower, upper, scale;
};

HugeInterval HugeFiniteInterval() {
	const Real unit = std::ldexp(Real(1), std::numeric_limits<Real>::max_exponent - 7);
	return {100 * unit, 110 * unit, 10 * unit};
}

[[maybe_unused]] inline Real Chebyshev14(Real x) {
	Real previous = 1, current = x;
	for (int n = 2; n <= 14; ++n) {
		const Real next = 2 * x * current - previous;
		previous = current;
		current = next;
	}
	return current;
}
}

TEST_CASE("Closed-loop work uses a dimensionally consistent cancellation scale", "[integration-contract]") {
	const ParametricCurve<3> circle([](Real t) { return VectorN<Real, 3>{std::cos(t), std::sin(t), 0}; });
	struct AnalyticCircle : ParametricCurve<3>, ICurveDerivativeProvider<3> {
		AnalyticCircle() : ParametricCurve<3>([](Real t) {
			return VectorN<Real, 3>{std::cos(t), std::sin(t), 0};
		}) {}
		VectorN<Real, 3> EvaluateDerivative(Real t, int order) const override {
			if (order == 1) return {-std::sin(t), std::cos(t), 0};
			if (order == 2) return {-std::cos(t), -std::sin(t), 0};
			if (order == 3) return {std::sin(t), -std::cos(t), 0};
			throw ArgumentError("Circle supports derivative orders 1 to 3");
		}
	} analyticCircle;
	struct Field : IVectorFunction<3> {
		Real scale;
		explicit Field(Real magnitude) : scale(magnitude) {}
		VectorN<Real, 3> operator()(const VectorN<Real, 3>& p) const override {
			return {2 * scale * p[0], 2 * scale * p[1], 0};
		}
	};
	const Real smallScale = std::is_same_v<Real, float> ? std::ldexp(Real(1), -100) : Real(1e-100);
	const Real largeScale = std::is_same_v<Real, float> ? std::ldexp(Real(1), 100) : Real(1e100);
	const Real derivativeTolerance = std::is_same_v<Real, float> ? Precision::EqualityTolerance : Real(1e-7);
	for (const Real scale : {smallScale, Real(1), largeScale}) {
		REQUIRE(std::isnormal(scale));
		const Field field(scale);
		const auto result = PathIntegration::LineIntegralDetailed(field, circle, 0, 2 * Constants::PI, Real(1e-5));
		CHECK(result.converged);
		CHECK_THAT(result.value / scale, WithinAbs(Real(0), derivativeTolerance));
		CHECK(result.error_estimate <= Real(1e-5) * 4 * Constants::PI * scale);
		const auto analytic = PathIntegration::LineIntegralDetailed(field, analyticCircle, 0, 2 * Constants::PI, Real(1e-5));
		CHECK(analytic.converged);
		CHECK_THAT(analytic.value / scale, WithinAbs(Real(0), Real(1e-7)));
		CHECK(analytic.error_estimate <= Real(1e-5) * 4 * Constants::PI * scale);
	}
}

TEST_CASE("Integration finite huge intervals retain representable answers", "[integration-contract]") {
	const auto interval = HugeFiniteInterval();
	const Real a = interval.lower, b = interval.upper;
	REQUIRE(std::isnormal(a));
	REQUIRE(std::isnormal(b));
	REQUIRE(std::isnormal(interval.scale));
	REQUIRE(b > a);
	REQUIRE(a > std::numeric_limits<Real>::max() - b);
	RealFunctionFromStdFunc function([&](Real x) { return x / a; });
	for (const auto method : {TRAP, SIMPSON, ROMBERG, GAUSS10, GAUSS10KRONROD21}) {
		const auto result = IntegrateDetailed(function, a, b, method);
		CHECK(result.IsUsable());
		CHECK(result.error_message.empty());
		CHECK(std::isfinite(result.value));
		CHECK_THAT(result.value / interval.scale, WithinRel(Real(1.05), Precision::EqualityTolerance));
		CHECK(result.function_evaluations > 0);
	}
	const auto kronrod = Integration::IntegrateGK15Detailed([&](Real x) { return function(x); }, a, b);
	CHECK_THAT(kronrod.value / interval.scale, WithinRel(Real(1.05), Precision::EqualityTolerance));
	CHECK(kronrod.status == AlgorithmStatus::Completed);
	CHECK_FALSE(kronrod.converged);
}

TEST_CASE("Integration finite checking and converted failures retain real work", "[integration-contract]") {
	IntegrationConfig config;
	config.exception_policy = EvaluationExceptionPolicy::ConvertToStatus;
	int calls = 0;
	RealFunctionFromStdFunc invalid([&](Real) {
		++calls;
		return std::numeric_limits<Real>::quiet_NaN();
	});
	auto result = IntegrateGauss10Detailed(invalid, 0, 1, config);
	CHECK_FALSE(result.IsUsable());
	CHECK(result.status == AlgorithmStatus::NumericalInstability);
	CHECK(result.function_evaluations == calls);
	CHECK(calls > 0);
	CHECK_FALSE(result.error_message.empty());

	calls = 0;
	RealFunctionFromStdFunc interrupted([&](Real x) {
		if (++calls == 10) throw NumericalMethodError("controlled integration interruption");
		return x;
	});
	result = IntegrateTrapDetailed(interrupted, 0, 1, config);
	CHECK(result.status == AlgorithmStatus::NumericalInstability);
	CHECK(result.function_evaluations == 10);
	CHECK_THAT(result.value, WithinAbs(Real(.5), Precision::EqualityTolerance));
	CHECK(result.elapsed_time_ms > 0);
	CHECK_FALSE(result.converged);
}

TEST_CASE("Nested integration preserves inner failure uncertainty and callback work", "[integration-contract]") {
	ContextRestore restore;
	AlgorithmContext::Get().trapezoidIntegrationMaxSteps = 7;
	int calls2 = 0, calls3 = 0;
	struct Function2 : IScalarFunction<2> {
		int& calls;
		explicit Function2(int& count) : calls(count) {}
		Real operator()(const VectorN<Real, 2>& p) const override { ++calls; return std::pow(p[1], MML::Real(14)); }
	} function2(calls2);
	struct Function3 : IScalarFunction<3> {
		int& calls;
		explicit Function3(int& count) : calls(count) {}
		Real operator()(const VectorN<Real, 3>& p) const override { ++calls; return std::pow(p[2], MML::Real(14)); }
	} function3(calls3);
	const auto result2 = Integrate2DDetailed(function2, TRAP, 0, 1, Lower, Upper);
	const auto result3 = Integrate3DDetailed(function3, TRAP, 0, 1, Lower, Upper, Lower3, Upper3);
	for (const auto* result : {&result2, &result3}) {
		CHECK_FALSE(result->IsUsable());
		CHECK_FALSE(result->converged);
		CHECK(result->error_estimate > Real(1e-4));
		CHECK(std::isfinite(result->value));
	}
	CHECK(result2.function_evaluations == calls2);
	CHECK(result3.function_evaluations == calls3);
}

TEST_CASE("Adaptive ND converges after subdivision and counts discarded parent work", "[integration-contract]") {
	int calls2 = 0, calls3 = 0;
	const Real tolerance = std::is_same_v<Real, float> ? 16 * Precision::RoundoffZeroThreshold : Real(1e-10);
	const auto integrand = [](Real x) {
		if constexpr (std::is_same_v<Real, float>)
			return Chebyshev14(2 * x - 1);
		else
			return std::pow(x, 14);
	};
	// Integral of T_14(2x-1) on [0,1] is -1/195, versus 1/15 for x^14.
	const Real expected = std::is_same_v<Real, float> ? Real(-1) / 195 : Real(1) / 15;
	auto function2 = [&](Real x, Real) { ++calls2; return integrand(x); };
	auto function3 = [&](Real x, Real, Real) { ++calls3; return integrand(x); };
	const auto parent2 = Integration::IntegrateAdaptive2DCore(function2, 0, 1, 0, 1, tolerance, 0, 225);
	const auto parent3 = Integration::IntegrateAdaptive3DCore(function3, 0, 1, 0, 1, 0, 1, tolerance, 0, 3375);
	CHECK_FALSE(parent2.converged);
	CHECK_FALSE(parent3.converged);
	CHECK(parent2.error_estimate > tolerance);
	CHECK(parent3.error_estimate > tolerance);
	CHECK(calls2 == 225);
	CHECK(calls3 == 3375);
	CHECK(parent2.function_evaluations == calls2);
	CHECK(parent3.function_evaluations == calls3);
	calls2 = calls3 = 0;
	const auto result2 = Integration::IntegrateAdaptive2DCore(function2, 0, 1, 0, 1, tolerance, 1, 1125);
	const auto result3 = Integration::IntegrateAdaptive3DCore(function3, 0, 1, 0, 1, 0, 1, tolerance, 1, 30375);
	CHECK(result2.converged);
	CHECK(result3.converged);
	CHECK(result2.status == AlgorithmStatus::Success);
	CHECK(result3.status == AlgorithmStatus::Success);
	CHECK(result2.function_evaluations == calls2);
	CHECK(result3.function_evaluations == calls3);
	CHECK(calls2 == 1125);
	CHECK(calls3 == 30375);
	CHECK(result2.error_estimate <= tolerance);
	CHECK(result3.error_estimate <= tolerance);
	CHECK_THAT(result2.value, WithinAbs(expected, Precision::EqualityTolerance));
	CHECK_THAT(result3.value, WithinAbs(expected, Precision::EqualityTolerance));
}

TEST_CASE("Adaptive ND reserves complete subdivisions before exhausting budget", "[integration-contract]") {
	int calls2 = 0, calls3 = 0;
	const auto result2 = Integration::IntegrateAdaptive2DCore([&](Real x, Real) { ++calls2; return std::pow(x, 14); },
		0, 1, 0, 1, Real(1e-10), 5, 226);
	const auto result3 = Integration::IntegrateAdaptive3DCore([&](Real x, Real, Real) { ++calls3; return std::pow(x, 14); },
		0, 1, 0, 1, 0, 1, Real(1e-10), 5, 3376);
	CHECK(calls2 <= 226);
	CHECK(calls3 <= 3376);
	CHECK(calls2 == 225);
	CHECK(calls3 == 3375);
	CHECK(result2.function_evaluations == calls2);
	CHECK(result3.function_evaluations == calls3);
	CHECK(result2.status == AlgorithmStatus::ResourceLimitExceeded);
	CHECK(result3.status == AlgorithmStatus::ResourceLimitExceeded);
	CHECK_FALSE(result2.converged);
	CHECK_FALSE(result3.converged);
	CHECK_THAT(result2.value, WithinAbs(Real(1) / 15, Precision::EqualityTolerance));
	CHECK_THAT(result3.value, WithinAbs(Real(1) / 15, Precision::EqualityTolerance));
	Integration::AdaptiveConfig3D config;
	config.rule = Integration::GKRule::GK31;
	CHECK_THROWS_AS(Integration::IntegrateAdaptive3DCore([](Real, Real, Real) { return Real(1); }, 0, 1, 0, 1, 0, 1, config),
		NotImplementedError);
}

TEST_CASE("Adaptive GK enforces the aggregate tolerance on cancelling integrals", "[integration-contract]") {
	int calls = 0;
	const Real tolerance = std::is_same_v<Real, float> ? 5 * Precision::RoundoffZeroThreshold : Real(1e-8);
	const auto integrand = [](Real x) {
		if constexpr (std::is_same_v<Real, float>)
			return Chebyshev14(x) + Real(1) / 195;
		else
			return std::exp(10 * x) - std::sinh(Real(10)) / 10;
	};
	if constexpr (std::is_same_v<Real, float>) {
		const auto left = Integration::IntegrateGK15Core(integrand, -1, 0);
		const auto right = Integration::IntegrateGK15Core(integrand, 0, 1);
		const Real leftError = std::max(left.error_estimate, left.roundoff_floor);
		const Real rightError = std::max(right.error_estimate, right.roundoff_floor);
		// Both panels meet a local budget, but their combined error must force more work.
		REQUIRE(leftError <= tolerance);
		REQUIRE(rightError <= tolerance);
		REQUIRE(leftError + rightError > tolerance);
	}
	const auto result = Integration::IntegrateGKAdaptiveCore([&](Real x) {
		if (++calls > 15000) throw ResourceLimitError("bounded cancellation witness");
		return integrand(x);
	}, -1, 1, tolerance, tolerance, 10);
	CHECK(result.function_evaluations == calls);
	CHECK(result.converged);
	CHECK(result.status == AlgorithmStatus::Success);
	CHECK(result.error_estimate <= std::max(tolerance, tolerance * std::abs(result.value)));
	CHECK(std::abs(result.value) <= tolerance);
	if constexpr (std::is_same_v<Real, float>) {
		CHECK(calls > 45);
		calls = 0;
		const Real unresolvedTolerance = Real(1e-8);
		const auto unresolved = Integration::IntegrateGKAdaptiveCore([&](Real x) {
			if (++calls > 15000) throw ResourceLimitError("bounded unresolved cancellation witness");
			return std::exp(10 * x) - std::sinh(Real(10)) / 10;
		}, -1, 1, unresolvedTolerance, unresolvedTolerance, 10);
		CHECK_FALSE(unresolved.converged);
		CHECK(unresolved.status == AlgorithmStatus::ToleranceUnachievable);
		CHECK(unresolved.error_estimate > std::max(unresolvedTolerance, unresolvedTolerance * std::abs(unresolved.value)));
		CHECK(std::isfinite(unresolved.value));
		CHECK(unresolved.function_evaluations == calls);
	}
}

TEST_CASE("Adaptive surface flux uses one global error budget", "[integration-contract]") {
	const RectSurface3D surface(Pnt3Cart(0, 0, 0), Pnt3Cart(1, 0, 0), Pnt3Cart(1, 1, 0), Pnt3Cart(0, 1, 0));
	VectorFunction<3> field([](const VectorN<Real, 3>& p) { return VectorN<Real, 3>{0, 0, p[0] * p[0]}; });
	const auto result = SurfaceIntegration::SurfaceIntegralDetailed(field, surface, Real(.001));
	CHECK(result.converged);
	CHECK(result.error_estimate <= Real(.001));
	CHECK_THAT(result.value, WithinAbs(Real(1) / 3, Real(.001)));
	CHECK(result.function_evaluations > 0);
	CHECK(result.iterations > 0);
	CHECK(result.iterations < 7);
}

TEST_CASE("Parallel integration converted exceptions count every executed callback", "[integration-contract]") {
	std::atomic<int> calls{0};
	ParallelIntegrationConfig config;
	config.thread_count = 2;
	config.partitions = 4;
	config.method = GAUSS10;
	config.exception_policy = EvaluationExceptionPolicy::ConvertToStatus;
	const auto result = IntegrateParallelDetailed([&](Real) -> Real {
		++calls;
		throw NumericalMethodError("controlled partition failure");
	}, 0, 1, config);
	CHECK(result.status == AlgorithmStatus::NumericalInstability);
	CHECK(calls == 4);
	CHECK(result.function_evaluations == calls.load());
	CHECK(result.elapsed_time_ms > 0);
	CHECK_FALSE(result.IsUsable());
}

TEST_CASE("QMC never reports Completed for a nonrepresentable final value", "[integration-contract]") {
	int calls = 0;
	struct Function : IScalarFunction<1> {
		int& calls;
		explicit Function(int& count) : calls(count) {}
		Real operator()(const VectorN<Real, 1>&) const override {
			++calls;
			return std::numeric_limits<Real>::max() / 2;
		}
	} function(calls);
	REQUIRE(std::isfinite(std::numeric_limits<Real>::max() / 2));
	REQUIRE(std::numeric_limits<Real>::max() / 2 > std::numeric_limits<Real>::max() / 10);
	MonteCarloIntegrator<1> integrator;
	CHECK_THROWS_AS(integrator.integrateQuasiRandom(function, VectorN<Real, 1>{0}, VectorN<Real, 1>{10}, 16),
		NumericalMethodError);
	CHECK(calls == 16);
}

TEST_CASE("Integration controls are checked before invoking callbacks", "[integration-contract]") {
	int calls = 0;
	IntegrationConfig config;
	config.exception_policy = EvaluationExceptionPolicy::ConvertToStatus;
	RealFunctionFromStdFunc function([&](Real x) { ++calls; return x; });
	const auto invalid = IntegrateTrapDetailed(function, 0, std::numeric_limits<Real>::infinity(), config);
	CHECK(invalid.status == AlgorithmStatus::InvalidInput);
	CHECK(calls == 0);
	CHECK_THROWS(IntegrateTrapCore(function, 0, 1, 0));
	CHECK(calls == 0);
	const auto insufficient = Integration::IntegrateAdaptive2DDetailed(
		[&](Real, Real) { ++calls; return Real(1); }, 0, 1, 0, 1, config, Real(1e-10), 2, 224);
	CHECK(insufficient.status == AlgorithmStatus::ResourceLimitExceeded);
	CHECK(insufficient.function_evaluations == 0);
	CHECK(std::isnan(insufficient.value));
	CHECK_FALSE(insufficient.error_estimated);
	CHECK(calls == 0);
}

TEST_CASE("Adaptive integration converted interruption retains whole-domain progress", "[integration-contract]") {
	IntegrationConfig config;
	config.exception_policy = EvaluationExceptionPolicy::ConvertToStatus;
	int calls = 0;
	const auto result = Integration::IntegrateAdaptive2DDetailed([&](Real x, Real) {
		if (++calls == 226) throw NumericalMethodError("controlled child failure");
		return std::pow(x, 14);
	}, 0, 1, 0, 1, config, Real(1e-10), 2, 1125);
	CHECK(result.status == AlgorithmStatus::NumericalInstability);
	CHECK(result.function_evaluations == calls);
	CHECK(calls == 226);
	CHECK_THAT(result.value, WithinAbs(Real(1) / 15, Precision::EqualityTolerance));
	CHECK(result.error_estimated);
	CHECK(std::isfinite(result.error_estimate));
	CHECK_FALSE(result.converged);
}

TEST_CASE("Adaptive GK diagnoses unresolved cancellation instead of false success", "[integration-contract]") {
	int calls = 0;
	const auto result = Integration::IntegrateGKAdaptiveCore([&](Real x) {
		if (++calls > 15000) throw ResourceLimitError("bounded precision-limited witness");
		return std::exp(30 * x) - std::sinh(Real(30)) / 30;
	}, -1, 1, Real(1e-8), Real(1e-8), 20);
	CHECK_FALSE(result.IsUsable());
	CHECK(result.status == AlgorithmStatus::ToleranceUnachievable);
	CHECK(result.function_evaluations == calls);
	CHECK(result.error_estimate > Real(1e-8));
	CHECK_FALSE(result.error_message.empty());
}

TEST_CASE("Integration configuration changes checking and presentation without inventing accuracy", "[integration-contract]") {
	int calls = 0;
	IntegrationConfig config;
	config.exception_policy = EvaluationExceptionPolicy::ConvertToStatus;
	config.check_finite = false;
	RealFunctionFromStdFunc invalid([&](Real) { ++calls; return std::numeric_limits<Real>::quiet_NaN(); });
	const auto failed = IntegrateGauss10Detailed(invalid, 0, 1, config);
	CHECK(failed.status == AlgorithmStatus::NumericalInstability);
	CHECK(calls == 10);
	CHECK(failed.function_evaluations == calls);
	config.estimate_error = false;
	const auto fixed = Integration::IntegrateGK21Detailed([](Real x) { return x; }, 0, 1, config);
	CHECK(fixed.status == AlgorithmStatus::Completed);
	CHECK_FALSE(fixed.converged);
	CHECK_FALSE(fixed.error_estimated);
	CHECK(std::isnan(fixed.error_estimate));
	const auto adaptive = Integration::IntegrateAdaptive2DDetailed([](Real x, Real y) { return x * y; }, 0, 1, 0, 1, config);
	CHECK(adaptive.converged);
	CHECK(adaptive.status == AlgorithmStatus::Success);
	CHECK_FALSE(adaptive.error_estimated);
}

TEST_CASE("Computed Gaussian rules validate storage and retain failed sample counts", "[integration-contract]") {
	const auto interval = HugeFiniteInterval();
	RealFunctionFromStdFunc linear([&](Real x) { return x / interval.lower; });
	const auto large = IntegrateGaussLegendreDetailed(linear, interval.lower, interval.upper, 20);
	CHECK_THAT(large.value / interval.scale, WithinRel(Real(1.05), Precision::EqualityTolerance));
	CHECK(large.function_evaluations == 20);
	CHECK_THROWS_AS(GaussQuadratureRule(0), ArgumentError);
	GaussQuadratureRule malformed(2);
	malformed.nodes.resize(1);
	CHECK_THROWS_AS(IntegrateWithRule(linear, malformed), ArgumentError);
	int calls = 0;
	RealFunctionFromStdFunc interrupted([&](Real x) {
		if (++calls == 3) throw NumericalMethodError("controlled Gaussian interruption");
		return x;
	});
	IntegrationConfig config;
	config.exception_policy = EvaluationExceptionPolicy::ConvertToStatus;
	const auto failed = IntegrateGaussLegendreDetailed(interrupted, 0, 1, 20, config);
	CHECK(failed.status == AlgorithmStatus::NumericalInstability);
	CHECK(failed.function_evaluations == 3);
}

TEST_CASE("Parallel Monte Carlo counts failed work rather than the configured sample budget", "[integration-contract]") {
	std::atomic<int> calls{0};
	ParallelMonteCarloConfig config;
	config.num_samples = 16;
	config.block_size = 4;
	config.min_parallel_samples = 0;
	config.thread_count = 2;
	config.exception_policy = EvaluationExceptionPolicy::ConvertToStatus;
	const auto result = IntegrateMonteCarloParallelDetailed<1>([&](const VectorN<Real, 1>&) -> Real {
		++calls;
		throw NumericalMethodError("controlled MC block failure");
	}, VectorN<Real, 1>{0}, VectorN<Real, 1>{1}, config);
	CHECK(result.status == AlgorithmStatus::NumericalInstability);
	CHECK(calls == 4);
	CHECK(result.function_evaluations == calls.load());
	CHECK(result.elapsed_time_ms > 0);
}

TEST_CASE("Open and improper integration retain work and truthful partial payloads", "[integration-contract]") {
	IntegrationConfig config;
	config.exception_policy = EvaluationExceptionPolicy::ConvertToStatus;
	int calls = 0;
	RealFunctionFromStdFunc interrupted([&](Real) {
		if (++calls == 10) throw NumericalMethodError("controlled open refinement failure");
		return Real(1);
	});
	const auto open = IntegrateOpenDetailed(interrupted, 0, 1, config);
	CHECK(open.status == AlgorithmStatus::NumericalInstability);
	CHECK(open.function_evaluations == calls);
	CHECK_THAT(open.value, WithinAbs(Real(1), Precision::EqualityTolerance));
	calls = 0;
	const auto improper = IntegrateUpperInfDetailed(interrupted, 0, config);
	CHECK(improper.status == AlgorithmStatus::NumericalInstability);
	CHECK(improper.function_evaluations == calls);
	CHECK(std::isnan(improper.value));
	CHECK(improper.elapsed_time_ms > 0);
	calls = 0;
	RealFunctionFromStdFunc exponential([&](Real x) { ++calls; return std::exp(-x); });
	const auto completed = IntegrateUpperInfDetailed(exponential, 0, config);
	CHECK(completed.IsUsable());
	CHECK(completed.function_evaluations == calls);
	CHECK_THAT(completed.value, WithinAbs(Real(1), Real(1e-6)));
}

TEST_CASE("Mesh flux aggregates one tolerance across every face and counts actual field calls", "[integration-contract]") {
	const RectSurface3D square(Pnt3Cart(0, 0, 0), Pnt3Cart(1, 0, 0), Pnt3Cart(1, 1, 0), Pnt3Cart(0, 1, 0));
	struct Faces : BodyWithRectSurfaces {
		explicit Faces(const RectSurface3D& face) { _surfaces = {face, face}; }
	} faces(square);
	int calls = 0;
	struct Field : IVectorFunction<3> {
		int& calls;
		explicit Field(int& counter) : calls(counter) {}
		VectorN<Real, 3> operator()(const VectorN<Real, 3>& p) const override {
			++calls;
			return {0, 0, p[0] * p[0]};
		}

	} field(calls);
	const auto result = SurfaceIntegration::SurfaceIntegralDetailed(field, faces, Real(.001));
	CHECK(result.converged);
	CHECK(result.error_estimate <= Real(.001));
	CHECK_THAT(result.value, WithinAbs(Real(2) / 3, Real(.001)));
	CHECK(result.function_evaluations == calls);
	struct Triangles : BodyWithTriangleSurfaces {
		Triangles() {
			_surfaces.emplace_back(Pnt3Cart(0, 0, 0), Pnt3Cart(1, 0, 0), Pnt3Cart(1, 1, 0));
			_surfaces.emplace_back(Pnt3Cart(0, 0, 0), Pnt3Cart(1, 1, 0), Pnt3Cart(0, 1, 0));
		}
	} triangles;
	calls = 0;
	const auto triangular = SurfaceIntegration::SurfaceIntegralDetailed(field, triangles, Real(.001));
	CHECK(triangular.converged);
	CHECK(triangular.error_estimate <= Real(.001));
	CHECK_THAT(triangular.value, WithinAbs(Real(1) / 3, Real(.001)));
	CHECK(triangular.function_evaluations == calls);
}

TEST_CASE("Quadrature public function pointer signatures remain compatible", "[integration-contract]") {
	using Function = Real(*)(Real);
	using Rule = Integration::GKResult(*)(Function, Real, Real);
	const Rule rules[] = {&Integration::IntegrateGK15Core<Function>, &Integration::IntegrateGK21Core<Function>,
		&Integration::IntegrateGK31Core<Function>};
	for (const auto rule : rules)
		CHECK_THAT(rule(+[](Real x) { return x; }, 0, 1).value, WithinAbs(Real(.5), Precision::EqualityTolerance));
	using ComputedRule = Real(*)(const IRealFunction&, const GaussQuadratureRule&);
	const ComputedRule computed = &IntegrateWithRule;
	RealFunction linear([](Real x) { return x; });
	CHECK_THAT(computed(linear, GaussLegendreRule(4, 0, 1)), WithinAbs(Real(.5), Precision::EqualityTolerance));
}

TEST_CASE("Short Romberg budgets retain partial estimates without inventing extrapolation", "[integration-contract]") {
	ContextRestore restore;
	AlgorithmContext::Get().rombergIntegrationMaxSteps = 2;
	AlgorithmContext::Get().rombergIntegrationUsedPnts = 5;
	RealFunction constant([](Real) { return Real(1); });
	const auto closed = IntegrateRombergDetailed(constant, 0, 1);
	const auto open = IntegrateOpenDetailed(constant, 0, 1);
	for (const auto* result : {&closed, &open}) {
		CHECK(result->status == AlgorithmStatus::MaxIterationsExceeded);
		CHECK(result->value == 1);
		CHECK_FALSE(result->error_estimated);
		CHECK(result->function_evaluations == 3);
	}
}
