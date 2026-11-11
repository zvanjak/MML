#include <catch2/catch_all.hpp>
#include "../../TestPrecision.h"
#ifndef MML_SINGLE_HEADER_INCLUDE
#include <mml/core/Derivation/DerivationRealFunction.h>
#include <mml/core/Derivation/DerivationParametricCurve.h>
#include <mml/core/Derivation/Hessians.h>
#endif

#include <array>
#include <vector>

using namespace MML;

namespace {
	using Point = std::array<Real, 2>;

	struct RecordedReal : IRealFunction {
		mutable std::vector<Real> samples;
		Real amplitude = 1;
		int degree = 0;
		bool bounded = false;
		Real operator()(Real x) const override {
			samples.push_back(x);
			if (bounded && (x < 0 || x > 1)) throw DomainError("Second derivative left [0,1]");
			return amplitude * (degree ? std::pow(x, MML::Real(degree)) : std::exp(x));
		}
	};

	struct RecordedCurve : IParametricCurve<3> {
		RecordedReal function;
		Real getMinT() const override { return function.bounded ? Real(0) : Constants::NegInf; }
		Real getMaxT() const override { return function.bounded ? Real(1) : Constants::PosInf; }
		VectorN<Real, 3> operator()(Real x) const override {
			const Real value = function(x);
			return {value, 2 * value, -2 * value};
		}
	};

	struct RecordedScalar : IScalarFunction<2> {
		mutable std::vector<Point> samples;
		Real amplitude = 1;
		int degree = 0;
		bool mixed = false, bounded = false;
		Real operator()(const VectorN<Real, 2>& x) const override {
			samples.push_back({x[0], x[1]});
			if (bounded && (x[0] < 0 || x[0] > 1 || x[1] < 0 || x[1] > 1))
				throw DomainError("Second partial left [0,1]^2");
			return amplitude * (degree ? (mixed ? std::pow(x[0], MML::Real(degree - 1)) * x[1] : std::pow(x[0], MML::Real(degree)))
				: std::exp(x[0] + x[1]));
		}
	};

	Real BaseStep(int order) {
		return order == 2 ? Precision::FourthRootEps : std::pow(Constants::Eps, Real(1) / 6);
	}
	Real ExponentialBudget(int order, Real h, bool mixed = false) {
		const Real balance = order == 2 ? std::sqrt(Constants::Eps) : std::pow(Constants::Eps, Real(2) / 3);
		// Retained coordinate scaling enlarges the known exponential Taylor term as h^p.
		const Real truncation = order == 2 ? h * h / (mixed ? 3 : 12)
			: h * h * h * h * (mixed ? Real(8) / 45 : Real(1) / 90);
		return 16 * balance + 2 * truncation;
	}
	Real RealSecond(int order, const IRealFunction& f, Real x, Real* error = nullptr) {
		return order == 2 ? Derivation::NSecDer2(f, x, error) : Derivation::NSecDer4(f, x, error);
	}
	Real RealExplicit(int order, const IRealFunction& f, Real x, Real h, Real* error = nullptr) {
		return order == 2 ? Derivation::NSecDer2(f, x, h, error) : Derivation::NSecDer4(f, x, h, error);
	}
	VectorN<Real, 3> CurveSecond(int order, const IParametricCurve<3>& f, Real x, Real* error = nullptr) {
		return order == 2 ? Derivation::NSecDer2(f, x, error) : Derivation::NSecDer4(f, x, error);
	}
	VectorN<Real, 3> CurveExplicit(int order, const IParametricCurve<3>& f, Real x, Real h, Real* error = nullptr) {
		return order == 2 ? Derivation::NSecDer2(f, x, h, error) : Derivation::NSecDer4(f, x, h, error);
	}
	Real ScalarSecond(int order, const IScalarFunction<2>& f, int i, int j, const VectorN<Real, 2>& p,
		Real* error = nullptr) {
		return order == 2 ? Derivation::NSecDer2Partial(f, i, j, p, error)
			: Derivation::NSecDer4Partial(f, i, j, p, error);
	}
	Real ScalarExplicit(int order, const IScalarFunction<2>& f, int i, int j, const VectorN<Real, 2>& p,
		Real h, Real* error = nullptr) {
		return order == 2 ? Derivation::NSecDer2Partial(f, i, j, p, h, error)
			: Derivation::NSecDer4Partial(f, i, j, p, h, error);
	}
	void CheckCoordinate(Real actual, Real expected) {
		CHECK_THAT(actual, Catch::Matchers::WithinAbs(expected,
			8 * Constants::Eps * std::max(Real(1), std::abs(expected))));
	}
	std::vector<int> Offsets(int order, bool error) {
		std::vector<int> offsets{0, 1, -1};
		if (order == 4 || error) { offsets.push_back(2); offsets.push_back(-2); }
		if (order == 4 && error) { offsets.push_back(3); offsets.push_back(-3); }
		return offsets;
	}
	void CheckLine(std::vector<Real> actual, Real x, Real h, int order, bool error) {
		std::vector<Real> expected;
		for (int offset : Offsets(order, error)) expected.push_back(x + offset * h);
		REQUIRE(actual.size() == expected.size());
		std::sort(actual.begin(), actual.end());
		std::sort(expected.begin(), expected.end());
		for (size_t index = 0; index < actual.size(); ++index) CheckCoordinate(actual[index], expected[index]);
	}
	std::vector<Point> ExpectedPoints(const VectorN<Real, 2>& p, int i, int j, Real h, int order, bool error) {
		std::vector<Point> expected;
		if (i == j) {
			for (int offset : Offsets(order, error)) {
				Point x{p[0], p[1]};
				x[i] += offset * h;
				expected.push_back(x);
			}
		} else {
			std::vector<Real> steps{h};
			if (order == 4) steps.push_back(2 * h);
			if (error) steps.push_back(h / 2);
			for (Real step : steps)
				for (int a : {-1, 1}) for (int b : {-1, 1})
					expected.push_back({p[0] + a * step, p[1] + b * step});
		}
		return expected;
	}
	void CheckPoints(std::vector<Point> actual, std::vector<Point> expected) {
		REQUIRE(actual.size() == expected.size());
		std::sort(actual.begin(), actual.end());
		std::sort(expected.begin(), expected.end());
		for (size_t index = 0; index < actual.size(); ++index)
			for (int axis = 0; axis < 2; ++axis) CheckCoordinate(actual[index][axis], expected[index][axis]);
	}
}

TEST_CASE("Real second derivative defaults account for derivative and accuracy orders", "[second-steps]") {
	for (int order : {2, 4}) for (Real x : {Real(0), Real(0.375), Real(-3)}) for (bool error : {false, true}) {
		CAPTURE(order, x, error);
		const Real h = Derivation::ScaleStep(BaseStep(order), x);
		RecordedReal f;
		Real estimate;
		RealSecond(order, f, x, error ? &estimate : nullptr);
		CheckLine(f.samples, x, h, order, error);
		f.samples.clear();
		DerivativeConfig config;
		config.estimate_error = error;
		const auto detailed = order == 2 ? Derivation::NSecDer2Detailed(f, x, config)
			: Derivation::NSecDer4Detailed(f, x, config);
		CHECK(detailed.status == AlgorithmStatus::Success);
		CHECK(detailed.function_evaluations == int(Offsets(order, error).size()));
		CheckCoordinate(detailed.step_used, h);
		CheckLine(f.samples, x, h, order, error);
		f.samples.clear();
		auto callable = [&](Real t) { return f(t); };
		if (order == 2) Derivation::NSecDer2(callable, x, error ? &estimate : nullptr);
		else Derivation::NSecDer4(callable, x, error ? &estimate : nullptr);
		CheckLine(f.samples, x, h, order, error);
		f.samples.clear();
		const auto callableDetailed = order == 2 ? Derivation::NSecDer2Detailed(callable, x, config)
			: Derivation::NSecDer4Detailed(callable, x, config);
		CHECK(callableDetailed.function_evaluations == detailed.function_evaluations);
		CheckCoordinate(callableDetailed.step_used, h);
		CheckLine(f.samples, x, h, order, error);
	}
}

TEST_CASE("Curve and scalar second defaults retain coordinate scaling and sample costs", "[second-steps]") {
	for (int order : {2, 4}) for (bool error : {false, true}) {
		CAPTURE(order, error);
		RecordedCurve curve;
		Real estimate;
		const Real x = -3, h = Derivation::ScaleStep(BaseStep(order), x);
		CurveSecond(order, curve, x, error ? &estimate : nullptr);
		CheckLine(curve.function.samples, x, h, order, error);
		for (const auto& indices : {std::pair{0, 0}, std::pair{1, 1}, std::pair{0, 1}, std::pair{1, 0}}) {
			RecordedScalar scalar;
			const VectorN<Real, 2> p{2, -3};
			const Real step = Derivation::ScaleStep(BaseStep(order),
				std::max(std::abs(p[indices.first]), std::abs(p[indices.second])));
			ScalarSecond(order, scalar, indices.first, indices.second, p, error ? &estimate : nullptr);
			CheckPoints(scalar.samples, ExpectedPoints(p, indices.first, indices.second, step, order, error));
		}
	}
}

TEST_CASE("First and third derivative defaults are unchanged by second-step repair", "[second-steps]") {
	CHECK(Derivation::NDer2_h == std::pow(3 * Constants::Eps, Real(1) / Real(3)));
	CHECK(Derivation::NDer4_h == std::pow(Real(11.25) * Constants::Eps, Real(1) / Real(5)));
	const Real x = Real(0.375);
	RecordedReal f;
	RecordedCurve curve;
	RecordedScalar scalar;
	const VectorN<Real, 2> p{x, Real(-2)};
	CHECK(Derivation::NDer2(f, x) == Derivation::NDer2(f, x, Derivation::ScaleStep(Derivation::NDer2_h, x)));
	CHECK(Derivation::NDer4(f, x) == Derivation::NDer4(f, x, Derivation::ScaleStep(Derivation::NDer4_h, x)));
	CHECK(Derivation::NDer2(curve, x) == Derivation::NDer2(curve, x, Derivation::ScaleStep(Derivation::NDer2_h, x)));
	CHECK(Derivation::NDer4(curve, x) == Derivation::NDer4(curve, x, Derivation::ScaleStep(Derivation::NDer4_h, x)));
	CHECK(Derivation::NDer2Partial(scalar, 0, p) == Derivation::NDer2Partial(scalar, 0, p, Derivation::ScaleStep(Derivation::NDer2_h, x)));
	CHECK(Derivation::NDer4Partial(scalar, 0, p) == Derivation::NDer4Partial(scalar, 0, p, Derivation::ScaleStep(Derivation::NDer4_h, x)));
	CHECK(Derivation::NThirdDer2(f, x) == Derivation::NThirdDer2(f, x, Derivation::ScaleStep(Derivation::NDer4_h, x)));
	CHECK(Derivation::NThirdDer4(curve, x) == Derivation::NThirdDer4(curve, x, Derivation::ScaleStep(Derivation::NDer4_h, x)));
}

TEST_CASE("Default non-surface second derivatives meet roundoff and truncation budgets", "[second-steps]") {
	for (int order : {2, 4}) for (Real amplitude : {Real(1e-8), Real(1), Real(1e8)})
		for (Real x : {Real(-0.7), Real(0.3), Real(1.4), Real(3.25)}) {
			CAPTURE(order, amplitude, x);
			const Real balance = order == 2 ? std::sqrt(Constants::Eps) : std::pow(Constants::Eps, Real(2) / 3);
			const Real truth = amplitude * std::exp(x);
			RecordedReal f;
			f.amplitude = amplitude;
			RecordedCurve curve;
			curve.function.amplitude = amplitude;
			RecordedScalar scalar;
			scalar.amplitude = amplitude;
			Real estimate;
			const Real h = Derivation::ScaleStep(BaseStep(order), x);
			CHECK(std::abs(RealSecond(order, f, x, &estimate) / truth - 1) <= ExponentialBudget(order, h));
			CHECK(estimate >= 0);
			CHECK(estimate / truth <= 32 * balance);
			const auto acceleration = CurveSecond(order, curve, x, &estimate);
			CHECK((acceleration / truth - VectorN<Real, 3>{1, 2, -2}).NormL2() / 3 <= ExponentialBudget(order, h));
			CHECK(estimate >= 0);
			CHECK(estimate / (3 * truth) <= 32 * balance);
			const VectorN<Real, 2> p{x, 0};
			for (const auto& indices : {std::pair{0, 0}, std::pair{1, 1}, std::pair{0, 1}}) {
				const Real step = Derivation::ScaleStep(BaseStep(order),
					std::max(std::abs(p[indices.first]), std::abs(p[indices.second])));
				CHECK(std::abs(ScalarSecond(order, scalar, indices.first, indices.second, p, &estimate) / truth - 1) <=
					ExponentialBudget(order, step, indices.first != indices.second));
				CHECK(estimate >= 0);
				CHECK(estimate / truth <= 32 * balance);
			}
		}
}

TEST_CASE("Explicit second steps preserve convergence estimates and returned values", "[second-steps]") {
	for (int order : {2, 4}) {
		RecordedReal f;
		f.degree = order + 2;
		RecordedCurve curve;
		curve.function.degree = f.degree;
		RecordedScalar scalar;
		scalar.degree = f.degree;
		std::array<Real, 4> previous{};
		for (int refinement = 0; refinement < 3; ++refinement) {
			const Real h = std::ldexp(Real(1) / 4, -refinement);
			CAPTURE(order, h);
			std::array<Real, 4> estimates;
			const Real real = RealExplicit(order, f, 0, h, &estimates[0]);
			const auto vector = CurveExplicit(order, curve, 0, h, &estimates[1]);
			const Real pure = ScalarExplicit(order, scalar, 0, 0, {0, 0}, h, &estimates[2]);
			scalar.mixed = true;
			const Real mixed = ScalarExplicit(order, scalar, 0, 1, {0, 0}, h, &estimates[3]);
			scalar.mixed = false;
			CHECK(real == RealExplicit(order, f, 0, h));
			CHECK(vector == CurveExplicit(order, curve, 0, h));
			CHECK(pure == ScalarExplicit(order, scalar, 0, 0, {0, 0}, h));
			scalar.mixed = true;
			CHECK(mixed == ScalarExplicit(order, scalar, 0, 1, {0, 0}, h));
			scalar.mixed = false;
			const std::array<Real, 4> actual{std::abs(real), vector.NormL2(), std::abs(pure), std::abs(mixed)};
			for (size_t method = 0; method < actual.size(); ++method) {
				CHECK(actual[method] > 0);
				CHECK(estimates[method] >= actual[method] / 2);
				CHECK(estimates[method] <= 2 * actual[method]);
				if (refinement != 0) {
					CHECK(previous[method] / actual[method] > (order == 2 ? 3 : 12));
					CHECK(previous[method] / actual[method] < (order == 2 ? 5 : 20));
				}
				previous[method] = actual[method];
			}
		}
	}
}

TEST_CASE("Detailed second-step overrides remain absolute and explicit steps take priority", "[second-steps]") {
	for (int order : {2, 4}) {
		RecordedReal f;
		DerivativeConfig config;
		config.step = Real(1) / 8;
		config.estimate_error = false;
		auto result = order == 2 ? Derivation::NSecDer2Detailed(f, Real(4), config)
			: Derivation::NSecDer4Detailed(f, Real(4), config);
		CHECK(result.step_used == config.step);
		CheckLine(f.samples, 4, config.step, order, false);
		f.samples.clear();
		const Real explicitStep = Real(1) / 16;
		config.estimate_error = true;
		result = order == 2 ? Derivation::NSecDer2Detailed(f, Real(4), explicitStep, config)
			: Derivation::NSecDer4Detailed(f, Real(4), explicitStep, config);
		CHECK(result.step_used == explicitStep);
		CheckLine(f.samples, 4, explicitStep, order, true);
	}
}

TEST_CASE("Default aliases and numerical Hessians inherit the second-derivative scale", "[second-steps]") {
	RecordedReal f;
	Real estimate;
	Derivation::DeriveSec(f, Real(2), &estimate);
	CheckLine(f.samples, 2, 2 * BaseStep(4), 4, true);
	RecordedCurve curve;
	Derivation::DeriveCurveSec<3>(curve, Real(2), &estimate);
	CheckLine(curve.function.samples, 2, 2 * BaseStep(4), 4, true);
	RecordedScalar scalar;
	const VectorN<Real, 2> p{2, -3};
	const auto hessian = Derivation::CalcHessian(scalar, p);
	auto expected = ExpectedPoints(p, 0, 0, 2 * BaseStep(4), 4, false);
	for (const auto& indices : {std::pair{0, 1}, std::pair{1, 1}}) {
		const auto samples = ExpectedPoints(p, indices.first, indices.second, 3 * BaseStep(4), 4, false);
		expected.insert(expected.end(), samples.begin(), samples.end());
	}
	CheckPoints(scalar.samples, expected);
	CHECK(hessian(0, 1) == hessian(1, 0));
	CHECK(hessian(0, 0) == Derivation::NSecDer4Partial(scalar, 0, 0, p));
	CHECK(hessian(1, 1) == Derivation::NSecDer4Partial(scalar, 1, 1, p));
	const Real step = Real(1) / 8;
	const auto explicitHessian = Derivation::CalcHessian(scalar, p, step);
	CHECK(explicitHessian(0, 1) == Derivation::NSecDer4Partial(scalar, 0, 1, p, step));
}

TEST_CASE("Second defaults retain centered boundaries and detailed exception policy", "[second-steps]") {
	for (int order : {2, 4}) {
		const Real h = BaseStep(order);
		RecordedReal f;
		f.bounded = true;
		RecordedCurve curve;
		curve.function.bounded = true;
		RecordedScalar scalar;
		scalar.bounded = true;
		CHECK_THROWS_AS(RealSecond(order, f, 0), DomainError);
		CHECK_THROWS_AS(CurveSecond(order, curve, 0), DomainError);
		CHECK_THROWS_AS(ScalarSecond(order, scalar, 0, 0, {0, Real(0.5)}), DomainError);
		const Real center = (order == 2 ? Real(1.5) : Real(2.5)) * h;
		CHECK_NOTHROW(RealSecond(order, f, center));
		CHECK_NOTHROW(CurveSecond(order, curve, center));
		Real error = -17;
		CHECK_THROWS_AS(RealSecond(order, f, center, &error), DomainError);
		CHECK(error == Real(-17));
		CHECK_THROWS_AS(CurveSecond(order, curve, center, &error), DomainError);
		CHECK(error == Real(-17));
		DerivativeConfig config;
		config.exception_policy = EvaluationExceptionPolicy::ConvertToStatus;
		const auto result = order == 2 ? Derivation::NSecDer2Detailed(f, Real(0), config)
			: Derivation::NSecDer4Detailed(f, Real(0), config);
		CHECK(result.status == AlgorithmStatus::AlgorithmSpecificFailure);
		CheckCoordinate(result.step_used, h);
	}
}
