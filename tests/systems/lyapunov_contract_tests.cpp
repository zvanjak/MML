#include "../TestPrecision.h"
#ifndef MML_SINGLE_HEADER_INCLUDE
#include <mml/systems/DiscreteMaps.h>
#include <mml/systems/DynamicalSystem/LyapunovAnalysis.h>
#include <mml/systems/DynamicalSystem/DynamicalSystemBase.h>
#endif
#include <catch2/catch_all.hpp>
#include <algorithm>
#include <cmath>
#include <functional>
#include <limits>
#include <stdexcept>
#include <type_traits>

using namespace MML;
using namespace MML::Systems;

namespace {
	struct LyapunovStretches {
		Real small, large;
	};
	LyapunovStretches ExtremeLyapunovStretches() {
		constexpr int exponent = 3 * std::numeric_limits<Real>::max_exponent / 4;
		return {std::ldexp(Real(1), -exponent), std::ldexp(Real(1), exponent)};
	}
	Real ResolvedStrongFlowStep() {
		if constexpr (std::is_same_v<Real, float>)
			return Real(1) / 2048;
		else
			return Real(0.01);
	}
	template<int N>
	class LyapunovLinearMap : public MML::IDiscreteMap<N> {
	public:
		Matrix<Real> matrix;
		mutable int mapCalls = 0;
		mutable int jacobianCalls = 0;
		std::function<void(Vector<Real>&)> mapOverride;
		std::function<void(Matrix<Real>&)> jacobianOverride;
		explicit LyapunovLinearMap(const Matrix<Real>& value) : matrix(value) {}
		int getNumParam() const override { return 0; }
		Real getParam(int) const override { throw ArgumentError("No parameters"); }
		void setParam(int, Real) override { throw ArgumentError("No parameters"); }
		void map(const Vector<Real>& x, Vector<Real>& next) const override {
			++mapCalls;
			if (mapOverride) mapOverride(next);
			else next = matrix * x;
		}
		void jacobian(const Vector<Real>&, Matrix<Real>& value) const override {
			++jacobianCalls;
			if (jacobianOverride) jacobianOverride(value);
			else value = matrix;
		}
	};
	template<int N>
	class LyapunovLinearFlow : public DynamicalSystemBase<N, 0> {
	public:
		Matrix<Real> matrix;
		int dimension = N;
		mutable int derivativeCalls = 0;
		mutable int jacobianCalls = 0;
		std::function<void(Real, const Vector<Real>&, Vector<Real>&)> derivativeOverride;
		std::function<void(Real, const Vector<Real>&, Matrix<Real>&)> jacobianOverride;
		explicit LyapunovLinearFlow(const Matrix<Real>& value) : matrix(value) {}
		int getDim() const override { return dimension; }
		void derivs(Real t, const Vector<Real>& state, Vector<Real>& derivative) const override {
			++derivativeCalls;
			if (derivativeOverride) derivativeOverride(t, state, derivative);
			else derivative = matrix * state;
		}
		void jacobian(Real t, const Vector<Real>& state, Matrix<Real>& value) const override {
			++jacobianCalls;
			if (jacobianOverride) jacobianOverride(t, state, value);
			else value = matrix;
		}
	};
}

TEST_CASE("Discrete Lyapunov spectra are sorted before maximum and chaos classification", "[LyapunovContracts]") {
	LyapunovLinearMap<2> map(Matrix<Real>(2, 2, {Real(0.5), 0, 0, 2}));
	const auto result = DiscreteMapLyapunov<2>::Compute(map, Vector<Real>{0, 0}, 8, 0);
	REQUIRE(result.exponents.size() == 2);
	CHECK(IsNearlyEqual(result.exponents[0], std::log(Real(2))));
	CHECK(IsNearlyEqual(result.exponents[1], -std::log(Real(2))));
	CHECK(result.maxExponent == result.exponents[0]);
	CHECK(result.isChaotic);
}

TEST_CASE("Discrete Lyapunov logarithms preserve tiny and huge nonzero stretches", "[LyapunovContracts]") {
	const auto extremes = ExtremeLyapunovStretches();
	for (Real stretch : {Real(1e-40), extremes.small, extremes.large}) {
		CAPTURE(stretch);
		REQUIRE(stretch > 0);
		REQUIRE(std::isfinite(stretch));
		LyapunovLinearMap<1> map(Matrix<Real>(1, 1, {stretch}));
		const auto result = DiscreteMapLyapunov<1>::Compute(map, Vector<Real>{0}, 8, 0);
		REQUIRE(result.exponents.size() == 1);
		CHECK(std::isfinite(result.maxExponent));
		CHECK(IsNearlyEqual(result.maxExponent, std::log(stretch)));
	}
}

TEST_CASE("Exactly zero discrete stretch is an explicit singular failure", "[LyapunovContracts]") {
	LogisticMap map(2);
	CHECK_THROWS_AS(DiscreteMapLyapunov<1>::Compute(map, Vector<Real>{Real(0.5)}, 8, 0), SingularMatrixError);
	const Real erased = std::numeric_limits<Real>::denorm_min() / 2;
	REQUIRE(erased == 0);
	LyapunovLinearMap<1> zero(Matrix<Real>(1, 1, {erased}));
	CHECK_THROWS_AS(DiscreteMapLyapunov<1>::Compute(zero, Vector<Real>{0}, 8, 0), SingularMatrixError);
	CHECK(zero.jacobianCalls == 1);
	CHECK(zero.mapCalls == 0);
	for (Real invalid : {std::numeric_limits<Real>::infinity(), std::numeric_limits<Real>::quiet_NaN()}) {
		LyapunovLinearMap<1> nonfinite(Matrix<Real>(1, 1, {invalid}));
		CHECK_THROWS_AS(DiscreteMapLyapunov<1>::Compute(nonfinite, Vector<Real>{0}, 8, 0), NumericalMethodError);
		CHECK(nonfinite.jacobianCalls == 1);
		CHECK(nonfinite.mapCalls == 0);
	}
}

TEST_CASE("Dependent discrete tangent columns cannot fabricate a neutral exponent", "[LyapunovContracts]") {
	LyapunovLinearMap<2> map(Matrix<Real>(2, 2, {1, 1, 0, 0}));
	CHECK_THROWS_AS(DiscreteMapLyapunov<2>::Compute(map, Vector<Real>{0, 0}, 8, 0), SingularMatrixError);
}

TEST_CASE("Continuous strong contraction is resolved with either requested QR interval", "[LyapunovContracts]") {
	if constexpr (std::is_same_v<Real, float>) {
		// The RK4 stability polynomial independently checks the fixture's rate-bias budget.
		const double step = ResolvedStrongFlowStep();
		for (double rate : {100.0, -100.0, -120.0}) {
			const double z = rate * step;
			const double increment = z + z * z / 2 + z * z * z / 6 + z * z * z * z / 24;
			const double bias = std::abs(std::log1p(increment) / step - rate);
			REQUIRE(bias < Precision::CoarseTolerance / 4);
		}
	}
	for (Real interval : {Real(1), Real(0.125)}) {
		CAPTURE(interval);
		LyapunovLinearFlow<1> flow(Matrix<Real>(1, 1, {-100}));
		const auto result = LyapunovAnalyzer::Compute(flow, Vector<Real>{0}, 1, interval, ResolvedStrongFlowStep());
		CHECK(IsWithinAbs(result.maxExponent, Real(-100), Precision::CoarseTolerance));
		CHECK(IsWithinAbs(result.sum, Real(-100), Precision::CoarseTolerance));
		CHECK_FALSE(result.isChaotic);
		CHECK(result.kaplanYorkeDimension == 0);
	}
}

TEST_CASE("QR logarithms do not require a representable unscaled Euclidean norm", "[LyapunovContracts]") {
	const Real scale = std::numeric_limits<Real>::max() * Real(0.75);
	LyapunovLinearMap<2> map(Matrix<Real>(2, 2, {scale, scale, scale, -scale}));
	const auto result = DiscreteMapLyapunov<2>::Compute(map, Vector<Real>{0, 0}, 1, 0);
	const Real expected = std::log(scale) + std::log(Real(2)) / 2;
	for (Real exponent : result.exponents) {
		CHECK(std::isfinite(exponent));
		CHECK(IsNearlyEqual(exponent, expected));
	}
}

TEST_CASE("Column-relative QR scaling preserves a mixed extreme discrete spectrum", "[LyapunovContracts]") {
	const auto extremes = ExtremeLyapunovStretches();
	REQUIRE(std::isnormal(extremes.small));
	REQUIRE(std::isnormal(extremes.large));
	REQUIRE(extremes.small * extremes.small == 0);
	REQUIRE_FALSE(std::isfinite(extremes.large * extremes.large));
	REQUIRE(extremes.small / extremes.large == 0);
	LyapunovLinearMap<2> map(Matrix<Real>(2, 2, {extremes.small, 0, 0, extremes.large}));
	const auto result = DiscreteMapLyapunov<2>::Compute(map, Vector<Real>{0, 0}, 8, 0);
	REQUIRE(result.exponents.size() == 2);
	CHECK(IsNearlyEqual(result.exponents[0], std::log(extremes.large)));
	CHECK(IsNearlyEqual(result.exponents[1], std::log(extremes.small)));
	CHECK(result.maxExponent == result.exponents[0]);
	CHECK(result.isChaotic);
}

TEST_CASE("Scalar discrete QR handles the finite representable endpoints", "[LyapunovContracts]") {
	for (Real scale : {std::numeric_limits<Real>::denorm_min(), std::numeric_limits<Real>::min(), std::numeric_limits<Real>::max()}) {
		REQUIRE(scale > 0);
		REQUIRE(std::isfinite(scale));
		LyapunovLinearMap<1> map(Matrix<Real>(1, 1, {scale}));
		const auto result = DiscreteMapLyapunov<1>::Compute(map, Vector<Real>{0}, 2, 0);
		CHECK(std::isfinite(result.maxExponent));
		CHECK(IsNearlyEqual(result.maxExponent, std::log(scale)));
	}
}

TEST_CASE("Discrete rank diagnostics are column relative rather than absolute", "[LyapunovContracts]") {
	const Real smallScale = std::is_same_v<Real, float> ? std::ldexp(Real(1), -100) : Real(1e-150);
	const Real largeScale = std::is_same_v<Real, float> ? std::ldexp(Real(1), 100) : Real(1e150);
	for (Real scale : {smallScale, largeScale}) {
		REQUIRE(std::isnormal(scale));
		REQUIRE(scale * Precision::Eps > 0);
		REQUIRE(std::isfinite(scale * Precision::Eps));
		LyapunovLinearMap<2> map(Matrix<Real>(2, 2, {scale, scale, 0, scale * Precision::Eps}));
		CHECK_THROWS_AS(DiscreteMapLyapunov<2>::Compute(map, Vector<Real>{0, 0}, 1, 0), SingularMatrixError);
		const Real perpendicular = scale * std::sqrt(Precision::Eps);
		LyapunovLinearMap<2> independent(Matrix<Real>(2, 2, {scale, scale, 0, perpendicular}));
		const auto result = DiscreteMapLyapunov<2>::Compute(independent, Vector<Real>{0, 0}, 1, 0);
		REQUIRE(result.exponents.size() == 2);
		CHECK(IsNearlyEqual(result.exponents[0], std::log(scale)));
		CHECK(IsNearlyEqual(result.exponents[1], std::log(perpendicular)));
		CHECK(independent.jacobianCalls == 1);
		CHECK(independent.mapCalls == 1);
	}
}

TEST_CASE("Late discrete rank loss invalidates the full spectrum", "[LyapunovContracts]") {
	LyapunovLinearMap<1> map(Matrix<Real>(1, 1, {2}));
	map.jacobianOverride = [&](Matrix<Real>& value) { value(0, 0) = map.jacobianCalls < 3 ? Real(2) : Real(0); };
	CHECK_THROWS_AS(DiscreteMapLyapunov<1>::Compute(map, Vector<Real>{0}, 8, 0), SingularMatrixError);
	CHECK(map.jacobianCalls == 3);
	CHECK(map.mapCalls == 2);
}

TEST_CASE("Discrete Lyapunov validates requests before invoking callbacks", "[LyapunovContracts]") {
	LyapunovLinearMap<1> map(Matrix<Real>(1, 1, {1}));
	CHECK_THROWS_AS(DiscreteMapLyapunov<1>::Compute(map, Vector<Real>{0}, 0, 0), ArgumentError);
	CHECK_THROWS_AS(DiscreteMapLyapunov<1>::Compute(map, Vector<Real>{0}, -1, 0), ArgumentError);
	CHECK_THROWS_AS(DiscreteMapLyapunov<1>::Compute(map, Vector<Real>{0}, 1, -1), ArgumentError);
	CHECK_THROWS_AS(DiscreteMapLyapunov<1>::Compute(map, Vector<Real>{}, 1, 0), VectorDimensionError);
	CHECK_THROWS_AS(DiscreteMapLyapunov<1>::Compute(map, Vector<Real>{0, 0}, 1, 0), VectorDimensionError);
	CHECK_THROWS_AS(DiscreteMapLyapunov<1>::Compute(map, Vector<Real>{std::numeric_limits<Real>::quiet_NaN()}, 1, 0), ArgumentError);
	CHECK(map.mapCalls == 0);
	CHECK(map.jacobianCalls == 0);
}

TEST_CASE("Discrete Lyapunov rejects malformed and unwritten callback outputs", "[LyapunovContracts]") {
	for (int transient : {0, 1}) {
		LyapunovLinearMap<1> map(Matrix<Real>(1, 1, {1}));
		map.mapOverride = [](Vector<Real>& next) { next.resize(0); };
		CHECK_THROWS_AS(DiscreteMapLyapunov<1>::Compute(map, Vector<Real>{0}, 1, transient), VectorDimensionError);
		map.mapOverride = [](Vector<Real>&) {};
		CHECK_THROWS_AS(DiscreteMapLyapunov<1>::Compute(map, Vector<Real>{0}, 1, transient), NumericalMethodError);
	}
	LyapunovLinearMap<1> map(Matrix<Real>(1, 1, {1}));
	map.jacobianOverride = [](Matrix<Real>& matrix) { matrix.resize(2, 1); };
	CHECK_THROWS_AS(DiscreteMapLyapunov<1>::Compute(map, Vector<Real>{0}, 1, 0), MatrixDimensionError);
	map.jacobianOverride = [](Matrix<Real>&) {};
	CHECK_THROWS_AS(DiscreteMapLyapunov<1>::Compute(map, Vector<Real>{0}, 1, 0), NumericalMethodError);
	map.jacobianOverride = [](Matrix<Real>& matrix) { matrix(0, 0) = std::numeric_limits<Real>::infinity(); };
	CHECK_THROWS_AS(DiscreteMapLyapunov<1>::Compute(map, Vector<Real>{0}, 1, 0), NumericalMethodError);
}

TEST_CASE("Continuous adaptive QR resolves expansion and mixed growth without changing physical time", "[LyapunovContracts]") {
	LyapunovLinearFlow<2> flow(Matrix<Real>(2, 2, {100, 0, 0, -100}));
	const auto result = LyapunovAnalyzer::Compute(flow, Vector<Real>{0, 0}, 1, 1, ResolvedStrongFlowStep());
	CHECK(IsWithinAbs(result.exponents[0], Real(100), Precision::CoarseTolerance));
	CHECK(IsWithinAbs(result.exponents[1], Real(-100), Precision::CoarseTolerance));
	CHECK(IsWithinAbs(result.sum, Real(0), Precision::CoarseTolerance));
	CHECK(IsWithinAbs(result.kaplanYorkeDimension, Real(2), Precision::CoarseTolerance));
	CHECK(result.isChaotic);
	CHECK(result.totalTime == 1);
	CHECK(result.qr_refinements > 0);
	CHECK(result.numOrthonormalizations > 1);
	CHECK(result.function_evaluations == flow.derivativeCalls);
	CHECK(flow.jacobianCalls == flow.derivativeCalls);
}

TEST_CASE("Continuous adaptive QR preserves the complete trajectory and nonautonomous clock", "[LyapunovContracts]") {
	LyapunovLinearFlow<1> flow(Matrix<Real>(1, 1, {-100}));
	Real maximumStateError = 0, maximumTime = 0;
	flow.derivativeOverride = [&](Real t, const Vector<Real>& state, Vector<Real>& derivative) {
		maximumTime = std::max(maximumTime, t);
		maximumStateError = std::max(maximumStateError, std::abs(state[0] - t));
		derivative[0] = -100 * (state[0] - t) + 1;
	};
	const auto trajectory = LyapunovAnalyzer::Compute(flow, Vector<Real>{0}, 1, 1, ResolvedStrongFlowStep());
	CHECK(IsWithinAbs(trajectory.maxExponent, Real(-100), Precision::CoarseTolerance));
	CHECK(IsWithinAbs(maximumTime, Real(1), Precision::EqualityTolerance));
	CHECK(maximumStateError < Precision::CoarseTolerance);
	CHECK(trajectory.qr_refinements > 0);
	CHECK(trajectory.totalTime == 1);

	flow.derivativeOverride = [](Real t, const Vector<Real>& x, Vector<Real>& derivative) { derivative[0] = -(100 + 20 * t) * x[0]; };
	flow.jacobianOverride = [](Real t, const Vector<Real>&, Matrix<Real>& matrix) { matrix(0, 0) = -(100 + 20 * t); };
	const auto clock = LyapunovAnalyzer::Compute(flow, Vector<Real>{0}, 1, 1, ResolvedStrongFlowStep());
	CHECK(IsWithinAbs(clock.maxExponent, Real(-110), Precision::CoarseTolerance));
	CHECK(clock.totalTime == 1);
}

TEST_CASE("Continuous QR retains ordinary decay rotation and shear controls", "[LyapunovContracts]") {
	for (const Matrix<Real>& matrix : {
		Matrix<Real>(2, 2, {-1, 0, 0, -2}),
		Matrix<Real>(2, 2, {-1, 1000, 0, -2})}) {
		LyapunovLinearFlow<2> flow(matrix);
		const auto result = LyapunovAnalyzer::Compute(flow, Vector<Real>{0, 0}, 2, Real(0.5));
		CHECK(IsWithinAbs(result.exponents[0], Real(-1), Precision::CoarseTolerance));
		CHECK(IsWithinAbs(result.exponents[1], Real(-2), Precision::CoarseTolerance));
		CHECK(result.qr_refinements == 0);
		CHECK(result.numOrthonormalizations == 4);
	}
	LyapunovLinearFlow<2> rotation(Matrix<Real>(2, 2, {0, -20, 20, 0}));
	const auto result = LyapunovAnalyzer::Compute(rotation, Vector<Real>{0, 0}, 2, Real(0.5));
	for (Real exponent : result.exponents)
		CHECK(IsWithinAbs(exponent, Real(0), Precision::CoarseTolerance));
}

TEST_CASE("Detailed Lyapunov reports actual RHS work including rejected QR trials", "[LyapunovContracts]") {
	LyapunovLinearFlow<1> flow(Matrix<Real>(1, 1, {-100}));
	const auto result = ComputeLyapunovDetailed(flow, Vector<Real>{0}, 1, 1, ResolvedStrongFlowStep());
	REQUIRE(result.IsSuccess());
	CHECK(result.error_message.empty());
	CHECK(result.function_evaluations == flow.derivativeCalls);
	CHECK(result.lyapunov.function_evaluations == result.function_evaluations);
	CHECK(result.function_evaluations > result.lyapunov.numOrthonormalizations);
	CHECK(result.lyapunov.qr_refinements > 0);
	CHECK(IsWithinAbs(result.lyapunov.maxExponent, Real(-100), Precision::CoarseTolerance));
}

TEST_CASE("Continuous Lyapunov validates finite controls and augmented dimensions before work", "[LyapunovContracts]") {
	DynSysConfig config;
	config.exception_policy = EvaluationExceptionPolicy::ConvertToStatus;
	LyapunovLinearFlow<1> flow(Matrix<Real>(1, 1, {-1}));
	for (Real invalid : {Real(0), Real(-1), std::numeric_limits<Real>::quiet_NaN(), std::numeric_limits<Real>::infinity()}) {
		for (const auto& result : {
			ComputeLyapunovDetailed(flow, Vector<Real>{0}, invalid, 1, Real(0.01), config),
			ComputeLyapunovDetailed(flow, Vector<Real>{0}, 1, invalid, Real(0.01), config),
			ComputeLyapunovDetailed(flow, Vector<Real>{0}, 1, 1, invalid, config)}) {
			CHECK(result.status == AlgorithmStatus::InvalidInput);
			CHECK_FALSE(result.error_message.empty());
			CHECK(result.function_evaluations == 0);
		}
	}
	for (const auto& initial : {Vector<Real>{}, Vector<Real>{0, 0}, Vector<Real>{std::numeric_limits<Real>::quiet_NaN()}}) {
		const auto result = ComputeLyapunovDetailed(flow, initial, 1, 1, Real(0.01), config);
		CHECK(result.status == AlgorithmStatus::InvalidInput);
	}
	for (int invalid : {0, -1, std::numeric_limits<int>::max()}) {
		flow.dimension = invalid;
		CHECK(ComputeLyapunovDetailed(flow, Vector<Real>{0}, 1, 1, Real(0.01), config).status == AlgorithmStatus::InvalidInput);
	}
	CHECK(flow.derivativeCalls == 0);
	CHECK(flow.jacobianCalls == 0);
}

TEST_CASE("Continuous Lyapunov fails safely on malformed and nonfinite callbacks", "[LyapunovContracts]") {
	DynSysConfig config;
	config.exception_policy = EvaluationExceptionPolicy::ConvertToStatus;
	LyapunovLinearFlow<1> flow(Matrix<Real>(1, 1, {-1}));
	flow.derivativeOverride = [](Real, const Vector<Real>&, Vector<Real>& derivative) { derivative.resize(0); };
	CHECK(ComputeLyapunovDetailed(flow, Vector<Real>{0}, 1, 1, Real(0.01), config).status == AlgorithmStatus::InvalidInput);
	flow.derivativeOverride = [](Real, const Vector<Real>&, Vector<Real>&) {};
	CHECK(ComputeLyapunovDetailed(flow, Vector<Real>{0}, 1, 1, Real(0.01), config).status == AlgorithmStatus::NumericalInstability);
	flow.derivativeOverride = {};
	flow.jacobianOverride = [](Real, const Vector<Real>&, Matrix<Real>& matrix) { matrix.resize(2, 1); };
	CHECK(ComputeLyapunovDetailed(flow, Vector<Real>{0}, 1, 1, Real(0.01), config).status == AlgorithmStatus::InvalidInput);
	flow.jacobianOverride = [](Real, const Vector<Real>&, Matrix<Real>&) {};
	const auto result = ComputeLyapunovDetailed(flow, Vector<Real>{0}, 1, 1, Real(0.01), config);
	CHECK(result.status == AlgorithmStatus::NumericalInstability);
	CHECK_FALSE(result.error_message.empty());
	CHECK(result.function_evaluations == 1);
	CHECK(result.lyapunov.exponents.empty());
}

TEST_CASE("Continuous QR refinement is bounded and cannot silently reset an unusable frame", "[LyapunovContracts]") {
	DynSysConfig config;
	config.exception_policy = EvaluationExceptionPolicy::ConvertToStatus;
	LyapunovLinearFlow<1> flow(Matrix<Real>(1, 1, {std::numeric_limits<Real>::max()}));
	const auto result = ComputeLyapunovDetailed(flow, Vector<Real>{0}, 1, 1, Real(0.01), config);
	CHECK(result.status == AlgorithmStatus::NumericalInstability);
	CHECK(result.error_message.find("budget") != std::string::npos);
	CHECK(result.function_evaluations == flow.derivativeCalls);
	CHECK(result.function_evaluations > 0);
	CHECK(result.lyapunov.exponents.empty());
}

TEST_CASE("User callback failures are propagated or converted without QR retries", "[LyapunovContracts]") {
	DynSysConfig config;
	config.exception_policy = EvaluationExceptionPolicy::ConvertToStatus;
	LyapunovLinearFlow<1> flow(Matrix<Real>(1, 1, {-1}));
	flow.derivativeOverride = [](Real, const Vector<Real>&, Vector<Real>&) { throw std::runtime_error("User RHS failure"); };
	const auto result = ComputeLyapunovDetailed(flow, Vector<Real>{0}, 1, 1, Real(0.01), config);
	CHECK(result.status == AlgorithmStatus::AlgorithmSpecificFailure);
	CHECK(result.error_message == "User RHS failure");
	CHECK(result.function_evaluations == 1);
	CHECK(flow.derivativeCalls == 1);
	CHECK_THROWS_WITH(ComputeLyapunovDetailed(flow, Vector<Real>{0}, 1), "User RHS failure");
	CHECK(flow.derivativeCalls == 2);
	flow.derivativeOverride = [](Real, const Vector<Real>&, Vector<Real>&) { throw SingularMatrixError("User singular failure"); };
	const auto singular = ComputeLyapunovDetailed(flow, Vector<Real>{0}, 1, 1, Real(0.01), config);
	CHECK(singular.status == AlgorithmStatus::SingularMatrix);
	CHECK(singular.error_message == "User singular failure");
	CHECK(singular.function_evaluations == 1);
	CHECK(flow.derivativeCalls == 3);
}
