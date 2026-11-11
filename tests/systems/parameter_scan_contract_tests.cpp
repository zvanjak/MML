#include "../TestPrecision.h"
#ifndef MML_SINGLE_HEADER_INCLUDE
#include <mml/systems/DynamicalSystem.h>
#endif
#include <catch2/catch_all.hpp>
#include <cmath>
#include <functional>
#include <limits>
#include <numeric>
#include <stdexcept>

using namespace MML;
using namespace MML::Systems;

namespace {
	class ParameterScanFixture : public DynamicalSystemBase<1, 1> {
	public:
		mutable int calls = 0;
		int dimension = 1;
		Real velocity = 0;
		std::vector<Real> writes;
		std::function<void(Real)> onWrite;
		std::function<void(const Vector<Real>&)> onDerivative;
		ParameterScanFixture() { _params[0] = 7; }
		int getDim() const override { return dimension; }
		void SetOriginal(Real value) { _params[0] = value; }
		void setParam(int index, Real value) override {
			writes.push_back(value);
			if (onWrite) onWrite(value);
			if (!std::isfinite(value)) throw ArgumentError("nonfinite parameter");
			DynamicalSystemBase<1, 1>::setParam(index, value);
		}
		void derivs(Real, const Vector<Real>& state, Vector<Real>& derivative) const override {
			++calls;
			if (onDerivative) onDerivative(state);
			derivative[0] = velocity;
		}
		void jacobian(Real, const Vector<Real>&, Matrix<Real>& matrix) const override {
			matrix = Matrix<Real>(1, 1, {0});
		}
	};

	class GridFixture : public DynamicalSystemBase<2, 0> {
	public:
		mutable int calls = 0;
		int dimension = 2;
		bool autonomous = true;
		int getDim() const override { return dimension; }
		bool isAutonomous() const override { return autonomous; }
		void derivs(Real, const Vector<Real>&, Vector<Real>& derivative) const override {
			++calls;
			derivative[0] = derivative[1] = 0;
		}
		void jacobian(Real, const Vector<Real>&, Matrix<Real>& matrix) const override {
			matrix = Matrix<Real>(2, 2, {0, 0, 0, 0});
		}
	};

	std::string ExceptionMessage(std::exception_ptr error) {
		REQUIRE(error);
		try { std::rethrow_exception(error); }
		catch (const std::runtime_error& exception) { return exception.what(); }
	}
}

TEST_CASE("Chaos parameter scans evaluate a single sample at the minimum", "[ParameterScanContracts][ScanBaseline]") {
	ParameterScanFixture flow;
	const auto result = DynamicalSystemAnalyzer<>(flow).ScanChaosVsParameter(0, 1, 2, 1, Vector<Real>{0}, Real(0.02));
	REQUIRE(result.size() == 1);
	CHECK(result[0].first == 1);
	CHECK(result[0].second == 0);
	CHECK(flow.getParam(0) == 7);
}

TEST_CASE("Chaos parameter scans restore after a later callback fails", "[ParameterScanContracts][ScanBaseline]") {
	ParameterScanFixture flow;
	flow.onDerivative = [&](const Vector<Real>&) {
		if (flow.getParam(0) == 2) throw std::runtime_error("operation sentinel");
	};
	CHECK_THROWS_WITH(DynamicalSystemAnalyzer<>(flow).ScanChaosVsParameter(0, 1, 2, 2, Vector<Real>{0}, Real(0.02)),
		"operation sentinel");
	CHECK(flow.getParam(0) == 7);
}

TEST_CASE("Bifurcation supports opposite extreme finite endpoints", "[ParameterScanContracts][ScanBaseline]") {
	ParameterScanFixture flow;
	const Real bound = std::numeric_limits<Real>::max();
	const auto result = BifurcationAnalyzer::Sweep(flow, 0, -bound, bound, 3, Vector<Real>{0}, 0, 0, 0);
	REQUIRE(result.parameterValues.size() == 3);
	CHECK(result.parameterValues == std::vector<Real>{-bound, 0, bound});
	CHECK(flow.calls == 0);
	CHECK(flow.getParam(0) == 7);
}

TEST_CASE("Box grids support opposite extreme finite bounds", "[ParameterScanContracts][ScanBaseline]") {
	ParameterScanFixture flow;
	const Real bound = std::numeric_limits<Real>::max();
	const auto result = DynamicalSystemAnalyzer<>(flow).FindFixedPointsInBox(Vector<Real>{-bound}, Vector<Real>{bound}, 3);
	REQUIRE(result.size() == 3);
	CHECK(result[0].location[0] == -bound);
	CHECK(result[1].location[0] == 0);
	CHECK(result[2].location[0] == bound);
}

TEST_CASE("Chaos scans reject invalid requests before mutation or RHS work", "[ParameterScanContracts]") {
	ParameterScanFixture flow;
	DynamicalSystemAnalyzer<> analyzer(flow);
	for (int count : {0, -1})
		CHECK_THROWS_AS(analyzer.ScanChaosVsParameter(0, 1, 2, count, Vector<Real>{0}), ArgumentError);
	for (int index : {-1, 1})
		CHECK_THROWS_AS(analyzer.ScanChaosVsParameter(index, 1, 2, 2, Vector<Real>{0}), ArgumentError);
	for (Real invalid : {std::numeric_limits<Real>::quiet_NaN(), std::numeric_limits<Real>::infinity()}) {
		CHECK_THROWS_AS(analyzer.ScanChaosVsParameter(0, invalid, 2, 2, Vector<Real>{0}), ArgumentError);
		CHECK_THROWS_AS(analyzer.ScanChaosVsParameter(0, 1, invalid, 2, Vector<Real>{0}), ArgumentError);
		CHECK_THROWS_AS(analyzer.ScanChaosVsParameter(0, 1, 2, 2, Vector<Real>{invalid}), ArgumentError);
	}
	for (Real time : {Real(0), Real(-1), std::numeric_limits<Real>::quiet_NaN(), std::numeric_limits<Real>::infinity()})
		CHECK_THROWS_AS(analyzer.ScanChaosVsParameter(0, 1, 2, 2, Vector<Real>{0}, time), ArgumentError);
	CHECK_THROWS_AS(analyzer.ScanChaosVsParameter(0, 1, 2, 2, Vector<Real>{}), VectorDimensionError);
	CHECK(flow.writes.empty());
	CHECK(flow.calls == 0);
	CHECK(flow.getParam(0) == 7);
}

TEST_CASE("Chaos scans preflight augmented-dimension capacity", "[ParameterScanContracts]") {
	ParameterScanFixture flow;
	flow.dimension = 50000;
	CHECK_THROWS_AS(DynamicalSystemAnalyzer<>(flow).ScanChaosVsParameter(0, 1, 2, 2, Vector<Real>(flow.dimension, Real(0))), ArgumentError);
	CHECK(flow.writes.empty());
	CHECK(flow.calls == 0);
}

TEST_CASE("Scans require a finite original parameter before mutation", "[ParameterScanContracts]") {
	ParameterScanFixture flow;
	flow.SetOriginal(std::numeric_limits<Real>::quiet_NaN());
	CHECK_THROWS_AS(DynamicalSystemAnalyzer<>(flow).ScanChaosVsParameter(0, 1, 2, 2, Vector<Real>{0}), ArgumentError);
	CHECK_THROWS_AS(BifurcationAnalyzer::Sweep(flow, 0, 1, 2, 2, Vector<Real>{0}, 0, 0, 0), ArgumentError);
	CHECK(flow.writes.empty());
	CHECK(flow.calls == 0);
}

TEST_CASE("Chaos scans preserve constant ranges and exact endpoint order", "[ParameterScanContracts]") {
	ParameterScanFixture flow;
	DynamicalSystemAnalyzer<> analyzer(flow);
	const Real bound = std::numeric_limits<Real>::max();
	for (const auto& range : {std::pair<Real, Real>{-bound, bound}, {bound, -bound},
		{bound / 2, bound}, {-bound / 2, -bound}, {Real(1), Real(1)}}) {
		const auto result = analyzer.ScanChaosVsParameter(0, range.first, range.second, 3, Vector<Real>{0}, Real(0.02));
		REQUIRE(result.size() == 3);
		CHECK(result.front().first == range.first);
		CHECK(result.back().first == range.second);
		CHECK(result[1].first == std::midpoint(range.first, range.second));
		for (const auto& sample : result) {
			CHECK(std::isfinite(sample.first));
			CHECK(sample.first >= std::min(range.first, range.second));
			CHECK(sample.first <= std::max(range.first, range.second));
			CHECK(sample.second == 0);
		}
		CHECK(flow.getParam(0) == 7);
	}
}

TEST_CASE("Bifurcation ranges preserve all finite endpoints and remain monotonic", "[ParameterScanContracts]") {
	ParameterScanFixture flow;
	const Real bound = std::numeric_limits<Real>::max();
	for (const auto& range : {std::pair<Real, Real>{-bound, bound}, {bound, -bound},
		{bound / 2, bound}, {-bound / 2, -bound}}) {
		const auto result = BifurcationAnalyzer::Sweep(flow, 0, range.first, range.second, 7, Vector<Real>{0}, 0, 0, 0);
		REQUIRE(result.parameterValues.size() == 7);
		CHECK(result.parameterValues.front() == range.first);
		CHECK(result.parameterValues.back() == range.second);
		for (std::size_t index = 0; index < result.parameterValues.size(); ++index) {
			CHECK(std::isfinite(result.parameterValues[index]));
			if (index > 0) CHECK((result.parameterValues[index] >= result.parameterValues[index - 1]) == (range.second > range.first));
		}
		CHECK(flow.getParam(0) == 7);
	}
}

TEST_CASE("Chaos samples restart from the supplied initial state", "[ParameterScanContracts]") {
	ParameterScanFixture flow;
	flow.velocity = 1;
	std::vector<Real> starts;
	Real previous = 7;
	flow.onDerivative = [&](const Vector<Real>& state) {
		if (flow.getParam(0) != previous) { starts.push_back(state[0]); previous = flow.getParam(0); }
	};
	DynamicalSystemAnalyzer<>(flow).ScanChaosVsParameter(0, 1, 3, 3, Vector<Real>{3}, Real(0.02));
	CHECK(starts == std::vector<Real>{3, 3, 3});
	CHECK(flow.getParam(0) == 7);
}

TEST_CASE("Bifurcation samples preserve continuation through both phases", "[ParameterScanContracts]") {
	ParameterScanFixture flow;
	flow.velocity = 1;
	std::vector<Real> starts;
	Real previous = 7;
	flow.onDerivative = [&](const Vector<Real>& state) {
		if (flow.getParam(0) != previous) { starts.push_back(state[0]); previous = flow.getParam(0); }
	};
	BifurcationAnalyzer::Sweep(flow, 0, 1, 3, 3, Vector<Real>{3}, 0, Real(0.02), Real(0.02));
	REQUIRE(starts.size() == 3);
	for (int index = 0; index < 3; ++index)
		CHECK(IsWithinAbs(starts[index], Real(3) + Real(0.04) * index, Precision::EqualityTolerance));
	CHECK(flow.getParam(0) == 7);
}

TEST_CASE("Parameter setters that mutate then throw are followed by restoration", "[ParameterScanContracts]") {
	ParameterScanFixture flow;
	flow.onWrite = [&](Real value) {
		if (value == 2) { flow.SetOriginal(value); throw std::runtime_error("setter sentinel"); }
	};
	SECTION("chaos") {
		CHECK_THROWS_WITH(DynamicalSystemAnalyzer<>(flow).ScanChaosVsParameter(0, 1, 2, 2, Vector<Real>{0}, Real(0.02)), "setter sentinel");
	}
	SECTION("bifurcation") {
		CHECK_THROWS_WITH(BifurcationAnalyzer::Sweep(flow, 0, 1, 2, 2, Vector<Real>{0}, 0, 0, 0), "setter sentinel");
	}
	CHECK(flow.getParam(0) == 7);
	REQUIRE(flow.writes.size() == 3);
	CHECK(flow.writes.back() == 7);
}

TEST_CASE("Restoration errors after successful scans are explicit and retain their cause", "[ParameterScanContracts]") {
	ParameterScanFixture flow;
	flow.onWrite = [](Real value) { if (value == 7) throw std::runtime_error("restoration sentinel"); };
	const auto check = [&](auto&& compute) {
		try { compute(); FAIL("restoration failure must not return success"); }
		catch (const ParameterRestorationError& error) {
			CHECK_FALSE(error.OperationException());
			CHECK(ExceptionMessage(error.RestorationException()) == "restoration sentinel");
			CHECK(std::string(error.what()).find("model state is unknown") != std::string::npos);
		}
	};
	SECTION("chaos") {
		check([&] { DynamicalSystemAnalyzer<>(flow).ScanChaosVsParameter(0, 1, 2, 2, Vector<Real>{0}, Real(0.02)); });
	}
	SECTION("bifurcation") {
		check([&] { BifurcationAnalyzer::Sweep(flow, 0, 1, 2, 2, Vector<Real>{0}, 0, 0, 0); });
	}
	CHECK(flow.getParam(0) == 2);
	CHECK(flow.writes.size() == 3);
}

TEST_CASE("Simultaneous operation and restoration errors retain both original exceptions", "[ParameterScanContracts]") {
	ParameterScanFixture flow;
	flow.onWrite = [](Real value) { if (value == 7) throw std::runtime_error("restoration sentinel"); };
	flow.onDerivative = [](const Vector<Real>&) { throw std::runtime_error("operation sentinel"); };
	const auto check = [&](auto&& compute) {
		try { compute(); FAIL("dual failure must not return success"); }
		catch (const ParameterRestorationError& error) {
			CHECK(ExceptionMessage(error.OperationException()) == "operation sentinel");
			CHECK(ExceptionMessage(error.RestorationException()) == "restoration sentinel");
			CHECK(std::string(error.what()).find("operation sentinel") != std::string::npos);
			CHECK(std::string(error.what()).find("restoration sentinel") != std::string::npos);
		}
	};
	SECTION("chaos") {
		check([&] { DynamicalSystemAnalyzer<>(flow).ScanChaosVsParameter(0, 1, 2, 2, Vector<Real>{0}, Real(0.02)); });
	}
	SECTION("bifurcation") {
		check([&] { BifurcationAnalyzer::Sweep(flow, 0, 1, 2, 2, Vector<Real>{0}, 0, Real(0.02), 0); });
	}
	CHECK(flow.calls == 1);
	CHECK(flow.writes == std::vector<Real>{1, 7});
}

TEST_CASE("Bifurcation Detailed reports restoration failures with actual work", "[ParameterScanContracts]") {
	ParameterScanFixture flow;
	flow.onWrite = [](Real value) { if (value == 7) throw std::runtime_error("restoration sentinel"); };
	flow.onDerivative = [](const Vector<Real>&) { throw std::runtime_error("operation sentinel"); };
	DynSysConfig config;
	config.exception_policy = EvaluationExceptionPolicy::ConvertToStatus;
	const auto result = SweepBifurcationDetailed(flow, 0, 1, 2, 2, Vector<Real>{0}, 0, Real(0.02), 0, Real(0.01), config);
	CHECK(result.status == AlgorithmStatus::AlgorithmSpecificFailure);
	CHECK(result.function_evaluations == 1);
	CHECK(result.function_evaluations == flow.calls);
	CHECK(result.error_message.find("operation sentinel") != std::string::npos);
	CHECK(result.error_message.find("restoration sentinel") != std::string::npos);
}

TEST_CASE("Non-standard callback exceptions are restored and rethrown unchanged", "[ParameterScanContracts]") {
	ParameterScanFixture flow;
	flow.onDerivative = [](const Vector<Real>&) { throw 42; };
	CHECK_THROWS_AS(DynamicalSystemAnalyzer<>(flow).ScanChaosVsParameter(0, 1, 2, 2, Vector<Real>{0}, Real(0.02)), int);
	CHECK(flow.getParam(0) == 7);
}

TEST_CASE("Nested scan restoration failures do not bypass the outer restoration", "[ParameterScanContracts]") {
	ParameterScanFixture outer, inner;
	inner.onWrite = [](Real value) { if (value == 7) throw std::runtime_error("inner restore"); };
	outer.onDerivative = [&](const Vector<Real>&) {
		BifurcationAnalyzer::Sweep(inner, 0, 1, 2, 2, Vector<Real>{0}, 0, 0, 0);
	};
	CHECK_THROWS_AS(DynamicalSystemAnalyzer<>(outer).ScanChaosVsParameter(0, 1, 2, 2, Vector<Real>{0}, Real(0.02)),
		ParameterRestorationError);
	CHECK(outer.getParam(0) == 7);
	CHECK(inner.getParam(0) == 2);
}

TEST_CASE("Dual non-standard failures remain inspectable without termination", "[ParameterScanContracts]") {
	ParameterScanFixture flow;
	flow.onWrite = [](Real value) { if (value == 7) throw 17; };
	flow.onDerivative = [](const Vector<Real>&) { throw 42; };
	try {
		DynamicalSystemAnalyzer<>(flow).ScanChaosVsParameter(0, 1, 2, 2, Vector<Real>{0}, Real(0.02));
		FAIL("dual failure must not return success");
	} catch (const ParameterRestorationError& error) {
		try { std::rethrow_exception(error.OperationException()); }
		catch (int value) { CHECK(value == 42); }
		try { std::rethrow_exception(error.RestorationException()); }
		catch (int value) { CHECK(value == 17); }
		CHECK(std::string(error.what()).find("non-standard exception") != std::string::npos);
	}
	CHECK(flow.writes == std::vector<Real>{1, 7});
}

TEST_CASE("Zero-work Detailed sweeps do not disguise restoration failure as success", "[ParameterScanContracts]") {
	ParameterScanFixture flow;
	flow.onWrite = [](Real value) { if (value == 7) throw std::runtime_error("restore refused"); };
	DynSysConfig config;
	config.exception_policy = EvaluationExceptionPolicy::ConvertToStatus;
	const auto result = SweepBifurcationDetailed(flow, 0, 1, 2, 2, Vector<Real>{0}, 0, 0, 0, Real(0.01), config);
	CHECK_FALSE(result.IsSuccess());
	CHECK(result.status == AlgorithmStatus::AlgorithmSpecificFailure);
	CHECK(result.function_evaluations == 0);
	CHECK(flow.calls == 0);
	CHECK(result.error_message.find("restore refused") != std::string::npos);
	CHECK(result.error_message.find("Operation error: none") != std::string::npos);
}

TEST_CASE("Box searches reject malformed bounds before callbacks", "[ParameterScanContracts]") {
	GridFixture flow;
	DynamicalSystemAnalyzer<> analyzer(flow);
	CHECK_THROWS_AS(analyzer.FindFixedPointsInBox(Vector<Real>{0}, Vector<Real>{1, 1}), VectorDimensionError);
	CHECK_THROWS_AS(analyzer.FindFixedPointsInBox(Vector<Real>{0, 0}, Vector<Real>{1}), VectorDimensionError);
	CHECK_THROWS_AS(analyzer.FindFixedPointsInBox(Vector<Real>{0, 2}, Vector<Real>{1, 1}), ArgumentError);
	for (Real invalid : {std::numeric_limits<Real>::quiet_NaN(), std::numeric_limits<Real>::infinity()}) {
		CHECK_THROWS_AS(analyzer.FindFixedPointsInBox(Vector<Real>{invalid, 0}, Vector<Real>{1, 1}), ArgumentError);
		CHECK_THROWS_AS(analyzer.FindFixedPointsInBox(Vector<Real>{0, 0}, Vector<Real>{invalid, 1}), ArgumentError);
	}
	CHECK(flow.calls == 0);
}

TEST_CASE("Box searches reject invalid grid sizes and budgets before work", "[ParameterScanContracts]") {
	GridFixture flow;
	DynamicalSystemAnalyzer<> analyzer(flow);
	for (int invalid : {0, -1}) {
		CHECK_THROWS_AS(analyzer.FindFixedPointsInBox(Vector<Real>{0, 0}, Vector<Real>{1, 1}, invalid), ArgumentError);
		CHECK_THROWS_AS(analyzer.FindFixedPointsInBox(Vector<Real>{0, 0}, Vector<Real>{1, 1}, 2, invalid), ArgumentError);
	}
	CHECK_THROWS_AS(analyzer.FindFixedPointsInBox(Vector<Real>{0, 0}, Vector<Real>{1, 1}, std::numeric_limits<int>::max()), ResourceLimitError);
	CHECK_THROWS_AS(analyzer.FindFixedPointsInBox(Vector<Real>{0, 0}, Vector<Real>{1, 1}, 317), ResourceLimitError);
	CHECK(flow.calls == 0);
}

TEST_CASE("Box grid cardinality checks exact caller budget boundaries", "[ParameterScanContracts]") {
	GridFixture flow;
	DynamicalSystemAnalyzer<> analyzer(flow);
	CHECK_THROWS_AS(analyzer.FindFixedPointsInBox(Vector<Real>{-1, -1}, Vector<Real>{1, 1}, 2, 3), ResourceLimitError);
	CHECK(flow.calls == 0);
	const auto result = analyzer.FindFixedPointsInBox(Vector<Real>{-1, -1}, Vector<Real>{1, 1}, 2, 4);
	CHECK(result.size() == 4);
	CHECK(flow.calls == 4);
}

TEST_CASE("One-point box grids retain the midpoint and allow degenerate axes", "[ParameterScanContracts]") {
	GridFixture flow;
	const Real bound = std::numeric_limits<Real>::max();
	const auto result = DynamicalSystemAnalyzer<>(flow).FindFixedPointsInBox(Vector<Real>{-bound, 2}, Vector<Real>{bound, 2}, 1, 1);
	REQUIRE(result.size() == 1);
	CHECK(result.front().location == Vector<Real>{0, 2});
	CHECK(flow.calls == 1);
}

TEST_CASE("Box search rejects nonautonomous and nonpositive-dimensional systems before work", "[ParameterScanContracts]") {
	GridFixture flow;
	flow.autonomous = false;
	CHECK_THROWS_AS(DynamicalSystemAnalyzer<>(flow).FindFixedPointsInBox(Vector<Real>{0, 0}, Vector<Real>{1, 1}), ArgumentError);
	flow.autonomous = true;
	for (int dimension : {0, -1}) {
		flow.dimension = dimension;
		CHECK_THROWS_AS(DynamicalSystemAnalyzer<>(flow).FindFixedPointsInBox(Vector<Real>{}, Vector<Real>{}), ArgumentError);
	}
	CHECK(flow.calls == 0);
}
