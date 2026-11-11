#include "../TestPrecision.h"
#ifndef MML_SINGLE_HEADER_INCLUDE
#include <mml/systems/DynamicalSystem.h>
#endif
#include <catch2/catch_all.hpp>
#include <cmath>
#include <functional>
#include <limits>
#include <stdexcept>

using namespace MML;
using namespace MML::Systems;

namespace {
	class ClockFlow : public DynamicalSystemBase<1, 1> {
		mutable bool _firstCallback = false;
	public:
		bool cosine = false;
		mutable int calls = 0;
		int writes = 0;
		mutable std::vector<std::pair<Real, Real>> starts;
		std::function<void(Real)> callback;
		ClockFlow() { _params[0] = 7; }
		bool isAutonomous() const override { return false; }
		void setParam(int index, Real value) override {
			++writes;
			_firstCallback = true;
			DynamicalSystemBase<1, 1>::setParam(index, value);
		}
		void derivs(Real time, const Vector<Real>& state, Vector<Real>& derivative) const override {
			++calls;
			if (_firstCallback) {
				starts.push_back({time, state[0]});
				_firstCallback = false;
			}
			if (callback) callback(time);
			derivative[0] = cosine ? std::cos(time) : time;
		}
	};

	void CheckPolynomialStarts(const ClockFlow& flow, Real sampleDuration) {
		REQUIRE(flow.starts.size() == 3);
		for (int index = 0; index < 3; ++index) {
			const Real time = index * sampleDuration;
			CHECK(flow.starts[index].first == time);
			CHECK(IsWithinAbs(flow.starts[index].second, time * time / 2, Precision::CoarseTolerance));
		}
		CHECK(flow.getParam(0) == 7);
	}
}

TEST_CASE("Bifurcation recording continues the transient forcing phase", "[BifurcationClockContracts][ClockBaseline]") {
	ClockFlow flow;
	flow.cosine = true;
	const auto diagram = BifurcationAnalyzer::Sweep(flow, 0, 1, 2, 1, Vector<Real>{0}, 0,
		Constants::PI / 2, Constants::PI, Real(0.02));
	REQUIRE(diagram.attractorValues.size() == 1);
	CHECK(diagram.attractorValues.front().empty());
	CHECK(flow.getParam(0) == 7);
}

TEST_CASE("Bifurcation continuation carries state and time through both phases", "[BifurcationClockContracts][ClockBaseline]") {
	ClockFlow flow;
	BifurcationAnalyzer::Sweep(flow, 0, 1, 3, 3, Vector<Real>{0}, 0, 1, 1, Real(0.1));
	CheckPolynomialStarts(flow, 2);
}

TEST_CASE("Zero-transient bifurcation sweeps still carry the recording clock", "[BifurcationClockContracts][ClockBaseline]") {
	ClockFlow flow;
	BifurcationAnalyzer::Sweep(flow, 0, 1, 3, 3, Vector<Real>{0}, 0, 0, 1, Real(0.1));
	CheckPolynomialStarts(flow, 1);
}

TEST_CASE("Zero-record bifurcation sweeps still carry the transient clock", "[BifurcationClockContracts][ClockBaseline]") {
	ClockFlow flow;
	BifurcationAnalyzer::Sweep(flow, 0, 1, 3, 3, Vector<Real>{0}, 0, 1, 0, Real(0.1));
	CheckPolynomialStarts(flow, 1);
}

TEST_CASE("Bifurcation continuation does not replay periodic maxima at each parameter", "[BifurcationClockContracts]") {
	ClockFlow flow;
	flow.cosine = true;
	const auto diagram = BifurcationAnalyzer::Sweep(flow, 0, 1, 3, 3, Vector<Real>{0}, 0,
		0, Constants::PI, Real(0.02));
	REQUIRE(diagram.attractorValues.size() == 3);
	REQUIRE(diagram.attractorValues[0].size() == 1);
	CHECK(diagram.attractorValues[1].empty());
	REQUIRE(diagram.attractorValues[2].size() == 1);
	CHECK(IsWithinAbs(diagram.attractorValues[0][0], Real(1), Precision::CoarseTolerance));
	CHECK(IsWithinAbs(diagram.attractorValues[2][0], Real(1), Precision::CoarseTolerance));
	CHECK(flow.getParam(0) == 7);
}

TEST_CASE("Bifurcation Detailed preserves physical time with optional checks enabled or disabled", "[BifurcationClockContracts]") {
	for (bool checkFinite : {true, false}) {
		ClockFlow flow;
		DynSysConfig config;
		config.check_finite = checkFinite;
		const auto result = SweepBifurcationDetailed(flow, 0, 1, 3, 3, Vector<Real>{0}, 0, 1, 1, Real(0.1), config);
		REQUIRE(result.IsSuccess());
		CHECK(result.function_evaluations == flow.calls);
		CHECK(result.function_evaluations > 0);
		CheckPolynomialStarts(flow, 2);
	}
}

TEST_CASE("Bifurcation RHS clocks remain inside each sample physical interval", "[BifurcationClockContracts]") {
	ClockFlow flow;
	int zeroTimes = 0;
	flow.callback = [&](Real time) {
		const Real start = 2 * (flow.getParam(0) - 1);
		CHECK(time >= start);
		CHECK(time <= start + 2);
		if (time == 0) ++zeroTimes;
	};
	BifurcationAnalyzer::Sweep(flow, 0, 1, 3, 3, Vector<Real>{0}, 0, 1, 1, Real(0.1));
	CHECK(zeroTimes == 1);
	CHECK(flow.getParam(0) == 7);
}

TEST_CASE("Overflowing bifurcation phase schedules fail before parameter writes or RHS work", "[BifurcationClockContracts]") {
	ClockFlow flow;
	DynSysConfig config;
	config.exception_policy = EvaluationExceptionPolicy::ConvertToStatus;
	const Real large = std::numeric_limits<Real>::max();
	SECTION("within the first sample") {
		const auto result = SweepBifurcationDetailed(flow, 0, 1, 2, 1, Vector<Real>{0}, 0, large, large, large, config);
		CHECK(result.status == AlgorithmStatus::InvalidInput);
		CHECK_FALSE(result.error_message.empty());
		CHECK(result.function_evaluations == 0);
	}
	SECTION("across samples") {
		CHECK_THROWS_AS(BifurcationAnalyzer::Sweep(flow, 0, 1, 2, 2, Vector<Real>{0}, 0,
			large * Real(0.75), 0, large), ArgumentError);
	}
	CHECK(flow.calls == 0);
	CHECK(flow.writes == 0);
	CHECK(flow.getParam(0) == 7);
}

TEST_CASE("Unresolvable positive phase durations are rejected before any sweep mutation", "[BifurcationClockContracts]") {
	ClockFlow flow;
	const Real tiny = std::numeric_limits<Real>::denorm_min();
	SECTION("recording cannot advance after the transient") {
		CHECK_THROWS_AS(BifurcationAnalyzer::Sweep(flow, 0, 1, 2, 2, Vector<Real>{0}, 0, 1, tiny, 1), ArgumentError);
	}
	SECTION("a later transient cannot advance the continued clock") {
		CHECK_THROWS_AS(BifurcationAnalyzer::Sweep(flow, 0, 1, 2, 2, Vector<Real>{0}, 0, tiny, 1, 1), ArgumentError);
	}
	CHECK(flow.calls == 0);
	CHECK(flow.writes == 0);
	CHECK(flow.getParam(0) == 7);
}

TEST_CASE("Later physical-time callback failures preserve actual work and restore the parameter", "[BifurcationClockContracts]") {
	ClockFlow flow;
	flow.callback = [&](Real time) {
		if (flow.getParam(0) == 2 && time >= 2) throw std::runtime_error("continued clock sentinel");
	};
	DynSysConfig config;
	config.exception_policy = EvaluationExceptionPolicy::ConvertToStatus;
	const auto result = SweepBifurcationDetailed(flow, 0, 1, 3, 3, Vector<Real>{0}, 0, 1, 1, Real(0.1), config);
	CHECK(result.status == AlgorithmStatus::AlgorithmSpecificFailure);
	CHECK(result.error_message == "continued clock sentinel");
	CHECK(result.function_evaluations == flow.calls);
	CHECK(flow.calls > 1);
	REQUIRE(flow.starts.size() == 2);
	CHECK(flow.starts.back().first == 2);
	CHECK(flow.getParam(0) == 7);
}

TEST_CASE("Zero-duration sweeps preserve zero work and empty samples", "[BifurcationClockContracts]") {
	ClockFlow flow;
	const auto result = SweepBifurcationDetailed(flow, 0, 1, 3, 3, Vector<Real>{0}, 0, 0, 0);
	REQUIRE(result.IsSuccess());
	REQUIRE(result.diagram.attractorValues.size() == 3);
	for (const auto& values : result.diagram.attractorValues) CHECK(values.empty());
	CHECK(result.function_evaluations == 0);
	CHECK(flow.calls == 0);
	CHECK(flow.writes == 4);
	CHECK(result.diagram.numTransientSteps == 0);
	CHECK(result.diagram.numRecordedPoints == 0);
	CHECK(flow.getParam(0) == 7);
}

TEST_CASE("Separate bifurcation invocations each start at physical time zero", "[BifurcationClockContracts]") {
	ClockFlow flow;
	BifurcationAnalyzer::Sweep(flow, 0, 1, 2, 1, Vector<Real>{0}, 0, 1, 0);
	BifurcationAnalyzer::Sweep(flow, 0, 1, 2, 1, Vector<Real>{0}, 0, 1, 0);
	REQUIRE(flow.starts.size() == 2);
	CHECK(flow.starts[0] == std::pair<Real, Real>{Real(0), Real(0)});
	CHECK(flow.starts[1] == std::pair<Real, Real>{Real(0), Real(0)});
	CHECK(flow.getParam(0) == 7);
}

TEST_CASE("Lifted Duffing systems remain accepted by bifurcation sweeps", "[BifurcationClockContracts]") {
	DuffingSystem flow;
	const Real original = flow.getParam(3);
	const auto result = SweepBifurcationDetailed(flow, 3, Real(0.1), Real(0.2), 2,
		Vector<Real>{0, 0, 0}, 0, Real(0.01), Real(0.02));
	CHECK(result.IsSuccess());
	CHECK(result.diagram.parameterValues.size() == 2);
	CHECK(result.function_evaluations > 0);
	CHECK(flow.getParam(3) == original);
}

TEST_CASE("The analyzer bifurcation facade shares the corrected forcing phase", "[BifurcationClockContracts]") {
	ClockFlow flow;
	flow.cosine = true;
	const auto result = DynamicalSystemAnalyzer<>(flow).ComputeBifurcation(
		0, 1, 2, 1, Vector<Real>{0}, 0, Constants::PI / 2, Constants::PI);
	REQUIRE(result.attractorValues.size() == 1);
	CHECK(result.attractorValues.front().empty());
	CHECK(flow.getParam(0) == 7);
}

TEST_CASE("Legacy nominal sweep metadata remains per phase rather than cumulative", "[BifurcationClockContracts]") {
	ClockFlow flow;
	const auto result = SweepBifurcationDetailed(flow, 0, 1, 3, 3, Vector<Real>{0}, 0, Real(0.5), Real(0.5), Real(0.1));
	REQUIRE(result.IsSuccess());
	CHECK(result.diagram.numTransientSteps == 5);
	CHECK(result.diagram.numRecordedPoints == 5);
	CheckPolynomialStarts(flow, 1);
}
